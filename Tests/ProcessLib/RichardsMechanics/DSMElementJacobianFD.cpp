/**
 * \file
 * \copyright
 * Copyright (c) 2012-2026, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Verification of the assembled Newton tangent of the RichardsMechanics DSM
 * element against a central finite difference of the assembled residual,
 * block by block (pp, pu, up, uu). Report:
 * docs/dsm_jacobian_fd_verification_2026-10.md.
 *
 * Harness. The element is the single axisymmetric 1 cm quad of
 * Tests/Data/RichardsMechanics/square_1e-2_quad_1e0.vtu, set up by reading the
 * deck ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj through the regular project
 * reader (ProjectData), so the local assembler, the material set and every
 * potential_exchange switch are those of the deck. The switch-set variants
 * (swelling_stress_form = level; + micro_potential_partner_z_weight = true)
 * and the strained variant (top u_z = -1e-4 m, eps_v = -0.01) are written as
 * patched copies of the deck into a temporary directory; the deck itself is
 * not touched. The time loop of the deck drives the element to the
 * representative states; at each state the local assembler of the element is
 * called directly with x_prev = the last accepted solution and x = the first
 * Newton iterate of the next step (x_prev with the Dirichlet values of the
 * pressure ramp at t + dt; and, for the "strain step" state, a 1 %
 * compressive volumetric strain added in x only).
 *
 * The integration-point states (current_states_, output_data_) are snapshot
 * before the base evaluation and restored before every perturbed evaluation,
 * so the FD is that of the residual R(x; S0) at the fixed lagged state S0 the
 * Newton iteration starts from (the same convention as the in-assembler
 * micro_ceiling_fd_check). A second "relag" FD evaluates the perturbed
 * residual n_relag + 1 times in sequence from S0 (the lagged values converge
 * to those of the perturbed point), to separate lag effects from tangent
 * terms.
 *
 * Access to the time loop solution, the next time increment, the local
 * assemblers and their integration-point states uses the explicit-
 * instantiation access rule ([temp.spec.general]/6) instead of changing the
 * library classes: no existing source file is modified.
 *
 * Environment: OGS_DSM_JACFD_CSV = path of the CSV written (default: none,
 * the table is printed to stdout); OGS_DSM_JACFD_MATRICES = directory to
 * write the analytic and FD matrices of every case (default: none).
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Applications/ApplicationsLib/ProjectData.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/Logging.h"
#include "BaseLib/PrjProcessing.h"
#include "InfoLib/TestInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/RichardsMechanics/LocalAssemblerInterface.h"
#include "ProcessLib/RichardsMechanics/RichardsMechanicsProcess.h"
#include "ProcessLib/TimeLoop.h"

namespace
{
namespace RM = ProcessLib::RichardsMechanics;

// ── access to private / protected members by explicit instantiation ─────────
template <typename Tag, typename Tag::type M>
struct Access
{
    friend typename Tag::type member(Tag) { return M; }
};

struct TimeLoopSolutions
{
    using type = std::vector<GlobalVector*> ProcessLib::TimeLoop::*;
    friend type member(TimeLoopSolutions);
};
template struct Access<TimeLoopSolutions,
                       &ProcessLib::TimeLoop::_process_solutions>;

struct TimeLoopDt
{
    using type = NumLib::TimeIncrement ProcessLib::TimeLoop::*;
    friend type member(TimeLoopDt);
};
template struct Access<TimeLoopDt, &ProcessLib::TimeLoop::_dt>;

struct RMLocalAssemblers
{
    using type = std::vector<std::unique_ptr<RM::LocalAssemblerInterface<2>>>
        RM::RichardsMechanicsProcess<2>::*;
    friend type member(RMLocalAssemblers);
};
template struct Access<RMLocalAssemblers,
                       &RM::RichardsMechanicsProcess<2>::local_assemblers_>;

struct LACurrentStates
{
    using type =
        std::vector<RM::StatefulData<2>> RM::LocalAssemblerInterface<2>::*;
    friend type member(LACurrentStates);
};
template struct Access<LACurrentStates,
                       &RM::LocalAssemblerInterface<2>::current_states_>;

struct LAOutputData
{
    using type =
        std::vector<RM::OutputData<2>> RM::LocalAssemblerInterface<2>::*;
    friend type member(LAOutputData);
};
template struct Access<LAOutputData,
                       &RM::LocalAssemblerInterface<2>::output_data_>;

// ── deck handling ───────────────────────────────────────────────────────────
std::string const deck_dir =
    TestInfoLib::TestInfo::data_path +
    "/RichardsMechanics/ANCHORS_MS33_ModelI";
std::string const deck_name = "ms33_modelI_dd1600.prj";

// Pressure Dirichlet ramp of the deck (curve suction_ramp, p_bc_scale = 1):
// -1e8 Pa at t = 0, 0 at t = 1728000 s and after.
double pressureBC(double const t)
{
    double constexpr t_ramp = 1728000.0;
    double constexpr p0 = -1e8;
    return t >= t_ramp ? 0.0 : p0 * (1.0 - t / t_ramp);
}

void replaceOnce(std::string& s, std::string const& what,
                 std::string const& with)
{
    auto const pos = s.find(what);
    ASSERT_NE(pos, std::string::npos) << "deck text not found: " << what;
    ASSERT_EQ(s.find(what, pos + 1), std::string::npos)
        << "deck text not unique: " << what;
    s.replace(pos, what.size(), with);
}

struct SwitchSet
{
    std::string name;
    bool level_form;
    bool z_weight;
};

// Writes the patched deck into dir and returns its path.
std::string writeDeck(std::filesystem::path const& dir, SwitchSet const& sw,
                      bool const strained)
{
    std::ifstream in(deck_dir + "/" + deck_name);
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string s = buffer.str();

    std::string const k_level =
        "<swelling_stress_k_level>true</swelling_stress_k_level>";
    std::string extra;
    if (sw.level_form)
    {
        extra += "\n<swelling_stress_form>level</swelling_stress_form>";
    }
    if (sw.z_weight)
    {
        extra +=
            "\n<micro_potential_partner_z_weight>true</"
            "micro_potential_partner_z_weight>";
    }
    replaceOnce(s, k_level, k_level + extra);
    if (strained)
    {
        replaceOnce(s,
                    "<boundary_condition><mesh>square_1e-2_quad_1e0_top</"
                    "mesh><type>Dirichlet</type><component>1</"
                    "component><parameter>zero</parameter></"
                    "boundary_condition>",
                    "<boundary_condition><mesh>square_1e-2_quad_1e0_top</"
                    "mesh><type>Dirichlet</type><component>1</"
                    "component><parameter>u_top_strained</parameter></"
                    "boundary_condition>");
        replaceOnce(s, "    </parameters>",
                    "<parameter><name>u_top_strained</name><type>Constant</"
                    "type><value>-1e-4</value></parameter>\n    </"
                    "parameters>");
    }
    auto const path = dir / (sw.name + (strained ? "_strained" : "") + ".prj");
    std::ofstream out(path);
    out << s;
    return path.string();
}

// ── one deck, driven by its time loop ───────────────────────────────────────
struct DeckRun
{
    std::unique_ptr<BaseLib::ConfigTree> config;
    std::unique_ptr<ProjectData> project;

    DeckRun(std::string const& prj, std::string const& outdir)
    {
        std::stringstream prj_stream;
        BaseLib::prepareProjectFile(prj_stream, prj, {}, false, outdir);
        config = std::make_unique<BaseLib::ConfigTree>(BaseLib::makeConfigTree(
            prj, true, "OpenGeoSysProject", prj_stream));
        config->ignoreConfigParameter("test_definition");
        config->ignoreConfigParameter("insitu");
        project = std::make_unique<ProjectData>(*config, outdir, deck_dir,
                                                deck_dir);
        for (auto& p : project->getProcesses())
        {
            p->initialize(project->getMedia());
        }
        BaseLib::checkAndInvalidate(*config);
        project->getTimeLoop().initialize();
    }

    ProcessLib::TimeLoop& timeLoop() { return project->getTimeLoop(); }
    double t() { return timeLoop().currentTime()(); }
    double dt() { return (timeLoop().*member(TimeLoopDt{}))(); }

    // Accepted steps until t >= t_target.
    bool runTo(double const t_target)
    {
        while (t() < t_target)
        {
            if (!timeLoop().executeTimeStep())
            {
                // a rejected step is repeated by the time stepper
            }
            if (!timeLoop().calculateNextTimeStep())
            {
                return t() >= t_target;
            }
        }
        return true;
    }

    ProcessLib::Process& process() { return *project->getProcesses()[0]; }

    RM::LocalAssemblerInterface<2>& localAssembler()
    {
        auto& rm = dynamic_cast<RM::RichardsMechanicsProcess<2>&>(process());
        auto& las = rm.*member(RMLocalAssemblers{});
        EXPECT_EQ(las.size(), 1u);
        return *las[0];
    }

    std::vector<double> localSolution()
    {
        auto const& x = *(timeLoop().*member(TimeLoopSolutions{}))[0];
        auto const indices =
            NumLib::getIndices(0, process().getDOFTable(0));
        return x.get(indices);
    }
};

// ── FD comparison ───────────────────────────────────────────────────────────
constexpr int n_p = 4;
constexpr int n_u = 8;
constexpr int n_dof = n_p + n_u;
constexpr double p_ref = 1e8;   // [Pa]  deck initial suction magnitude
constexpr double u_ref = 1e-4;  // [m]   1 % of the element height

struct BlockError
{
    double max_abs = 0;        // max |J_an - J_fd|
    double max_rel = 0;        // max_abs / max(max|J_an|, max|J_fd|)
    double max_entry_rel = 0;  // max_ij |d_ij| / max(|a_ij|,|f_ij|), entries
                               // above 1e-8 of the block maximum
    double max_an = 0;
    double max_fd = 0;
};

BlockError compareBlock(Eigen::MatrixXd const& a, Eigen::MatrixXd const& f)
{
    BlockError e;
    Eigen::MatrixXd const d = a - f;
    e.max_abs = d.cwiseAbs().maxCoeff();
    e.max_an = a.cwiseAbs().maxCoeff();
    e.max_fd = f.cwiseAbs().maxCoeff();
    double const den = std::max(e.max_an, e.max_fd);
    e.max_rel = den > 0 ? e.max_abs / den : 0.0;
    for (int i = 0; i < a.rows(); ++i)
    {
        for (int j = 0; j < a.cols(); ++j)
        {
            double const m = std::max(std::abs(a(i, j)), std::abs(f(i, j)));
            if (m > 1e-8 * den && m > 0)
            {
                e.max_entry_rel =
                    std::max(e.max_entry_rel, std::abs(d(i, j)) / m);
            }
        }
    }
    return e;
}

struct HMode
{
    std::string name;
    // returns h for dof j
    std::function<double(std::vector<double> const&, int)> h;
};

std::vector<HMode> hModes()
{
    auto rel = [](double const r)
    {
        return [r](std::vector<double> const& x, int const j)
        {
            double const ref = j < n_p ? p_ref : u_ref;
            return r * (std::abs(x[j]) > 0 ? std::abs(x[j]) : ref);
        };
    };
    return {
        {"abs", [](std::vector<double> const&, int const j)
         { return j < n_p ? 1.0 : 1e-9; }},
        {"rel1e-2", rel(1e-2)},
        {"rel1e-7", rel(1e-7)},
    };
}

struct CaseResult
{
    std::string switch_set;
    std::string state;
    std::string h_mode;
    double t;
    double dt;
    int kkt_active;
    double n_l_mean;
    double eps_v_mean;
    std::map<std::string, BlockError> blocks;
    double fro_rel;  // ||J_an - J_fd||_F / ||J_fd||_F, rows scaled by
                     // their block-row maximum of |J_fd|
};

struct Evaluator
{
    RM::LocalAssemblerInterface<2>& la;
    double t;   // time of the evaluation (t_n + dt)
    double dt;
    std::vector<double> x;
    std::vector<double> x_prev;
    std::vector<RM::StatefulData<2>> states0;
    std::vector<RM::OutputData<2>> output0;

    Evaluator(RM::LocalAssemblerInterface<2>& la_, double t_, double dt_,
              std::vector<double> x_, std::vector<double> x_prev_)
        : la(la_), t(t_), dt(dt_), x(std::move(x_)), x_prev(std::move(x_prev_))
    {
        states0 = la.*member(LACurrentStates{});
        output0 = la.*member(LAOutputData{});
    }

    void restore()
    {
        la.*member(LACurrentStates{}) = states0;
        la.*member(LAOutputData{}) = output0;
    }

    // residual R = -local_rhs, with the state restored to S0 first and
    // n_relag lag-converging evaluations before the one read.
    Eigen::VectorXd residual(std::vector<double> const& xx, int n_relag,
                             Eigen::MatrixXd* jac = nullptr)
    {
        restore();
        std::vector<double> rhs;
        std::vector<double> J;
        for (int r = 0; r <= n_relag; ++r)
        {
            rhs.clear();
            J.clear();
            la.assembleWithJacobian(t, dt, xx, x_prev, rhs, J);
        }
        if (jac)
        {
            *jac = Eigen::Map<Eigen::Matrix<double, n_dof, n_dof,
                                            Eigen::RowMajor>>(J.data());
        }
        Eigen::VectorXd R = -Eigen::Map<Eigen::VectorXd>(rhs.data(), n_dof);
        restore();
        return R;
    }

    Eigen::MatrixXd analytic()
    {
        Eigen::MatrixXd J;
        residual(x, 0, &J);
        return J;
    }

    Eigen::MatrixXd fd(HMode const& mode, int n_relag)
    {
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n_dof, n_dof);
        for (int j = 0; j < n_dof; ++j)
        {
            double const h = mode.h(x, j);
            auto xp = x;
            auto xm = x;
            xp[j] += h;
            xm[j] -= h;
            double const h_total = (xp[j] - x[j]) + (x[j] - xm[j]);
            J.col(j) = (residual(xp, n_relag) - residual(xm, n_relag)) / h_total;
        }
        return J;
    }
};

std::map<std::string, BlockError> compareAll(Eigen::MatrixXd const& a,
                                             Eigen::MatrixXd const& f,
                                             double& fro_rel)
{
    std::map<std::string, BlockError> b;
    b["pp"] = compareBlock(a.block(0, 0, n_p, n_p), f.block(0, 0, n_p, n_p));
    b["pu"] = compareBlock(a.block(0, n_p, n_p, n_u), f.block(0, n_p, n_p, n_u));
    b["up"] = compareBlock(a.block(n_p, 0, n_u, n_p), f.block(n_p, 0, n_u, n_p));
    b["uu"] = compareBlock(a.block(n_p, n_p, n_u, n_u),
                           f.block(n_p, n_p, n_u, n_u));
    // row scaling by the equation's (p or u rows) maximum |J_fd| entry, so
    // that both equations count (their units differ).
    double const sp = std::max(f.topRows(n_p).cwiseAbs().maxCoeff(), 1e-300);
    double const su =
        std::max(f.bottomRows(n_u).cwiseAbs().maxCoeff(), 1e-300);
    Eigen::MatrixXd d = a - f;
    Eigen::MatrixXd fs = f;
    d.topRows(n_p) /= sp;
    fs.topRows(n_p) /= sp;
    d.bottomRows(n_u) /= su;
    fs.bottomRows(n_u) /= su;
    fro_rel = d.norm() / fs.norm();
    return b;
}

void writeMatrix(std::ostream& os, Eigen::MatrixXd const& m)
{
    os << std::setprecision(17);
    for (int i = 0; i < m.rows(); ++i)
    {
        for (int j = 0; j < m.cols(); ++j)
        {
            os << m(i, j) << (j + 1 < m.cols() ? "," : "\n");
        }
    }
}

std::string fmt3(double v)
{
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << v;
    return o.str();
}

struct StateSummary
{
    int kkt_active = 0;
    double n_l_mean = 0;
    double eps_v_mean = 0;
};

StateSummary summarise(RM::LocalAssemblerInterface<2>& la)
{
    StateSummary s;
    auto const& states = la.*member(LACurrentStates{});
    for (auto const& st : states)
    {
        // MicroCeilingKktStatus::Active == 1 (RichardsMechanicsFEM-impl.h)
        s.kkt_active += *std::get<RM::MicroCeilingStatus>(st) == 1.0 ? 1 : 0;
        s.n_l_mean += *std::get<RM::MicroWaterContent>(st) / states.size();
        auto const& eps = std::get<RM::StrainData<2>>(st).eps;
        s.eps_v_mean += (eps[0] + eps[1] + eps[2]) / states.size();
    }
    return s;
}

}  // namespace

TEST(RichardsMechanics, DSMElementJacobianFD)
{
    auto const tmp = std::filesystem::temp_directory_path() /
                     "ogs_dsm_element_jacobian_fd";
    std::filesystem::create_directories(tmp);

    std::vector<SwitchSet> const switch_sets = {
        {"deck_step_klevel", false, false},
        {"level", true, false},
        {"level_z", true, true},
    };
    // states: name, deck strained?, time to run to, add strain in x?
    struct State
    {
        std::string name;
        bool strained_deck;
        double t_target;
        bool strain_step;
    };
    std::vector<State> const states = {
        {"dry", false, 0.0, false},
        {"mid", false, 864000.0, false},
        {"mid_eps_step", false, 864000.0, true},
        {"mid_eps_hist", true, 864000.0, false},
        {"wet_ceiling", false, 8640000.0, false},
    };

    char const* const csv_env = std::getenv("OGS_DSM_JACFD_CSV");
    char const* const mat_env = std::getenv("OGS_DSM_JACFD_MATRICES");
    std::ofstream csv;
    if (csv_env)
    {
        csv.open(csv_env);
        csv << "switch_set,state,h_mode,relag,t,dt,kkt_active_ip,n_l_mean,"
               "eps_v_mean,block,max_abs_err,max_rel_err,max_entry_rel_err,"
               "max_abs_analytic,max_abs_fd,fro_rel\n";
    }

    std::vector<CaseResult> results;

    for (auto const& sw : switch_sets)
    {
        for (bool const strained : {false, true})
        {
            auto const prj = writeDeck(tmp, sw, strained);
            auto const outdir =
                tmp / ("out_" + sw.name + (strained ? "_strained" : ""));
            std::filesystem::create_directories(outdir);
            DeckRun run(prj, outdir.string());

            for (auto const& st : states)
            {
                if (st.strained_deck != strained)
                {
                    continue;
                }
                ASSERT_TRUE(run.runTo(st.t_target))
                    << sw.name << " " << st.name;
                auto& la = run.localAssembler();
                double const t_n = run.t();
                double const dt = run.dt();
                auto const x_prev = run.localSolution();
                auto x = x_prev;
                for (int k = 0; k < n_p; ++k)
                {
                    x[k] = pressureBC(t_n + dt);  // every node is Dirichlet
                }
                if (st.strain_step)
                {
                    // eps_zz = -0.01: u_z = -0.01 (z - z_min); u_z dofs are
                    // local dofs n_p + 4 .. n_p + 7 (component-wise order).
                    auto const& e = *run.process().getMesh().getElement(0);
                    double z_min = 1e300;
                    for (unsigned k = 0; k < 4; ++k)
                    {
                        z_min = std::min(z_min, (*e.getNode(k))[1]);
                    }
                    for (unsigned k = 0; k < 4; ++k)
                    {
                        x[n_p + 4 + k] += -0.01 * ((*e.getNode(k))[1] - z_min);
                    }
                }
                auto const summary = summarise(la);

                Evaluator ev(la, t_n + dt, dt, x, x_prev);
                Eigen::MatrixXd const J_an = ev.analytic();
                auto const after = summarise(la);  // restored = S0
                (void)after;

                if (mat_env)
                {
                    std::ofstream m(std::string(mat_env) + "/" + sw.name +
                                    "_" + st.name + "_analytic.csv");
                    writeMatrix(m, J_an);
                }

                for (auto const& mode : hModes())
                {
                    for (int const relag : {0, 5})
                    {
                        if (relag > 0 && mode.name != "abs")
                        {
                            continue;
                        }
                        Eigen::MatrixXd const J_fd = ev.fd(mode, relag);
                        CaseResult r;
                        r.switch_set = sw.name;
                        r.state = st.name;
                        r.h_mode = mode.name + (relag ? "+relag5" : "");
                        r.t = t_n + dt;
                        r.dt = dt;
                        r.kkt_active = summary.kkt_active;
                        r.n_l_mean = summary.n_l_mean;
                        r.eps_v_mean = summary.eps_v_mean;
                        r.blocks = compareAll(J_an, J_fd, r.fro_rel);
                        results.push_back(r);
                        if (mat_env)
                        {
                            std::ofstream m(std::string(mat_env) + "/" +
                                            sw.name + "_" + st.name + "_fd_" +
                                            r.h_mode + ".csv");
                            writeMatrix(m, J_fd);
                        }
                        if (csv)
                        {
                            csv << std::setprecision(10);
                            for (auto const& [name, b] : r.blocks)
                            {
                                csv << r.switch_set << "," << r.state << ","
                                    << mode.name << "," << relag << "," << r.t
                                    << "," << r.dt << "," << r.kkt_active
                                    << "," << r.n_l_mean << "," << r.eps_v_mean
                                    << "," << name << "," << b.max_abs << ","
                                    << b.max_rel << "," << b.max_entry_rel
                                    << "," << b.max_an << "," << b.max_fd << ","
                                    << r.fro_rel << "\n";
                            }
                        }
                    }
                }
            }
        }
    }

    // stdout table: block max relative error (block-normalised)
    std::cout << "\nDSM element Jacobian vs central FD, max_rel = max|J_an - "
                 "J_fd| / max(max|J_an|, max|J_fd|) per block\n";
    std::cout << std::left << std::setw(18) << "switch_set" << std::setw(14)
              << "state" << std::setw(14) << "h" << std::setw(5) << "kkt"
              << std::setw(11) << "pp" << std::setw(11) << "pu"
              << std::setw(11) << "up" << std::setw(11) << "uu"
              << std::setw(11) << "fro_rel" << "\n";
    for (auto const& r : results)
    {
        std::cout << std::left << std::setw(18) << r.switch_set
                  << std::setw(14) << r.state << std::setw(14) << r.h_mode
                  << std::setw(5) << r.kkt_active << std::setw(11)
                  << fmt3(r.blocks.at("pp").max_rel) << std::setw(11)
                  << fmt3(r.blocks.at("pu").max_rel) << std::setw(11)
                  << fmt3(r.blocks.at("up").max_rel) << std::setw(11)
                  << fmt3(r.blocks.at("uu").max_rel) << std::setw(11)
                  << fmt3(r.fro_rel) << "\n";
    }

    // Regression guard: the blocks measured exact to round-off in the report
    // (abs h, fixed lag) must stay so.
    for (auto const& r : results)
    {
        if (r.h_mode != "abs")
        {
            continue;
        }
        EXPECT_LT(r.blocks.at("uu").max_rel, 1e-6)
            << r.switch_set << " " << r.state;
    }
}
