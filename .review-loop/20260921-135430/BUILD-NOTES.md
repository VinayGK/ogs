# Compilation notes — `dsm_upstream_mr_2026-09-20` on Apple silicon / macOS 27

Companion to `HANDOFF.md`. Written 2026-09-21 on the laptop; everything below was
run and observed on that host, not inferred.

---

## 1. The short version

The branch **does not compile unaided** on this host. One build-dir-only flag fixes it:

```
-DCMAKE_CXX_FLAGS="-D_LIBCPP_TEMPLATE_VIS="
```

No source file is modified, so **the fix does not travel with the branch** — it lives
in the build directory's `CMakeCache.txt` and must be reapplied wherever you configure
next. The cause is a pre-existing dependency/toolchain incompatibility, not anything
the branch introduced.

## 2. Host this was established on

| | |
|---|---|
| OS | macOS **27.0** (build 26A428), arm64 |
| Compiler | Apple clang **21.0.0** (clang-2100.3.34.2), target `arm64-apple-darwin27.0.0` |
| Toolchain | `/usr/bin/c++`, Xcode at `/Applications/Xcode.app/Contents/Developer` |
| SDK | **27.0** (`MacOSX.sdk`) |
| CMake / Ninja | 4.4.3 / 1.13.2 |
| Cores | 18 physical, 18 logical |

Note for anyone reading older notes of mine: I first wrote "macOS 26" — that was wrong.
It came from misreading `CMAKE_OSX_DEPLOYMENT_TARGET:STRING=26.0` in an older build
directory as an OS version. `sw_vers` says 27.0.

## 3. The failure

Configure succeeds. The build then fails across many targets with, e.g.:

```
/Users/vinaykumar/git/build/shared_cpm_cache/range-v3/c704/include/meta/meta.hpp:3793:7:
  note: forward declaration of 'std::_LIBCPP_TEMPLATE_VIS'
 3793 | class META_TEMPLATE_VIS allocator;
      |       ^
.../meta/meta_fwd.hpp:221:27: note: expanded from macro 'META_TEMPLATE_VIS'
  221 | #define META_TEMPLATE_VIS _LIBCPP_TEMPLATE_VIS
fatal error: too many errors emitted, stopping now [-ferror-limit=]
20 errors generated.
```

Failing targets (from `ninja -k 0`):

```
BaseLib/CMakeFiles/BaseLib.dir/cmake_pch.hxx.pch
GeoLib/CMakeFiles/GeoLib.dir/cmake_pch.hxx.pch
MathLib/CMakeFiles/MathLib.dir/cmake_pch.hxx.pch
NumLib/CMakeFiles/NumLib.dir/cmake_pch.hxx.pch
MeshLib, MeshGeoToolsLib, MeshToolsLib, ParameterLib  (same, PCH)
MaterialLib/Utils/CMakeFiles/MaterialLib_Utils.dir/Unity/unity_0_cxx.cxx.o
MaterialLib/FractureModels/.../Unity/unity_0_cxx.cxx.o
```

## 4. Why it is not the branch's fault

`range-v3` is pinned at **c704** (`ogs.minimum_version.range-v3`, consumed in
`scripts/cmake/Dependencies.cmake:405-411`). Its `meta_fwd.hpp` defines
`META_TEMPLATE_VIS` as libc++'s `_LIBCPP_TEMPLATE_VIS`, a visibility macro the
libc++ shipped with SDK 27 **no longer defines**. Every TU that pulls range-v3 in
therefore parses `_LIBCPP_TEMPLATE_VIS` as a type name and fails.

Two independent checks that this is pre-existing, not branch-introduced:

1. **The failing targets are not in the branch diff.** The diff touches only
   `ProcessLib/RichardsMechanics`, `Tests`, `Tests/Data` and `web`. `BaseLib`,
   `GeoLib`, `MathLib`, `NumLib` and `MaterialLib` are untouched.
2. **The pin is identical on the pre-branch tree.** `dsm_native_maxwell_conjugate`
   resolves the same `range-v3` version, so upstream master has the same problem on
   this toolchain.

The older build directories (`ogs-dsm-active`, `maxwell_rebased_2026-08-26`, …) still
appear to work only because their objects predate the SDK bump. Reconfiguring any of
them from scratch on this host would hit the same wall.

## 5. The fix, and how it was narrowed

**Use exactly this:**

```bash
cmake -S <worktree> -B <builddir> -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DOGS_USE_MFRONT=ON \
  -DOGS_BUILD_TESTING=ON \
  -DOGS_BUILD_UTILS=ON \
  -DCPM_SOURCE_CACHE=<cpm cache dir> \
  -DCMAKE_CXX_FLAGS="-D_LIBCPP_TEMPLATE_VIS="
```

On the laptop `<builddir>` was `/Users/vinaykumar/git/build/dsm_upstream_mr_2026-09-20`
and `<cpm cache dir>` was `/Users/vinaykumar/git/build/shared_cpm_cache`.
`OGS_BUILD_TESTING=ON` is **required** — the branch adds unit tests, and the older DSM
build dirs have it OFF, so they cannot verify this branch.

**How the flag set was narrowed.** My first attempt passed six defines
(`_LIBCPP_TEMPLATE_VIS`, `_LIBCPP_ENUM_VIS`, `_LIBCPP_EXCEPTION_ABI`, `_LIBCPP_HIDDEN`,
`_LIBCPP_FUNC_VIS`, `_LIBCPP_TYPE_VIS`) on the assumption that the whole visibility
macro family had gone. That built, but it was over-broad and one entry was wrong:

- Grepping the SDK's `__config` shows **`_LIBCPP_HIDDEN` is still defined**
  (`__config:219` and `:231`). Passing `-D_LIBCPP_HIDDEN=` overrode a live macro,
  emptying a visibility attribute libc++ still applies to its own internals — a
  deviation with no justification behind it.
- The other five are genuinely absent from this libc++.
- Compiling `#include <range/v3/all.hpp>` directly: no flags → 20 errors;
  `-D_LIBCPP_TEMPLATE_VIS=` alone → 0 errors.
- The **full** `ogs` + `testrunner` build with `-D_LIBCPP_TEMPLATE_VIS=` alone →
  330/330 targets, exit 0. So no other dependency needed the rest.

The build dir was left reconfigured on the single flag, and the unit suite re-run on
that binary gives the same 44 pass / 2 skip. **Use the one flag, not the six.**

## 6. Verification command and expected result

```bash
cd <builddir> && ninja -j <cores> ogs testrunner \
  && ./bin/testrunner --gtest_filter='RichardsMechanics*'
```

Expected, and what the loop saw in every round:

```
46 tests from 4 test suites ran.
[  PASSED  ] 44 tests.
[  SKIPPED ] 2 tests:
             RichardsMechanicsExactFilmPair.LiquidCarrierEnergyPressureConsistency
             RichardsMechanicsExactFilmPair.ExpulsionProbeDrainedRamp
```

Suites present: `RichardsMechanics`, `RichardsMechanicsExactFilmPair`,
`RichardsMechanicsStrainedFilm`, `RichardsMechanicsLiveKOfRhoD`.

Both skips are pre-existing and carry stated reasons tied to named open modelling
decisions — they are not a symptom of the build workaround.

**Integration ctests were never run.** The loop's verification is unit tests only, so
nothing in this run speaks to benchmark or reference-VTU status. `ctest -LE large
--output-on-failure` from a release build is the documented command
(`AGENTS.md:37-38`) if you want that coverage on the mini.

## 7. What to do on the mac mini

1. Check the toolchain first — the whole problem is SDK-version-dependent:
   ```bash
   sw_vers; /usr/bin/c++ --version; xcrun --show-sdk-version
   ```
2. **Try configuring without `CMAKE_CXX_FLAGS` first.** If the mini is on an older SDK
   whose libc++ still defines `_LIBCPP_TEMPLATE_VIS`, it will build unaided and you
   should not carry the flag over — it would be an unnecessary deviation there.
   Confirm either way with:
   ```bash
   grep -n "define _LIBCPP_TEMPLATE_VIS" "$(xcrun --show-sdk-path)/usr/include/c++/v1/__config"
   ```
   Output → macro present, no flag needed. No output → apply the flag.
3. If it does need the flag, use `-D_LIBCPP_TEMPLATE_VIS=` only.
4. Build dirs on this machine map to sources as follows, if you need a reference
   configuration:

   | build dir | source worktree |
   |---|---|
   | `dsm_upstream_mr_2026-09-20` | `dsm_upstream_mr_2026-09-20_wt` (this branch) |
   | `ogs-dsm-active`, `maxwell_*` | `dsm_native_maxwell_conjugate_wt` |
   | `dsm_exact_weight_2026-09-02` | `dsm_exact_weight_2026-09-02_wt` |
   | `ogs_6.5.8_local` | `tag_6.5.8_local` |

## 8. Standing caveats

- The workaround changes a visibility macro for every TU in the build. It is the
  minimum needed to compile, but the binary under test is not bit-identical to what an
  unaffected toolchain would produce. Worth stating in any report that quotes results
  built this way.
- **The real fix is upstream**, not here: bumping the `range-v3` pin to a version that
  does not reference the removed libc++ macros. That is a dependency-version decision
  on OGS master and is outside this MR's scope — deliberately not done.
- If someone deletes or reconfigures the build dir from scratch without the flag, the
  branch will appear to be broken. It is not; check this file first.
