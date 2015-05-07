# Source files
set(SOURCES
	mainwindow.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../Utils/OGSFileConverter/OGSFileConverter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../Utils/OGSFileConverter/FileListDialog.cpp
)

# Moc Header files
set(MOC_HEADERS
	mainwindow.h
	${CMAKE_CURRENT_SOURCE_DIR}/../Utils/OGSFileConverter/OGSFileConverter.h
	${CMAKE_CURRENT_SOURCE_DIR}/../Utils/OGSFileConverter/FileListDialog.h
)

# UI files
set(UIS
	mainwindow.ui
	../Utils/OGSFileConverter/OGSFileConverter.ui
	../Utils/OGSFileConverter/FileList.ui
)


# Run Qts user interface compiler uic on .ui files
qt4_wrap_ui(UI_HEADERS ${UIS} )

qt4_add_resources(QTRESOURCES ./Img/icons.qrc )

# Run Qts meta object compiler moc on header files
qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS} )

# Include the headers which are generated by uic and moc
# and include additional header
set(SOURCE_DIR_REL ${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(
	${SOURCE_DIR_REL}/BaseLib
	${SOURCE_DIR_REL}/MathLib
	${SOURCE_DIR_REL}/GeoLib
	${SOURCE_DIR_REL}/FileIO
	${SOURCE_DIR_REL}/MeshLib
	${SOURCE_DIR_REL}/MeshLibGEOTOOLS
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_BINARY_DIR}/Base
	${CMAKE_CURRENT_BINARY_DIR}/DataView
	${CMAKE_CURRENT_BINARY_DIR}/DataView/StratView
	${CMAKE_CURRENT_BINARY_DIR}/DataView/DiagramView
	${CMAKE_CURRENT_BINARY_DIR}/VtkVis
	${CMAKE_CURRENT_BINARY_DIR}/VtkAct
	${CMAKE_CURRENT_BINARY_DIR}/Applications/Utils/OGSFileConverter
	${CMAKE_CURRENT_SOURCE_DIR}/Base
	${CMAKE_CURRENT_SOURCE_DIR}/DataView
	${CMAKE_CURRENT_SOURCE_DIR}/DataView/StratView
	${CMAKE_CURRENT_SOURCE_DIR}/DataView/DiagramView
	${CMAKE_CURRENT_SOURCE_DIR}/VtkVis
	${CMAKE_CURRENT_SOURCE_DIR}/VtkAct
)

# Put moc files in a project folder
source_group("UI Files" REGULAR_EXPRESSION "\\w*\\.ui")
source_group("Moc Files" REGULAR_EXPRESSION "moc_.*")

# Application icon
set(APP_ICON ${SOURCE_DIR_REL}/scripts/packaging/ogs-de-icon.icns)

# Create the executable
add_executable(DataExplorer MACOSX_BUNDLE
	main.cpp
	${SOURCES}
	${MOC_HEADERS}
	${MOC_SOURCES}
	${UIS}
	${QTRESOURCES}
	${APP_ICON}
	exe-icon.rc
)

target_link_libraries(DataExplorer
	${QT_LIBRARIES}
	ApplicationsLib
	BaseLib
	GeoLib
	FileIO
	InSituLib
	MeshLib
	#MSHGEOTOOLS
	QtBase
	QtDataView
	QtStratView
	VtkVis
	VtkAct
	${Boost_LIBRARIES}
	${CATALYST_LIBRARIES}
	zlib
	shp
)

if(VTK_NETCDF_FOUND)
	target_link_libraries(DataExplorer vtkNetCDF vtkNetCDF_cxx )
else()
	target_link_libraries(DataExplorer ${Shapelib_LIBRARIES} )
endif () # Shapelib_FOUND

if (GEOTIFF_FOUND)
	target_link_libraries(DataExplorer ${GEOTIFF_LIBRARIES} )
endif () # GEOTIFF_FOUND

add_dependencies (DataExplorer VtkVis)

if(MSVC)
	# Set linker flags
	set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:MSVCRT /IGNORE:4099")
	target_link_libraries(DataExplorer winmm)
endif()

### OpenSG support ###
if (VTKOSGCONVERTER_FOUND)
	USE_OPENSG(DataExplorer)
	include_directories(${VTKOSGCONVERTER_INCLUDE_DIRS})
	target_link_libraries(DataExplorer ${VTKOSGCONVERTER_LIBRARIES})
endif ()

if(VTKFBXCONVERTER_FOUND)
	target_link_libraries(DataExplorer ${VTKFBXCONVERTER_LIBRARIES})
endif()

include(AddCatalystDependency)
ADD_CATALYST_DEPENDENCY(DataExplorer)

set_property(TARGET DataExplorer PROPERTY FOLDER "DataExplorer")


####################
### Installation ###
####################
if(APPLE)
	include(packaging/PackagingMacros)
	ConfigureMacOSXBundle(DataExplorer ${APP_ICON})

	install(TARGETS DataExplorer DESTINATION .)
	set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION .)
	include(InstallRequiredSystemLibraries)
	include(DeployQt4)
	INSTALL_QT4_EXECUTABLE(DataExplorer.app "" "" "" "" "" ogs_gui)
else()
	install(TARGETS DataExplorer RUNTIME DESTINATION bin COMPONENT ogs_gui)
endif()

cpack_add_component(ogs_gui
	DISPLAY_NAME "OGS Data Explorer"
	DESCRIPTION "The graphical user interface for OpenGeoSys."
	GROUP Applications
)
set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "DataExplorer" "OGS Data Explorer" PARENT_SCOPE)
set(CPACK_NSIS_MENU_LINKS ${CPACK_NSIS_MENU_LINKS} "bin/DataExplorer.exe" "Data Explorer" PARENT_SCOPE)
if(APPLE)
	return()
endif()

include(packaging/InstallDependencies)
InstallDependencies(DataExplorer ogs_gui)
