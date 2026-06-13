cmake_minimum_required(VERSION 3.21)
project(OpenFDTD_X
    VERSION 1.0.0
    DESCRIPTION "Multi-Domain FDTD GUI (EM / Optical / Acoustic / Underwater)"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# Qt6 Widgets is the only hard requirement (classic C++ Qt, NOT QML/Quick).
find_package(Qt6 REQUIRED COMPONENTS Widgets)

option(USE_HDF5   "Enable HDF5 time-series export"            OFF)
option(USE_LIBIGL "Use libigl for STL import + voxelization"  OFF)
option(BUILD_TESTS "Build the .ofd round-trip self test"      ON)

if(USE_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS C)
endif()

if(USE_LIBIGL)
    include(FetchContent)
    set(LIBIGL_GLFW OFF CACHE BOOL "")
    set(LIBIGL_IMGUI OFF CACHE BOOL "")
    set(LIBIGL_OPENGL OFF CACHE BOOL "")
    set(LIBIGL_PNG OFF CACHE BOOL "")
    FetchContent_Declare(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.5.0
    )
    FetchContent_MakeAvailable(libigl)
endif()

set(CORE_SOURCES
    src/core/Project.cpp
    src/core/Project.h
    src/core/Domain.h
    src/core/Material.h
    src/core/Geometry.h
    src/core/Source.h
    src/core/MeshAxis.h
    src/core/PostOpts.h

    src/io/OfdIO.cpp
    src/io/OfdIO.h
    src/io/Touchstone.cpp
    src/io/Touchstone.h
    src/io/Tidy3dExporter.cpp
    src/io/Tidy3dExporter.h
    src/io/H5Writer.cpp
    src/io/H5Writer.h
    src/io/EvReader.h
    src/io/StlImporter.cpp
    src/io/StlImporter.h

    src/kernel/Runner.cpp
    src/kernel/Runner.h
)

set(GUI_SOURCES
    src/main.cpp
    src/MainWindow.cpp
    src/MainWindow.h
    src/DomainBar.cpp
    src/DomainBar.h
    src/RightDock.cpp
    src/RightDock.h
    src/I18n.cpp
    src/I18n.h

    src/widgets/SectionBox.cpp
    src/widgets/SectionBox.h
    src/widgets/UnitNav.cpp
    src/widgets/UnitNav.h
    src/widgets/Viewport3D.cpp
    src/widgets/Viewport3D.h
    src/widgets/PlotPanel.cpp
    src/widgets/PlotPanel.h
    src/widgets/LogConsole.cpp
    src/widgets/LogConsole.h
    src/widgets/EvViewer.cpp
    src/widgets/EvViewer.h

    src/tabs/GeneralTab.cpp
    src/tabs/GeneralTab.h
    src/tabs/MeshTab.cpp
    src/tabs/MeshTab.h
    src/tabs/MaterialTab.cpp
    src/tabs/MaterialTab.h
    src/tabs/GeometryTab.cpp
    src/tabs/GeometryTab.h
    src/tabs/SourceTab.cpp
    src/tabs/SourceTab.h
    src/tabs/Post1Tab.cpp
    src/tabs/Post1Tab.h
    src/tabs/Post2Tab.cpp
    src/tabs/Post2Tab.h
    src/tabs/OpticalTab.cpp
    src/tabs/OpticalTab.h
    src/tabs/AcousticTab.cpp
    src/tabs/AcousticTab.h
    src/tabs/UnderwaterTab.cpp
    src/tabs/UnderwaterTab.h
    src/tabs/Tidy3dTab.cpp
    src/tabs/Tidy3dTab.h

    resources/openfdtd.qrc
)

qt_add_executable(openfdtd_x ${CORE_SOURCES} ${GUI_SOURCES})

target_include_directories(openfdtd_x PRIVATE src)
target_link_libraries(openfdtd_x PRIVATE Qt6::Widgets)

if(USE_HDF5)
    target_link_libraries(openfdtd_x PRIVATE ${HDF5_C_LIBRARIES})
    target_include_directories(openfdtd_x PRIVATE ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(openfdtd_x PRIVATE OFD_USE_HDF5)
endif()

if(USE_LIBIGL)
    target_link_libraries(openfdtd_x PRIVATE igl::core)
    target_compile_definitions(openfdtd_x PRIVATE OFD_USE_LIBIGL)
endif()

if(MSVC)
    target_compile_options(openfdtd_x PRIVATE /utf-8 /W3)
else()
    target_compile_options(openfdtd_x PRIVATE -Wall -Wextra)
endif()

set_target_properties(openfdtd_x PROPERTIES
    WIN32_EXECUTABLE ON
    MACOSX_BUNDLE    ON
)

# ── Self test: .ofd round-trip against ../data/sample ───────────────────────
if(BUILD_TESTS)
    qt_add_executable(ofdx_selftest ${CORE_SOURCES} tests/selftest.cpp)
    target_include_directories(ofdx_selftest PRIVATE src)
    target_link_libraries(ofdx_selftest PRIVATE Qt6::Widgets)
    set_target_properties(ofdx_selftest PROPERTIES
        WIN32_EXECUTABLE OFF
        MACOSX_BUNDLE    OFF
    )
endif()
