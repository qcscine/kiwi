find_package(Eigen3 REQUIRED)

# Dependencies
include(ImportCore)
import_core()
include(ImportUtilsOS)
import_utils_os()
include(ImportIntegralEvaluator)
import_integral_evaluator()


# Module libraries
#add_subdirectory(HartreeFock)
#add_subdirectory(Module)

# Only look for OpenMP if parallelization is sought
if(SCINE_PARALLELIZE)
    find_package(OpenMP REQUIRED)
endif()


#
# KiwiUtils Target
#
include(${CMAKE_CURRENT_SOURCE_DIR}/Kiwi/KiwiUtils/Files.cmake)

add_library(KiwiUtils STATIC ${KIWIUTILS_HEADERS} ${KIWIUTILS_SOURCES})
if(NOT BUILD_SHARED_LIBS)
    set_target_properties(KiwiUtils PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

set_target_properties(KiwiUtils PROPERTIES
        OUTPUT_NAME kiwiutils
        WINDOWS_EXPORT_ALL_SYMBOLS ON
        )

target_link_libraries(KiwiUtils
        PUBLIC
        Scine::UtilsOS
        Scine::LibintIntegrals
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

target_include_directories(KiwiUtils PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/Scine/>
        )


# -- Define Install
# Headers
install(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/KiwiUtils/
        DESTINATION include/Scine/KiwiUtils/
        FILES_MATCHING PATTERN "*.h"
)

# Targets
install(
        TARGETS KiwiUtils
        EXPORT kiwiUtilsTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

# Component-level cmake package files
scine_install_component_cmake_files(
        COMPONENT ScineKiwiUtils
        EXPORT_NAME kiwiUtilsTargets
)


#
# HartreeFock Target
#
include(${CMAKE_CURRENT_SOURCE_DIR}/Kiwi/HartreeFock/Files.cmake)

add_library(HartreeFock STATIC ${HARTREEFOCK_HEADERS} ${HARTREEFOCK_SOURCES})
if(NOT BUILD_SHARED_LIBS)
    set_target_properties(HartreeFock PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

set_target_properties(HartreeFock PROPERTIES
        OUTPUT_NAME hartreefock
        WINDOWS_EXPORT_ALL_SYMBOLS ON
        )

target_link_libraries(HartreeFock
        PUBLIC
        Scine::UtilsOS
        KiwiUtils
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

target_include_directories(HartreeFock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/Scine/>
        )

# -- Define Install
# Headers
install(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/HartreeFock/
        DESTINATION include/Scine/HartreeFock/
        FILES_MATCHING PATTERN "*.h"
)

# Targets
install(
        TARGETS HartreeFock KiwiUtils
        EXPORT hartreeFockTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)

# Component-level cmake package files
scine_install_component_cmake_files(
        COMPONENT ScineHartreeFock
        EXPORT_NAME hartreeFockTargets
)

#
# KiwiOpt Target
#
include(${CMAKE_CURRENT_SOURCE_DIR}/Kiwi/KiwiOpt/Files.cmake)

# add_library(KiwiOpt STATIC ${KIWIOPT_HEADERS} ${KIWIOPT_SOURCES})
# if(NOT BUILD_SHARED_LIBS)
#     set_target_properties(KiwiOpt PROPERTIES POSITION_INDEPENDENT_CODE ON)
# endif()
#
# set_target_properties(KiwiOpt PROPERTIES
#         OUTPUT_NAME kiwiopt
#         WINDOWS_EXPORT_ALL_SYMBOLS ON
#         )
#
# target_link_libraries(KiwiOpt
#         PUBLIC
#         Scine::OptOS
#         Scine::LibintIntegrals
#         $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
#
# target_include_directories(KiwiOpt PUBLIC
#         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#         $<INSTALL_INTERFACE:include/Scine/>
#         )


# -- Define Install
# Headers
install(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/KiwiOpt/
        DESTINATION include/Scine/KiwiOpt/
        FILES_MATCHING PATTERN "*.h"
)

# # Targets
# install(
#         TARGETS KiwiOpt
#         EXPORT kiwiOptTargets
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib
# )

# # Component-level cmake package files
# scine_install_component_cmake_files(
#         COMPONENT ScineKiwi
#         EXPORT_NAME kiwiOptTargets
# )




# App
add_subdirectory(App)

if (SCINE_BUILD_TESTS)
    add_subdirectory(Tests)
endif()


if(WIN32)
    set_target_properties(KiwiUtils PROPERTIES PREFIX "" IMPORT_PREFIX "")
    if(MINGW)
        target_link_libraries(KiwiUtils PRIVATE ws2_32)
    endif()
endif()
