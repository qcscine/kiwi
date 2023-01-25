project(KiwiApp)

set(headers
        CalculationManager.h
        ReferenceCalculationTask.h
        ParticleDensityTask.h
        NaturalOrbitalsTask.h
        Ao2MoTask.h
        TaskFactory.h
        Task.h
        Keywords.h
        )
set(sources
        main.cpp
        ReferenceCalculationTask.cpp
        ParticleDensityTask.cpp
        NaturalOrbitalsTask.cpp
        Ao2MoTask.cpp
        CalculationManager.cpp
        )

# Executable
find_package(Boost COMPONENTS program_options system REQUIRED)

# Set the RPATH to be used when installing. $ORIGIN refers to the installation directory of the binary.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

add_executable(KiwiApp ${headers} ${sources})

if (NOT APPLE)
    set_target_properties(KiwiApp PROPERTIES
            OUTPUT_NAME kiwi
            )
endif()


target_compile_definitions(KiwiApp
        PUBLIC
        BOOST_ALL_DYN_LINK
        )

target_compile_options(KiwiApp
        PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
        )


target_link_libraries(KiwiApp
        PUBLIC
        KiwiUtils
        HartreeFock
        Scine::Core
        Eigen3::Eigen
        Scine::UtilsOS
        Boost::program_options
        yaml-cpp
        )

# Add link library specific for MinGW
if(WIN32)
    target_link_libraries(KiwiApp PRIVATE $<$<BOOL:MINGW>:ws2_32>)
endif()

set_target_properties(KiwiApp PROPERTIES OUTPUT_NAME "kiwiBinary")

install(TARGETS KiwiApp RUNTIME DESTINATION bin)

