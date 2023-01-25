project(KiwiTests)

find_package(Boost REQUIRED)

set(TEST_LIBRARIES_HF
        PRIVATE
        gtest_main
        gmock
        KiwiUtils
        HartreeFock
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
        Libint2::cxx
        Boost::boost)

add_executable(kiwiOptTests KiwiOptTests.cpp)
target_link_libraries(kiwiOptTests ${TEST_LIBRARIES_HF})

add_test(NAME KiwiOptTests COMMAND kiwiOptTests)
target_compile_options(kiwiOptTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(kiwiUtilsTests KiwiUtilsTests.cpp)
target_link_libraries(kiwiUtilsTests ${TEST_LIBRARIES_HF})

add_test(NAME KiwiUtilsTests COMMAND kiwiUtilsTests)
target_compile_options(kiwiUtilsTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(bfgsHartreeFockTests BFGSHartreeFockTests.cpp)
target_link_libraries(bfgsHartreeFockTests ${TEST_LIBRARIES_HF})

add_test(NAME BFGSHartreeFockTests COMMAND bfgsHartreeFockTests)
target_compile_options(bfgsHartreeFockTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(sadTests SadTests.cpp)
target_link_libraries(sadTests ${TEST_LIBRARIES_HF})

add_test(NAME SadTests COMMAND sadTests)
target_compile_options(sadTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(preBoAugmentedRoothaanHallTests PreBoAugmentedRoothaanHallTests.cpp)
target_link_libraries(preBoAugmentedRoothaanHallTests ${TEST_LIBRARIES_HF})

add_test(NAME PreBoAugmentedRoothaanHallTests COMMAND preBoAugmentedRoothaanHallTests)
target_compile_options(preBoAugmentedRoothaanHallTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(boAugmentedRoothaanHallTests BoAugmentedRoothaanHallTests.cpp)
target_link_libraries(boAugmentedRoothaanHallTests ${TEST_LIBRARIES_HF})

add_test(NAME BoAugmentedRoothaanHallTests COMMAND boAugmentedRoothaanHallTests)
target_compile_options(boAugmentedRoothaanHallTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(preBoDirectHartreeFockTests PreBoDirectHartreeFockTests.cpp)
target_link_libraries(preBoDirectHartreeFockTests ${TEST_LIBRARIES_HF})

add_test(NAME PreBoDirectHartreeFockTests COMMAND preBoDirectHartreeFockTests)
target_compile_options(preBoDirectHartreeFockTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(boDirectHartreeFockTests BoDirectHartreeFockTests.cpp)
target_link_libraries(boDirectHartreeFockTests ${TEST_LIBRARIES_HF})

add_test(NAME BoDirectHartreeFockTests COMMAND boDirectHartreeFockTests)
target_compile_options(boDirectHartreeFockTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(boRoothaanHallTests BoRoothaanHallTests.cpp)
target_link_libraries(boRoothaanHallTests ${TEST_LIBRARIES_HF})

add_test(NAME BoRoothaanHallTests COMMAND boRoothaanHallTests)
target_compile_options(boRoothaanHallTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(preBoRoothaanHallTests PreBoRoothaanHallTests.cpp)
target_link_libraries(preBoRoothaanHallTests ${TEST_LIBRARIES_HF})

add_test(NAME PreBoRoothaanHallTests COMMAND preBoRoothaanHallTests)
target_compile_options(preBoRoothaanHallTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

add_executable(iOTests IOTests.cpp)
target_link_libraries(iOTests ${TEST_LIBRARIES_HF})

add_test(NAME IOTests COMMAND iOTests)
target_compile_options(iOTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

set(TEST_INTEGRALS_LIBRARIES
        PRIVATE
        gtest_main
        gmock
        Scine::LibintIntegrals
        $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
        Libint2::cxx
        Boost::boost)

add_executable(integralsTests IntegralsTests.cpp)
target_link_libraries(integralsTests ${TEST_INTEGRALS_LIBRARIES})

add_test(NAME IntegralsTests COMMAND integralsTests)
target_compile_options(integralsTests PUBLIC
        $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)
