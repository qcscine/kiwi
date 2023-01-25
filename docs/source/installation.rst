Installation
------------

Dependencies
............

Required software, minimum required versions in parantheses, for this SCINE project are:

- A C++ compiler supporting the C++17 standard (at the moment only GCC is supported)
- Libint (2.7.2, versions < 2.7 are not compatible and versions > 2.7 are not tested)
- CMake (3.9)
- Boost (1.65.0)
- Eigen3 (3.3.2)

Installation
............

In order to compile this as a SCINE developer, execute the following
commands::

    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst ..
    make -j 4
    make test
    make install
