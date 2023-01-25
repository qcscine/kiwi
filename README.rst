SCINE - Kiwi
=============================

Introduction
------------

The Kiwi program is an electronic and nuclear-electronic Hartree--Fock code
with advaned SCF algorithms with some additional functionalities:

- Born--Oppenheimer (BO) and nuclear-electronic (NE) first order and, exact and approximate, second order SCF algorithms.
- Storing of an FCIDUMP file on disk.
- Calculation of molecular orbitals (MOs) and particle densities.

License and Copyright Information
---------------------------------

For license and copyright information, see the file ``LICENSE.txt`` in this
directory.

Dependencies
----------------------

Required software, minimum required versions in brackets, for this SCINE project are:

- A C++ compiler supporting the C++17 standard (at the moment only GCC is supported)
- Libint 2.7.2 (other versions are not tested) [Will be automatically downloaded if not found]
- CMake (3.9)
- Boost (1.65.0)
- Eigen3 (3.3.2)

Installation
----------------------

In order to compile Kiwi, execute the following
commands:

.. code-block:: bash

    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst ..
    make -j 4

To complete setting up Kiwi, the following bash variable has to be set such that
libint can find the correct basis sets:

.. code-block:: bash

    export LIBINT_DATA_PATH=/path/to/kiwi/data

Now, the tests can be executed and kiwi can be installed:

.. code-block:: bash

    make test
    make install


Manual
------

The manual can be compiled with `Sphinx <https://www.sphinx-doc.org>`_ according to:

.. code-block:: bash

    cd docs
    make html
    firefox build/html/index.html

Alternatively, a PDF can be obtained:

.. code-block:: bash

    cd docs
    make latexpdf


Known Issues
------------

Intel compilers do not work due to libint at the moment.

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message
to scine@phys.chem.ethz.ch.

Third-Party Libraries Used
--------------------------

SCINE Kiwi makes use of the following third-party libraries:

- `Boost <https://www.boost.org/>`_
- `Eigen <http://eigen.tuxfamily.org>`_
- `Google Test <https://github.com/google/googletest>`_
- `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_
- `libint <https://github.com/evaleev/libint>`_

