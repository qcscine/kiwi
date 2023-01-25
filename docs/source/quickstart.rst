Quickstart
==================

Kiwi executable
-----------------

The Kiwi executable is called ``kiwiBinary`` and is located in the build
directory at: ``src/Kiwi/App/kiwiBinary``. If Kiwi is just compiled and not
installed, the program is executed according to

.. code-block:: bash

  $ /path/to/build/directory/src/Kiwi/App/kiwiBinary


YAML input
--------------

The Kiwi program relies entirely on YAML input files with the file extension ``.yml``.
The input file is provided as a second argument following the Kiwi binary, 
``kiwiBinary``:

.. code-block:: bash

  $ kiwiBinary <input.yml>

A minimal ``.yml`` input file requires to have several blocks where the ordering does
not play a role. 


XYZ file
.............

The first block is the XYZ block:

.. code-block:: yaml

  xyz file: /path/to/molecule.xyz

The XYZ files of Kiwi can have an additional ``Q`` or ``q`` at the end of the
line that gives the coordinates of an atom. This signals whether the nucleus
should be treated quantum mechanically. Below we give the example of the XYZ
file for the HeHHe\ :sup:`+`\  system where the proton is treated quantum
mechanically:


.. code-block:: 

  3
  // comment line  
  H 0.0 0.0 0.0 Q
  He 0.9236 0.0 0.0
  He -0.9236 0.0 0.0


Molecule
.............

The second block is the molecule block. It is, in principle, optional but most
of the time the user wants to specify some parameters that can be set here. The
keywords and their options are given in the table below:

+----------------------+-------------+--------------------------------------------+
| Keyword              | Value  type | Explanation                                |
+======================+=============+============================================+
| ``mult``             | ``int``     | multiplicity of the electrons              |
+----------------------+-------------+--------------------------------------------+
| ``charge``           | ``int``     | charge of the molecule                     |
+----------------------+-------------+--------------------------------------------+
| ``restricted``       | ``bool``    | ``true`` for restricted HF                 |
+----------------------+-------------+--------------------------------------------+
| ``high spin``        | ``bool``    | ``true`` for high-spin approx. for nuclei  |
+----------------------+-------------+--------------------------------------------+
| ``pure spherical``   | ``bool``    | ``true`` for spherical harmoic Gaussian    |
+----------------------+-------------+--------------------------------------------+
| ``uncontract``       | ``bool``    | ``true`` for uncontracting the basis set   |
+----------------------+-------------+--------------------------------------------+

The high-spin approximation here means that all nuclei are treated as alpha
particles. If this key is set to ``false`` the spin-multiplicity is
automatically minimized for all nuclei. At the moment it is not possible to set
a specific multiplicity for the nuclei.

Basis
.............

This block selects the basis set for the electrons.
It can either be specified in a single line as

.. code-block:: yaml

  basis: cc-pVTZ

which selects the same basis set for all element types. 
Alternatively, the user can provide the symbol for the element types contained
in the molecule and separately specify the basis set.

.. code-block:: yaml

  basis: 
      H: cc-pVQZ
      He: cc-pVDZ


**Note:** Check the directory ``data/basis`` for available basis sets.
The user can manually add basis sets in the ``g94`` format in this directory and
they will be automatically available.

Nuclear basis
...............

This block selects the basis set for the nuclei.
In this case, the user must provide the symbol for the quantum mechanical
nucleus and select a basis.
At the moment, only basis sets for protons are available.

.. code-block:: yaml

  nuclear basis: 
      H: PB4-D 


**Note:** Check the directory ``data/basis`` for available basis sets.
The user can manually add basis sets in the ``g94`` format in this directory and
they will be automatically available.

Tasks
.............

The ``tasks`` block controls the workflow of the program. In the quickstart-section, we will familiarize the user with the Hartree--Fock task and refer to
the detailed documentation to the task-specific section.

The tasks will be executed in the order provided in the input file. Hence, the
user should always start with the Hartree--Fock task, otherwise no orbitals will
be constructed.

In the most simple example, we can just specify the Hartree--Fock task and
choose the default options. The task that we want to execute is always specified
by the ``type`` keyword.


.. code-block:: yaml

  tasks: 
      - type: hf 


Complete input file example
..............................

Now, the user just has to assemble all blocks into a complete input file and
execute the program.

.. code-block:: yaml

  molecule:
    charge: 1
    mult: 1
  tasks:
    - type: hf
  xyz file: hehhe+.xyz
  basis: def2-SVP
  nuclear basis:
    H: PB4-D

After the program executes the SCF optimization, it will write orbital files for
each particle type to the disk. Those orbitals can be used by specifying the 
``guess`` block in the ``hf`` task with the parameter ``read``:


.. code-block:: yaml

  tasks: 
    - type: hf 
      guess: read

Note that there is no name specified for the orbital file. Instead, it uses the file name of the input file. 


Additional nuclear center
..........................

Additionally, we provide the option to add more centers at which nuclear basis
functions can be added, e.g., for the malonaldehyde molecule, the user might
want to place basis function in both wells. This can be achieved by specifying
the ``xyz nuclear centers`` keyword, followed by another XYZ file. Note that the
``Q`` at the end of the line is not necessary here. 

