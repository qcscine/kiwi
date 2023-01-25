Tasks
==================

Hartree--Fock
--------------

.. code-block:: yaml

  tasks:
      - type: hf

The Hartree--Fock task is chosen by entering ``hf`` as the ``type``. 

+-----------------------------+----------------------------------------------------+------------+
| Keywords                    | Value type                                         | Default    | 
+=============================+====================================================+============+
| ``print keywords``          | ``bool``                                           | ``false``  |
+-----------------------------+----------------------------------------------------+------------+
| ``scf type``                | ``serial``, ``bfgs``, ``alternating``, ``trah``    | ``trah``   |
+-----------------------------+----------------------------------------------------+------------+
| ``scf``                     | ``bool``                                           | ``true``   |
+-----------------------------+----------------------------------------------------+------------+
| ``trah settings``           | --> see according section                          |            |
+-----------------------------+----------------------------------------------------+------------+
| ``max iter``                | ``int``                                            | ``500``    |
+-----------------------------+----------------------------------------------------+------------+
| ``guess``                   | ``sad``, ``hueckel``, ``core``, ``read``, ``bo``   | ``sad``    |
+-----------------------------+----------------------------------------------------+------------+
| ``nuclear guess``           | ``snd``, ``core``                                  | ``snd``    |
+-----------------------------+----------------------------------------------------+------------+
| ``integral direct``         | ``bool``                                           | ``true``   |
+-----------------------------+----------------------------------------------------+------------+
| ``reset incremental``       | ``int``                                            | ``10``     |
+-----------------------------+----------------------------------------------------+------------+
| ``disable incremental``     | ``bool``                                           | ``true``   |
+-----------------------------+----------------------------------------------------+------------+
| ``max iter nested``         | ``int``                                            | ``20``     |
+-----------------------------+----------------------------------------------------+------------+
| ``thresh``                  | ``float``                                          | ``1e-9``   |
+-----------------------------+----------------------------------------------------+------------+
| ``loewdin thresh``          | ``float``                                          | ``1e-9``   |
+-----------------------------+----------------------------------------------------+------------+
| ``thresh fock``             | ``float``                                          | ``1e-5``   |
+-----------------------------+----------------------------------------------------+------------+
| ``mix angle``               | ``float``                                          | ``0.2``    |
+-----------------------------+----------------------------------------------------+------------+
| ``mixes``                   | ``int``                                            | ``10``     |
+-----------------------------+----------------------------------------------------+------------+
| ``perturbations``           | ``int``                                            | ``0``      |
+-----------------------------+----------------------------------------------------+------------+
| ``accelerator``             | ``ediis diis``, ``ediis``, ``diis``, ``none``      | ``diis``   |
+-----------------------------+----------------------------------------------------+------------+
| ``diis thresh``             | ``float``                                          | ``1.0``    |
+-----------------------------+----------------------------------------------------+------------+
| ``nuclear diis thresh``     | ``float``                                          | ``1.0``    |
+-----------------------------+----------------------------------------------------+------------+
| ``ediis thresh``            | ``float``                                          | ``20.0``   |
+-----------------------------+----------------------------------------------------+------------+
| ``nuclear ediis thresh``    | ``float``                                          | ``20.0``   |
+-----------------------------+----------------------------------------------------+------------+
| ``stability analysis``      | ``bool``                                           | ``false``  |
+-----------------------------+----------------------------------------------------+------------+

In the following, we will explain some of the meanings of the keywords.


SCF type
...............

The SCF type keyword describes which optimization strategy is used. By default,
the trust region augmented Hessian, ``trah``, method is chosen since it is a
black-box method and the settings rarely have to be changed to guarantee
convergence. The TRAH method can be used with different optimizers,
which are the exact Newton method and the Augmented Roothaan--Hall (ARH) method.
The details are explained below. 
The TRAH method is in our experience the most reliable strategy for electronic
and nuclear-electronic calculations.

There is one exception, in the nuclear-electronic case, for very small systems,
e.g, up to 4 or 5 atoms, we recommend to choose the ``serial`` SCF type which is
implemented after Aodong Liu et al.

  Liu, Aodong, et al. "Simultaneous Optimization of Nuclear–Electronic Orbitals." The Journal of Physical Chemistry A 126.39 (2022): 7033-7039.

The ``bfgs`` strategy is usually not recommended. It follows the strategy of
Van Voorhis and Head-Gordon:

  Van Voorhis, Troy, and Martin Head-Gordon. "A geometric approach to direct minimization." Molecular Physics 100.11 (2002): 1713-1721.

The ``alternating`` SCF strategy is an outdated strategy for optimizing
nuclear-electronic Hartree--Fock wave functions and is also never recommended.


Do SCF?
..........

The ``scf`` key can be set to ``false`` such that no SCF optimization is performed.
It is used, in case there are orbitals stored on disk that should be used
exactly as they are. This enables reproducible FCIDUMP files or can be used in
conjunction with natural orbitals that are constructed from an external RDM,
stored on disk.


Guess
...........

The ``guess`` is, in principle, self-explanatory. ``sad``, ``hueckel``, and
``core``, are standard guesses for the electronic part and ``bo`` is used for
nuclear-electronic calculations. In this case, first a BO-SCF is performed and
the resulting electronic density matrix is used as a guess for the
nuclear-electronic optimization. This is always recommended because it can
drastically reduce the overall computation time.
``read`` can be selected if orbitals that are stored on disk should be used.
In this case, the base-name of the orbital files must be the
same as the input file name minus the YAML file extension.

The ``nuclear guess`` should always be set to ``snd`` (superposition of
nuclear densities). The ``core`` guess works extremely badly for nuclei,
Details are given here:

  arXiv preprint arXiv:2210.10170 (accepted in JCTC)


Integral settings
...................

``integral direct``: if this keyword is set to ``false``, all integrals will be
stored in memory. This gives a huge speed-up for small molecules, but because of
the O(N\ :sup:`4`\) scaling of the memory requirement, it becomes unfeasible
very quickly. Also, the advantage is lost for larger systems, because no integral
screening can be used.

``disable incremental``: if this keyword is set to ``true``, the Fock matrix is
constructed incrementally and is re-build from scratch all ``reset incremental``
iterations.


Thresholds
............

``thresh``: energy convergence threshold.

``thresh fock``: Fock gradient convergence threshold.

``loewdin thresh``: threshold for discarding nearly linear independent orbitals
in the Loewdin orthogonalization.

``diis thresh``, ``ediis thresh``: thresholds when DIIS or EDIIS should be
activated, depending on the Fock gradient. 

**Important note:** they are only used in the electronic ``serial`` and the
nuclear-electronic ``alternating`` cases.

``nuclear diis thresh``, ``nuclear ediis thresh``: thresholds when DIIS or EDIIS should be
activated, depending on the Fock gradient. 

**Important note:** they are only used in the nuclear-electronic ``serial``.


Orbital steering
...................

Orbital steering can be used in cases where the user suspects that not the
global, but a local minimum is found. The details are given in

  Vaucher, Alain C., and Markus Reiher. "Steering orbital optimization out of local minima and saddle points toward lower energy." Journal of Chemical Theory and Computation 13.3 (2017): 1219-1228.

Associated keywords are ``perturbations``, ``mixes``, ``mix angle``. 
We recommend only changing the parameter ``perturbations``. E.g., by setting it
to ``1`` the orbitals are mixed once at the beginning of the SCF. This is
recommended for optimizing unrestricted electronic Hartree--Fock wave functions.


Stability analysis
....................

If the ``stability analysis`` is enabled, the exact Hessian is diagonalized
after the optimization with the Davidson method, and the lowest eigenvalue is
printed. We note here that this can take many iterations in the
nuclear-electronic case and may not converge at all. If the lowest eigenvalue is
less than 0, the optimization converged to a saddle point. In this case, try the
ARH method or even the Newton method instead.


The trust region augmented Hessian (TRAH) method
...................................................

All details regarding the TRAH method can be found in

  arXiv preprint arXiv:2210.10170 (accpeted in JCTC)

We provide all keywords below, although we rarely recommend to change them, with
a few exceptions.

+-----------------------------------+-----------------------+--------------+
| Keywords                          | Value types           | Default      |
+===================================+=======================+==============+
| ``print keywords``                | ``bool``              | ``false``    |
+-----------------------------------+-----------------------+--------------+
| ``optimizer``                     | ``ARH``, ``Newton``   | ``ARH``      |
+-----------------------------------+-----------------------+--------------+
| ``initial trust radius``          | ``float``             | ``0.5``      |
+-----------------------------------+-----------------------+--------------+
| ``max trust radius``              | ``float``             | ``1.0``      |
+-----------------------------------+-----------------------+--------------+
| ``max davidson iterations``       | ``int``               | ``32``       |
+-----------------------------------+-----------------------+--------------+
| ``max davidson subspace dim``     | ``int``               | ``70``       |
+-----------------------------------+-----------------------+--------------+
| ``grad scaling``                  | ``float``             | ``0.1``      |
+-----------------------------------+-----------------------+--------------+
| ``min thresh``                    | ``float``             | ``0.01``     |
+-----------------------------------+-----------------------+--------------+
| ``local thresh``                  | ``float``             | ``0.001``    |
+-----------------------------------+-----------------------+--------------+
| ``max arh dim``                   | ``int``               | ``25``       |
+-----------------------------------+-----------------------+--------------+
| ``2nd start vector``              | ``bool``              | ``true``     |
+-----------------------------------+-----------------------+--------------+
| ``2nd start vector noise``        | ``bool``              | ``false``    |
+-----------------------------------+-----------------------+--------------+
| ``3rd start vector``              | ``bool``              | ``true``     |
+-----------------------------------+-----------------------+--------------+
| ``3rd start vector noise``        | ``bool``              | ``false``    |
+-----------------------------------+-----------------------+--------------+

Keywords we recommend adjusting occasionally:

``optimizer``: in cases that are extremely difficult to converge, the ``Newton``
optimizer may be more reliable than the ``ARH`` optimizer. 

**Important note:** the nuclear-electronic integral-direct version is not optimized for
speed.

``initial trust radius``, ``max trust radius``: for small systems, e.g., up to
4 or 5 atoms, those two values can be decreased.

``max davidson iterations``: for strongly correlated electronic systems when the
Newton method is used this number should be decreased to 16.

``2nd start vector noise``, ``3rds start vector noise``: in the case of an
unrestricted electronic optimization, and if no orbital steering is used, this
can help to break the restricted symmetry of the initial guess and steer the
optimization to a lower minimum.


AO to MO transformation
-------------------------

The atomic orbital (AO) to molecular orbital (MO) transformation is activated
with the ``ao2mo`` value for the ``type`` keyword:

.. code-block:: yaml

  tasks:
      - type: hf


Available options are further:

``write``: if set to ``true``, the integrals are written on disk.

``thresh``: threshold for screening out integrals, by default it is set to ``1e-16``.

By default, in the electronic case, the FCIDUMP file is written in the
conventional electronic format.
In the nuclear-electronic case, we adapt the QCMaquis format, as explained in
the QCMaquis manual (https://github.com/qcscine/qcmaquis).

**Important note:** at the moment, only restricted electrons with nuclei in the high-spin
approximation are supported.

In short: the first line, with only a single float, contains a constant shift in
energy, e.g., point-charge repulsion of classical nuclei. The following lines
contain the terms in the Hamiltonian where pairs of numbers ``m-n`` describe a
single second-quantization operator. The first number denotes the particle type
(starting with 0) and the second number is the orbital index associated with the
given particle type (also starting with 0 for all types).
We assume physics notation and normal-ordering


Natural orbitals
-----------------

With this task, the user can read RDM files from disk that come, e.g., from a
correlated calculation (DMRG, CI, ...) and the corresponding natural orbitals
are constructed. 
It is enabled with the ``natural orbitals`` option for the ``type`` keyword.
The ``rdm files`` key must be followed by a particle type symbol, e.g., ``e``,
``H`` followed by the rdm file names. The files must be given in the csv format.

One can provide one rdm file for each spin and list them separated by a white
space in the same line (starting with alpha spin).

.. code-block:: yaml

  tasks:
      - type: natural orbitals
        rdm files:
          e: natural_orbitals.e.csv
          H: natural_orbitals.alpha.H.csv natural_orbitals.beta.H.csv


Particle densities and molecular orbitals
------------------------------------------

The ``density`` task enables the evaluation of molecular orbitals or particle
densities on a grid. The output format, ``format``, can be a cube file,
``cube``,
(http://paulbourke.net/dataformats/cube/) or a KIWI-specific format that
contains the three XYZ coordinates of the grid point followed by the value (of
the MO or the density) on a line, ``grid``.

The possible options are


+-----------------------+-----------------------+--------------+
| Keywords              | Value types           | Default      |
+=======================+=======================+==============+
| ``mode``              | ``density``, ``MO``   | REQUIRED     |
+-----------------------+-----------------------+--------------+
| ``format``            | ``grid``, ``cube``    | REQUIRED     |
+-----------------------+-----------------------+--------------+
| ``print keywords``    | ``bool``              | ``false``    |
+-----------------------+-----------------------+--------------+
| ``particle type``     | Particle type symbol  | REQUIRED     |
+-----------------------+-----------------------+--------------+
| ``xmin``              | ``float``             | ``0.0``      |
+-----------------------+-----------------------+--------------+
| ``xmax``              | ``float``             | ``0.0``      |
+-----------------------+-----------------------+--------------+
| ``ymin``              | ``float``             | ``0.0``      |
+-----------------------+-----------------------+--------------+
| ``ymax``              | ``float``             | ``0.0``      |
+-----------------------+-----------------------+--------------+
| ``zmin``              | ``float``             | ``0.0``      |
+-----------------------+-----------------------+--------------+
| ``zmax``              | ``float``             | ``0.0``      |
+-----------------------+-----------------------+--------------+
| ``N``                 | ``int``               | ``100``      |
+-----------------------+-----------------------+--------------+
| ``NX``                | ``int``               | ``N``        |
+-----------------------+-----------------------+--------------+
| ``NY``                | ``int``               | ``N``        |
+-----------------------+-----------------------+--------------+
| ``NZ``                | ``int``               | ``N``        |
+-----------------------+-----------------------+--------------+
| ``rdm file``          | csv file name         | OPTIONAL     |
+-----------------------+-----------------------+--------------+
| ``rdm file alpha``    | csv file name         | OPTIONAL     |
+-----------------------+-----------------------+--------------+
| ``rdm file beta``     | csv file name         | OPTIONAL     |  
+-----------------------+-----------------------+--------------+
| ``mo index``          | ``int``               | ``0``        |
+-----------------------+-----------------------+--------------+
| ``ms index``          | ``int``               | ``0``        |
+-----------------------+-----------------------+--------------+

The ``mode`` keyword signals whether a molecular orbital, ``MO``, or the
particle density, ``density``, is to be evaluated.

Both modes require a ``particle type`` to be specified. In the MO mode, further
the ``mo index`` (starting from ``0``) should be given and also the ``ms index``,
where ``0`` stands for alpha and ``1`` for beta.

In addition, one can provide an external density matrix (RDM). Then, instead of
the Hartree--Fock density matrix, the particle density can be evaluated with a
correlated denstiy, e.g., from a DMRG calculation. The RDM should again be given
in the csv format.

**Important note:** ``rdm file`` is used for restricted particle types only.

The keyword ``N`` is used to specify the size of the grid in either all three
dimensions or the size can also be specified for each direction separately. 
The dimensions of the box can be given with ``(i)min``, and ``(i)max``.

**Important note:** the input and output units are Angstrom. 
