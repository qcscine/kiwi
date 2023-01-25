How to use Kiwi 
==================

Advanced Usage
--------------

This section explains the advanced usage of the Kiwi program, which boils down
to the individual tasks that can be executed in the ordering specified in the
input file. For example, we might want to do a Hartree--Fock optimization and
afterwards print the FCIDUMP file. The ``tasks`` block would then look like:

.. code-block:: yaml
  
  tasks:
    - type: hf
      guess: read
    - type: ao2mo
      write: true

In the following sections, we will explain all available tasks in detail.

The available tasks are:

+--------------------------+---------------------------------------------------------------+
| Keyword                  | Task                                                          |
+==========================+===============================================================+
| ``hf``                   | Perform a Hartree--Fock orbital optimization                  |  
+--------------------------+---------------------------------------------------------------+
| ``ao2mo``                | Perform the AO to MO transformation                           |    
+--------------------------+---------------------------------------------------------------+
| ``density``              | Evaluate the particle density or molecular orbitals on a grid |
+--------------------------+---------------------------------------------------------------+
| ``natural orbitals``     | Calculate natural orbitals from RDMs that are read from disk  |
+--------------------------+---------------------------------------------------------------+

