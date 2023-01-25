/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_MOLECULARORBITALSIO_H
#define KIWI_MOLECULARORBITALSIO_H

#include <string>

namespace Scine {

namespace Utils {
class MolecularOrbitals;
} // namespace Utils
namespace Kiwi {

/*!
 * @class MolecularOrbitalsIO @file MolecularOrbitalsIO.h
 * @brief to write molecular orbitals on disk or read them from disk, in binary format.
 * TODO: save space and only write and read triangular matrix (since symmetric)?
 * TODO: Make this class employ EigenMatrixIO? NB: then the memory layout would change
 */
class MolecularOrbitalsIO {
 public:
  /*! Write the density matrix to disk. */
  static void write(const std::string& filename, const Utils::MolecularOrbitals& m);

  /*! Read a density matrix from the disk. */
  static Utils::MolecularOrbitals read(const std::string& filename);
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_MOLECULARORBITALSIO_H
