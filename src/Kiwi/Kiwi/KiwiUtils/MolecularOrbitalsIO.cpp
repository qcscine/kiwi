/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularOrbitalsIO.h"
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <fstream>

namespace Scine {
namespace Kiwi {

using namespace std;

void MolecularOrbitalsIO::write(const std::string& filename, const Utils::MolecularOrbitals& orbitals) {
  /*
   * Format : - restricted / unrestricted
   *          - number atomic orbitals
   *          - number molecular orbitals
   *          - one or two matrices, depending on restricted or unrestricted
   */

  ofstream fout;
  fout.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);

  auto unrestricted = static_cast<int8_t>(orbitals.isUnrestricted());
  fout.write(reinterpret_cast<char*>(&unrestricted), sizeof(int8_t)); // NOLINT

  if (unrestricted != 0) {
    auto nAOs = static_cast<int32_t>(orbitals.alphaMatrix().rows());
    fout.write(reinterpret_cast<char*>(&nAOs), sizeof(int32_t)); // NOLINT
    auto nMOs = static_cast<int32_t>(orbitals.alphaMatrix().cols());
    fout.write(reinterpret_cast<char*>(&nAOs), sizeof(int32_t)); // NOLINT
    const auto& alpha = orbitals.alphaMatrix();
    const auto& beta = orbitals.betaMatrix();
    fout.write(reinterpret_cast<const char*>(alpha.data()), nAOs * nMOs * sizeof(double)); // NOLINT
    fout.write(reinterpret_cast<const char*>(beta.data()), nAOs * nMOs * sizeof(double));  // NOLINT
  }
  else {
    auto nAOs = static_cast<int32_t>(orbitals.restrictedMatrix().rows());
    fout.write(reinterpret_cast<char*>(&nAOs), sizeof(int32_t)); // NOLINT
    auto nMOs = static_cast<int32_t>(orbitals.restrictedMatrix().cols());
    fout.write(reinterpret_cast<char*>(&nAOs), sizeof(int32_t)); // NOLINT
    const auto& matrix = orbitals.restrictedMatrix();
    fout.write(reinterpret_cast<const char*>(matrix.data()), nAOs * nMOs * sizeof(double)); // NOLINT
  }
}

Utils::MolecularOrbitals MolecularOrbitalsIO::read(const std::string& filename) {
  /*
   * Format : - restricted / unrestricted
   *          - number atomic orbitals
   *          - number molecular orbitals
   *          - one or two matrices, depending on restricted or unrestricted
   */

  ifstream fin;
  fin.open(filename, ios_base::in | ios_base::binary);
  int8_t unrestricted8 = 0;
  int32_t nAOs = 0;
  int32_t nMOs = 0;
  fin.read(reinterpret_cast<char*>(&unrestricted8), sizeof(int8_t)); // NOLINT
  auto unrestricted = unrestricted8 != 0;
  fin.read(reinterpret_cast<char*>(&nAOs), sizeof(int32_t)); // NOLINT
  fin.read(reinterpret_cast<char*>(&nMOs), sizeof(int32_t)); // NOLINT

  Utils::MolecularOrbitals orbitals;
  if (unrestricted) {
    Eigen::MatrixXd alpha(nAOs, nMOs);
    Eigen::MatrixXd beta(nAOs, nMOs);
    fin.read(reinterpret_cast<char*>(alpha.data()), nAOs * nMOs * sizeof(double)); // NOLINT
    fin.read(reinterpret_cast<char*>(beta.data()), nAOs * nMOs * sizeof(double));  // NOLINT
    orbitals = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients(std::move(alpha), std::move(beta));
  }
  else {
    Eigen::MatrixXd matrix(nAOs, nAOs);
    fin.read(reinterpret_cast<char*>(matrix.data()), nAOs * nMOs * sizeof(double)); // NOLINT
    orbitals = Utils::MolecularOrbitals::createFromRestrictedCoefficients(std::move(matrix));
  }

  return orbitals;
}

} // namespace Kiwi
} // namespace Scine
