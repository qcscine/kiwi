/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_SND_H
#define KIWI_SND_H
#include <Utils/DataStructures/BasisSet.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Geometry/AtomCollection.h>
#include <map>

namespace Scine {
namespace Kiwi {

/**
 * @class SNDGuess
 * @brief This class performs the superposition of nuclear densities (SND) initial guess.
 *        This is the pre-BO analogue to the SAD guess with the differnce that we do not have to take care of
 *        the spherical symmetry.
 */
class SNDGuess {
 private:
  Utils::DensityMatrix D_;
  Utils::MolecularOrbitals C_;
  bool verbose_ = true;

 public:
  int getSize() const {
    return size_;
  }

 private:
  Utils::AtomCollection positions_;
  Utils::AtomCollection atom_;
  Utils::Integrals::BasisSet nuclearBasis_;
  Utils::Integrals::BasisSet electronicBasis_;
  int dim_;
  std::size_t size_;
  Utils::ElementType type_;
  int n_alpha;
  int n_beta;

  std::vector<int> offest_;

  auto makeBasis(const Utils::Integrals::BasisSet& nuclearBasis, const Utils::Integrals::BasisSet& electronicBasis) -> void;

 public:
  SNDGuess(Utils::AtomCollection positions, const Utils::Integrals::BasisSet& nuclearBasis,
           const Utils::Integrals::BasisSet& electronicBasis, int n_alpha, int n_beta);

  [[nodiscard]] auto getDensity() const -> Utils::DensityMatrix;

  [[nodiscard]] auto getMOs() const -> Utils::MolecularOrbitals;

  auto setVerbose(bool verbose) -> void;

  auto runScf() -> void;

 private:
  auto finish() -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_SND_H
