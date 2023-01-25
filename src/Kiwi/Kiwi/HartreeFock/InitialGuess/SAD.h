/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_SAD_H
#  define KIWI_SAD_H

#  include <Utils/DataStructures/BasisSet.h>
#  include <Utils/DataStructures/DensityMatrix.h>
#  include <Utils/DataStructures/MolecularOrbitals.h>
#  include <Utils/Geometry/AtomCollection.h>
#  include <map>

namespace Scine {
namespace Kiwi {

/**
 * @class SADGuess
 * @brief This class performs a superposition of atomic densities guess.
 * In my experience, this is the best available initial guess. It is implemented after:
 * https://pubs.rsc.org/en/content/articlehtml/2009/cp/b901987a
 */
class SADGuess {
 private:
  Utils::DensityMatrix D_;
  Utils::MolecularOrbitals C_;
  Eigen::VectorXd energies_;
  bool verbose_ = true;

  std::size_t size_;

 public:
  int getSize() const {
    return size_;
  }

 private:
  Utils::AtomCollection geometry_;
  bool isRestricted_;
  const Utils::Integrals::BasisSet& basis_;
  int dim_;

  std::vector<int> offest_;
  std::vector<std::vector<long>> shellIndices_;
  std::map<Utils::ElementType, int> orbitalsPerAtomMap_;
  std::map<Utils::ElementType, std::vector<int>> atomIndexMap_;
  std::map<Utils::ElementType, std::vector<double>> occNumMap_;
  std::map<Utils::ElementType, Eigen::MatrixXd> occNumMat_;
  std::map<Utils::ElementType, Utils::DensityMatrix> densityMatrixMap_;
  std::map<Utils::ElementType, Utils::MolecularOrbitals> mosMap_;
  std::map<Utils::ElementType, Eigen::MatrixXd> coreMatrixMap_;
  std::map<Utils::ElementType, Eigen::MatrixXd> overlapMap_;
  std::map<Utils::ElementType, Utils::Integrals::BasisSet> basisMap_;

  inline auto atomIndexToParticleType(int atomIndex) -> Utils::ElementType {
    return geometry_.at(atomIndex).getElementType();
  }

 public:
  SADGuess(Utils::AtomCollection geometry, bool isRestricted, const Utils::Integrals::BasisSet& basis);

  [[nodiscard]] auto getDensity() const -> Utils::DensityMatrix;

  [[nodiscard]] auto getMOs() const -> Utils::MolecularOrbitals;

  [[nodiscard]] auto getEnergies() const -> Eigen::VectorXd;

  auto setVerbose(bool verbose) -> void;

  auto computeInitialDensityMatrices() -> void;

  auto runScf() -> void;

 private:
  auto finish() -> void;

  auto numElectronsToOccupation(Utils::ElementType type) -> void;

  auto computeOccupationNumberVector() -> void;

  auto computeInitialDensityMatrixAt(Utils::ElementType type) -> void;

  auto evaluateProton(Utils::ElementType type) -> void;

  auto atomicScfAt(Utils::ElementType type) -> void;

  auto makeCoulomb(Utils::ElementType type) -> Eigen::MatrixXd;

  static auto makeExchange(Eigen::MatrixXd& coulomb) -> Eigen::MatrixXd;

  static auto buildFockMatrix(Eigen::MatrixXd& F, Eigen::MatrixXd& J, Eigen::MatrixXd& K, const Eigen::MatrixXd& coulomb,
                              const Eigen::MatrixXd& exchange, const Eigen::MatrixXd& H, Eigen::MatrixXd& D) -> void;

  static auto computeEnergy(const Eigen::MatrixXd& H, const Eigen::MatrixXd& F, const Eigen::MatrixXd& D) -> double;

  static auto updateDensity(const Eigen::MatrixXd& F, const Eigen::MatrixXd& S, const Eigen::MatrixXd& occupation,
                            Eigen::MatrixXd& D, Eigen::MatrixXd& C, Eigen::VectorXd& energies) -> void;

  [[nodiscard]] const std::map<Utils::ElementType, std::vector<double>>& getOccNumVec() const {
    return occNumMap_;
  }
  [[nodiscard]] const std::map<Utils::ElementType, Eigen::MatrixXd>& getOccNumMat() const {
    return occNumMat_;
  }
  [[nodiscard]] const std::vector<int>& getOffest() const {
    return offest_;
  }
  [[nodiscard]] const std::vector<std::vector<long>>& getShellIndices() const {
    return shellIndices_;
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_SAD_H

/**

*/
