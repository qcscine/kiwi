/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_SIMULTANEOUSDIIS_H
#define KIWI_SIMULTANEOUSDIIS_H

#include <Kiwi/KiwiUtils/Data.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Scf/ConvergenceAccelerators/DiisError.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Kiwi {

/*!
 * Class performing the calculation of a Simultaneous matrix based on the direct inversion of the iterative subspace
 * (DIIS) algorithm.
 */
class SimultaneousDiis {
 public:
  SimultaneousDiis();

  SimultaneousDiis(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data)
    : molecule_(std::move(molecule)), data_(std::move(data)) {
  }

  void setSubspaceSize(int n);

  void setOrthonormal(bool useOaoBasis) {
    useOrthonormalBasis = useOaoBasis;
  }

  void addMatrices();
  void restart();

  auto setMixedFockMatrix() -> void;

 private:
  bool useOrthonormalBasis = true;

  std::shared_ptr<Molecule> molecule_;
  std::shared_ptr<Data> data_;

  void resizeMembers();
  void updateBMatrix();

  int subspaceSize_ = 5;

  int index_;
  int lastAdded_;
  int iterationNo_;

  std::vector<std::map<Utils::ElementType, Utils::SpinAdaptedMatrix>> fockMatrices;
  std::vector<std::map<Utils::ElementType, Utils::SpinAdaptedMatrix>> errorMatrices;

  Eigen::MatrixXd B;
  Eigen::VectorXd rhs;
  Eigen::VectorXd C;

  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> calculateLinearCombination();

  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> getMixedFockMatrix();

  Utils::SpinAdaptedMatrix calculateLinearCombination(const Utils::ElementType type);

  auto evalauteBMatrixElement(int i, int j) -> double;

  auto evalauteError(int idx) -> void;

  static Eigen::MatrixXd calculateErrorMatrixOrthonormal(const Eigen::MatrixXd& D, const Eigen::MatrixXd& F);

  static Eigen::MatrixXd calculateErrorMatrix(const Eigen::MatrixXd& D, const Eigen::MatrixXd& F, const Eigen::MatrixXd& S);
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_SIMULTANEOUSDIIS_H
