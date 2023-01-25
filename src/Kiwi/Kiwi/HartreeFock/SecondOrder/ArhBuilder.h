/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_ARHBUILDER_H
#define KIWI_ARHBUILDER_H

#include <Kiwi/KiwiUtils/SymmetryEnums.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/DataStructures/SpinAdaptedVector.h>
#include <Utils/Geometry/ElementTypes.h>
#include <deque>
#include <map>
#include <memory>
#include <vector>

namespace Scine {
namespace Kiwi {

class Molecule;
class Data;

struct ArhData {
  int size;
  int maxSize;

  std::shared_ptr<Molecule> molecule;
  std::shared_ptr<Data> data;

  // Updated in every evaluation:
  std::map<Utils::ElementType, Utils::DensityMatrix> Xi;
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> HXi;
  // Gradient
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> gradient;
  // The fock matrix is required for Gradient and Roothaan-Hall part.
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> F_OAO;

  // Intermediate1_j := Tr ( [Dn, X] D_jn )
  std::map<Utils::ElementType, Utils::SpinAdaptedVector> Intermediate1;

  // Intermediate2_i = Sum_j [T^-1]_ij I1_i
  std::map<Utils::ElementType, Utils::SpinAdaptedVector> Intermediate2;

  // T_ij := Tr ( D_in D_jn)
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> T;

  /*
   * The G and I matrices are required for the augmented Roothaan-Hall part.
   * Note that in the BO procedure, instead for G and I, F is used, but
   * the more general pre-BO procedure differs slightly in this respect.
   */

  std::map<Utils::ElementType, std::deque<Utils::SpinAdaptedMatrix>> G;
  std::map<Utils::ElementType, std::deque<Utils::SpinAdaptedMatrix>> J;
  std::map<Utils::ElementType, std::deque<std::map<Utils::ElementType, Eigen::MatrixXd>>> L;
  std::map<Utils::ElementType, std::deque<Utils::DensityMatrix>> D_OAO;
  // D_in := D_i - D_n
  std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>> D_in;

  // Containers that are re-evaluated at every iteration:
  //// G_in := G_i - G_n
  // std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>> G_in;
  //// G_ov_vo_in := P_o G_in P_v - P_v G_in P_o
  // std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>> G_ov_vo_in;
  //// J_in := J_i - J_n
  // std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>> J_in;
  //// J_ov_vo_in := P_o J_in P_v - P_v J_in P_o
  // std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>> J_ov_vo_in;
  //// L_in := L_i - L_n
  // std::map<Utils::ElementType, std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>>> L_in;
  //// L_ov_vo_in := P_o L_in P_v - P_v L_in P_o
  // std::map<Utils::ElementType, std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>>> L_ov_vo_in;
};

class ArhBuilder {
 public:
  void setUseExactHessian(bool useExactHessian);

 private:
  bool useExactHessian = false;

  bool debug = false;

  const int dim;

  const SpinFunction spin;

  const Utils::ElementType type;

  const std::shared_ptr<ArhData> arhData;

  // Difference between virt-virt and occ-occ contribution: F_vv - F_oo
  Utils::SpinAdaptedMatrix F_vv_oo;
  // G_in := G_i - G_n
  std::vector<Utils::SpinAdaptedMatrix> G_in;
  // G_ov_vo_in := P_o G_in P_v - P_v G_in P_o
  std::vector<Utils::SpinAdaptedMatrix> G_ov_vo_in;
  // J_in := J_i - J_n
  std::vector<Utils::SpinAdaptedMatrix> J_in;
  // J_ov_vo_in := P_o J_in P_v - P_v J_in P_o
  std::vector<Utils::SpinAdaptedMatrix> J_ov_vo_in;
  // L_in := L_i - L_n
  std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>> L_in;
  // L_ov_vo_in := P_o L_in P_v - P_v L_in P_o
  std::map<Utils::ElementType, std::vector<Utils::SpinAdaptedMatrix>> L_ov_vo_in;

  Eigen::MatrixXd kappa_;

 public:
  ArhBuilder(int dimension, Utils::ElementType tp, SpinFunction sp, std::shared_ptr<ArhData> dat);

  auto evaluateGradient() -> void;

  auto getTransformationMatrix() -> Eigen::MatrixXd;

  auto evaluate() -> Eigen::VectorXd;

  auto addMOs() -> void;

  [[nodiscard]] auto getGradient() const -> const Eigen::MatrixXd&;

  [[nodiscard]] auto getRhDiagonal(bool print = false) const -> Eigen::DiagonalMatrix<double, Eigen::Dynamic>;

  auto getGradientVector() -> Eigen::VectorXd;

  [[nodiscard]] auto getThirdTrialVector() const -> Eigen::VectorXd;

  auto makeAllContributions() -> void;

  auto makeXindependentContributions() -> void;

  auto makeXdependentContributions() -> void;

  auto updateX(Eigen::VectorXd& vectorX) -> void;

  void setMu(double m);

  void addLevelShift(bool addLS);

 private:
  bool addLevelShift_ = false;

  double mu = 0;

  Eigen::MatrixXd HXi;

  Eigen::MatrixXd G;

  int occ;
  int virt;

  Eigen::MatrixXd Theta_virt;
  Eigen::MatrixXd Theta_occ;
  Eigen::MatrixXd Theta;

  auto makeDin() -> void;
  auto makeIntermediate1() -> void;
  auto makeIntermediate2() -> void;
  //!!! Attention: addRH contribution 1st
  auto addRhContribution() -> void;
  //!!! Attention: add TRAH specific contribution after RH contribution
  auto addBoArhSpecificContribution() -> void;
  auto addUnrestrictedArhSpecificContribution() -> void;
  auto addPreBoArhSpecificContribution() -> void;

  //!!! Attention: add TRAH specific contribution after RH contribution
  auto addBo2ndOrderSpecificContribution() -> void;
  auto addPreBO2ndOrderSpecificContribution() -> void;

  auto makeFvvoo() -> void;
  auto makeT() -> void;
  auto makeGinVector() -> void;
  auto makeJinVector() -> void;
  auto makeLinVector() -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_ARHBUILDER_H
