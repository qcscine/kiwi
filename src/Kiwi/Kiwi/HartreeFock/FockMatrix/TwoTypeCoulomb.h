/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_TWOTYPECOULOMB_H
#define KIWI_TWOTYPECOULOMB_H

#include <Kiwi/KiwiUtils/Data.h>
#include <LibintIntegrals/LibintIntegrals.h>

namespace Scine {
namespace Kiwi {
namespace FockMatrix {

//! In-memory evaluation of the pre-BO contribution to the Fock matrix.
inline auto buildLmatrices(const std::map<Utils::ElementType, Utils::DensityMatrix>& D,
                           const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data,
                           Utils::ElementType type1, Utils::ElementType type2) -> std::array<Eigen::MatrixXd, 2> {
  //
  // Note: The contraction of the two-body tensor for different particle type_s, depends on the ordering
  //       of the particles: V_ijIJ != V_IJij.
  //       Hence: we must keep track of the ordering.
  //

  Eigen::MatrixXd L1;
  Eigen::MatrixXd L2;

  auto typeCorrectOrder1 = type1;
  auto typeCorrectOrder2 = type2;
  if (type1 > type2) {
    typeCorrectOrder1 = type2;
    typeCorrectOrder2 = type1;
  }
  {
    // TODO check the transpose
    Eigen::MatrixXd Dt = D.at(typeCorrectOrder1).restrictedMatrix().transpose();
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
    Eigen::VectorXd tmpVec = data->Coulomb[Data::getUnique(typeCorrectOrder1, typeCorrectOrder2)] * coefficientVector;
    L2 = Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(typeCorrectOrder2).LAO, molecule->at(typeCorrectOrder2).LAO);
  }
  // TODO check the transpose
  {
    Eigen::MatrixXd Dt = D.at(typeCorrectOrder2).restrictedMatrix().transpose();
    Dt.transposeInPlace();
    Eigen::RowVectorXd coefficientRowVector(Eigen::Map<Eigen::RowVectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
    Eigen::RowVectorXd tmpVec = coefficientRowVector * data->Coulomb[Data::getUnique(typeCorrectOrder1, typeCorrectOrder2)];
    Dt.transposeInPlace();
    L1 = Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(typeCorrectOrder1).LAO, molecule->at(typeCorrectOrder1).LAO);
  }
  if (type1 == typeCorrectOrder1) {
    return {L1, L2};
  }
  return {L2, L1};
}

//! Integral-direct evaluation of the pre-BO contribution to the Fock matrix.
inline auto buildLmatricesDirect(const std::map<Utils::ElementType, Utils::DensityMatrix>& D,
                                 const std::shared_ptr<Molecule>& molecule, Utils::ElementType type1,
                                 Utils::ElementType type2) -> std::array<Eigen::MatrixXd, 2> {
  Utils::Integrals::IntegralSpecifier specifier;

  auto tp1 = molecule->at(type1);
  const auto& densityMatrix1 = D.at(type1);
  auto& basis1 = molecule->at(type1).basisSet;

  auto tp2 = molecule->at(type2);
  const auto& densityMatrix2 = D.at(type2);
  auto& basis2 = molecule->at(type2).basisSet;
  specifier.typeVector = {tp1.typeInfo, tp2.typeInfo};
  specifier.op = Utils::Integrals::Operator::Coulomb;

  auto result =
      Integrals::LibintIntegrals::evaluateTwoBodyDirectPreBo(specifier, basis1, basis2, densityMatrix1, densityMatrix2);

  return {result.first, result.second};
}

} // namespace FockMatrix
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_TWOTYPECOULOMB_H
