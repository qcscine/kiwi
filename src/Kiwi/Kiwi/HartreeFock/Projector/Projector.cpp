/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/HartreeFock/Projector/Projector.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <iomanip>

namespace Scine {
namespace Kiwi {

auto Projector::updateDensity(Utils::ElementType type) -> void {
  auto const& X = data_->X[type];
  Utils::SpinAdaptedMatrix& Fortho = data_->hartreeFockData->F_OAO[type];

  if (molecule_->at(type).isRestricted) {
    Utils::MolecularOrbitals C;
    Utils::MolecularOrbitals C_OAO;
    es_.compute(Fortho.restrictedMatrix());
    Eigen::MatrixXd tmp = X * es_.eigenvectors();
    Eigen::MatrixXd tmp_oao = es_.eigenvectors();
    C = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp));
    C_OAO = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp_oao));
    data_->C.at(type) = C;
    data_->C_OAO.at(type) = C_OAO;
    data_->D.at(type) = HartreeFockUtils::makeDensity(type, molecule_, C);
    data_->D_OAO.at(type) = HartreeFockUtils::makeDensity(type, molecule_, C_OAO);
  }
  else {
    Eigen::MatrixXd C_alpha = Eigen::MatrixXd::Zero(data_->X[type].rows(), data_->X[type].cols());
    Eigen::MatrixXd C_beta = Eigen::MatrixXd::Zero(data_->X[type].rows(), data_->X[type].cols());
    Eigen::MatrixXd C_alpha_OAO = Eigen::MatrixXd::Zero(data_->X[type].cols(), data_->X[type].cols());
    Eigen::MatrixXd C_beta_OAO = Eigen::MatrixXd::Zero(data_->X[type].cols(), data_->X[type].cols());
    if (molecule_->at(type).msVector[0] > 0) {
      es_.compute(Fortho.alphaMatrix());
      C_alpha = X * es_.eigenvectors();
      C_alpha_OAO = es_.eigenvectors();
    }
    if (molecule_->at(type).msVector[1] > 0) {
      es_.compute(Fortho.betaMatrix());
      C_beta = X * es_.eigenvectors();
      C_beta_OAO = es_.eigenvectors();
    }
    auto C = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Eigen::MatrixXd>(std::move(C_alpha),
                                                                                           std::move(C_beta));
    auto C_OAO = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Eigen::MatrixXd>(std::move(C_alpha_OAO),
                                                                                               std::move(C_beta_OAO));
    data_->C.at(type) = C;
    data_->C_OAO.at(type) = C_OAO;
    data_->D.at(type) = HartreeFockUtils::makeDensity(type, molecule_, C);
    data_->D_OAO.at(type) = HartreeFockUtils::makeDensity(type, molecule_, C_OAO);
  }
}

auto Projector::computeNaturalOrbitals(Utils::ElementType type) -> void {
  if (molecule_->at(type).isRestricted) {
    es_.compute(data_->D_OAO[type].restrictedMatrix());
    data_->O_NO.at(type).restrictedMatrix() = es_.eigenvectors().rowwise().reverse();
  }
  else {
    es_.compute(data_->D_OAO[type].alphaMatrix());
    data_->O_NO.at(type).alphaMatrix() = es_.eigenvectors().rowwise().reverse();
    if (molecule_->at(type).msVector[1] > 0) {
      es_.compute(data_->D_OAO[type].betaMatrix());
      data_->O_NO.at(type).betaMatrix() = es_.eigenvectors().rowwise().reverse();
    }
  }
}

auto Projector::updateDensityFromRotation(Utils::ElementType type, const Utils::SpinAdaptedMatrix& orbitalRotation) -> void {
  if (molecule_->at(type).isRestricted) {
    auto const& X = data_->X[type];

    Eigen::MatrixXd tmp_oao = orbitalRotation.restrictedMatrix() * data_->C_OAO[type].restrictedMatrix();
    Eigen::MatrixXd tmp = X * tmp_oao;

    Utils::MolecularOrbitals C;
    Utils::MolecularOrbitals C_OAO;
    C = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp));
    C_OAO = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp_oao));
    data_->C.at(type) = C;
    data_->C_OAO.at(type) = C_OAO;
    data_->D.at(type) = Kiwi::HartreeFockUtils::makeDensity(type, molecule_, C);
    data_->D_OAO.at(type) = Kiwi::HartreeFockUtils::makeDensity(type, molecule_, C_OAO);
  }
  else {
    const auto& X = data_->X.at(type);
    Eigen::MatrixXd tmp_oao_alpha = orbitalRotation.alphaMatrix() * data_->C_OAO[type].alphaMatrix();
    Eigen::MatrixXd tmp_alpha = X * tmp_oao_alpha;

    Eigen::MatrixXd tmp_oao_beta;
    tmp_oao_beta.resizeLike(tmp_oao_alpha);
    Eigen::MatrixXd tmp_beta;
    tmp_beta.resizeLike(tmp_alpha);
    if (molecule_->at(type).msVector[1] > 0) {
      tmp_oao_beta = orbitalRotation.betaMatrix() * data_->C_OAO[type].betaMatrix();
      tmp_beta = X * tmp_oao_beta;
    }
    Utils::MolecularOrbitals C;
    Utils::MolecularOrbitals C_OAO;
    C = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp_alpha),
                                                                                      std::move(tmp_beta));
    C_OAO = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp_oao_alpha),
                                                                                          std::move(tmp_oao_beta));
    data_->C.at(type) = C;
    data_->C_OAO.at(type) = C_OAO;
    data_->D.at(type) = Kiwi::HartreeFockUtils::makeDensity(type, molecule_, C);
    data_->D_OAO.at(type) = Kiwi::HartreeFockUtils::makeDensity(type, molecule_, C_OAO);
  }
}

std::chrono::duration<double, std::milli> Projector::getFockMatrixFormationTime() {
  return fockMatrixFormationTime_;
}

auto Projector::evaluateDensity(bool incremental) -> void {
  if (data_->integralDirect && incremental) {
    for (auto& elem : *molecule_) {
      incrementalDensityMatrix_[elem.first] = data_->D.at(elem.first);
      incrementalDensityMatrix_[elem.first].addDensity(lastDensityMatrix_.at(elem.first), -1);
    }
  }
  else {
    for (auto& elem : *molecule_) {
      incrementalDensityMatrix_[elem.first] = data_->D.at(elem.first);
    }
  }
}

auto Projector::setLastDensity(const std::map<Utils::ElementType, Utils::DensityMatrix>& lastDensity) -> void {
  lastDensityMatrix_ = lastDensity;
}

auto Projector::resetDensity(const std::map<Utils::ElementType, Utils::DensityMatrix>& D,
                             const std::map<Utils::ElementType, Utils::DensityMatrix>& D_OAO,
                             const std::map<Utils::ElementType, Utils::MolecularOrbitals>& C,
                             const std::map<Utils::ElementType, Utils::MolecularOrbitals>& C_OAO) -> void {
  data_->D = D;
  data_->D_OAO = D_OAO;
  data_->C = C;
  data_->C_OAO = C_OAO;
}

} // namespace Kiwi
} // namespace Scine