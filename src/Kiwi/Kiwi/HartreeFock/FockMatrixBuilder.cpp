/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/FockMatrixBuilder.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <iostream>
#include <utility>

#define UNUSED(expr) \
  do {               \
    (void)(expr);    \
  } while (0)

namespace Scine {
namespace Kiwi {

FockMatrixBuilder::FockMatrixBuilder(Utils::ElementType type, std::shared_ptr<Data> data,
                                     std::shared_ptr<Molecule> molecule, std::shared_ptr<HartreeFockData> hartreeFockData)
  : type_(type), data_(std::move(data)), molecule_(std::move(molecule)), hartreeFockData_(std::move(hartreeFockData)) {
  elapsedTime_ = std::chrono::duration<double, std::milli>::zero();
}

auto FockMatrixBuilder::init() -> void {
  if (hartreeFockData_->F.find(type_) == hartreeFockData_->F.end())
    hartreeFockData_->F.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->F_OAO.find(type_) == hartreeFockData_->F_OAO.end())
    hartreeFockData_->F_OAO.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->G.find(type_) == hartreeFockData_->G.end())
    hartreeFockData_->G.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->I.find(type_) == hartreeFockData_->I.end())
    hartreeFockData_->I.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->J.find(type_) == hartreeFockData_->J.end())
    hartreeFockData_->J.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->K.find(type_) == hartreeFockData_->K.end())
    hartreeFockData_->K.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->L.find(type_) == hartreeFockData_->L.end()) {
    std::map<Utils::ElementType, Eigen::MatrixXd> tmp;
    for (auto const& elem : *molecule_) {
      if (elem.first != type_) {
        tmp[elem.first] = Eigen::MatrixXd(molecule_->at(type_).LAO, molecule_->at(type_).LAO);
      }
    }
    hartreeFockData_->L.insert({type_, tmp});
  }
  hartreeFockData_->F.at(type_).resize(molecule_->at(type_).LAO);
  hartreeFockData_->F_OAO.at(type_).resize(molecule_->at(type_).LMO);
  hartreeFockData_->G.at(type_).resize(molecule_->at(type_).LAO);
  hartreeFockData_->I.at(type_).resize(molecule_->at(type_).LAO);
  hartreeFockData_->J.at(type_).resize(molecule_->at(type_).LAO);
  hartreeFockData_->K.at(type_).resize(molecule_->at(type_).LAO);
}

auto FockMatrixBuilder::initGeneralized() -> void {
  if (hartreeFockData_->F.find(type_) == hartreeFockData_->F.end())
    hartreeFockData_->F.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->F_OAO.find(type_) == hartreeFockData_->F_OAO.end())
    hartreeFockData_->F_OAO.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->G.find(type_) == hartreeFockData_->G.end())
    hartreeFockData_->G.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->I.find(type_) == hartreeFockData_->I.end())
    hartreeFockData_->I.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->J.find(type_) == hartreeFockData_->J.end())
    hartreeFockData_->J.insert({type_, Utils::SpinAdaptedMatrix()});
  if (hartreeFockData_->K.find(type_) == hartreeFockData_->K.end())
    hartreeFockData_->K.insert({type_, Utils::SpinAdaptedMatrix()});
  hartreeFockData_->F.at(type_).resize(2 * molecule_->at(type_).LAO);
  hartreeFockData_->F_OAO.at(type_).resize(molecule_->at(type_).LMO);
  hartreeFockData_->G.at(type_).resize(2 * molecule_->at(type_).LAO);
  hartreeFockData_->I.at(type_).resize(2 * molecule_->at(type_).LAO);
  hartreeFockData_->J.at(type_).resize(2 * molecule_->at(type_).LAO);
  hartreeFockData_->K.at(type_).resize(2 * molecule_->at(type_).LAO);
}

template<>
auto FockMatrixBuilder::buildJKmatrices<SpinSymmetry::Generalized>(std::map<Utils::ElementType, Utils::DensityMatrix>& D)
    -> void {
  if (molecule_->at(type_).N < 2) {
    return;
  }
  hartreeFockData_->J.at(type_).restrictedMatrix().setZero();
  hartreeFockData_->K.at(type_).restrictedMatrix().setZero();

  int L = molecule_->at(type_).LAO;

  Eigen::MatrixXd Dalpha = D.at(type_).restrictedMatrix().block(0, 0, L, L).transpose();
  Eigen::MatrixXd Dbeta = D.at(type_).restrictedMatrix().block(L, L, L, L).transpose();
  Eigen::MatrixXd Dtot = Dalpha + Dbeta;

  // Coulomb
  {
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dtot.data(), Dtot.rows() * Dtot.cols()));
    Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, type_)] * coefficientVector;
    hartreeFockData_->J.at(type_).restrictedMatrix().block(0, 0, L, L) =
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
    hartreeFockData_->J.at(type_).restrictedMatrix().block(L, L, L, L) =
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
  }
  // Alpha exchange
  if (molecule_->at(type_).msVector[0] > 0) {
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dalpha.data(), Dalpha.rows() * Dalpha.cols()));
    Eigen::VectorXd tmpVec = data_->Exchange[type_] * coefficientVector;
    hartreeFockData_->K.at(type_).restrictedMatrix().block(0, 0, L, L) =
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
  }
  // Beta exchange
  if (molecule_->at(type_).msVector[1] > 0) {
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dbeta.data(), Dbeta.rows() * Dbeta.cols()));
    Eigen::VectorXd tmpVec = data_->Exchange[type_] * coefficientVector;
    hartreeFockData_->K.at(type_).restrictedMatrix().block(L, L, L, L) =
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
  }
  //// Generalized alpha-beta and beta-alpha exchange
  if (molecule_->at(type_).msVector[0] > 0 && molecule_->at(type_).msVector[1] > 0) {
    // alpha-beta
    {
      Eigen::MatrixXd D_alpha_beta = D.at(type_).restrictedMatrix().block(0, L, L, L).transpose();
      Eigen::VectorXd coefficientVector(
          Eigen::Map<Eigen::VectorXd>(D_alpha_beta.data(), D_alpha_beta.rows() * D_alpha_beta.cols()));
      Eigen::VectorXd tmpVec = data_->Exchange[type_] * coefficientVector;
      hartreeFockData_->K.at(type_).restrictedMatrix().block(0, L, L, L) =
          Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
    }
    // beta-alpha
    {
      Eigen::MatrixXd D_beta_alpha = D.at(type_).restrictedMatrix().block(L, 0, L, L).transpose();
      Eigen::VectorXd coefficientVector(
          Eigen::Map<Eigen::VectorXd>(D_beta_alpha.data(), D_beta_alpha.rows() * D_beta_alpha.cols()));
      Eigen::VectorXd tmpVec = data_->Exchange[type_] * coefficientVector;
      hartreeFockData_->K.at(type_).restrictedMatrix().block(L, 0, L, L) =
          Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
    }
  }
}

template<>
auto FockMatrixBuilder::buildJKmatrices<SpinSymmetry::Restricted>(std::map<Utils::ElementType, Utils::DensityMatrix>& D)
    -> void {
  Eigen::MatrixXd Dt = D.at(type_).restrictedMatrix().transpose();

  Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
  {
    Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, type_)] * coefficientVector;
    hartreeFockData_->J.at(type_).setRestrictedMatrix(
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO));
  }
  {
    Eigen::VectorXd tmpVec = data_->Exchange[type_] * coefficientVector;
    hartreeFockData_->K.at(type_).setRestrictedMatrix(
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO));
  }
}

template<>
auto FockMatrixBuilder::buildJKmatrices<SpinSymmetry::Unrestricted>(std::map<Utils::ElementType, Utils::DensityMatrix>& D)
    -> void {
  if (molecule_->at(type_).N < 2) {
    return;
  }

  Eigen::MatrixXd Dalpha = D.at(type_).alphaMatrix().transpose();
  Eigen::MatrixXd Dbeta = D.at(type_).betaMatrix().transpose();

  // Coulomb alpha
  {
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dalpha.data(), Dalpha.rows() * Dalpha.cols()));
    Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, type_)] * coefficientVector;
    hartreeFockData_->J.at(type_).setAlphaMatrix(
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO));
  }
  // Coulomb beta
  if (molecule_->at(type_).msVector[1] > 0) {
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dbeta.data(), Dbeta.rows() * Dbeta.cols()));
    Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, type_)] * coefficientVector;
    hartreeFockData_->J.at(type_).setBetaMatrix(
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO));
  }
  // Alpha exchange
  if (molecule_->at(type_).msVector[0] > 0) {
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dalpha.data(), Dalpha.rows() * Dalpha.cols()));
    Eigen::VectorXd tmpVec = data_->Exchange[type_] * coefficientVector;
    hartreeFockData_->K.at(type_).setAlphaMatrix(
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO));
  }
  // Beta exchange
  if (molecule_->at(type_).msVector[1] > 0) {
    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dbeta.data(), Dbeta.rows() * Dbeta.cols()));
    Eigen::VectorXd tmpVec = data_->Exchange[type_] * coefficientVector;
    hartreeFockData_->K.at(type_).setBetaMatrix(
        Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO));
  }
}

template<>
auto FockMatrixBuilder::buildFockmatrices<SpinSymmetry::Restricted>(std::map<Utils::ElementType, Utils::DensityMatrix>& D)
    -> void {
  UNUSED(D);

  hartreeFockData_->I.at(type_).restrictedMatrix().setZero();
  auto& Itype_ = hartreeFockData_->I.at(type_).restrictedMatrix();
  for (auto const& elem : *molecule_) {
    if (elem.first == type_) {
      continue;
    }
    auto otherType = elem.first;
    // Unrestricted
    hartreeFockData_->I.at(type_).restrictedMatrix() += hartreeFockData_->L.at(type_).at(otherType);
  } // molecule_
  hartreeFockData_->G.at(type_).restrictedMatrix() =
      hartreeFockData_->J.at(type_).restrictedMatrix() - 0.5 * hartreeFockData_->K.at(type_).restrictedMatrix();
  hartreeFockData_->F.at(type_).restrictedMatrix() = data_->H[type_] + hartreeFockData_->G.at(type_).restrictedMatrix() + Itype_;
  hartreeFockData_->F_OAO.at(type_).restrictedMatrix() =
      data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).restrictedMatrix() * data_->X.at(type_);
}

template<>
auto FockMatrixBuilder::buildFockmatrices<SpinSymmetry::Unrestricted>(std::map<Utils::ElementType, Utils::DensityMatrix>& D)
    -> void {
  UNUSED(D);

  hartreeFockData_->I.at(type_).restrictedMatrix().setZero();
  if (molecule_->at(type_).msVector[0] > 0) {
    hartreeFockData_->G.at(type_).alphaMatrix().setZero();
    if (molecule_->at(type_).N > 1) {
      hartreeFockData_->G.at(type_).alphaMatrix() +=
          hartreeFockData_->J.at(type_).alphaMatrix() - hartreeFockData_->K.at(type_).alphaMatrix();
    }
  }
  if (molecule_->at(type_).msVector[1] > 0) {
    hartreeFockData_->G.at(type_).alphaMatrix() += hartreeFockData_->J.at(type_).betaMatrix();
    hartreeFockData_->G.at(type_).betaMatrix().setZero();
    if (molecule_->at(type_).N > 1) {
      hartreeFockData_->G.at(type_).betaMatrix() +=
          hartreeFockData_->J.at(type_).betaMatrix() - hartreeFockData_->K.at(type_).betaMatrix();
      hartreeFockData_->G.at(type_).betaMatrix() += hartreeFockData_->J.at(type_).alphaMatrix();
    }
  }
  for (auto const& elem : *molecule_) {
    if (elem.first == type_) {
      continue;
    }
    auto otherType = elem.first;
    hartreeFockData_->I.at(type_).restrictedMatrix() += hartreeFockData_->L.at(type_).at(otherType);
  } // molecule_

  if (molecule_->at(type_).msVector[0] > 0) {
    hartreeFockData_->F.at(type_).alphaMatrix() =
        data_->H[type_] + hartreeFockData_->G.at(type_).alphaMatrix() + hartreeFockData_->I.at(type_).restrictedMatrix();
    hartreeFockData_->F_OAO.at(type_).alphaMatrix() =
        data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).alphaMatrix() * data_->X.at(type_);
  }
  if (molecule_->at(type_).msVector[1] > 0) {
    hartreeFockData_->F.at(type_).betaMatrix() =
        data_->H[type_] + hartreeFockData_->G.at(type_).betaMatrix() + hartreeFockData_->I.at(type_).restrictedMatrix();
    hartreeFockData_->F_OAO.at(type_).betaMatrix() =
        data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).betaMatrix() * data_->X.at(type_);
  }
}

// TODO adapt generalized builder to the new procedure with the "L" matrices.
template<>
auto FockMatrixBuilder::buildFockmatrices<SpinSymmetry::Generalized>(std::map<Utils::ElementType, Utils::DensityMatrix>& D)
    -> void {
  int L = molecule_->at(type_).LMO / 2;
  hartreeFockData_->I.at(type_).restrictedMatrix().setZero();
  hartreeFockData_->G.at(type_).restrictedMatrix().setZero();
  hartreeFockData_->F.at(type_).restrictedMatrix().setZero();
  hartreeFockData_->F_OAO.at(type_).restrictedMatrix().setZero();

  hartreeFockData_->G.at(type_).restrictedMatrix() +=
      hartreeFockData_->J.at(type_).restrictedMatrix() - hartreeFockData_->K.at(type_).restrictedMatrix();

  for (auto const& elem : *molecule_) {
    if (elem.first == type_) {
      continue;
    }
    int L_loop = molecule_->at(elem.first).LMO / 2;
    // Only the unrestricted-unrestricted case is possible.
    Eigen::MatrixXd Dalpha = D.at(elem.first).restrictedMatrix().block(0, 0, L_loop, L_loop).transpose();
    Eigen::MatrixXd Dbeta = D.at(elem.first).restrictedMatrix().block(L_loop, L_loop, L_loop, L_loop).transpose();
    Eigen::MatrixXd Dt = Dalpha + Dbeta;
    //    if (elem.second.isRestricted) {
    if (type_ < elem.first) {
      Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
      Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, elem.first)] * coefficientVector;
      hartreeFockData_->I.at(type_).restrictedMatrix().block(0, 0, L, L) +=
          Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
      hartreeFockData_->I.at(type_).restrictedMatrix().block(L, L, L, L) +=
          hartreeFockData_->I.at(type_).restrictedMatrix().block(0, 0, L, L);
    }
    else {
      Dt.transposeInPlace(); // Projected density matrix is not symmetric!
      Eigen::RowVectorXd coefficientRowVector(Eigen::Map<Eigen::RowVectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
      Eigen::RowVectorXd tmpVec = coefficientRowVector * data_->Coulomb[Data::getUnique(type_, elem.first)];
      Dt.transposeInPlace(); // Projected density matrix is not symmetric!
      hartreeFockData_->I.at(type_).restrictedMatrix().block(0, 0, L, L) +=
          Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
      hartreeFockData_->I.at(type_).restrictedMatrix().block(L, L, L, L) +=
          hartreeFockData_->I.at(type_).restrictedMatrix().block(0, 0, L, L);
    }
    //    } // Unrestricted-unrestriced
  } // molecule_
  if (molecule_->at(type_).msVector[0] > 0) {
    hartreeFockData_->F.at(type_).restrictedMatrix().block(0, 0, L, L) +=
        data_->H[type_] + hartreeFockData_->G.at(type_).restrictedMatrix().block(0, 0, L, L) +
        hartreeFockData_->I.at(type_).restrictedMatrix().block(0, 0, L, L);
    hartreeFockData_->F_OAO.at(type_).restrictedMatrix().block(0, 0, L, L) =
        data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).restrictedMatrix().block(0, 0, L, L) *
        data_->X.at(type_);
  }
  if (molecule_->at(type_).msVector[1] > 0) {
    hartreeFockData_->F.at(type_).restrictedMatrix().block(L, L, L, L) +=
        data_->H[type_] + hartreeFockData_->G.at(type_).restrictedMatrix().block(L, L, L, L) +
        hartreeFockData_->I.at(type_).restrictedMatrix().block(L, L, L, L);
    hartreeFockData_->F_OAO.at(type_).restrictedMatrix().block(L, L, L, L) =
        data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).restrictedMatrix().block(L, L, L, L) *
        data_->X.at(type_);
  }
  // Generalized alpha-beta and beta-alpha exchange
  if (molecule_->at(type_).msVector[0] > 0 && molecule_->at(type_).msVector[1] > 0) {
    // alpha-beta
    hartreeFockData_->F.at(type_).restrictedMatrix().block(0, L, L, L) +=
        hartreeFockData_->G.at(type_).restrictedMatrix().block(0, L, L, L);
    hartreeFockData_->F_OAO.at(type_).restrictedMatrix().block(0, L, L, L) =
        data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).restrictedMatrix().block(0, L, L, L) *
        data_->X.at(type_);
    // beta-alpha
    hartreeFockData_->F.at(type_).restrictedMatrix().block(L, 0, L, L) +=
        hartreeFockData_->G.at(type_).restrictedMatrix().block(L, 0, L, L);
    hartreeFockData_->F_OAO.at(type_).restrictedMatrix().block(L, 0, L, L) =
        data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).restrictedMatrix().block(L, 0, L, L) *
        data_->X.at(type_);
  }
}

// auto FockMatrixBuilder::update(std::map<Utils::ElementType, Utils::DensityMatrix>& D, bool incremental,
//                               std::size_t numRebuilds) -> void {
//  auto start = std::chrono::high_resolution_clock::now();
//  ++fockMatrixFunctionCalls_;
//
//  formIncremental_ = incremental;
//  numRebuilds_ = numRebuilds;
//
//  if (!data_->integralDirect) {
//    buildLmatrices(D);
//    if (molecule_->at(type_).isRestricted) {
//      buildJKmatrices<SpinSymmetry::Restricted>(D);
//    }
//    else {
//      buildJKmatrices<SpinSymmetry::Unrestricted>(D);
//    }
//  }
//  else {
//    buildJKmatricesDirect(D);
//    // direct evaluation is of no use, yet.
//    buildLmatricesDirect(data_->D);
//  }
//
//  if (molecule_->at(type_).isRestricted) {
//    buildFockmatrices<SpinSymmetry::Restricted>(D);
//  }
//  else {
//    buildFockmatrices<SpinSymmetry::Unrestricted>(D);
//  }
//
//  auto stop = std::chrono::high_resolution_clock::now();
//
//  elapsedTime_ += stop - start;
//  std::chrono::duration<double, std::milli> duration = stop - start;
//}

auto FockMatrixBuilder::updateGeneralized(std::map<Utils::ElementType, Utils::DensityMatrix>& D, bool verbose) -> void {
  if (verbose) {
    std::cout << "Updating Fock matrices... " << std::endl;
    Clock& clock = Clock::getInstance();
    clock.time("Fock");
  }
  buildJKmatrices<SpinSymmetry::Generalized>(D);
  buildFockmatrices<SpinSymmetry::Generalized>(D);
  if (verbose) {
    Clock& clock = Clock::getInstance();
    clock.time("Fock");
  }
}

const std::shared_ptr<HartreeFockData>& FockMatrixBuilder::getHartreeFockData() const {
  return hartreeFockData_;
}

auto FockMatrixBuilder::buildLmatrices(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void {
  //
  // Note: The contraction of the two-body tensor for different particle type_s, depends on the ordering
  //       of the particles: V_ijIJ != V_IJij.
  //       Hence: we must keep track of the ordering.
  //

  for (auto const& elem : *molecule_) {
    if (elem.first == type_) {
      continue;
    }
    auto otherType = elem.first;
    // Unrestricted
    if (!elem.second.isRestricted) {
      hartreeFockData_->L.at(type_).at(otherType).setZero();
      if (type_ < otherType) {
        // TODO use total density matrix instead of alpha and beta separately
        // alpha
        {
          // TODO check the transpose
          Eigen::MatrixXd Dt_alpha = D[otherType].alphaMatrix().transpose();
          Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dt_alpha.data(), Dt_alpha.rows() * Dt_alpha.cols()));
          Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, otherType)] * coefficientVector;
          hartreeFockData_->L.at(type_).at(otherType) +=
              Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
        }
        // beta
        if (elem.second.msVector[1] > 0) {
          // TODO check the transpose
          Eigen::MatrixXd Dt_beta = D[otherType].betaMatrix().transpose();
          Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dt_beta.data(), Dt_beta.rows() * Dt_beta.cols()));
          Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, otherType)] * coefficientVector;
          hartreeFockData_->L.at(type_).at(otherType) +=
              Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
        }
      }
      else {
        // alpha
        {
          // TODO check the transpose
          Eigen::MatrixXd Dt_alpha = D[otherType].alphaMatrix().transpose();
          Dt_alpha.transposeInPlace();
          Eigen::RowVectorXd coefficientRowVector(
              Eigen::Map<Eigen::RowVectorXd>(Dt_alpha.data(), Dt_alpha.rows() * Dt_alpha.cols()));
          Eigen::RowVectorXd tmpVec = coefficientRowVector * data_->Coulomb[Data::getUnique(type_, otherType)];
          Dt_alpha.transposeInPlace();
          hartreeFockData_->L.at(type_).at(otherType) +=
              Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
        }
        // beta
        if (elem.second.msVector[1] > 0) {
          // TODO check the transpose
          Eigen::MatrixXd Dt_beta = D[otherType].betaMatrix().transpose();
          Dt_beta.transposeInPlace();
          Eigen::RowVectorXd coefficientRowVector(
              Eigen::Map<Eigen::RowVectorXd>(Dt_beta.data(), Dt_beta.rows() * Dt_beta.cols()));
          Eigen::RowVectorXd tmpVec = coefficientRowVector * data_->Coulomb[Data::getUnique(type_, otherType)];
          Dt_beta.transposeInPlace();
          hartreeFockData_->L.at(type_).at(otherType) +=
              Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
        }
      }
    } // Unrestricted
    // Restricted
    else if (elem.second.isRestricted) {
      hartreeFockData_->L.at(type_).at(otherType).setZero();
      if (type_ < otherType) {
        // TODO check the transpose
        Eigen::MatrixXd Dt = D[otherType].restrictedMatrix().transpose();
        Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
        Eigen::VectorXd tmpVec = data_->Coulomb[Data::getUnique(type_, otherType)] * coefficientVector;
        hartreeFockData_->L.at(type_).at(otherType) +=
            Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
      }
      else {
        // TODO check the transpose
        Eigen::MatrixXd Dt = D[otherType].restrictedMatrix().transpose();
        Dt.transposeInPlace();
        Eigen::RowVectorXd coefficientRowVector(Eigen::Map<Eigen::RowVectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
        Eigen::RowVectorXd tmpVec = coefficientRowVector * data_->Coulomb[Data::getUnique(type_, otherType)];
        Dt.transposeInPlace();
        hartreeFockData_->L.at(type_).at(otherType) +=
            Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule_->at(type_).LAO, molecule_->at(type_).LAO);
      }
    } // Restriced
  }   // molecule_
}

// auto FockMatrixBuilder::buildJKmatricesDirect(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void {
//  const auto& densityMatrix = D.at(type_);
//  auto& basis = molecule_->at(type_).basisSet;
//
//  Utils::Integrals::IntegralSpecifier specifier;
//  specifier.op = Utils::Integrals::Operator::Coulomb;
//
//  double prescreeningThreshold = 1e-8 / (3.0 * double(basis.nbf()));
//
//  if (formIncremental_) {
//    for (auto i = 0UL; i < numRebuilds_; ++i) {
//      prescreeningThreshold /= 100;
//    }
//  }
//
//  auto prescreener = Integrals::TwoBody::CauchySchwarzDensityPrescreener(basis, densityMatrix, prescreeningThreshold);
//
//  auto saver = Integrals::TwoBody::CoulombExchangeDigester(basis, basis, specifier, densityMatrix);
//  auto evaluator =
//      Integrals::TwoBody::Evaluator<Integrals::TwoBody::CoulombExchangeDigester,
//      Integrals::TwoBody::CauchySchwarzDensityPrescreener>(
//          basis, basis, specifier, std::move(saver), std::move(prescreener));
//  evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
//
//  if (formIncremental_) {
//    if (molecule_->at(type_).isRestricted) {
//      hartreeFockData_->J.at(type_).restrictedMatrix() += evaluator.getResult().first.restrictedMatrix();
//      hartreeFockData_->K.at(type_).restrictedMatrix() += evaluator.getResult().second.restrictedMatrix();
//    }
//    else {
//      hartreeFockData_->J.at(type_).alphaMatrix() += evaluator.getResult().first.alphaMatrix();
//      hartreeFockData_->K.at(type_).alphaMatrix() += evaluator.getResult().second.alphaMatrix();
//      if (molecule_->at(type_).msVector[1] > 0) {
//        hartreeFockData_->J.at(type_).betaMatrix() += evaluator.getResult().first.betaMatrix();
//        hartreeFockData_->K.at(type_).betaMatrix() += evaluator.getResult().second.betaMatrix();
//      }
//    }
//  }
//  else {
//    hartreeFockData_->J.at(type_) = evaluator.getResult().first;
//    hartreeFockData_->K.at(type_) = evaluator.getResult().second;
//  }
//}
//
// auto FockMatrixBuilder::buildLmatricesDirect(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void {
//  Utils::Integrals::IntegralSpecifier specifier;
//
//  auto tp1 = molecule_->at(type_);
//  const auto& densityMatrix1 = D.at(type_);
//  auto& basis1 = molecule_->at(type_).basisSet;
//
//  for (auto const& elem : *molecule_) {
//    if (elem.first == type_) {
//      continue;
//    }
//    auto otherType = elem.first;
//    auto tp2 = molecule_->at(otherType);
//    const auto& densityMatrix2 = D.at(otherType);
//    auto& basis2 = molecule_->at(otherType).basisSet;
//    specifier.typeVector = {tp1.typeInfo, tp2.typeInfo};
//    specifier.op = Utils::Integrals::Operator::Coulomb;
//
//    auto prescreener = Integrals::TwoBody::VoidPrescreener();
//
//    auto saver = Integrals::TwoBody::TwoTypeCoulombDigester(basis1, basis2, specifier, densityMatrix1,
//    densityMatrix2); auto evaluator = Integrals::TwoBody::Evaluator<Integrals::TwoBody::TwoTypeCoulombDigester>(
//        basis1, basis2, specifier, std::move(saver), std::move(prescreener));
//    evaluator.evaluateTwoBodyIntegrals<libint2::Operator::coulomb>();
//
//    hartreeFockData_->L.at(type_).at(otherType) = evaluator.getResult().first;
//  }
//}

} // namespace Kiwi
} // namespace Scine