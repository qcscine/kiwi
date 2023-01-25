/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_BOHELPER_H
#define KIWI_BOHELPER_H

#include <Kiwi/HartreeFock/FockMatrix/CoulombExchangeMatrices.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <chrono>

namespace Scine {
namespace Kiwi {
namespace FockMatrix {

/**
 * @class BoHelper @file BoHelper.h
 * @brief This class handles the calulcation of the BO terms in the Fock matrix.
 * It can handle the in-memory calculation and the integral-direct evaluation.
 * Also supported is the incremental Fock building.
 */
class BoHelper {
 public:
  BoHelper() = default;

  BoHelper(Utils::ElementType type, std::shared_ptr<Data> data, std::shared_ptr<Molecule> molecule,
           std::shared_ptr<HartreeFockData> hartreeFockData)
    : type_(type), data_(std::move(data)), molecule_(std::move(molecule)), hartreeFockData_(std::move(hartreeFockData)) {
    elapsedTime_ = std::chrono::duration<double, std::milli>::zero();

    init();
  }

  auto update(const std::map<Utils::ElementType, Utils::DensityMatrix>& D, bool incremental, std::size_t numRebuilds) -> void {
    updateJK(D, data_->integralDirect, incremental, numRebuilds);
    updateFockMatrix();
  }

  auto updateJK(const std::map<Utils::ElementType, Utils::DensityMatrix>& D, bool integralDirect,
                bool incremental = false, std::size_t numRebuilds = 0) -> void {
    auto start = std::chrono::high_resolution_clock::now();
    ++fockMatrixFunctionCalls_;

    if (integralDirect) {
      auto JK = buildJKmatricesDirect(D.at(type_), molecule_, type_, incremental, numRebuilds);
      if (incremental) {
        if (molecule_->at(type_).isRestricted) {
          hartreeFockData_->J.at(type_).restrictedMatrix() += JK[0].restrictedMatrix();
          hartreeFockData_->K.at(type_).restrictedMatrix() += JK[1].restrictedMatrix();
        }
        else {
          hartreeFockData_->J.at(type_).alphaMatrix() += JK[0].alphaMatrix();
          hartreeFockData_->K.at(type_).alphaMatrix() += JK[1].alphaMatrix();
          if (molecule_->at(type_).msVector[1] > 0) {
            hartreeFockData_->J.at(type_).betaMatrix() += JK[0].betaMatrix();
            hartreeFockData_->K.at(type_).betaMatrix() += JK[1].betaMatrix();
          }
        }
      }
      else {
        hartreeFockData_->J.at(type_) = std::move(JK[0]);
        hartreeFockData_->K.at(type_) = std::move(JK[1]);
      }
    }
    else {
      auto JK = buildJKmatrices(D.at(type_), molecule_, data_, type_);
      hartreeFockData_->J.at(type_) = std::move(JK[0]);
      hartreeFockData_->K.at(type_) = std::move(JK[1]);
    }

    auto stop = std::chrono::high_resolution_clock::now();

    elapsedTime_ += stop - start;
  }

  auto updateFockMatrix() -> void {
    auto start = std::chrono::high_resolution_clock::now();

    if (molecule_->at(type_).isRestricted) {
      hartreeFockData_->I.at(type_).restrictedMatrix().setZero();

      auto& Itype_ = hartreeFockData_->I.at(type_).restrictedMatrix();
      for (auto const& elem : *molecule_) {
        if (elem.first == type_) {
          continue;
        }
        auto otherType = elem.first;
        // Unrestricted
        hartreeFockData_->I.at(type_).restrictedMatrix() += hartreeFockData_->L.at(type_).at(otherType);
      } // molecule

      hartreeFockData_->G.at(type_).restrictedMatrix() =
          hartreeFockData_->J.at(type_).restrictedMatrix() - 0.5 * hartreeFockData_->K.at(type_).restrictedMatrix();
      hartreeFockData_->F.at(type_).restrictedMatrix() =
          data_->H[type_] + hartreeFockData_->G.at(type_).restrictedMatrix() + Itype_;
      hartreeFockData_->F_OAO.at(type_).restrictedMatrix() =
          data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).restrictedMatrix() * data_->X.at(type_);
    }
    else {
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
        hartreeFockData_->F.at(type_).alphaMatrix() = data_->H[type_] + hartreeFockData_->G.at(type_).alphaMatrix() +
                                                      hartreeFockData_->I.at(type_).restrictedMatrix();
        hartreeFockData_->F_OAO.at(type_).alphaMatrix() =
            data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).alphaMatrix() * data_->X.at(type_);
      }
      if (molecule_->at(type_).msVector[1] > 0) {
        hartreeFockData_->F.at(type_).betaMatrix() = data_->H[type_] + hartreeFockData_->G.at(type_).betaMatrix() +
                                                     hartreeFockData_->I.at(type_).restrictedMatrix();
        hartreeFockData_->F_OAO.at(type_).betaMatrix() =
            data_->X.at(type_).transpose() * hartreeFockData_->F.at(type_).betaMatrix() * data_->X.at(type_);
      }
    }

    auto stop = std::chrono::high_resolution_clock::now();

    elapsedTime_ += stop - start;
  }

  auto getElapsedTime() -> std::chrono::duration<double, std::milli> {
    return elapsedTime_ / fockMatrixFunctionCalls_;
  }

 private:
  Utils::ElementType type_;
  std::shared_ptr<Data> data_;
  std::shared_ptr<Molecule> molecule_;
  std::shared_ptr<HartreeFockData> hartreeFockData_;

  std::chrono::duration<double, std::milli> elapsedTime_;

  int fockMatrixFunctionCalls_ = 0;

  auto init() -> void {
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
    hartreeFockData_->F.at(type_).resize(molecule_->at(type_).LAO);
    hartreeFockData_->F_OAO.at(type_).resize(molecule_->at(type_).LMO);
    hartreeFockData_->G.at(type_).resize(molecule_->at(type_).LAO);
    hartreeFockData_->I.at(type_).resize(molecule_->at(type_).LAO);
    hartreeFockData_->J.at(type_).resize(molecule_->at(type_).LAO);
    hartreeFockData_->K.at(type_).resize(molecule_->at(type_).LAO);
  }
};

} // namespace FockMatrix
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_BOHELPER_H
