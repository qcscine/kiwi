/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_PREBOHELPER_H
#define KIWI_PREBOHELPER_H

#include <Kiwi/HartreeFock/FockMatrix/TwoTypeCoulomb.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <chrono>

namespace Scine {
namespace Kiwi {
namespace FockMatrix {

/**
 * @class PreBoHelper @file PreBoHelper.h
 * @brief This class handles the calulcation of the pre-BO terms in the Fock matrix.
 * It can handle the in-memory calculation and the integral-direct evaluation.
 */
class PreBoHelper {
 public:
  PreBoHelper() = default;

  PreBoHelper(std::shared_ptr<Data> data, std::shared_ptr<Molecule> molecule, std::shared_ptr<HartreeFockData> hartreeFockData)
    : data_(std::move(data)), molecule_(std::move(molecule)), hartreeFockData_(std::move(hartreeFockData)) {
    elapsedTime_ = std::chrono::duration<double, std::milli>::zero();

    init();
  }

  auto updateL(std::map<Utils::ElementType, Utils::DensityMatrix>& D, bool integralDirect) -> void {
    auto start = std::chrono::high_resolution_clock::now();
    ++updateFunctionCalls_;

    for (const auto& uniquePair : data_->uniquePairs) {
      if (uniquePair.first == uniquePair.second) {
        continue;
      }
      std::array<Eigen::MatrixXd, 2> L12;
      if (integralDirect) {
        L12 = buildLmatricesDirect(D, molecule_, uniquePair.first, uniquePair.second);
      }
      else {
        L12 = buildLmatrices(D, molecule_, data_, uniquePair.first, uniquePair.second);
      }
      hartreeFockData_->L.at(uniquePair.first).at(uniquePair.second) = std::move(L12[0]);
      hartreeFockData_->L.at(uniquePair.second).at(uniquePair.first) = std::move(L12[1]);
    }

    auto stop = std::chrono::high_resolution_clock::now();

    elapsedTime_ += stop - start;
  }

  auto getElapsedTime() -> std::chrono::duration<double, std::milli> {
    return elapsedTime_ / updateFunctionCalls_;
  }

 private:
  std::shared_ptr<Data> data_;
  std::shared_ptr<Molecule> molecule_;
  std::shared_ptr<HartreeFockData> hartreeFockData_;

  std::chrono::duration<double, std::milli> elapsedTime_;

  int updateFunctionCalls_ = 0;

  auto init() -> void {
    for (auto const& elem1 : *molecule_) {
      if (hartreeFockData_->L.find(elem1.first) == hartreeFockData_->L.end()) {
        std::map<Utils::ElementType, Eigen::MatrixXd> tmp;
        for (auto const& elem2 : *molecule_) {
          if (elem2.first != elem1.first) {
            tmp[elem2.first] = Eigen::MatrixXd(molecule_->at(elem1.first).LAO, molecule_->at(elem1.first).LAO);
          }
        }
        hartreeFockData_->L.insert({elem1.first, tmp});
      }
    }
  }
};

} // namespace FockMatrix
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_PREBOHELPER_H
