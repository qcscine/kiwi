/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_INTERFACE_H
#define KIWI_INTERFACE_H

#include <Kiwi/HartreeFock/BFGS/Builder.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/KiwiOpt/Optimizer.h>
#include <numeric>

namespace Scine {
namespace Kiwi {
namespace BFGS {

template<SymmetryType Symmetry>
class Interface {
 private:
  std::shared_ptr<Molecule> molecule_;
  std::shared_ptr<Kiwi::Data> data_;
  std::shared_ptr<BFGS::Data<Symmetry>> firstOrderData_;

  std::map<std::pair<Utils::ElementType, SpinFunction>, Builder<Symmetry>> mapBuilder_;

  std::map<int, std::pair<Utils::ElementType, SpinFunction>> typeMapping_;

  std::vector<int> dim;

  std::vector<int> offset;

  int numberOfEvaluations_ = 0;

  int size;

  double E_alpha;

  double E_0;

 public:
  Interface() = delete;

  Interface(std::shared_ptr<Molecule> molecule, std::shared_ptr<Kiwi::Data> data,
            ProjectionParameters projectionParameters = ProjectionParameters())
    : molecule_(std::move(molecule)),
      data_(std::move(data)),
      firstOrderData_(std::make_shared<BFGS::Data<Symmetry>>(molecule_, data_, projectionParameters)) {
    // 1. Get the mapping for the solution vector of the Newton-equations:
    //
    //    H x = -h
    //    x collects in a column vector the orbital rotation parameters for all paricle
    //    types and spins.
    //    x = ( e-_up e-_up e-_down e-_down p+_up p+_up p+_down p+_down )^T
    //    We can access the correct segment of the vector with `dim` and `offset`.
    //    On the other hand, with type mapping, we can access the correct type and spin,
    //    given the index n, which corresponds to the n-th distinguishable particle type.
    //
    int key = 0;
    for (auto const& elem : *molecule_) {
      if (elem.second.isRestricted) {
        typeMapping_[key] = {elem.first, SpinFunction::Restricted};
        ++key;
        dim.push_back(elem.second.occ.restricted * elem.second.virt.restricted);
      }
      else {
        dim.push_back(elem.second.occ.alpha * elem.second.virt.alpha);
        typeMapping_[key] = {elem.first, SpinFunction::Alpha};
        ++key;
        if (elem.second.msVector[1] > 0) {
          dim.push_back(elem.second.occ.beta * elem.second.virt.beta);
          typeMapping_[key] = {elem.first, SpinFunction::Beta};
          ++key;
        }
      }
    }
    size = dim.size();
    offset.resize(dim.size());
    offset[0] = 0;
    for (auto i = 1UL; i < dim.size(); ++i) {
      offset[i] = dim[i - 1] + offset[i - 1];
    }

    for (auto const& elem : typeMapping_) {
      mapBuilder_.insert({elem.second, BFGS::Builder(molecule_->at(elem.second.first).LMO, elem.second.first,
                                                     elem.second.second, firstOrderData_)});
    }

    for (auto& builder : mapBuilder_) {
      builder.second.setUseDiis(false);
    }

    // 3.
    initData();
  }

  // **********************************************************************
  // *  The functions below are to be used in the sigma vector evaluator  *
  // **********************************************************************

  auto evaluateGradient() -> void {
    for (auto& builder : mapBuilder_) {
      builder.second.evaluateGradient();
    }
  }

  auto evaluateDirection(bool BFGS = true) -> void {
    for (auto& builder : mapBuilder_) {
      builder.second.evaluateDirection(BFGS);
      builder.second.storeGradient();
    }
  }

  auto evaluate(double alpha) -> void {
    firstOrderData_->projector.resetDensity(firstOrderData_->D_0, firstOrderData_->D_OAO_0, firstOrderData_->C_0,
                                            firstOrderData_->C_OAO_0);

    for (auto const& elem : *molecule_) {
      auto orbitalRotation = getOrbitalRotationMatrix(elem.first, alpha);
      // std::cout << "R\n" << orbitalRotation.restrictedMatrix() << std::endl;
      firstOrderData_->projector.updateDensityFromRotation(elem.first, orbitalRotation);
    }

    // firstOrderData_->projector.evaluateEnergy();
    // std::cout << "E before\n" << firstOrderData_->projector.getEnergy() << std::endl;

    firstOrderData_->projector.evaluateDensity();
    firstOrderData_->projector.updateFockMatrices();
    firstOrderData_->projector.evaluateEnergy();

    E_alpha = firstOrderData_->projector.getEnergy();

    // std::cout << "E after\n" << firstOrderData_->projector.getEnergy() << std::endl;

    evaluateGradient();
  }

  auto applyUpdate(bool BFGS = true, double alpha = 1) -> void {
    firstOrderData_->D_0 = data_->D;
    firstOrderData_->D_OAO_0 = data_->D_OAO;
    firstOrderData_->C_0 = data_->C;
    firstOrderData_->C_OAO_0 = data_->C_OAO;
    for (auto& builder : mapBuilder_) {
      builder.second.addMOs();
    }

    if (BFGS) {
      for (auto& builder : mapBuilder_) {
        builder.second.updateBFGS(alpha);
      }
    }

    evaluateDirection(BFGS);
  }

  auto getRestrictionAdaptedGradient() -> Eigen::VectorXd {
    int fullDim = std::accumulate(dim.begin(), dim.end(), 0);

    Eigen::VectorXd ret(fullDim);
    for (auto i = 0UL; i < dim.size(); ++i) {
      double factor = 1;
      if (molecule_->at(typeMapping_.at(i).first).isRestricted) {
        factor = 2;
      }
      ret.segment(offset[i], dim[i]) = factor * getGradient(i);
    }

    return ret;
  }

  // auto evaluate(int at) -> Eigen::VectorXd {
  //  return mapBuilder_.at(typeMapping_[at]).evaluate();
  //}

  auto getDirection(int at) -> Eigen::VectorXd {
    return mapBuilder_.at(typeMapping_[at]).getDirectionVector();
  }

  auto getGradient(int at) -> Eigen::VectorXd {
    return mapBuilder_.at(typeMapping_[at]).getGradientVector();
  }

  auto getOrbitalRotationMatrix(Utils::ElementType type, double alpha) -> Utils::SpinAdaptedMatrix {
    Utils::SpinAdaptedMatrix ret;

    if (molecule_->at(type).isRestricted) {
      ret.restrictedMatrix() = mapBuilder_.at({type, SpinFunction::Restricted}).getRetractionMatrix(alpha);
    }
    else {
      ret.alphaMatrix() = mapBuilder_.at({type, SpinFunction::Alpha}).getRetractionMatrix(alpha);
      if (molecule_->at(type).msVector[1] > 0) {
        ret.betaMatrix() = mapBuilder_.at({type, SpinFunction::Beta}).getRetractionMatrix(alpha);
      }
    }

    return ret;
  }

  // **************************************************************************
  // *  The functions below are to be used in the macro optimization routine  *
  // **************************************************************************

  // *************************
  // *  Getters and Setters  *
  // *************************

  [[nodiscard]] double getEnergyAlpha() const {
    return E_alpha;
  }

  [[nodiscard]] int getSize() const {
    return size;
  }

  [[nodiscard]] const std::vector<int>& getDim() const {
    return dim;
  }
  [[nodiscard]] const std::vector<int>& getOffset() const {
    return offset;
  }
  [[nodiscard]] const std::map<std::pair<Utils::ElementType, SpinFunction>, BFGS::Builder<Symmetry>>& getMapArhBuilder() const {
    return mapBuilder_;
  }
  const std::shared_ptr<BFGS::Data<Symmetry>>& getPtrData() const {
    return firstOrderData_;
  }
  [[nodiscard]] const std::map<int, std::pair<Utils::ElementType, SpinFunction>>& getTypeMapping() const {
    return typeMapping_;
  }

 private:
  auto initData() -> void {
    firstOrderData_->projector.init(true);

    E_alpha = firstOrderData_->projector.getEnergy();

    applyUpdate(false);

    for (auto& builder : mapBuilder_) {
      builder.second.initBFGS();
    }
  }
};

} // namespace BFGS
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_INTERFACE_H
