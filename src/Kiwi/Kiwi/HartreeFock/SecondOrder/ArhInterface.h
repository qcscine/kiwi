/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_ARHINTERFACE_H
#define KIWI_ARHINTERFACE_H

#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/HartreeFock/SecondOrder/ArhBuilder.h>
#include <numeric>

namespace Scine {
namespace Kiwi {

class ArhInterface {
 private:
  std::shared_ptr<Molecule> molecule_;
  std::shared_ptr<Data> data_;
  std::shared_ptr<ArhData> ptrArhData_;

  std::map<std::pair<Utils::ElementType, SpinFunction>, ArhBuilder> mapArhBuilder_;

  std::map<int, std::pair<Utils::ElementType, SpinFunction>> typeMapping_;

  std::vector<int> dim;

  std::vector<int> offset;

 private:
  int numberOfEvaluations_ = 0;

  int size;

 public:
  ArhInterface(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data)
    : molecule_(std::move(molecule)), data_(std::move(data)) {
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

    // 2.
    ptrArhData_ = std::make_shared<ArhData>();
    ptrArhData_->size = 0;
    ptrArhData_->molecule = molecule_;
    ptrArhData_->data = data_;

    for (auto const& elem : typeMapping_) {
      mapArhBuilder_.insert({elem.second, ArhBuilder(molecule_->at(elem.second.first).LMO, elem.second.first,
                                                     elem.second.second, ptrArhData_)});
    }

    // 3.
    initAdata();
  }

  // **********************************************************************
  // *  The functions below are to be used in the sigma vector evaluator  *
  // **********************************************************************

  auto updateX(const Eigen::VectorXd& vectorX) -> void {
    for (auto i = 0UL; i < dim.size(); ++i) {
      Eigen::VectorXd seg = vectorX.segment(offset[i], dim[i]);
      mapArhBuilder_.at(typeMapping_[i]).updateX(seg);
    }
  }

  auto evaluateGradient() -> void {
    for (auto& elem : mapArhBuilder_) {
      elem.second.addMOs();
      elem.second.evaluateGradient();
    }
  }

  auto getGradient() -> Eigen::VectorXd {
    int fullDim = std::accumulate(dim.begin(), dim.end(), 0);

    Eigen::VectorXd ret(fullDim);
    for (auto i = 0UL; i < dim.size(); ++i) {
      ret.segment(offset[i], dim[i]) = getGradient(i);
    }

    return ret;
  }

  auto evaluate(int at) -> Eigen::VectorXd {
    return mapArhBuilder_.at(typeMapping_[at]).evaluate();
  }

  auto getRhDiagonal(int at) -> Eigen::DiagonalMatrix<double, Eigen::Dynamic> {
    return mapArhBuilder_.at(typeMapping_[at]).getRhDiagonal();
  }

  auto getThirdTrialVector(int at) -> Eigen::VectorXd {
    return mapArhBuilder_.at(typeMapping_[at]).getThirdTrialVector();
  }

  void setMaxSize(int maxSize) {
    ptrArhData_->maxSize = maxSize;
  }

  void setUseExactHessian(bool useExactHessian) {
    for (auto& builder : mapArhBuilder_) {
      builder.second.setUseExactHessian(useExactHessian);
    }
  }

  auto getGradient(int at) -> Eigen::VectorXd {
    return mapArhBuilder_.at(typeMapping_[at]).getGradientVector();
  }

  auto makeAllContributions() -> void {
    ++numberOfEvaluations_;
    for (auto& builder : mapArhBuilder_) {
      builder.second.makeAllContributions();
    }
  }
  auto makeContributionsIndepentOfX() -> void {
    for (auto& builder : mapArhBuilder_) {
      builder.second.makeXindependentContributions();
    }
  }

 private:
 public:
  auto makeContributionsDependentOnX() -> void {
    ++numberOfEvaluations_;
    for (auto& builder : mapArhBuilder_) {
      builder.second.makeXdependentContributions();
    }
  }

  auto setLevelShift(double levelShift) -> void {
    for (auto& builder : mapArhBuilder_) {
      builder.second.setMu(levelShift);
      builder.second.addLevelShift(true);
    }
  }

  auto removeLevelShift() -> void {
    for (auto& builder : mapArhBuilder_) {
      builder.second.addLevelShift(false);
    }
  }

  // **************************************************************************
  // *  The functions below are to be used in the macro optimization routine  *
  // **************************************************************************

  auto addIteration() -> void {
    ptrArhData_->F_OAO = data_->hartreeFockData->F_OAO;

    if (ptrArhData_->size < ptrArhData_->maxSize) {
      ++ptrArhData_->size;
      for (auto const& elem : *molecule_) {
        addDataAt(elem.first);
      }
    }
    else {
      for (auto const& elem : *molecule_) {
        auto type = elem.first;
        ptrArhData_->D_OAO.at(type).pop_front();
        ptrArhData_->G.at(type).pop_front();
        if (!elem.second.isRestricted && elem.second.msVector[1] > 0) {
          ptrArhData_->J.at(type).pop_front();
        }
        if (molecule_->size() > 1) {
          ptrArhData_->L.at(type).pop_front();
        }
        addDataAt(type);
      }
    }

    makeContributionsIndepentOfX();
  }

  auto removeLastIteration() -> void {
    --ptrArhData_->size;
    for (auto const& elem : *molecule_) {
      auto type = elem.first;
      ptrArhData_->D_OAO.at(type).pop_back();
      ptrArhData_->G.at(type).pop_back();
      if (!elem.second.isRestricted && elem.second.msVector[1] > 0) {
        ptrArhData_->J.at(type).pop_back();
      }
      if (molecule_->size() > 1) {
        ptrArhData_->L.at(type).pop_back();
      }
    }
  }

  auto reset() -> void {
    for (auto i = 0; i < ptrArhData_->size; ++i) {
      removeLastIteration();
    }
  }

  auto getOrbitalRotationMatrix(Utils::ElementType type) -> Utils::SpinAdaptedMatrix {
    Utils::SpinAdaptedMatrix ret;

    if (molecule_->at(type).isRestricted) {
      ret.restrictedMatrix() = mapArhBuilder_.at({type, SpinFunction::Restricted}).getTransformationMatrix();
    }
    else {
      ret.alphaMatrix() = mapArhBuilder_.at({type, SpinFunction::Alpha}).getTransformationMatrix();
      if (molecule_->at(type).msVector[1] > 0) {
        ret.betaMatrix() = mapArhBuilder_.at({type, SpinFunction::Beta}).getTransformationMatrix();
      }
    }

    return ret;
  }

  // *************************
  // *  Getters and Setters  *
  // *************************

  int getSize() const {
    return size;
  }

  const std::vector<int>& getDim() const {
    return dim;
  }
  const std::vector<int>& getOffset() const {
    return offset;
  }
  const std::map<std::pair<Utils::ElementType, SpinFunction>, ArhBuilder>& getMapArhBuilder() const {
    return mapArhBuilder_;
  }
  const std::shared_ptr<ArhData>& getPtrArhData() const {
    return ptrArhData_;
  }
  const std::map<int, std::pair<Utils::ElementType, SpinFunction>>& getTypeMapping() const {
    return typeMapping_;
  }

 private:
  auto initAdata() -> void {
    for (auto const& elem : *molecule_) {
      const auto& type = elem.first;
      // This is here such that Xi is once initialized as a proper density matrix (needed for the integral routines)
      ptrArhData_->Xi[type] = data_->D_OAO.at(type);
      ptrArhData_->HXi[type] = Utils::SpinAdaptedMatrix();
      ptrArhData_->gradient[type] = Utils::SpinAdaptedMatrix();
      ptrArhData_->F_OAO[type] = Utils::SpinAdaptedMatrix();
      ptrArhData_->T[type] = Utils::SpinAdaptedMatrix();
      ptrArhData_->Intermediate1[type] = Utils::SpinAdaptedVector();
      ptrArhData_->Intermediate2[type] = Utils::SpinAdaptedVector();
      ptrArhData_->G[type] = std::deque<Utils::SpinAdaptedMatrix>();
      ptrArhData_->J[type] = std::deque<Utils::SpinAdaptedMatrix>();
      ptrArhData_->L[type] = std::deque<std::map<Utils::ElementType, Eigen::MatrixXd>>();
      ptrArhData_->D_OAO[type] = std::deque<Utils::DensityMatrix>();
      ptrArhData_->D_in[type] = std::vector<Utils::SpinAdaptedMatrix>();
    }
  }

  auto addDataAt(const Utils::ElementType& type) -> void {
    ptrArhData_->D_OAO.at(type).push_back(data_->D_OAO.at(type));
    const auto& AO2OAO = data_->X.at(type);

    // Restricted
    if (data_->molecule->at(type).isRestricted) {
      ptrArhData_->D_OAO.at(type).back().getRestrictedMatrix() *= 0.5;
      ptrArhData_->G.at(type).push_back(Utils::SpinAdaptedMatrix());
      ptrArhData_->G.at(type).back().restrictedMatrix() =
          AO2OAO.transpose() * data_->hartreeFockData->G.at(type).restrictedMatrix() * AO2OAO;
    }
    // Unrestricted
    else {
      ptrArhData_->G.at(type).push_back(Utils::SpinAdaptedMatrix());
      Eigen::MatrixXd G_alpha =
          data_->hartreeFockData->J.at(type).alphaMatrix() - data_->hartreeFockData->K.at(type).alphaMatrix();
      ptrArhData_->G.at(type).back().alphaMatrix() = AO2OAO.transpose() * G_alpha * AO2OAO;

      if (data_->molecule->at(type).msVector[1] > 0) {
        // Alpha part:
        ptrArhData_->J.at(type).push_back(Utils::SpinAdaptedMatrix());
        ptrArhData_->J.at(type).back().alphaMatrix() =
            AO2OAO.transpose() * data_->hartreeFockData->J.at(type).alphaMatrix() * AO2OAO;

        // Beta part:
        Eigen::MatrixXd G_beta =
            data_->hartreeFockData->J.at(type).betaMatrix() - data_->hartreeFockData->K.at(type).betaMatrix();
        ptrArhData_->G.at(type).back().betaMatrix() = AO2OAO.transpose() * G_beta * AO2OAO;
        ptrArhData_->J.at(type).back().betaMatrix() =
            AO2OAO.transpose() * data_->hartreeFockData->J.at(type).betaMatrix() * AO2OAO;
      }
    }

    // Pre-BO specific part:
    if (molecule_->size() > 1) {
      std::map<Utils::ElementType, Eigen::MatrixXd> tmp;
      for (auto const& elem : *molecule_) {
        if (elem.first != type) {
          auto otherType = elem.first;
          tmp[otherType] = Eigen::MatrixXd();
          // Restricted
          if (elem.second.isRestricted) {
            tmp[elem.first] = AO2OAO.transpose() * data_->hartreeFockData->L.at(type).at(otherType) * AO2OAO;
          }
          // Unrestricted
          else {
            tmp[elem.first] = AO2OAO.transpose() * data_->hartreeFockData->L.at(type).at(otherType) * AO2OAO;
          }
        }
      }
      ptrArhData_->L.at(type).push_back(tmp);
    } // preBO
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_ARHINTERFACE_H
