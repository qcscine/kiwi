/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_BUILDER_H
#define KIWI_BUILDER_H

#include <Kiwi/HartreeFock/Projector/ProjectorTraits.h>
#include <Kiwi/KiwiOpt/Davidson.h>
#include <Kiwi/KiwiOpt/GDIIS.h>
#include <Kiwi/KiwiUtils/SymmetryEnums.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/DataStructures/SpinAdaptedVector.h>
#include <Utils/Geometry/ElementTypes.h>
#include <map>
#include <memory>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

namespace Scine {
namespace Kiwi {

class Molecule;
class Data;

namespace BFGS {

template<SymmetryType Symmetry>
struct Data {
  using ProjectorType = typename ProjectorTrait<Symmetry>::ProjectorType;

  int size;
  int maxSize;

  double E_0;

  double gradientDotKappa;

  std::shared_ptr<Molecule> molecule;
  std::shared_ptr<Kiwi::Data> data;
  ProjectorType projector;

  //// The fock matrix is required for Gradient and Roothaan-Hall part.
  // std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> F_OAO;

  std::map<Utils::ElementType, Utils::DensityMatrix> D_0;
  std::map<Utils::ElementType, Utils::DensityMatrix> D_OAO_0;
  std::map<Utils::ElementType, Utils::MolecularOrbitals> C_0;
  std::map<Utils::ElementType, Utils::MolecularOrbitals> C_OAO_0;

  Data() = delete;

  Data(std::shared_ptr<Molecule> mol, std::shared_ptr<Kiwi::Data> dat,
       ProjectionParameters projectionParameters = ProjectionParameters())
    : molecule(std::move(mol)), data(std::move(dat)), projector(ProjectorType(molecule, data, projectionParameters)) {
  }
};

template<SymmetryType Symmetry>
class Builder {
  bool init = true;

  Eigen::MatrixXd gradient;
  Eigen::MatrixXd direction;

  Eigen::MatrixXd B;
  std::shared_ptr<Eigen::MatrixXd> ptr_B;

  bool useDiis = false;
  bool useDamping = false;

 public:
  void setUseDiis(bool useDiis);

 private:
  Optimization::GDIIS gdiis;

  Eigen::MatrixXd gradient_old;
  Eigen::VectorXd s;
  Eigen::VectorXd y;
  Eigen::VectorXd s_By;

  Eigen::VectorXd B_diag;

  const int dim;

  const SpinFunction spin;

  const Utils::ElementType type;

  const std::shared_ptr<Data<Symmetry>> firstOrderData;

 public:
  Builder(int dimension, Utils::ElementType tp, SpinFunction sp, std::shared_ptr<BFGS::Data<Symmetry>> dat)
    : dim(dimension), spin(sp), type(tp), firstOrderData(std::move(dat)) {
    switch (spin) {
      case SpinFunction::Restricted: {
        virt = firstOrderData->molecule->at(type).virt.restricted;
        occ = firstOrderData->molecule->at(type).occ.restricted;
        break;
      }
      case SpinFunction::Alpha: {
        virt = firstOrderData->molecule->at(type).virt.alpha;
        occ = firstOrderData->molecule->at(type).occ.alpha;
        break;
      }
      case SpinFunction::Beta: {
        virt = firstOrderData->molecule->at(type).virt.beta;
        occ = firstOrderData->molecule->at(type).occ.beta;
        break;
      }
    }
  }

  auto addMOs() -> void;

  [[nodiscard]] auto getRetractionMatrix(double alpha = 1.) -> Eigen::MatrixXd;

  [[nodiscard]] auto getGradientVector() const -> Eigen::VectorXd;

  [[nodiscard]] auto getDirectionVector() const -> Eigen::VectorXd;

  auto evaluateDirection(bool BFGS = false) -> void;

  auto evaluateGradient() -> void;

  auto initBFGS() -> void;

  auto updateBFGS(double alpha) -> void;

  auto storeGradient() -> void;

 private:
  int occ;
  int virt;

  Eigen::MatrixXd Theta_0_virt;
  Eigen::MatrixXd Theta_0_occ;
  Eigen::MatrixXd Theta_0;
};

template<SymmetryType Symmetry>
auto Builder<Symmetry>::getRetractionMatrix(double alpha) -> Eigen::MatrixXd {
  // MO
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(occ + virt, occ + virt);
  Eigen::MatrixXd X(occ + virt, occ + virt);
  X.setZero();
  X.block(occ, 0, virt, occ) = alpha * direction;
  X.block(0, occ, occ, virt) = -alpha * direction.transpose();
  Eigen::MatrixXd mat = Id - X;
  Eigen::MatrixXd factor_occ = Eigen::MatrixXd::Identity(occ, occ);
  Eigen::MatrixXd factor_virt = Eigen::MatrixXd::Identity(virt, virt);
  Eigen::MatrixXd Z = mat.block(0, occ, occ, virt);
  Eigen::MatrixXd ZT = Z.transpose();

  factor_occ += Z * ZT;
  factor_virt += ZT * Z;

  Eigen::MatrixXd R11 = factor_occ.llt().matrixU();

  Eigen::MatrixXd R11_inv = R11.inverse();
  Eigen::MatrixXd R22 = factor_virt.llt().matrixU();
  Eigen::MatrixXd R22_inv = R22.inverse();

  Eigen::MatrixXd Q(occ + virt, occ + virt);

  Q.block(0, 0, occ, occ) = R11_inv;
  Q.block(occ, occ, virt, virt) = R22_inv;

  Q.block(0, occ, occ, virt) = Z * R22_inv;
  Q.block(occ, 0, virt, occ) = -ZT * R11_inv;

  return Optimization::Davidson::modifiedGramSchmidt(Theta_0 * Q * Theta_0.transpose());
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::addMOs() -> void {
  // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

  // switch (spin) {
  //  case SpinFunction::Restricted: {
  //    es.compute(firstOrderData->D_OAO.at(type).restrictedMatrix());
  //    break;
  //  }
  //  case SpinFunction::Alpha: {
  //    es.compute(firstOrderData->D_OAO.at(type).alphaMatrix());
  //    break;
  //  }
  //  case SpinFunction::Beta: {
  //    es.compute(firstOrderData->D_OAO.at(type).alphaMatrix());
  //    break;
  //  }
  //}
  ////Theta = es.eigenvectors().rowwise().reverse();
  // Theta = es.eigenvectors().rowwise().reverse();
  // Theta_virt = Theta.block(0, occ, dim, virt);
  // Theta_occ = Theta.block(0, 0, dim, occ);

  // switch (spin) {
  //  case SpinFunction::Restricted: {
  //    Eigen::MatrixXd& Fortho = firstOrderData->data->hartreeFockData->F_OAO.at(type).restrictedMatrix();
  //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  //    Eigen::MatrixXd F_MO = firstOrderData->C_OAO_0.at(type).restrictedMatrix().transpose() * Fortho *
  //                           firstOrderData->C_OAO_0.at(type).restrictedMatrix();

  //    Eigen::MatrixXd transformation;
  //    transformation.resizeLike(F_MO);
  //    transformation.setIdentity();
  //    es.compute(F_MO.block(0, 0, occ, occ));
  //    transformation.block(0, 0, occ, occ) = es.eigenvectors();
  //    es.compute(F_MO.block(occ, occ, virt, virt));
  //    transformation.block(occ, occ, virt, virt) = es.eigenvectors();
  //    Theta_0 = firstOrderData->C_OAO_0.at(type).restrictedMatrix() * transformation;
  //    Theta_0_virt = Theta_0.block(0, occ, dim, virt);
  //    Theta_0_occ = Theta_0.block(0, 0, dim, occ);
  //    break;
  //  }
  //  case SpinFunction::Alpha: {
  //    Eigen::MatrixXd& Fortho = firstOrderData->data->hartreeFockData->F_OAO.at(type).alphaMatrix();
  //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  //    Eigen::MatrixXd F_MO = firstOrderData->C_OAO_0.at(type).alphaMatrix().transpose() * Fortho *
  //                           firstOrderData->C_OAO_0.at(type).alphaMatrix();

  //    Eigen::MatrixXd transformation;
  //    transformation.resizeLike(F_MO);
  //    transformation.setIdentity();
  //    es.compute(F_MO.block(0, 0, occ, occ));
  //    transformation.block(0, 0, occ, occ) = es.eigenvectors();
  //    es.compute(F_MO.block(occ, occ, virt, virt));
  //    transformation.block(occ, occ, virt, virt) = es.eigenvectors();
  //    Theta_0 = firstOrderData->C_OAO_0.at(type).alphaMatrix() * transformation;
  //    Theta_0_virt = Theta_0.block(0, occ, dim, virt);
  //    Theta_0_occ = Theta_0.block(0, 0, dim, occ);
  //    break;
  //  }
  //  case SpinFunction::Beta: {
  //    Eigen::MatrixXd& Fortho = firstOrderData->data->hartreeFockData->F_OAO.at(type).betaMatrix();
  //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  //    Eigen::MatrixXd F_MO = firstOrderData->C_OAO_0.at(type).betaMatrix().transpose() * Fortho *
  //                           firstOrderData->C_OAO_0.at(type).betaMatrix();

  //    Eigen::MatrixXd transformation;
  //    transformation.resizeLike(F_MO);
  //    transformation.setIdentity();
  //    es.compute(F_MO.block(0, 0, occ, occ));
  //    transformation.block(0, 0, occ, occ) = es.eigenvectors();
  //    es.compute(F_MO.block(occ, occ, virt, virt));
  //    transformation.block(occ, occ, virt, virt) = es.eigenvectors();
  //    Theta_0 = firstOrderData->C_OAO_0.at(type).betaMatrix() * transformation;
  //    Theta_0_virt = Theta_0.block(0, occ, dim, virt);
  //    Theta_0_occ = Theta_0.block(0, 0, dim, occ);
  //    break;
  //  }
  //}

  switch (spin) {
    case SpinFunction::Restricted: {
      Theta_0 = firstOrderData->data->C_OAO.at(type).restrictedMatrix();
      break;
    }
    case SpinFunction::Alpha: {
      Theta_0 = firstOrderData->data->C_OAO.at(type).alphaMatrix();
      break;
    }
    case SpinFunction::Beta: {
      Theta_0 = firstOrderData->data->C_OAO.at(type).betaMatrix();
      break;
    }
  }

  Theta_0_virt = Theta_0.block(0, occ, dim, virt);
  Theta_0_occ = Theta_0.block(0, 0, dim, occ);
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::getGradientVector() const -> Eigen::VectorXd {
  Eigen::VectorXd ret = Eigen::Map<const Eigen::VectorXd>(gradient.data(), occ * virt);
  return ret;
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::getDirectionVector() const -> Eigen::VectorXd {
  Eigen::VectorXd ret = Eigen::Map<const Eigen::VectorXd>(direction.data(), occ * virt);
  return ret;
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::evaluateDirection(bool BFGS) -> void {
  if (BFGS) {
    // std::cout << "bfgs;\n";
    Eigen::VectorXd tmp = -B * Eigen::Map<Eigen::VectorXd>(gradient.data(), occ * virt);
    direction = Eigen::Map<Eigen::MatrixXd>(tmp.data(), virt, occ);
  }
  else {
    // std::cout << "no bfgs;\n";
    direction = -gradient;
  }
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::evaluateGradient() -> void {
  Eigen::MatrixXd Theta_current;
  Eigen::MatrixXd Theta_current_virt;
  Eigen::MatrixXd Theta_current_occ;

  switch (spin) {
    case SpinFunction::Restricted: {
      gradient.noalias() = 0.5 * (firstOrderData->data->D_OAO.at(type).restrictedMatrix() *
                                      firstOrderData->data->hartreeFockData->F_OAO.at(type).restrictedMatrix() -
                                  firstOrderData->data->hartreeFockData->F_OAO.at(type).restrictedMatrix() *
                                      firstOrderData->data->D_OAO.at(type).restrictedMatrix());
      Theta_current = firstOrderData->data->C_OAO.at(type).restrictedMatrix();
      break;
    }
    case SpinFunction::Alpha: {
      gradient.noalias() = (firstOrderData->data->D_OAO.at(type).alphaMatrix() *
                                firstOrderData->data->hartreeFockData->F_OAO.at(type).alphaMatrix() -
                            firstOrderData->data->hartreeFockData->F_OAO.at(type).alphaMatrix() *
                                firstOrderData->data->D_OAO.at(type).alphaMatrix());
      Theta_current = firstOrderData->data->C_OAO.at(type).alphaMatrix();
      break;
    }
    case SpinFunction::Beta: {
      gradient.noalias() = (firstOrderData->data->D_OAO.at(type).betaMatrix() *
                                firstOrderData->data->hartreeFockData->F_OAO.at(type).betaMatrix() -
                            firstOrderData->data->hartreeFockData->F_OAO.at(type).betaMatrix() *
                                firstOrderData->data->D_OAO.at(type).betaMatrix());
      Theta_current = firstOrderData->data->C_OAO.at(type).betaMatrix();
      break;
    }
  }

  // switch (spin) {
  //  case SpinFunction::Restricted: {
  //    Eigen::MatrixXd& Fortho = firstOrderData->data->hartreeFockData->F_OAO.at(type).restrictedMatrix();
  //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  //    Eigen::MatrixXd F_MO = firstOrderData->data->C_OAO.at(type).restrictedMatrix().transpose() * Fortho *
  //                           firstOrderData->data->C_OAO.at(type).restrictedMatrix();

  //    Eigen::MatrixXd transformation;
  //    transformation.resizeLike(F_MO);
  //    transformation.setIdentity();
  //    es.compute(F_MO.block(0, 0, occ, occ));
  //    transformation.block(0, 0, occ, occ) = es.eigenvectors();
  //    es.compute(F_MO.block(occ, occ, virt, virt));
  //    transformation.block(occ, occ, virt, virt) = es.eigenvectors();
  //    Theta_current = firstOrderData->data->C_OAO.at(type).restrictedMatrix() * transformation;
  //    Theta_current_virt = Theta_current.block(0, occ, dim, virt);
  //    Theta_current_occ = Theta_current.block(0, 0, dim, occ);
  //    break;
  //  }
  //  case SpinFunction::Alpha: {
  //    Eigen::MatrixXd& Fortho = firstOrderData->data->hartreeFockData->F_OAO.at(type).alphaMatrix();
  //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  //    Eigen::MatrixXd F_MO = firstOrderData->data->C_OAO.at(type).alphaMatrix().transpose() * Fortho *
  //                           firstOrderData->data->C_OAO.at(type).alphaMatrix();

  //    Eigen::MatrixXd transformation;
  //    transformation.resizeLike(F_MO);
  //    transformation.setIdentity();
  //    es.compute(F_MO.block(0, 0, occ, occ));
  //    transformation.block(0, 0, occ, occ) = es.eigenvectors();
  //    es.compute(F_MO.block(occ, occ, virt, virt));
  //    transformation.block(occ, occ, virt, virt) = es.eigenvectors();
  //    Theta_current = firstOrderData->data->C_OAO.at(type).alphaMatrix() * transformation;
  //    Theta_current_virt = Theta_current.block(0, occ, dim, virt);
  //    Theta_current_occ = Theta_current.block(0, 0, dim, occ);
  //    break;
  //  }
  //  case SpinFunction::Beta: {
  //    Eigen::MatrixXd& Fortho = firstOrderData->data->hartreeFockData->F_OAO.at(type).betaMatrix();
  //    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  //    Eigen::MatrixXd F_MO = firstOrderData->data->C_OAO.at(type).betaMatrix().transpose() * Fortho *
  //                           firstOrderData->data->C_OAO.at(type).betaMatrix();

  //    Eigen::MatrixXd transformation;
  //    transformation.resizeLike(F_MO);
  //    transformation.setIdentity();
  //    es.compute(F_MO.block(0, 0, occ, occ));
  //    transformation.block(0, 0, occ, occ) = es.eigenvectors();
  //    es.compute(F_MO.block(occ, occ, virt, virt));
  //    transformation.block(occ, occ, virt, virt) = es.eigenvectors();
  //    Theta_current = firstOrderData->data->C_OAO.at(type).betaMatrix() * transformation;
  //    Theta_current_virt = Theta_current.block(0, occ, dim, virt);
  //    Theta_current_occ = Theta_current.block(0, 0, dim, occ);
  //    break;
  //  }
  //}

  Theta_current_virt = Theta_current.block(0, occ, dim, virt);
  Theta_current_occ = Theta_current.block(0, 0, dim, occ);
  gradient = Theta_current_virt.transpose() * gradient * Theta_current_occ;
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::initBFGS() -> void {
  B = Eigen::MatrixXd::Identity(occ * virt, occ * virt);

  // preconditioning
  bool precond = false;
  if (precond) {
    Eigen::MatrixXd F_vv_oo_MO;

    switch (spin) {
      case SpinFunction::Restricted: {
        F_vv_oo_MO.noalias() = firstOrderData->data->hartreeFockData->F_OAO.at(type).restrictedMatrix() -
                               0.5 * firstOrderData->data->hartreeFockData->F_OAO.at(type).restrictedMatrix() *
                                   firstOrderData->data->D_OAO.at(type).restrictedMatrix() -
                               0.5 * firstOrderData->data->D_OAO.at(type).restrictedMatrix() *
                                   firstOrderData->data->hartreeFockData->F_OAO.at(type).restrictedMatrix();
        break;
      }
      case SpinFunction::Alpha: {
        F_vv_oo_MO.noalias() = firstOrderData->data->hartreeFockData->F_OAO.at(type).alphaMatrix() -
                               firstOrderData->data->hartreeFockData->F_OAO.at(type).alphaMatrix() *
                                   firstOrderData->data->D_OAO.at(type).alphaMatrix() -
                               firstOrderData->data->D_OAO.at(type).alphaMatrix() *
                                   firstOrderData->data->hartreeFockData->F_OAO.at(type).alphaMatrix();
        break;
      }
      case SpinFunction::Beta: {
        F_vv_oo_MO.noalias() = firstOrderData->data->hartreeFockData->F_OAO.at(type).betaMatrix() -
                               firstOrderData->data->hartreeFockData->F_OAO.at(type).betaMatrix() *
                                   firstOrderData->data->D_OAO.at(type).betaMatrix() -
                               firstOrderData->data->D_OAO.at(type).betaMatrix() *
                                   firstOrderData->data->hartreeFockData->F_OAO.at(type).betaMatrix();
        break;
      }
    }

    F_vv_oo_MO = Theta_0.transpose() * F_vv_oo_MO * Theta_0;

    Eigen::MatrixXd diag(virt, occ);

    diag.setZero();

    for (int a = occ; a < occ + virt; ++a) {
      for (int i = 0; i < occ; ++i) {
        diag(a - occ, i) = F_vv_oo_MO(a, a) + F_vv_oo_MO(i, i);
      }
    }

    B_diag = Eigen::Map<Eigen::VectorXd>(diag.data(), virt * occ);
    B_diag = B_diag.cwiseInverse();

    B.diagonal() = B_diag;
  }
  else {
    B_diag = Eigen::VectorXd::Ones(virt * occ);
  }

  if (useDiis) {
    ptr_B = std::make_shared<Eigen::MatrixXd>(B);
    gdiis = Optimization::GDIIS(ptr_B, 10);
  }

  s_By.resize(occ * virt);
  y.resize(occ * virt);
  s.resize(occ * virt);

  gradient_old = gradient;
  init = true;
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::updateBFGS(double alpha) -> void {
  s = alpha * Eigen::Map<const Eigen::VectorXd>(direction.data(), occ * virt);

  y = Eigen::Map<const Eigen::VectorXd>(gradient.data(), occ * virt) -
      Eigen::Map<const Eigen::VectorXd>(gradient_old.data(), occ * virt);

  double yt_s = y.dot(s);

  Eigen::VectorXd By = B * y;

  double ytBy = y.transpose() * By;

  if (useDamping) {
    /* Powell Update from Al-Baali Gandetti */
    const double sigma2 = 0.9;
    const double sigma3 = 9.0;
    double delta = 1.0;
    if (fabs(yt_s) < fabs((1.0 - sigma2) * ytBy)) {
      delta = sigma2 * ytBy / (ytBy - yt_s);
    }
    else if (fabs(yt_s) > fabs((1.0 + sigma3) * ytBy)) {
      delta = -sigma3 * ytBy / (ytBy - yt_s);
    }
    /* Update s by s = delta * s + (1.0 - delta) * invH * dg */
    if (fabs(delta - 1.0) > 1e-12) {
      s.noalias() = delta * s + (1.0 - delta) * By;
      yt_s = s.dot(y);
    }
  }

  // if (!init) {
  if (yt_s > 0) {
    // double yt_s_inv = 1 / yt_s;
    // double yt_s_inv2 = 1 / (yt_s * yt_s);
    // double s_By_t_y = s_By.dot(y);
    // s_By.noalias() = s - B * y;

    // B += yt_s_inv * s_By  * s.transpose();
    // B += yt_s_inv * s  * s_By.transpose();
    // B -= yt_s_inv2 * s_By_t_y * s * s.transpose();

    // B += (yt_s_inv * (s_By * s.transpose() + s * s_By.transpose()) - yt_s_inv2 * s_By_t_y * s *
    //      s.transpose()).eval();

    const double beta1 = (yt_s + ytBy) / (yt_s * yt_s);
    const double beta2 = 1.0 / yt_s;
    B += beta1 * (s * s.transpose()) - beta2 * (By * s.transpose() + s * y.transpose() * B);

    // B += yt_s_inv * (s_By * s.transpose() + s * s_By.transpose()) - yt_s_inv2 * s_By_t_y * s * s.transpose();
  }
  // Reset inverse Hessian if it is not positive definite
  else {
    std::cout << "reset;\n";
    if (useDiis) {
      gdiis.flush();
    }
    // B.setIdentity();
    // B *= 0.5;
    B.setZero();
    B.diagonal() = B_diag;
  }

  if (useDiis) {
    Eigen::VectorXd tmp = Eigen::Map<Eigen::VectorXd>(gradient.data(), occ * virt);
    gdiis.update(s, tmp);
  }
  //}
  // else {
  //  double yt_y = y.dot(y);
  //  B.diagonal() *= 2.0 * yt_s/yt_y;
  //  init = false;
  //}
}

template<SymmetryType Symmetry>
auto Builder<Symmetry>::storeGradient() -> void {
  gradient_old = gradient;
}

template<SymmetryType Symmetry>
void Builder<Symmetry>::setUseDiis(bool diis) {
  Builder::useDiis = diis;
}

} // namespace BFGS
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_BUILDER_H
