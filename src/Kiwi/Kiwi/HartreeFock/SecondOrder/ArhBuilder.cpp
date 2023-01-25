/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/FockMatrix/CoulombExchangeMatrices.h>
#include <Kiwi/HartreeFock/FockMatrix/TwoTypeCoulomb.h>
#include <Kiwi/HartreeFock/SecondOrder/ArhBuilder.h>
#include <Kiwi/KiwiOpt/Davidson.h>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

namespace Scine {
namespace Kiwi {

ArhBuilder::ArhBuilder(int dimension, Utils::ElementType tp, SpinFunction sp, std::shared_ptr<ArhData> dat)
  : dim(dimension), spin(sp), type(tp), arhData(std::move(dat)) {
  switch (spin) {
    case SpinFunction::Restricted: {
      virt = arhData->molecule->at(type).virt.restricted;
      occ = arhData->molecule->at(type).occ.restricted;
      break;
    }
    case SpinFunction::Alpha: {
      virt = arhData->molecule->at(type).virt.alpha;
      occ = arhData->molecule->at(type).occ.alpha;
      break;
    }
    case SpinFunction::Beta: {
      virt = arhData->molecule->at(type).virt.beta;
      occ = arhData->molecule->at(type).occ.beta;
      break;
    }
  }

  // pre-BO
  if (arhData->molecule->size() > 1) {
    for (auto const& elem : *arhData->molecule) {
      if (elem.first != type) {
        L_in[elem.first] = std::vector<Utils::SpinAdaptedMatrix>();
        L_ov_vo_in[elem.first] = std::vector<Utils::SpinAdaptedMatrix>();
      }
    }
  }
}

auto ArhBuilder::addMOs() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      Eigen::MatrixXd& Fortho = arhData->F_OAO.at(type).restrictedMatrix();
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
      Eigen::MatrixXd F_MO = arhData->data->C_OAO.at(type).restrictedMatrix().transpose() * Fortho *
                             arhData->data->C_OAO.at(type).restrictedMatrix();

      Eigen::MatrixXd transformation;
      transformation.resizeLike(F_MO);
      transformation.setIdentity();
      es.compute(F_MO.block(0, 0, occ, occ));
      transformation.block(0, 0, occ, occ) = es.eigenvectors();
      es.compute(F_MO.block(occ, occ, virt, virt));
      transformation.block(occ, occ, virt, virt) = es.eigenvectors();
      Theta = arhData->data->C_OAO.at(type).restrictedMatrix() * transformation;
      Theta_virt = Theta.block(0, occ, dim, virt);
      Theta_occ = Theta.block(0, 0, dim, occ);
      break;
    }
    case SpinFunction::Alpha: {
      Eigen::MatrixXd& Fortho = arhData->F_OAO.at(type).alphaMatrix();
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
      Eigen::MatrixXd F_MO =
          arhData->data->C_OAO.at(type).alphaMatrix().transpose() * Fortho * arhData->data->C_OAO.at(type).alphaMatrix();

      Eigen::MatrixXd transformation;
      transformation.resizeLike(F_MO);
      transformation.setIdentity();
      es.compute(F_MO.block(0, 0, occ, occ));
      transformation.block(0, 0, occ, occ) = es.eigenvectors();
      es.compute(F_MO.block(occ, occ, virt, virt));
      transformation.block(occ, occ, virt, virt) = es.eigenvectors();
      Theta = arhData->data->C_OAO.at(type).alphaMatrix() * transformation;
      Theta_virt = Theta.block(0, occ, dim, virt);
      Theta_occ = Theta.block(0, 0, dim, occ);
      break;
    }
    case SpinFunction::Beta: {
      Eigen::MatrixXd& Fortho = arhData->F_OAO.at(type).betaMatrix();
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
      Eigen::MatrixXd F_MO =
          arhData->data->C_OAO.at(type).betaMatrix().transpose() * Fortho * arhData->data->C_OAO.at(type).betaMatrix();

      Eigen::MatrixXd transformation;
      transformation.resizeLike(F_MO);
      transformation.setIdentity();
      es.compute(F_MO.block(0, 0, occ, occ));
      transformation.block(0, 0, occ, occ) = es.eigenvectors();
      es.compute(F_MO.block(occ, occ, virt, virt));
      transformation.block(occ, occ, virt, virt) = es.eigenvectors();
      Theta = arhData->data->C_OAO.at(type).betaMatrix() * transformation;
      Theta_virt = Theta.block(0, occ, dim, virt);
      Theta_occ = Theta.block(0, 0, dim, occ);
      break;
    }
  }
}

auto ArhBuilder::evaluateGradient() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      G = (arhData->D_OAO.at(type).back().restrictedMatrix() * arhData->F_OAO.at(type).restrictedMatrix() -
           arhData->F_OAO.at(type).restrictedMatrix() * arhData->D_OAO.at(type).back().restrictedMatrix());
      // This guarantees the symmetry of the Hessian.
      G *= 2;
      break;
    }
    case SpinFunction::Alpha: {
      G = (arhData->D_OAO.at(type).back().alphaMatrix() * arhData->F_OAO.at(type).alphaMatrix() -
           arhData->F_OAO.at(type).alphaMatrix() * arhData->D_OAO.at(type).back().alphaMatrix());
      break;
    }
    case SpinFunction::Beta: {
      G = (arhData->D_OAO.at(type).back().betaMatrix() * arhData->F_OAO.at(type).betaMatrix() -
           arhData->F_OAO.at(type).betaMatrix() * arhData->D_OAO.at(type).back().betaMatrix());
      break;
    }
  }

  // Virt-Occ
  G = Theta_virt.transpose() * G * Theta_occ;
}

auto ArhBuilder::addRhContribution() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      const auto& Xi = arhData->Xi.at(type).restrictedMatrix();
      const auto& F = arhData->F_OAO.at(type).restrictedMatrix();
      const auto& D = arhData->D_OAO.at(type).back().restrictedMatrix();

      Eigen::MatrixXd FD = F * D;
      Eigen::MatrixXd DF = FD.transpose();
      Eigen::MatrixXd DFD = D * FD;
      Eigen::MatrixXd Foovv = F - DF - FD + 2 * DFD;
      Eigen::MatrixXd comm = Foovv * Xi;

      HXi = -(comm - comm.transpose()).eval();
      break;
    }
    case SpinFunction::Alpha: {
      const auto& Xi = arhData->Xi.at(type).alphaMatrix();
      const auto& F = arhData->F_OAO.at(type).alphaMatrix();
      const auto& D = arhData->D_OAO.at(type).back().alphaMatrix();

      Eigen::MatrixXd FD = F * D;
      Eigen::MatrixXd DF = FD.transpose();
      Eigen::MatrixXd DFD = D * FD;
      Eigen::MatrixXd Foovv = F - DF - FD + 2 * DFD;
      Eigen::MatrixXd comm = Foovv * Xi;

      HXi = -(comm - comm.transpose()).eval();
      break;
    }
    case SpinFunction::Beta: {
      const auto& Xi = arhData->Xi.at(type).betaMatrix();
      const auto& F = arhData->F_OAO.at(type).betaMatrix();
      const auto& D = arhData->D_OAO.at(type).back().betaMatrix();

      Eigen::MatrixXd FD = F * D;
      Eigen::MatrixXd DF = FD.transpose();
      Eigen::MatrixXd DFD = D * FD;
      Eigen::MatrixXd Foovv = F - DF - FD + 2 * DFD;
      Eigen::MatrixXd comm = Foovv * Xi;

      HXi = -(comm - comm.transpose()).eval();
      break;
    }
  }
}

auto ArhBuilder::addBo2ndOrderSpecificContribution() -> void {
  auto const& AO2OAO = arhData->data->X.at(type);

  Utils::DensityMatrix DXXD = arhData->Xi.at(type);
  // switch (spin) {
  if (spin == SpinFunction::Restricted) {
    Eigen::MatrixXd G_exact;
    Eigen::MatrixXd J, K;
    const auto& D = arhData->D_OAO.at(type).back().restrictedMatrix();
    // auto const& DXXD=arhData->Xi.at(type);
    // Utils::DensityMatrix DXXD = arhData->D_OAO.at(type).back();
    // DXXD.getRestrictedMatrix() = D * arhData->X.at(type).restrictedMatrix() - arhData->X.at(type).restrictedMatrix()
    // * D; std::cout << "DXXD:\n"; std::cout << std::fixed << std::setprecision(4) << DXXD.restrictedMatrix() <<
    // std::endl; std::cout << "xi:\n"; std::cout << std::fixed << std::setprecision(4) <<
    // arhData->Xi.at(type).restrictedMatrix() << std::endl;
    // AO -> OAO
    // The trick here is to implicitly transform two indices of the 2-body tensors.
    DXXD.getRestrictedMatrix() = AO2OAO * DXXD.restrictedMatrix() * AO2OAO.transpose();

    if (arhData->data->integralDirect) {
      auto ret = FockMatrix::buildJKmatricesDirect(DXXD, arhData->molecule, type);
      J = ret[0].restrictedMatrix();
      K = ret[1].restrictedMatrix();
    }
    else {
      auto ret = FockMatrix::buildJKmatrices(DXXD, arhData->molecule, arhData->data, type);
      J = ret[0].restrictedMatrix();
      K = ret[1].restrictedMatrix();
    }

    // The remaining two indices are transformed here:
    // To see where the factor `2` comes from, do the derivative of the Fock matrix with respect
    // to the density matrix that was not scaled by 2.
    G_exact = AO2OAO.transpose() * 2 * (J - 0.5 * K) * AO2OAO;
    // We reduced the number of matrix multiplications according to
    //       D * A * (1-D) - (1-D) * A * D
    //           = DA - DAD - AD + DAD
    //           = DA - AD
    Eigen::MatrixXd G_exact_ov = D * G_exact;
    Eigen::MatrixXd G_exact_vo = G_exact * D;
    HXi += (G_exact_ov - G_exact_vo);
  }
  // This is very inefficient for the unrestricted case, because the integrals are calculated twice.
  else {
    const auto& D_alpha = arhData->D_OAO.at(type).back().alphaMatrix();
    Eigen::MatrixXd J_alpha, K_alpha, J_beta, K_beta;

    // AO -> OAO
    // The trick here is to implicitly transform two indices of the 2-body tensors.
    DXXD.getAlphaMatrix() = AO2OAO * DXXD.alphaMatrix() * AO2OAO.transpose();

    if (arhData->molecule->at(type).msVector[1] > 0) {
      DXXD.getBetaMatrix() = AO2OAO * DXXD.betaMatrix() * AO2OAO.transpose();
      DXXD.getRestrictedMatrix() = DXXD.alphaMatrix() + DXXD.betaMatrix();
    }
    else {
      DXXD.getRestrictedMatrix() = DXXD.alphaMatrix();
    }

    if (arhData->data->integralDirect) {
      auto ret = FockMatrix::buildJKmatricesDirect(DXXD, arhData->molecule, type);
      J_alpha = ret[0].alphaMatrix();
      K_alpha = ret[1].alphaMatrix();
      if (arhData->molecule->at(type).msVector[1] > 0) {
        J_beta = ret[0].betaMatrix();
        K_beta = ret[1].betaMatrix();
      }
    }
    else {
      auto ret = FockMatrix::buildJKmatrices(DXXD, arhData->molecule, arhData->data, type);
      J_alpha = ret[0].alphaMatrix();
      K_alpha = ret[1].alphaMatrix();
      if (arhData->molecule->at(type).msVector[1] > 0) {
        J_beta = ret[0].betaMatrix();
        K_beta = ret[1].betaMatrix();
      }
    }

    // The remaining two indices are transformed here:
    Eigen::MatrixXd G_exact_alpha = AO2OAO.transpose() * (J_alpha - K_alpha) * AO2OAO;
    // We reduced the number of matrix multiplications according to
    //       D * A * (1-D) - (1-D) * A * D
    //           = DA - DAD - AD + DAD
    //           = DA - AD

    if (spin == SpinFunction::Alpha) {
      Eigen::MatrixXd G_exact_ov_alpha = D_alpha * G_exact_alpha;
      Eigen::MatrixXd G_exact_vo_alpha = G_exact_alpha * D_alpha;
      HXi += G_exact_ov_alpha - G_exact_vo_alpha;

      if (arhData->molecule->at(type).msVector[1] > 0) {
        Eigen::MatrixXd G_exact_alpha_beta = AO2OAO.transpose() * J_beta * AO2OAO;
        Eigen::MatrixXd G_exact_ov_alpha_beta = D_alpha * G_exact_alpha_beta;
        Eigen::MatrixXd G_exact_vo_alpha_beta = G_exact_alpha_beta * D_alpha;
        HXi += G_exact_ov_alpha_beta - G_exact_vo_alpha_beta;
      }
    }
    else if (spin == SpinFunction::Beta) {
      const auto& D_beta = arhData->D_OAO.at(type).back().betaMatrix();

      Eigen::MatrixXd G_exact_beta = AO2OAO.transpose() * (J_beta - K_beta) * AO2OAO;
      Eigen::MatrixXd G_exact_ov_beta = D_beta * G_exact_beta;
      Eigen::MatrixXd G_exact_vo_beta = G_exact_beta * D_beta;
      HXi += G_exact_ov_beta - G_exact_vo_beta;

      Eigen::MatrixXd G_exact_beta_alpha = AO2OAO.transpose() * J_alpha * AO2OAO;
      Eigen::MatrixXd G_exact_ov_beta_alpha = D_beta * G_exact_beta_alpha;
      Eigen::MatrixXd G_exact_vo_beta_alpha = G_exact_beta_alpha * D_beta;
      HXi += G_exact_ov_beta_alpha - G_exact_vo_beta_alpha;
    }
  }
}

auto ArhBuilder::addPreBO2ndOrderSpecificContribution() -> void {
  // We reduced the number of matrix multiplications according to
  //       D * A * (1-D) - (1-D) * A * D
  //           = DA - DAD - AD + DAD
  //           = DA - AD

  Eigen::MatrixXd D;
  Eigen::MatrixXd J;
  Eigen::MatrixXd G_exact;
  J.resizeLike(HXi);

  auto const& AO2OAO = arhData->data->X.at(type);

  switch (spin) {
    case SpinFunction::Restricted: {
      D = arhData->D_OAO.at(type).back().restrictedMatrix();
      break;
    }
    case SpinFunction::Alpha: {
      D = arhData->D_OAO.at(type).back().alphaMatrix();
      break;
    }
    case SpinFunction::Beta: {
      D = arhData->D_OAO.at(type).back().betaMatrix();
      break;
    }
  }

  // Important note: this works only with symmetric density matrices
  for (auto const& elem : *arhData->molecule) {
    if (elem.first != type) {
      auto otherType = elem.first;
      Utils::DensityMatrix DXXD = arhData->Xi.at(otherType);
      auto& AO2OAO_other = arhData->data->X.at(otherType);
      if (elem.second.isRestricted) {
        // AO -> OAO
        // The trick here is to implicitly transform two indices of the 2-body tensors.

        // The factor of 2 that destroys the symmetry of the Hessian is here!!!!!!!!!!!!!!!!

        DXXD.getRestrictedMatrix() = 2 * AO2OAO_other * DXXD.restrictedMatrix() * AO2OAO_other.transpose();
        std::map<Utils::ElementType, Utils::DensityMatrix> D_map;
        D_map[type] = arhData->data->D.at(type);
        D_map[otherType] = DXXD;
        if (arhData->data->integralDirect) {
          auto result = FockMatrix::buildLmatricesDirect(D_map, arhData->molecule, otherType, type);
          J = std::move(result[1]);
        }
        else {
          auto result = FockMatrix::buildLmatrices(D_map, arhData->molecule, arhData->data, otherType, type);
          J = std::move(result[1]);
        }
      }
      else {
        // AO -> OAO
        // The trick here is to implicitly transform two indices of the 2-body tensors.
        DXXD.getAlphaMatrix() = AO2OAO_other * DXXD.alphaMatrix() * AO2OAO_other.transpose();
        if (elem.second.msVector[1] > 0) {
          DXXD.getBetaMatrix() = AO2OAO_other * DXXD.betaMatrix() * AO2OAO_other.transpose();
        }
        else {
          DXXD.getBetaMatrix() = Eigen::MatrixXd::Zero(DXXD.alphaMatrix().rows(), DXXD.alphaMatrix().cols());
        }
        DXXD.getRestrictedMatrix() = DXXD.alphaMatrix() + DXXD.betaMatrix();
        std::map<Utils::ElementType, Utils::DensityMatrix> D_map;
        D_map[type] = arhData->data->D.at(type);
        D_map[otherType] = DXXD;
        if (arhData->data->integralDirect) {
          auto result = FockMatrix::buildLmatricesDirect(D_map, arhData->molecule, otherType, type);
          J = std::move(result[1]);
        }
        else {
          auto result = FockMatrix::buildLmatrices(D_map, arhData->molecule, arhData->data, otherType, type);
          J = std::move(result[1]);
        }
      }
      // The remaining two indices are transformed here:
      G_exact = AO2OAO.transpose() * J * AO2OAO;
      // We reduced the number of matrix multiplications according to
      //       D * A * (1-D) - (1-D) * A * D
      //           = DA - DAD - AD + DAD
      //           = DA - AD
      Eigen::MatrixXd G_exact_ov = D * G_exact;
      Eigen::MatrixXd G_exact_vo = G_exact * D;
      HXi += G_exact_ov - G_exact_vo;
    }
  }
}

auto ArhBuilder::makeXindependentContributions() -> void {
  // RH part:
  // F_vv_oo X + X * F_vv_oo
  makeFvvoo();
  // ARH part:
  if (!useExactHessian && arhData->size > 1) {
    makeDin();
    makeT();
    makeGinVector();
    if (!arhData->molecule->at(type).isRestricted && arhData->molecule->at(type).msVector[1] > 0) {
      makeJinVector();
    }
    if (arhData->molecule->size() > 1) {
      makeLinVector();
    }
  }
}

auto ArhBuilder::makeXdependentContributions() -> void {
  if (!useExactHessian && arhData->size > 1) {
    makeIntermediate1();
    makeIntermediate2();
  }
}

auto ArhBuilder::makeAllContributions() -> void {
  makeXindependentContributions();
  makeXdependentContributions();
}

auto ArhBuilder::makeDin() -> void {
  arhData->D_in.at(type).resize(arhData->size - 1);

  switch (spin) {
    case SpinFunction::Restricted: {
      auto& Dn = arhData->D_OAO.at(type).back().restrictedMatrix();
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_OAOi = arhData->D_OAO.at(type).at(i).restrictedMatrix();
        arhData->D_in.at(type).at(i).restrictedMatrix() = D_OAOi - Dn;
      }
      break;
    }
    case SpinFunction::Alpha: {
      auto& Dn = arhData->D_OAO.at(type).back().alphaMatrix();
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_OAOi = arhData->D_OAO.at(type).at(i).alphaMatrix();
        arhData->D_in.at(type).at(i).alphaMatrix() = D_OAOi - Dn;
      }
      break;
    }
    case SpinFunction::Beta: {
      auto& Dn = arhData->D_OAO.at(type).back().betaMatrix();
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_OAOi = arhData->D_OAO.at(type).at(i).betaMatrix();
        arhData->D_in.at(type).at(i).betaMatrix() = D_OAOi - Dn;
      }
      break;
    }
  }
}

auto ArhBuilder::makeIntermediate1() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      const auto& Xi = arhData->Xi.at(type).restrictedMatrix();
      auto& int1 = arhData->Intermediate1.at(type).restrictedVector();
      int1.resize(arhData->size - 1);
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_in = arhData->D_in.at(type).at(i).restrictedMatrix();
        int1(i) = (Xi * D_in).trace();
      }
      break;
    }
    case SpinFunction::Alpha: {
      const auto& Xi = arhData->Xi.at(type).alphaMatrix();
      auto& int1 = arhData->Intermediate1.at(type).alphaVector();
      int1.resize(arhData->size - 1);
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_in = arhData->D_in.at(type).at(i).alphaMatrix();
        int1(i) = (Xi * D_in).trace();
      }
      break;
    }
    case SpinFunction::Beta: {
      const auto& Xi = arhData->Xi.at(type).betaMatrix();
      auto& int1 = arhData->Intermediate1.at(type).betaVector();
      int1.resize(arhData->size - 1);
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_in = arhData->D_in.at(type).at(i).betaMatrix();
        int1(i) = (Xi * D_in).trace();
      }
      break;
    }
  }
}

auto ArhBuilder::makeIntermediate2() -> void {
  // I2 = T_inv * I1;
  //  --> I2 = T_inv * I1;
  //  <=> T * I2 = I1 -> T x = I1, x => I2

  switch (spin) {
    case SpinFunction::Restricted: {
      arhData->Intermediate2.at(type).restrictedVector() =
          arhData->T.at(type).restrictedMatrix().fullPivLu().solve(arhData->Intermediate1.at(type).restrictedVector());
      break;
    }
    case SpinFunction::Alpha: {
      arhData->Intermediate2.at(type).alphaVector() =
          arhData->T.at(type).alphaMatrix().fullPivLu().solve(arhData->Intermediate1.at(type).alphaVector());
      break;
    }
    case SpinFunction::Beta: {
      arhData->Intermediate2.at(type).betaVector() =
          arhData->T.at(type).betaMatrix().fullPivLu().solve(arhData->Intermediate1.at(type).betaVector());
      break;
    }
  }
}

auto ArhBuilder::makeT() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      auto& T = arhData->T.at(type).restrictedMatrix();
      T.resize(arhData->size - 1, arhData->size - 1);
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_in = arhData->D_in.at(type).at(i).restrictedMatrix();
        for (auto j = 0; j <= i; ++j) {
          auto& D_jn = arhData->D_in.at(type).at(j).restrictedMatrix();
          T(i, j) = (D_in * D_jn).trace();
          if (i != j) {
            T(j, i) = T(i, j);
          }
        }
      }
      break;
    }
    case SpinFunction::Alpha: {
      auto& T = arhData->T.at(type).alphaMatrix();
      T.resize(arhData->size - 1, arhData->size - 1);
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_in = arhData->D_in.at(type).at(i).alphaMatrix();
        for (auto j = 0; j <= i; ++j) {
          auto& D_jn = arhData->D_in.at(type).at(j).alphaMatrix();
          T(i, j) = (D_in * D_jn).trace();
          if (i != j) {
            T(j, i) = T(i, j);
          }
        }
      }
      break;
    }
    case SpinFunction::Beta: {
      auto& T = arhData->T.at(type).betaMatrix();
      T.resize(arhData->size - 1, arhData->size - 1);
      for (auto i = 0; i < arhData->size - 1; ++i) {
        auto& D_in = arhData->D_in.at(type).at(i).betaMatrix();
        for (auto j = 0; j <= i; ++j) {
          auto& D_jn = arhData->D_in.at(type).at(j).betaMatrix();
          T(i, j) = (D_in * D_jn).trace();
          if (i != j) {
            T(j, i) = T(i, j);
          }
        }
      }
      break;
    }
  }
}

auto ArhBuilder::makeFvvoo() -> void {
  // Here we can reduce number of matrix multiplications:
  //  F_vv - F_oo
  //   = (1-D) F (1-D)      - DFD
  //   = F - FD - DF  + DFD - DFD
  //   = F - FD - DF
  switch (spin) {
    case SpinFunction::Restricted: {
      const auto& F = arhData->F_OAO.at(type).restrictedMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().restrictedMatrix();
      F_vv_oo.restrictedMatrix() = F - F * P_o - P_o * F;
      break;
    }
    case SpinFunction::Alpha: {
      const auto& F = arhData->F_OAO.at(type).alphaMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().alphaMatrix();
      F_vv_oo.alphaMatrix() = F - F * P_o - P_o * F;
      break;
    }
    case SpinFunction::Beta: {
      const auto& F = arhData->F_OAO.at(type).betaMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().betaMatrix();
      F_vv_oo.betaMatrix() = F - F * P_o - P_o * F;
      break;
    }
  }
}

auto ArhBuilder::makeGinVector() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      G_in.resize(arhData->size - 1);
      G_ov_vo_in.resize(arhData->size - 1);

      const auto& Gn = arhData->G.at(type).back().restrictedMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().restrictedMatrix();
      Eigen::MatrixXd G_in_ov;
      G_in_ov.resizeLike(P_o);

      for (auto i = 0; i < arhData->size - 1; ++i) {
        const auto& Gi = arhData->G.at(type).at(i).restrictedMatrix();
        G_in[i].restrictedMatrix() = Gi - Gn;
      }

      for (auto i = 0; i < arhData->size - 1; ++i) {
        // We reduced the number of matrix multiplications according to
        //       D * A * (1-D) - (1-D) * A * D
        //           = DA - DAD - AD + DAD
        //           = DA - AD
        // Exploit symmetry of A and D: AD=(DA)^T
        G_in_ov = P_o * G_in[i].restrictedMatrix();
        G_ov_vo_in[i].restrictedMatrix() = G_in_ov - G_in_ov.transpose();
      }

      break;
    }
    case SpinFunction::Alpha: {
      G_in.resize(arhData->size - 1);
      G_ov_vo_in.resize(arhData->size - 1);

      const auto& Gn = arhData->G.at(type).back().alphaMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().alphaMatrix();
      Eigen::MatrixXd G_in_ov;
      G_in_ov.resizeLike(P_o);

      for (auto i = 0; i < arhData->size - 1; ++i) {
        const auto& Gi = arhData->G.at(type).at(i).alphaMatrix();
        G_in[i].alphaMatrix() = Gi - Gn;
      }

      for (auto i = 0; i < arhData->size - 1; ++i) {
        // We reduced the number of matrix multiplications according to
        //       D * A * (1-D) - (1-D) * A * D
        //           = DA - DAD - AD + DAD
        //           = DA - AD
        // Exploit symmetry of A and D: DA=(AD)^T
        G_in_ov = P_o * G_in[i].alphaMatrix();
        G_ov_vo_in[i].alphaMatrix() = G_in_ov - G_in_ov.transpose();
      }

      break;
    }
    case SpinFunction::Beta: {
      G_in.resize(arhData->size - 1);
      G_ov_vo_in.resize(arhData->size - 1);

      const auto& Gn = arhData->G.at(type).back().betaMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().betaMatrix();
      Eigen::MatrixXd G_in_ov;
      G_in_ov.resizeLike(P_o);

      for (auto i = 0; i < arhData->size - 1; ++i) {
        const auto& Gi = arhData->G.at(type).at(i).betaMatrix();
        G_in[i].betaMatrix() = Gi - Gn;
      }

      for (auto i = 0; i < arhData->size - 1; ++i) {
        // We reduced the number of matrix multiplications according to
        //       D * A * (1-D) - (1-D) * A * D
        //           = DA - DAD - AD + DAD
        //           = DA - AD
        G_in_ov = P_o * G_in[i].betaMatrix();
        //  Exploit symmetry of A and D: DA=(AD)^T
        G_ov_vo_in[i].betaMatrix() = G_in_ov - G_in_ov.transpose();
      }

      break;
    }
  }
}

auto ArhBuilder::makeJinVector() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      throw std::runtime_error("How did we get here? --> No alpha-beta interaction in restricted calculation!");
    }
    case SpinFunction::Alpha: {
      J_in.resize(arhData->size - 1);
      J_ov_vo_in.resize(arhData->size - 1);

      const auto& Jn = arhData->J.at(type).back().betaMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().alphaMatrix();
      Eigen::MatrixXd J_in_ov;
      J_in_ov.resizeLike(P_o);

      for (auto i = 0; i < arhData->size - 1; ++i) {
        const auto& Ji = arhData->J.at(type).at(i).betaMatrix();
        J_in[i].alphaMatrix() = Ji - Jn;
      }

      for (auto i = 0; i < arhData->size - 1; ++i) {
        // We reduced the number of matrix multiplications according to
        //       D * A * (1-D) - (1-D) * A * D
        //           = DA - DAD - AD + DAD
        //           = DA - AD
        J_in_ov = P_o * J_in[i].alphaMatrix();
        J_ov_vo_in[i].alphaMatrix() = J_in_ov - J_in_ov.transpose();
      }

      break;
    }
    case SpinFunction::Beta: {
      J_in.resize(arhData->size - 1);
      J_ov_vo_in.resize(arhData->size - 1);

      const auto& Jn = arhData->J.at(type).back().alphaMatrix();
      const auto& P_o = arhData->D_OAO.at(type).back().betaMatrix();
      Eigen::MatrixXd J_in_ov;
      J_in_ov.resizeLike(P_o);

      for (auto i = 0; i < arhData->size - 1; ++i) {
        const auto& Ji = arhData->J.at(type).at(i).alphaMatrix();
        J_in[i].betaMatrix() = Ji - Jn;
      }

      for (auto i = 0; i < arhData->size - 1; ++i) {
        // We reduced the number of matrix multiplications according to
        //       D * A * (1-D) - (1-D) * A * D
        //           = DA - DAD - AD + DAD
        //           = DA - AD
        J_in_ov = P_o * J_in[i].betaMatrix();
        J_ov_vo_in[i].betaMatrix() = J_in_ov - J_in_ov.transpose();
      }

      break;
    }
  }
}

auto ArhBuilder::makeLinVector() -> void {
  // We reduced the number of matrix multiplications according to
  //       D * A * (1-D) - (1-D) * A * D
  //           = DA - DAD - AD + DAD
  //           = DA - AD

  Eigen::MatrixXd P_o;

  switch (spin) {
    case SpinFunction::Restricted: {
      P_o = arhData->D_OAO.at(type).back().restrictedMatrix();
      break;
    }
    case SpinFunction::Alpha: {
      P_o = arhData->D_OAO.at(type).back().alphaMatrix();
      break;
    }
    case SpinFunction::Beta: {
      P_o = arhData->D_OAO.at(type).back().betaMatrix();
      break;
    }
  }

  for (auto const& elem : *arhData->molecule) {
    if (elem.first != type) {
      auto otherType = elem.first;
      L_in.at(otherType).resize(arhData->size - 1);
      L_ov_vo_in.at(otherType).resize(arhData->size - 1);
      auto& L_in_other = L_in.at(otherType);
      auto& L_ov_vo_in_other = L_ov_vo_in.at(otherType);
      Eigen::MatrixXd L_in_other_ov;
      L_in_other_ov.resizeLike(P_o);

      if (arhData->molecule->at(otherType).isRestricted) {
        auto& Ln = arhData->L.at(type).back().at(otherType);
        for (auto i = 0; i < arhData->size - 1; ++i) {
          auto& Li = arhData->L.at(type).at(i).at(otherType);
          L_in_other[i].restrictedMatrix() = Li - Ln;
        }
        for (auto i = 0; i < arhData->size - 1; ++i) {
          L_in_other_ov = P_o * L_in_other[i].restrictedMatrix();
          L_ov_vo_in_other[i].restrictedMatrix() = L_in_other_ov - L_in_other_ov.transpose();
        }
      }
      else {
        // alpha
        {
          auto& Ln = arhData->L.at(type).back().at(otherType);
          for (auto i = 0; i < arhData->size - 1; ++i) {
            auto& Li = arhData->L.at(type).at(i).at(otherType);
            L_in_other[i].alphaMatrix() = Li - Ln;
          }
          for (auto i = 0; i < arhData->size - 1; ++i) {
            L_in_other_ov = P_o * L_in_other[i].alphaMatrix();
            L_ov_vo_in_other[i].alphaMatrix() = L_in_other_ov - L_in_other_ov.transpose();
          }
        }
        // beta
        if (elem.second.msVector[1] > 0) {
          auto& Ln = arhData->L.at(type).back().at(otherType);
          for (auto i = 0; i < arhData->size - 1; ++i) {
            auto& Li = arhData->L.at(type).at(i).at(otherType);
            L_in_other[i].betaMatrix() = Li - Ln;
          }
          for (auto i = 0; i < arhData->size - 1; ++i) {
            L_in_other_ov = P_o * L_in_other[i].betaMatrix();
            L_ov_vo_in_other[i].betaMatrix() = L_in_other_ov - L_in_other_ov.transpose();
          }
        }
      }
    }
  }
}

auto ArhBuilder::addBoArhSpecificContribution() -> void {
  switch (spin) {
    case SpinFunction::Restricted:
      for (auto i = 0; i < arhData->size - 1; ++i) {
        HXi += G_ov_vo_in.at(i).restrictedMatrix() * arhData->Intermediate2.at(type).restrictedVector()(i);
      }
      break;
    case SpinFunction::Alpha:
      for (auto i = 0; i < arhData->size - 1; ++i) {
        HXi += G_ov_vo_in.at(i).alphaMatrix() * arhData->Intermediate2.at(type).alphaVector()(i);
      }
      break;
    case SpinFunction::Beta:
      for (auto i = 0; i < arhData->size - 1; ++i) {
        HXi += G_ov_vo_in.at(i).betaMatrix() * arhData->Intermediate2.at(type).betaVector()(i);
      }
      break;
  }
}

auto ArhBuilder::addUnrestrictedArhSpecificContribution() -> void {
  switch (spin) {
    case SpinFunction::Restricted: {
      throw std::runtime_error("How did we get here? --> No unrestricted contribution in restricted calculation!");
    }
    case SpinFunction::Alpha: {
      for (auto i = 0; i < arhData->size - 1; ++i) {
        HXi += J_ov_vo_in.at(i).alphaMatrix() * arhData->Intermediate2.at(type).betaVector()(i);
      }
      break;
    }
    case SpinFunction::Beta: {
      for (auto i = 0; i < arhData->size - 1; ++i) {
        HXi += J_ov_vo_in.at(i).betaMatrix() * arhData->Intermediate2.at(type).alphaVector()(i);
      }
      break;
    }
  }
}

auto ArhBuilder::addPreBoArhSpecificContribution() -> void {
  for (auto const& elem : *arhData->molecule) {
    const auto& otherType = elem.first;
    if (otherType == type) {
      continue;
    }
    if (elem.second.isRestricted) {
      for (auto i = 0; i < arhData->size - 1; ++i) {
        HXi += L_ov_vo_in.at(otherType).at(i).restrictedMatrix() * arhData->Intermediate2.at(otherType).restrictedVector()(i);
      }
    }
    else {
      for (auto i = 0; i < arhData->size - 1; ++i) {
        HXi += L_ov_vo_in.at(otherType).at(i).alphaMatrix() * arhData->Intermediate2.at(otherType).alphaVector()(i);
      }
      if (elem.second.msVector[1] > 1) {
        for (auto i = 0; i < arhData->size - 1; ++i) {
          HXi += L_ov_vo_in.at(otherType).at(i).betaMatrix() * arhData->Intermediate2.at(otherType).betaVector()(i);
        }
      }
    }
  }
}

auto ArhBuilder::evaluate() -> Eigen::VectorXd {
  if (arhData->size == 0) {
    throw std::runtime_error("TRAH evaluator was called, but the size is Zero.");
  }

  addRhContribution();

  if (useExactHessian) {
    if (arhData->molecule->at(type).N > 1) {
      addBo2ndOrderSpecificContribution();
    }
    if (arhData->molecule->size() > 1) {
      addPreBO2ndOrderSpecificContribution();
    }
  }
  else if (!useExactHessian && arhData->size > 1) {
    // Only required if there is more than one particle of the given type.
    if (arhData->molecule->at(type).N > 1) {
      addBoArhSpecificContribution();
      if (!arhData->molecule->at(type).isRestricted && arhData->molecule->at(type).msVector[1] > 0) {
        addUnrestrictedArhSpecificContribution();
      }
    }
    if (arhData->molecule->size() > 1) {
      addPreBoArhSpecificContribution();
    }
  }

  if (addLevelShift_) {
    switch (spin) {
      case SpinFunction::Restricted:
        HXi -= mu * arhData->Xi.at(type).restrictedMatrix();
        break;
      case SpinFunction::Alpha:
        HXi -= mu * arhData->Xi.at(type).alphaMatrix();
        break;
      case SpinFunction::Beta:
        HXi -= mu * arhData->Xi.at(type).betaMatrix();
        break;
    }
  }

  // Virt-Occ
  HXi = Theta_virt.transpose() * HXi * Theta_occ;

  // This guarantees the symmetry of the Hessian
  if (spin == SpinFunction::Restricted) {
    HXi *= 2;
  }

  return Eigen::Map<Eigen::VectorXd>(HXi.data(), virt * occ, 1);
}

const Eigen::MatrixXd& ArhBuilder::getGradient() const {
  return G;
}

Eigen::VectorXd ArhBuilder::getGradientVector() {
  return Eigen::Map<Eigen::VectorXd>(G.data(), virt * occ, 1);
}

auto ArhBuilder::getTransformationMatrix() -> Eigen::MatrixXd {
  // MO
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(occ + virt, occ + virt);
  Eigen::MatrixXd mat = Id - kappa_;
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

  return Optimization::Davidson::modifiedGramSchmidt(Theta * Q * Theta.transpose());
}

void ArhBuilder::setMu(double m) {
  mu = m;
}

void ArhBuilder::addLevelShift(bool addLs) {
  addLevelShift_ = addLs;
}

auto ArhBuilder::getRhDiagonal(bool print) const -> Eigen::DiagonalMatrix<double, Eigen::Dynamic> {
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> ret(virt * occ);
  Eigen::MatrixXd diag(virt, occ);

  diag.setZero();
  Eigen::MatrixXd F_MO;

  switch (spin) {
    case SpinFunction::Restricted:
      F_MO = 2 * Theta.transpose() * F_vv_oo.restrictedMatrix() * Theta;
      break;
    case SpinFunction::Alpha:
      F_MO = Theta.transpose() * F_vv_oo.alphaMatrix() * Theta;
      break;
    case SpinFunction::Beta:
      F_MO = Theta.transpose() * F_vv_oo.betaMatrix() * Theta;
      break;
  }

  for (int a = occ; a < occ + virt; ++a) {
    for (int i = 0; i < occ; ++i) {
      diag(a - occ, i) = F_MO(a, a) + F_MO(i, i);
    }
  }

  ret.diagonal() = Eigen::Map<Eigen::VectorXd>(diag.data(), virt * occ);

  int negativeThetaCounter = 0;
  for (auto it = 0; it < ret.diagonal().size(); ++it) {
    if (ret.diagonal()(it) < 0) {
      ++negativeThetaCounter;
    }
  }
  if (print && negativeThetaCounter > 0) {
    std::cout << "WARNING: " << negativeThetaCounter << " diagonal Hessian elements are negative!" << std::endl;
  }
  return ret;
}

auto ArhBuilder::getThirdTrialVector() const -> Eigen::VectorXd {
  Eigen::VectorXd diag = getRhDiagonal(false).diagonal();
  Eigen::MatrixXd::Index indexMin = 0;
  diag.minCoeff(&indexMin);

  Eigen::VectorXd b3(virt * occ);
  b3.setZero();
  b3(indexMin) = diag(indexMin);
  return b3;
}

auto ArhBuilder::updateX(Eigen::VectorXd& vectorX) -> void {
  // MO 2 OAO transformation.
  Eigen::Ref<Eigen::MatrixXd> kappa_MO = Eigen::Map<Eigen::MatrixXd>(vectorX.data(), virt, occ);
  Eigen::MatrixXd xi(dim, dim);
  xi.setZero();
  xi.block(occ, 0, virt, occ) = kappa_MO;
  xi.block(0, occ, occ, virt) = kappa_MO.transpose();
  Eigen::MatrixXd X(dim, dim);
  X.setZero();
  X.block(occ, 0, virt, occ) = kappa_MO;
  X.block(0, occ, occ, virt) = -kappa_MO.transpose();

  kappa_ = X;

  switch (spin) {
    case SpinFunction::Restricted:
      arhData->Xi.at(type).getRestrictedMatrix() = -Theta * xi * Theta.transpose();
      // Related to symmetry of the Hessian
      // kappa_ *= 0.5;
      break;
    case SpinFunction::Alpha:
      arhData->Xi.at(type).getAlphaMatrix() = -Theta * xi * Theta.transpose();
      break;
    case SpinFunction::Beta:
      arhData->Xi.at(type).getBetaMatrix() = -Theta * xi * Theta.transpose();
      break;
  }
}

void ArhBuilder::setUseExactHessian(bool useExactHessian) {
  ArhBuilder::useExactHessian = useExactHessian;
}

} // namespace Kiwi
} // namespace Scine
