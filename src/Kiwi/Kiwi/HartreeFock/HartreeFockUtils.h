/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_HARTREEFOCKUTILS_H
#define KIWI_HARTREEFOCKUTILS_H

#include <Kiwi/HartreeFock/FockMatrixBuilder.h>
#include <Kiwi/KiwiUtils/Molecule.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Scf/LcaoUtils/DensityMatrixBuilder.h>
#include <iostream>

namespace Scine {
namespace Kiwi {

class HartreeFockUtils {
 public:
  struct SpinAdaptedDouble {
    double alpha;
    double beta;
    double restricted;
  };

  HartreeFockUtils() = default;
  ~HartreeFockUtils() = default;

  inline static auto makeDensity(Utils::ElementType type, const std::shared_ptr<Molecule>& molecule,
                                 const Utils::MolecularOrbitals& C) -> Utils::DensityMatrix {
    auto DBuilder = Utils::LcaoUtils::DensityMatrixBuilder(C);
    Utils::DensityMatrix D;
    if (molecule->at(type).isRestricted) {
      D = DBuilder.generateRestrictedForNumberElectrons(int(molecule->at(type).N));
    }
    else {
      D = DBuilder.generateUnrestrictedForNumberAlphaAndBetaElectrons(molecule->at(type).msVector[0],
                                                                      molecule->at(type).msVector[1]);
    }
    return D;
  }

  inline static auto getPointChargeRepulsion(const Utils::AtomCollection& atoms) -> double {
    double ret = 0;
    // prepare charges
    std::vector<double> charges(atoms.size());
    for (auto i = 0; i < atoms.size(); ++i) {
      charges[i] = Utils::ElementInfo::Z(atoms[i].getElementType());
    }
    for (auto i = 0; i < atoms.size(); ++i) {
      for (auto j = 0; j < i; ++j) {
        auto R = atoms[i].getPosition() - atoms[j].getPosition();
        ret += charges[i] * charges[j] / std::sqrt(R.dot(R));
      }
    }

    return ret;
  }

  inline static auto getEnergy(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data) -> double {
    double energy = 0;

    if (molecule->hasPointCharges()) {
      energy += getPointChargeRepulsion(molecule->getPointCharges());
    }
    for (auto const& elem : *molecule) {
      auto const& type = elem.first;
      if (elem.second.isRestricted) {
        auto const& D = data->D[type].restrictedMatrix();
        auto const& H = data->H[type];
        auto const& G = data->hartreeFockData->G[type].restrictedMatrix();
        auto const& I = data->hartreeFockData->I[type].restrictedMatrix();
        energy += (D * H).trace();
        energy += 0.5 * (D * G).trace();
        energy += 0.5 * D.cwiseProduct(I).sum();
      }
      else {
        auto const& D = data->D[type];
        auto const& H = data->H[type];
        auto const& G = data->hartreeFockData->G[type];
        auto const& I = data->hartreeFockData->I[type];
        energy += (D.restrictedMatrix() * H).trace();
        if (elem.second.msVector[0] > 0) {
          energy += 0.5 * (G.alphaMatrix() * D.alphaMatrix()).trace();
          energy += 0.5 * (I.restrictedMatrix() * D.alphaMatrix()).trace();
        }
        if (elem.second.msVector[1] > 0) {
          energy += 0.5 * (G.betaMatrix() * D.betaMatrix()).trace();
          energy += 0.5 * (I.restrictedMatrix() * D.betaMatrix()).trace();
        }
      }
    }
    return energy;
  }

  inline static auto getEnergy(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data,
                               const std::map<Utils::ElementType, Utils::DensityMatrix>& Dmap,
                               const std::shared_ptr<HartreeFockData>& hartreeFockData) -> double {
    double energy = 0;

    if (molecule->hasPointCharges()) {
      energy += getPointChargeRepulsion(molecule->getPointCharges());
    }
    for (auto const& elem : *molecule) {
      auto const& type = elem.first;
      if (elem.second.isRestricted) {
        // auto const& D = Dmap.at(type).restrictedMatrix();
        // auto const& H = data->H[type];
        // auto const& F = hartreeFockData->F[type].restrictedMatrix();
        // auto const& I = hartreeFockData->I[type].restrictedMatrix();
        // energy += 0.5 * (H * D).trace();
        // energy += 0.5 * (F * D).trace();
        // energy += 0.5 * D.cwiseProduct(I).sum();
        auto const& D = Dmap.at(type).restrictedMatrix();
        auto const& H = data->H[type];
        auto const& G = hartreeFockData->G[type].restrictedMatrix();
        auto const& I = hartreeFockData->I[type].restrictedMatrix();
        energy += (H * D).trace();
        energy += 0.5 * (G * D).trace();
        energy += D.cwiseProduct(I).sum();
      }
      else {
        auto const& D = Dmap.at(type);
        auto const& H = data->H[type];
        auto const& G = hartreeFockData->G[type];
        auto const& I = hartreeFockData->I[type];
        energy += (H * D.restrictedMatrix()).trace();
        if (elem.second.msVector[0] > 0) {
          energy += 0.5 * (G.alphaMatrix() * D.alphaMatrix()).trace();
          // energy += 0.5 * G.alphaMatrix().cwiseProduct(D.alphaMatrix()).sum();
          energy += 0.5 * (I.restrictedMatrix() * D.alphaMatrix()).trace();
          // energy += 0.5 * I.restrictedMatrix().cwiseProduct(D.alphaMatrix()).sum();
        }
        if (elem.second.msVector[1] > 0) {
          energy += 0.5 * (G.betaMatrix() * D.betaMatrix()).trace();
          // energy += 0.5 * G.betaMatrix().cwiseProduct(D.betaMatrix()).sum();
          energy += 0.5 * (I.restrictedMatrix() * D.betaMatrix()).trace();
          // energy += 0.5 * I.restrictedMatrix().cwiseProduct(D.betaMatrix()).sum();
        }
      }
    }

    return energy;
  }

  inline static auto getEnergyGeneral(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data,
                                      const std::map<Utils::ElementType, Utils::DensityMatrix>& Dmap,
                                      const std::shared_ptr<HartreeFockData>& hartreeFockData) -> double {
    double energy = 0;

    if (molecule->hasPointCharges()) {
      energy += getPointChargeRepulsion(molecule->getPointCharges());
    }
    for (auto const& elem : *molecule) {
      auto L = int(elem.second.LAO);
      auto const& type = elem.first;
      auto const& D = Dmap.at(type).restrictedMatrix();
      auto const& h = data->H[type];
      Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2 * L, 2 * L);
      H.block(0, 0, L, L) = h;
      H.block(L, L, L, L) = h;
      auto const& F = hartreeFockData->F[type].restrictedMatrix();
      //      auto const& I = hartreeFockData->I[type].restrictedMatrix();
      energy += 0.5 * ((H + F) * D).trace();
      //      energy += 0.5 * D.cwiseProduct(I).sum();
    }

    return energy;
  }

  inline static auto getDensityGradientError(const Eigen::MatrixXd& X, const Eigen::MatrixXd& S,
                                             const Eigen::MatrixXd& D, const Eigen::MatrixXd& F) -> double {
    auto e = X.transpose() * (F * D * S - S * D * F) * X;

    return std::sqrt((e * e.transpose()).trace());
  }

  inline static auto getDensityGradientError(const Eigen::MatrixXd& D, const Eigen::MatrixXd& F) -> double {
    auto e = F * D - D * F;

    return std::sqrt((e * e.transpose()).trace());
  }

  inline static auto getDensityGradientErrorMap(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data)
      -> std::map<Utils::ElementType, SpinAdaptedDouble> {
    auto ret = std::map<Utils::ElementType, SpinAdaptedDouble>();

    for (auto const& elem : *molecule) {
      SpinAdaptedDouble res;
      auto const& type = elem.first;
      auto const& X = data->X[type];
      auto const& D = data->D[type];
      auto const& S = data->S[type];
      auto const& F = data->hartreeFockData->F[type];
      if (elem.second.isRestricted) {
        res.restricted = getDensityGradientError(X, S, D.restrictedMatrix(), F.restrictedMatrix());
      }
      else {
        res.alpha = getDensityGradientError(X, S, D.alphaMatrix(), F.alphaMatrix());
        res.beta = getDensityGradientError(X, S, D.betaMatrix(), F.betaMatrix());
      }
      ret[type] = res;
    }
    return ret;
  }

  inline static auto getDensityGradientErrorMapOrthoTest(const std::shared_ptr<Molecule>& molecule,
                                                         const std::shared_ptr<Data>& data)
      -> std::map<Utils::ElementType, SpinAdaptedDouble> {
    auto ret = std::map<Utils::ElementType, SpinAdaptedDouble>();

    for (auto const& elem : *molecule) {
      SpinAdaptedDouble res;
      auto const& type = elem.first;
      auto const& D = data->D_OAO[type];
      Eigen::MatrixXd F(data->hartreeFockData->F_OAO[type].restrictedMatrix().rows(),
                        data->hartreeFockData->F_OAO[type].restrictedMatrix().cols());
      F.setZero();
      auto L = int(molecule->at(type).LMO);
      auto n = molecule->at(type).msVector[0] + molecule->at(type).msVector[1];
      F.block(0, n, n, L - n) += data->hartreeFockData->F_OAO[type].restrictedMatrix().block(0, n, n, L - n);
      F.block(n, 0, L - n, n) += data->hartreeFockData->F_OAO[type].restrictedMatrix().block(n, 0, L - n, n);
      res.restricted = getDensityGradientError(D.restrictedMatrix(), F);
      ret[type] = res;
    }
    return ret;
  }

  inline static auto getDensityGradientErrorMapOrtho(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data)
      -> std::map<Utils::ElementType, SpinAdaptedDouble> {
    auto ret = std::map<Utils::ElementType, SpinAdaptedDouble>();

    for (auto const& elem : *molecule) {
      SpinAdaptedDouble res;
      auto const& type = elem.first;
      auto const& D = data->D_OAO[type];
      auto const& F = data->hartreeFockData->F_OAO[type];
      if (elem.second.isRestricted) {
        res.restricted = getDensityGradientError(D.restrictedMatrix(), F.restrictedMatrix());
      }
      else {
        res.alpha = getDensityGradientError(D.alphaMatrix(), F.alphaMatrix());
        res.beta = getDensityGradientError(D.betaMatrix(), F.betaMatrix());
      }
      ret[type] = res;
    }
    return ret;
  }

  inline static auto S2(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data,
                        const Utils::ElementType type) -> std::pair<double, double> {
    if (molecule->at(type).isRestricted) {
      throw std::runtime_error("Tried to calculate S^2 for restricted particle type");
    }

    auto nAlpha = molecule->at(type).msVector[0];
    auto nBeta = molecule->at(type).msVector[1];

    auto const& C = data->C[type];
    auto const& overlap = data->S[type];

    auto Sab = C.alphaMatrix().block(0, 0, overlap.rows(), nAlpha).transpose() * overlap *
               C.betaMatrix().block(0, 0, overlap.rows(), nBeta);

    double exact = (nAlpha - nBeta) / 2. * ((nAlpha - nBeta) / 2. + 1);

    return {exact + nBeta - Sab.array().square().sum(), exact};
  }

  inline static auto S2(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data,
                        const Utils::ElementType type, const Eigen::MatrixXd& overlap) -> std::pair<double, double> {
    if (molecule->at(type).isRestricted) {
      throw std::runtime_error("Tried to calculate S^2 for restricted particle type");
    }

    auto nAlpha = molecule->at(type).msVector[0];
    auto nBeta = molecule->at(type).msVector[1];

    auto const& C = data->C[type];

    auto Sab = C.alphaMatrix().block(0, 0, overlap.rows(), nAlpha).transpose() * overlap *
               C.betaMatrix().block(0, 0, overlap.rows(), nBeta);

    double exact = (nAlpha - nBeta) / 2. * ((nAlpha - nBeta) / 2. + 1);

    return {exact + nBeta - Sab.array().square().sum(), exact};
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_HARTREEFOCKUTILS_H
