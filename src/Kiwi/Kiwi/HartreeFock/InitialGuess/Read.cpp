/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/HartreeFock/InitialGuess/Read.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/MolecularOrbitalsIO.h>
#include <Utils/Geometry/ElementInfo.h>
#include <cstdio>

namespace Scine {
namespace Kiwi {

auto Read::writeOrbitals() -> void {
  for (auto const& elem : *molecule_) {
    std::string ofname = "";
    ofname += molecule_->getInputFileName();
    ofname += ".";
    ofname += Utils::ElementInfo::symbol(elem.first);
    ofname += extension_;

    Kiwi::MolecularOrbitalsIO::write(ofname, data_->C_OAO.at(elem.first));
  }
}

auto Read::readOrbitals() -> void {
  for (auto const& elem : *molecule_) {
    std::string ifname = "";
    ifname += molecule_->getInputFileName();
    ifname += ".";
    ifname += Utils::ElementInfo::symbol(elem.first);
    ifname += extension_;

    data_->C_OAO[elem.first] = Kiwi::MolecularOrbitalsIO::read(ifname);

    if (elem.second.isRestricted != data_->C_OAO.at(elem.first).isRestricted()) {
      throw std::runtime_error(
          "Mismatch between spin-symmetry for particle type: " + Utils::ElementInfo::symbol(elem.first) + "!");
    }

    bool mismatch = false;

    if (elem.second.isRestricted) {
      if ((data_->C_OAO.at(elem.first).restrictedMatrix().rows() != int(elem.second.LAO)) ||
          (data_->C_OAO.at(elem.first).restrictedMatrix().cols() != int(elem.second.LMO))) {
        mismatch = true;
      }
    }
    else {
      if ((data_->C_OAO.at(elem.first).alphaMatrix().rows() != int(elem.second.LAO)) ||
          (data_->C_OAO.at(elem.first).alphaMatrix().cols() != int(elem.second.LMO))) {
        mismatch = true;
      }
      if ((data_->C_OAO.at(elem.first).betaMatrix().rows() != int(elem.second.LAO)) ||
          (data_->C_OAO.at(elem.first).betaMatrix().cols() != int(elem.second.LMO))) {
        mismatch = true;
      }
    }
    if (mismatch) {
      throw std::runtime_error("Mismatch between orbital number of type: " + Utils::ElementInfo::symbol(elem.first) + "!");
    }
    auto const& X = data_->X[elem.first];
    if (elem.second.isRestricted) {
      Eigen::MatrixXd tmp = X * data_->C_OAO.at(elem.first).restrictedMatrix();
      Utils::MolecularOrbitals C;
      C = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp));
      data_->C[elem.first] = C;
      data_->D[elem.first] = Kiwi::HartreeFockUtils::makeDensity(elem.first, molecule_, C);
      data_->D_OAO[elem.first] = Kiwi::HartreeFockUtils::makeDensity(elem.first, molecule_, data_->C_OAO.at(elem.first));
    }
    else {
      const auto& tmp_oao_alpha = data_->C_OAO.at(elem.first).alphaMatrix();

      Eigen::MatrixXd tmp_alpha = X * tmp_oao_alpha;
      Eigen::MatrixXd tmp_beta;
      tmp_beta.resizeLike(tmp_alpha);

      if (molecule_->at(elem.first).msVector[1] > 0) {
        const auto& tmp_oao_beta = data_->C_OAO.at(elem.first).betaMatrix();
        tmp_beta = X * tmp_oao_beta;
      }

      Utils::MolecularOrbitals C;
      C = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Eigen::MatrixXd>(std::move(tmp_alpha),
                                                                                        std::move(tmp_beta));
      data_->C[elem.first] = C;
      data_->D[elem.first] = Kiwi::HartreeFockUtils::makeDensity(elem.first, molecule_, C);
      data_->D_OAO[elem.first] = Kiwi::HartreeFockUtils::makeDensity(elem.first, molecule_, data_->C_OAO.at(elem.first));
    }
  }
}

auto Read::cleanUp() -> void {
  for (auto const& elem : *molecule_) {
    std::string ofname = "";
    ofname += molecule_->getInputFileName();
    ofname += ".";
    ofname += Utils::ElementInfo::symbol(elem.first);
    ofname += extension_;
    const char* c = ofname.c_str();
    std::remove(c);
  }
}

} // namespace Kiwi
} // namespace Scine
