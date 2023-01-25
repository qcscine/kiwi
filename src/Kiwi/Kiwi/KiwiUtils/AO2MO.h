/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_AO2MO_H
#define KIWI_AO2MO_H

#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/KiwiUtils/AO2MOHelper.h>
#include <Kiwi/KiwiUtils/BoFciDumper.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <fstream>
#include <utility>

namespace Scine {
namespace Kiwi {

enum class FciDumpFormat { PreBO, Standard };

/**
 * @class AO2MO @file AO2MO
 * @brief Class that performs the AO to MO transformation. Note that this class relies on integrals that are stored
 * in memory and is, therefore, only suitable for small basis sets up to ~ 100 basis functions. The class can also
 * write the integrals into an FCIDUMP file that is formatted according to the following layout:
 *
 * IMPORTANT: Note that the integrals are handled internally in the chemist's notation, but the FCIDUMP is printed in
 * the physicists' notation!
 *
 * Particle type - Basis Function ... integral
 * One body terms
 * i-j     i-l     double                    --> Same particle type, same spin,
 *                                               but different basis
 *                                               functions are possible.
 * Two body terms
 * i-j     m-n     m-p     i-l     double    --> The ordering must follow
 *                                               the rule: a(i)+ a(m)+ a(m)- a(i)-
 *                                               The matrix element would be
 *                                               zero otherwise.
 *
 */
class AO2MO {
 private:
  std::shared_ptr<Data> data_;
  std::shared_ptr<Molecule> molecule_;
  bool _print;
  FciDumpFormat format_;

  /* Two-body super-matrix */
  std::unique_ptr<std::map<ElementPair, Eigen::MatrixXd>> _twoBody;

 public:
  const std::unique_ptr<std::map<ElementPair, Eigen::MatrixXd>>& getTwoBody() const {
    return _twoBody;
  }
  const std::unique_ptr<std::map<Utils::ElementType, Eigen::MatrixXd>>& getOneBody() const {
    return _oneBody;
  }

 private:
  std::unique_ptr<std::map<Utils::ElementType, Eigen::MatrixXd>> _oneBody;

 public:
  AO2MO(std::shared_ptr<Data> data, bool print = true)
    : data_(std::move(data)), molecule_(data_->molecule), _print(print) {
    if (!molecule_->isRestricted()) {
      std::cout << "Unrestricted calculation: only alpha orbitals will be written in FCIDUMP." << std::endl;
    }
    if (!molecule_->useHighSpinApprox()) {
      throw std::runtime_error("AO2MO is only available with the high spin approximation.");
    }
    if (data_->integralDirect) {
      std::cout << "\n\nAO2MO was called with integral-direct mode.\nIntegrals will be calculated and stored in memory "
                   "now.\n\n";
      data_->integralDirect = false;
      data_->twoBodyIntegrals(true);
      data_->makeExchange(true);
    }
    if (molecule_->size() > 1) {
      format_ = FciDumpFormat::PreBO;
    }
    else {
      format_ = FciDumpFormat::Standard;
    }

    // This is for safety reasons:
    for (auto& elem : *molecule_) {
      if (elem.second.LMO == 0) {
        elem.second.LMO = data->X.at(elem.first).cols();
      }
    }
  }

  auto perform() -> void {
    _oneBody = std::make_unique<std::map<Utils::ElementType, Eigen::MatrixXd>>();
    _twoBody = std::make_unique<std::map<ElementPair, Eigen::MatrixXd>>();

    if (_print) {
      auto& clock = Clock::getInstance();
      std::cout << std::endl << std::endl;
      std::cout << "Contract one-body tensors";
      std::cout << std::string(12, ' ');
      clock.time("1body");
    }

    for (auto const& elem : *molecule_) {
      _oneBody->insert(std::pair<Utils::ElementType, Eigen::MatrixXd>(
          elem.first, Eigen::MatrixXd::Zero(elem.second.LMO, elem.second.LMO)));
      if (elem.second.isRestricted) {
        _oneBody->at(elem.first) = data_->C[elem.first].restrictedMatrix().transpose() * data_->H[elem.first] *
                                   data_->C[elem.first].restrictedMatrix();
      }
      else {
        _oneBody->at(elem.first) =
            data_->C[elem.first].alphaMatrix().transpose() * data_->H[elem.first] * data_->C[elem.first].alphaMatrix();
      }
    }
    if (_print) {
      auto& clock = Clock::getInstance();
      clock.time("1body");
      std::cout << std::endl;
    }

    if (_print) {
      auto& clock = Clock::getInstance();
      std::cout << "Contract two-body tensors";
      std::cout << std::string(12, ' ');
      clock.time("2body");
    }
    for (auto const& pair : data_->uniquePairs) {
      auto tp1 = molecule_->at(pair.first);
      auto tp2 = molecule_->at(pair.second);
      auto& AO = data_->Coulomb.at(pair);
      Eigen::MatrixXd C1;
      Eigen::MatrixXd C2;
      // type 1
      if (tp1.isRestricted) {
        C1 = data_->C.at(pair.first).restrictedMatrix();
      }
      else {
        C1 = data_->C.at(pair.first).alphaMatrix();
      }
      // type 2
      if (tp2.isRestricted) {
        C2 = data_->C.at(pair.second).restrictedMatrix();
      }
      else {
        C2 = data_->C.at(pair.second).alphaMatrix();
      }
      AO2MOHelper helper(AO, C1, C2);
      helper.perform();
      _twoBody->insert(std::pair<ElementPair, Eigen::MatrixXd>(pair, Eigen::MatrixXd()));
      _twoBody->at(pair) = helper.getResult();
    }
    if (_print) {
      auto& clock = Clock::getInstance();
      clock.time("2body");
      std::cout << std::endl;
    }
  }

  auto writePreBOFormat(std::string fileName, double integralThresh) const -> void {
    if (_print) {
      std::cout << "\n\n- writing the orbital integrals on disk..." << std::endl;
      auto& clock = Clock::getInstance();
      clock.time("fcidump");
    }
    std::ofstream outf(fileName.c_str(), std::ios::out);
    outf << std::setprecision(16);

    if (molecule_->hasPointCharges()) {
      outf << std::setw(38) << HartreeFockUtils::getPointChargeRepulsion(molecule_->getPointCharges()) << "\n";
    }

    for (auto const& elem : *molecule_) {
      for (auto i = 0UL; i < elem.second.LMO; ++i) {
        for (auto j = 0UL; j < elem.second.LMO; ++j) {
          auto tij = _oneBody->at(elem.first)(i, j);
          if (fabs(tij) > integralThresh) {
            outf << std::setw(10) << elem.second.index << "-" << i << std::setw(10) << elem.second.index << "-" << j
                 << std::setw(38) << tij << "\n";
          }
        }
      }
    }

    for (auto const& pair : data_->uniquePairs) {
      const auto& elem1 = molecule_->at(pair.first);
      const auto& elem2 = molecule_->at(pair.second);
      auto const& V = _twoBody->at(pair);
      for (auto i = 0UL; i < elem1.LMO; ++i) {
        auto index1 = i * elem1.LMO;
        for (auto k = 0UL; k < elem2.LMO; ++k) {
          auto index2 = k * elem2.LMO;
          for (auto l = 0UL; l < elem2.LMO; ++l) {
            for (auto j = 0UL; j < elem1.LMO; ++j) {
              auto Vijkl = V(index1 + j, index2 + l);
              if (elem1.index == elem2.index) {
                Vijkl *= 0.5;
              }
              // Chemist's to physicist's notation:
              // 1 1 2 2
              // i j k l
              //
              // 1 2 1 2
              // i k j l
              //
              // Permute positions 3 and 4, for 2nd quantization operator string:
              // 1 2 2 1
              // i k l j
              //
              if (fabs(Vijkl) > integralThresh) {
                outf << std::setw(10) << elem1.index << "-" << i << std::setw(10) << elem2.index << "-" << k
                     << std::setw(10) << elem2.index << "-" << l << std::setw(10) << elem1.index << "-" << j
                     << std::setw(40) << Vijkl << "\n";
              }
            }
          }
        }
      }
    }

    outf << std::endl;
    outf.close();
    if (_print) {
      std::cout << "    --> integrals are written successfully ";
      auto& clock = Clock::getInstance();
      std::cout << std::string(12, ' ');
      clock.time("fcidump");
      std::cout << std::endl;
    }
  }

  auto write(double integralThresh) const -> void {
    std::string fileName = "FCIDUMP." + molecule_->getInputFileName();

    switch (format_) {
      case FciDumpFormat::PreBO: {
        writePreBOFormat(fileName, integralThresh);
        break;
      }
      case FciDumpFormat::Standard: {
        BoFciDumpData boFciDumpData;
        boFciDumpData.coreEnergy = HartreeFockUtils::getPointChargeRepulsion(molecule_->getPointCharges());
        boFciDumpData.nElectrons = molecule_->at(Utils::ElementType::E).N;
        boFciDumpData.nOrbitals = molecule_->at(Utils::ElementType::E).LMO;
        boFciDumpData.spinPolarization =
            molecule_->at(Utils::ElementType::E).msVector[0] - molecule_->at(Utils::ElementType::E).msVector[1];
        boFciDumpData.orbitalSymmetry = std::vector<int>(boFciDumpData.nOrbitals, 1);
        BoFciDumper::integralThreshold = integralThresh;
        BoFciDumper::write(fileName, _oneBody->at(Utils::ElementType::E),
                           _twoBody->at({Utils::ElementType::E, Utils::ElementType::E}), boFciDumpData);
        break;
      }
    }
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_AO2MO_H
