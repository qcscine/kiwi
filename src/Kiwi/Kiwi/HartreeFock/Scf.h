/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_SCF_H
#define KIWI_SCF_H

#include <Core/Log.h>
#include <Kiwi/HartreeFock/ConvergenceAccelerator.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/Projector/ProjectorTraits.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Utils/Scf/OrbitalPerturbation/RandomOrbitalMixer.h>
#include <Eigen/Eigenvalues>
#include <iomanip>
#include <utility>

namespace Scine {
namespace Kiwi {

/**
 * @class Scf
 * @brief This is the base class for SCF-type algorithms.
 * @tparam Symmetry
 */
template<SymmetryType Symmetry>
class Scf {
  using ProjectorType = typename ProjectorTrait<Symmetry>::ProjectorType;

 public:
  Scf(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data, const HartreeFockSettings& settings,
      ProjectionParameters projectionParameters = ProjectionParameters(), bool verbose = true)
    : molecule_(std::move(molecule)),
      data_(std::move(data)),
      settings_(settings),
      projector_(ProjectorType(molecule_, data_, projectionParameters)),
      verbose_(verbose),
      convergenceAccelerator_(ConvergenceAccelerator(molecule_, data_, settings_.accelerator, verbose_)) {
    convergenceAccelerator_.setElectronicFockDiisThresh(settings_.electronicDiisThresh);
    convergenceAccelerator_.setNuclearFockDiisThresh(settings_.nuclearDiisThresh);
    convergenceAccelerator_.setElectronicEDiisThresh(settings_.electronicEDiisThresh);
    convergenceAccelerator_.setNuclearEDiisThresh(settings_.nuclearEDiisThresh);
    // Default settings. Can be overwritten later with `setTrahSettings`
    trahSettings_ = std::make_shared<TRAHSettings>();
  }

  virtual ~Scf() = default;

  virtual auto run() -> void final;

  virtual auto singleIteration() -> void final;

  auto getEnergy() -> double;

  auto setTrahSettings(const std::shared_ptr<TRAHSettings>& trahSettings) -> void {
    trahSettings_ = trahSettings;
  };

 protected:
  std::shared_ptr<Molecule> molecule_;
  std::shared_ptr<Data> data_;
  HartreeFockSettings settings_;
  ProjectorType projector_;
  bool verbose_ = true;
  ConvergenceAccelerator convergenceAccelerator_;

 public:
  void setVerbose(bool verbose);

 protected:
  std::shared_ptr<TRAHSettings> trahSettings_;

  double lastEnergy_;
  double currentEnergy_;
  std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble> errorMap_;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_;

  bool last_;
  bool first_;
  std::map<Utils::ElementType, bool> firstMap_;
  std::size_t restartAt_ = 15;
  bool converged_ = false;

  // variables for incremental Fock matrix formation
  std::size_t resetIncrementalCounter_ = 0;
  std::size_t numFockRebuilds_ = 0;
  bool incremental_ = false;
  bool direct_ = false;

  bool acceptSteer_;
  double lastEnergySteer_;
  double currentEnergySteer_;
  std::map<Utils::ElementType, Utils::MolecularOrbitals> lastOrbitals_;

  auto checkOrbitals() -> void;
  auto mixOrbitals() -> void;
  auto resetOrbitals() -> void;
  auto acceptOrbitals() -> void;

  auto constructFromCOAO(Utils::ElementType type) -> void;
  auto computeNoTransformation(Utils::ElementType type) -> void;

  virtual auto runScf() -> void = 0;
  virtual auto convergenceCheck() -> bool;

 private:
  virtual auto printHeader() -> void = 0;
  virtual auto printIteration() -> void = 0;
  virtual auto printAnalysis() -> void;
};

template<SymmetryType Symmetry>
auto Scf<Symmetry>::mixOrbitals() -> void {
  auto log = Core::Log();
  for (auto const& elem : *molecule_) {
    auto type = elem.first;

    if (molecule_->at(type).isRestricted) {
      int N = elem.second.N;
      auto mixer = Utils::OrbitalPerturbation::RandomOrbitalMixer(data_->C_OAO[type], N);
      mixer.considerAllOrbitals();
      mixer.setMaximalMixAngle(settings_.maxAngle);
      mixer.setNumberMixes(settings_.numberOfMixes);
      mixer.mix(log);
    }
    else {
      int beta = elem.second.msVector[1];
      auto mixer = Utils::OrbitalPerturbation::RandomOrbitalMixer(data_->C_OAO[type], elem.second.msVector[0], beta);
      mixer.considerAllOrbitals();
      mixer.setMaximalMixAngle(settings_.maxAngle);
      mixer.setNumberMixes(settings_.numberOfMixes);
      mixer.mix(log);
    }

    constructFromCOAO(type);
  }
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::resetOrbitals() -> void {
  data_->C_OAO = lastOrbitals_;

  for (auto const& elem : *molecule_) {
    auto type = elem.first;
    constructFromCOAO(type);
    if (Symmetry != SymmetryType::None) {
      computeNoTransformation(type);
    }
  }
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::acceptOrbitals() -> void {
  lastOrbitals_ = data_->C_OAO;
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::constructFromCOAO(Utils::ElementType type) -> void {
  if (molecule_->at(type).isRestricted) {
    data_->C[type].restrictedMatrix() = data_->X[type] * data_->C_OAO[type].restrictedMatrix();
  }
  else {
    data_->C[type].alphaMatrix() = data_->X[type] * data_->C_OAO[type].alphaMatrix();
    data_->C[type].betaMatrix() = data_->X[type] * data_->C_OAO[type].betaMatrix();
  }
  auto const& C = data_->C.at(type);
  auto const& C_OAO = data_->C_OAO.at(type);
  data_->D.at(type) = HartreeFockUtils::makeDensity(type, molecule_, C);
  data_->D_OAO.at(type) = HartreeFockUtils::makeDensity(type, molecule_, C_OAO);
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::computeNoTransformation(Utils::ElementType type) -> void {
  es_.compute(data_->D_OAO[type].restrictedMatrix());
  data_->O_NO.at(type).restrictedMatrix() = es_.eigenvectors().rowwise().reverse();
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::getEnergy() -> double {
  return lastEnergy_;
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::singleIteration() -> void {
  projector_.init(false);
  projector_.evaluateDensity(0);
  projector_.updateFockMatrices();
  for (auto const& elem : *molecule_) {
    projector_.updateDensity(elem.first);
  }
  projector_.finalize();
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::run() -> void {
  if (verbose_) {
    auto& clock = Clock::getInstance();
    clock.time("scf", true);
  }

  currentEnergySteer_ = -1;
  lastEnergySteer_ = 0;
  // if (Symmetry==SymmetryType::ParityAm) {
  if (molecule_->size() > 1) {
    convergenceAccelerator_.setSubspaceSize(16);
  }
  else {
    convergenceAccelerator_.setSubspaceSize(8);
  }

  if (settings_.maxIterations > 0) {
    if (settings_.numberOfPerturbations > 0) {
      projector_.init(verbose_);
      projector_.finalize();
      lastEnergy_ = projector_.getEnergy();

      auto perturb = 0;
      acceptSteer_ = true;

      // Set `lastOrbitals_`
      acceptOrbitals();

      while (true) {
        ++perturb;
        if (verbose_) {
          std::cout << "   *******************************\n";
          std::cout << "     Macro iteration " << perturb << std::endl;
          std::cout << "   *******************************\n";
        }
        mixOrbitals();

        runScf();
        currentEnergySteer_ = lastEnergy_;

        checkOrbitals();

        if (perturb == settings_.numberOfPerturbations) {
          break;
        }
      }

      if (!acceptSteer_ || !converged_) {
        if (verbose_) {
          std::cout << "\nLast orbitals were not accepted or last SCF did not converge\n    --> finalization with best "
                       "orbitals.\n\n";
        }
        resetOrbitals();
        runScf();
      }
    }
    else {
      runScf();
    }
  }
  else {
    projector_.init(verbose_);
    projector_.finalize();
    lastEnergy_ = projector_.getEnergy();
  }

  if (verbose_) {
    // TODO fix
    // std::cout << std::endl << "Avrg. time of Fock matrix formation:    ";
    // std::cout << std::fixed << std::setprecision(3) << std::setw(16)
    //          << projector_.getFockMatrixFormationTime().count() * 1.0e-3 << " sec" << std::endl;
    std::cout << "Time spent in SCF procedure:            ";
    auto& clock = Clock::getInstance();
    clock.time("scf");

    std::cout << "\n";
    printAnalysis();

    std::cout << "\n-------------------------   --------------------" << std::endl;
    std::cout << "FINAL SINGLE POINT ENERGY      " << std::setprecision(10) << lastEnergy_ << std::endl;
    std::cout << "-------------------------   --------------------" << std::endl;
  }
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::printAnalysis() -> void {
  using std::cout;
  using std::endl;
  using std::setw;

  bool foundUnrestricted = false;
  for (auto const& elem : *molecule_) {
    if (!elem.second.isRestricted) {
      foundUnrestricted = true;
    }
  }
  if (!foundUnrestricted) {
    return;
  }
  cout << std::left << setw(5) << "Type" << std::string(5, ' ') << setw(15) << "<S^2>" << std::string(5, ' ')
       << setw(15) << "S(S+1)" << std::string(5, ' ') << setw(15) << "<S^2> - S(S+1)" << endl;
  cout << std::string(70, '-') << endl;
  cout << std::fixed << std::right;
  for (auto const& elem : *molecule_) {
    if (!elem.second.isRestricted) {
      cout << setw(5) << Utils::ElementInfo::symbol(elem.first) << std::string(5, ' ');
      cout << setw(15) << std::setprecision(10) << HartreeFockUtils::S2(molecule_, data_, elem.first).first
           << std::string(5, ' ');
      cout << setw(15) << std::setprecision(10) << HartreeFockUtils::S2(molecule_, data_, elem.first).second
           << std::string(5, ' ');
      cout << setw(15) << std::setprecision(10)
           << HartreeFockUtils::S2(molecule_, data_, elem.first).first -
                  HartreeFockUtils::S2(molecule_, data_, elem.first).second;
      cout << endl;
    }
  }
  std::cout << "\n";
}

template<SymmetryType Symmetry>
auto Scf<Symmetry>::convergenceCheck() -> bool {
  /*
   * TODO: insert extensive criterium: e.g., -> Error / squareRoot( N )
   */
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      if (errorMap_[elem.first].restricted > settings_.fockMatrixNormThresh) {
        return false;
      }
    }
    else {
      if (errorMap_[elem.first].alpha > settings_.fockMatrixNormThresh) {
        return false;
      }
      if (molecule_->at(elem.first).msVector[1] > 0) {
        if (errorMap_[elem.first].beta > settings_.fockMatrixNormThresh) {
          return false;
        }
      }
    }
  }
  return std::abs(lastEnergy_ - currentEnergy_) <= settings_.energyThreshold;
}
template<SymmetryType Symmetry>
auto Scf<Symmetry>::checkOrbitals() -> void {
  if (currentEnergySteer_ < lastEnergySteer_) {
    std::cout << "\n    --> Orbitals accepted\n\n";
    acceptOrbitals();
    lastEnergySteer_ = currentEnergySteer_;
    acceptSteer_ = true;
  }
  else {
    std::cout << "\n    --> Orbitals rejected\n\n";
    resetOrbitals();
    acceptSteer_ = false;
  }
}

template<SymmetryType Symmetry>
void Scf<Symmetry>::setVerbose(bool verbose) {
  Scf::verbose_ = verbose;
}

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_SCF_H
