/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_SECONDORDER_H
#define KIWI_SECONDORDER_H

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/Scf.h>
#include <Kiwi/HartreeFock/SecondOrder/ArhInterface.h>
#include <Kiwi/HartreeFock/SecondOrder/Driver.h>
#include <Kiwi/HartreeFock/SecondOrder/SigmaVectorEvaluator.h>
#include <iomanip>

namespace Scine {
namespace Kiwi {

template<SymmetryType Symmetry>
class SecondOrder : public Scf<Symmetry> {
  using Scf<Symmetry>::Scf;
  using Scf<Symmetry>::Scf::molecule_;
  using Scf<Symmetry>::Scf::data_;
  using Scf<Symmetry>::Scf::projector_;
  using Scf<Symmetry>::Scf::currentEnergy_;
  using Scf<Symmetry>::Scf::errorMap_;
  using Scf<Symmetry>::Scf::lastEnergy_;
  using Scf<Symmetry>::Scf::settings_;
  using Scf<Symmetry>::Scf::trahSettings_;
  using Scf<Symmetry>::Scf::last_;
  using Scf<Symmetry>::Scf::convergenceAccelerator_;
  using Scf<Symmetry>::Scf::restartAt_;
  using Scf<Symmetry>::Scf::converged_;
  using Scf<Symmetry>::Scf::convergenceCheck;
  using Scf<Symmetry>::Scf::verbose_;

 private:
  int _macroIter;
  int numberOfRejections_;
  int _microIter;

  Eigen::VectorXd gradient_;
  // Eigen::VectorXd gradient_;
  Eigen::VectorXd solution_;

  std::string str_optimizer_;

  std::string str_method_;

  double gTx;
  double gTHx;
  double deltaEpredicted;
  double deltaE;

  double r_;

  double gradientNorm_;
  double thresh_;
  double davidsonError_;
  int davidsonMicroIterations_;
  bool localRegion = false;
  double trustRadius_;

  std::shared_ptr<ArhInterface> ptrArhInterface_;
  std::shared_ptr<Arh::SigmaVectorEvaluator> ptrSigmaVectorEvaluator_;

  auto runScf() -> void override;
  auto printHeader() -> void override;
  auto printIteration() -> void override;
};

template<SymmetryType Symmetry>
auto SecondOrder<Symmetry>::runScf() -> void {
  projector_.init(verbose_);
  lastEnergy_ = projector_.getEnergy();

  convergenceAccelerator_.init();
  convergenceAccelerator_.reset();

  const auto& maxIter = settings_.maxIterations;

  if (verbose_)
    printHeader();
  converged_ = false;
  _macroIter = 0;
  _microIter = 0;

  Data backUp;

  localRegion = false;
  bool lastStepRejected = false;
  TRAHSettings::TrahOptimizer opt = trahSettings_->opt;

  trustRadius_ = trahSettings_->initialTrustRadius;

  int numRejectionsInARow = 0;

  ptrArhInterface_ = std::make_shared<ArhInterface>(molecule_, data_);

  ptrArhInterface_->setMaxSize(trahSettings_->maxArhDim);

  if (opt == TRAHSettings::TrahOptimizer::ARH) {
    ptrArhInterface_->setUseExactHessian(false);
  }
  else if (opt == TRAHSettings::TrahOptimizer::Newton) {
    ptrArhInterface_->setUseExactHessian(true);
  }

  // variables for incremental Fock matrix formation
  std::size_t resetIncrementalCounter = 0;
  // We start with `1` in the ARH method, such that the prescreening threshold is tighter by a factor of 1e-2.
  std::size_t numFockRebuilds = 1;
  bool incremental = false;
  bool direct = data_->integralDirect;
  bool disableIncremental = settings_.disableIncrementalFock;
  last_ = false;
  bool close = false;

  // if (!disableIncremental && direct) {
  //  settings_.resetIncremental=5;
  //  std::cout << "Settings incremental build frequency to 5 to avoid numerical instabilities.\n";
  //}

  if (direct && !disableIncremental) {
    if (verbose_)
      std::cout << "                      *** Starting incremental Fock-matrix formation. ***" << std::endl;
    // if (settings_.resetIncremental>5) {
    //  settings_.resetIncremental = 5;
    //  std::cout << "                          Setting incremental build frequency to 5\n"
    //               "                          to avoid numerical instabilities.\n";
    //}
  }

  while (true) {
    // Back up:
    Kiwi::Data backUpData;
    Kiwi::HartreeFockData backupHFdata;
    backUpData.D = data_->D;
    backUpData.D_OAO = data_->D_OAO;
    backUpData.C = data_->C;
    backUpData.C_OAO = data_->C_OAO;
    backupHFdata = *data_->hartreeFockData;

    // Store the last density for incremental Fock matrix building
    projector_.setLastDensity(data_->D);

    // Hartree--Fock gradient
    errorMap_ = projector_.getErrorMap();
    //
    ptrArhInterface_->addIteration();
    ptrArhInterface_->evaluateGradient();
    ptrSigmaVectorEvaluator_ = std::make_shared<Arh::SigmaVectorEvaluator>(
        ptrArhInterface_, trahSettings_->maxDavidsonSubspaceDimension, trahSettings_->useSecondStartVector,
        trahSettings_->useWhiteNoiseSecondStartVec, trahSettings_->useThirdStartVector,
        trahSettings_->useWhiteNoiseThirdStartVec);

    gradient_ = ptrSigmaVectorEvaluator_->getGradient();
    gradientNorm_ = gradient_.norm();
    localRegion = gradientNorm_ < trahSettings_->localRegionThresh;

    if (lastStepRejected) {
      localRegion = false;
    }

    if (opt == TRAHSettings::TrahOptimizer::Newton) {
      str_method_ = "Newton";
    }
    else {
      str_method_ = "ARH";
    }

    thresh_ = gradientNorm_ * trahSettings_->microThreshGradScaling;

    if (thresh_ > trahSettings_->microMinThresh) {
      thresh_ = trahSettings_->microMinThresh;
    }

    Arh::Driver driver(ptrSigmaVectorEvaluator_, thresh_);

    driver.setTrustRadius(trustRadius_);

    if (!localRegion) {
      driver.evaluateTRAH(trahSettings_->maxDavidsonIterations);

      davidsonMicroIterations_ = driver.getIterations();
      solution_ = driver.getSolution();
      davidsonError_ = driver.getError();
      _microIter += davidsonMicroIterations_;
    }
    else {
      str_optimizer_ = "DAVID";
      driver.evaluateNR(trahSettings_->maxDavidsonIterations);
      // driver.evaluatePCG(trahSettings_->maxDavidsonIterations);
      davidsonMicroIterations_ = driver.getIterations();
      solution_ = driver.getSolution();
      trustRadius_ = solution_.norm();
      davidsonError_ = driver.getError();
      _microIter += davidsonMicroIterations_;
    }

    gradient_ = ptrArhInterface_->getGradient();
    gTx = gradient_.dot(solution_);
    gTHx = gradient_.dot(ptrSigmaVectorEvaluator_->evaluate(solution_));

    deltaEpredicted = gTx + 0.5 * gTHx;

    ptrArhInterface_->updateX(solution_);

    for (auto const& elem : *molecule_) {
      auto orbitalRotation = ptrArhInterface_->getOrbitalRotationMatrix(elem.first);
      projector_.updateDensityFromRotation(elem.first, orbitalRotation);
    }

    // TODO in cases that are prone to numerical instabilities, close to convergence the incremental formation should be
    // disabled
    if (!close && direct && !disableIncremental && gradientNorm_ < 1e-4) {
      // disableIncremental=true;
      incremental = false;
      ++numFockRebuilds;
      resetIncrementalCounter = 0;
      close = true;
      if (verbose_)
        std::cout << "                      *** Close to convergence: resetting incremental Fock-matrix formation. ***"
                  << std::endl;
    }
    else if (direct && !disableIncremental) {
      // if (direct && !disableIncremental) {
      if (!last_ && resetIncrementalCounter < settings_.resetIncremental) {
        incremental = true;
        ++resetIncrementalCounter;
      }
      else {
        incremental = false;
        resetIncrementalCounter = 0;
        ++numFockRebuilds;
        if (verbose_)
          std::cout << "                      *** Resetting incremental Fock-matrix formation. ***\n";
      }
    }

    // No incremental Fock-matrix formation in first iteration. The guess could be too far away.
    if (_macroIter == 0) {
      incremental = false;
    }

    projector_.evaluateDensity(incremental);

    projector_.updateFockMatrices(incremental, numFockRebuilds);

    projector_.evaluateEnergy();

    currentEnergy_ = projector_.getEnergy();

    deltaE = currentEnergy_ - lastEnergy_;

    ++_macroIter;

    if (verbose_)
      printIteration();

    /*
     *  SCF failed
     */
    if (_macroIter == int(maxIter)) {
      converged_ = false;
      if (verbose_)
        std::cout << "SCF did not converge after " << maxIter << " iterations." << std::endl;
      lastEnergy_ = currentEnergy_;
      break;
    }

    if (direct && !disableIncremental) {
      if (!last_) {
        if (convergenceCheck()) {
          if (verbose_)
            std::cout << "                      *** Energy converged: resetting incremental Fock matrix ***" << std::endl;
          last_ = true;
        }
      }
      else {
        if (convergenceCheck()) {
          if (verbose_)
            std::cout << "                      *** SCF converged ***" << std::endl;
          converged_ = true;
          lastEnergy_ = currentEnergy_;
          break;
        }
        if (verbose_) {
          std::cout << "                      *** Convergence criterium was violated in last iteration ***" << std::endl;
          std::cout << "                      *** SCF will continue ***" << std::endl;
        }
        convergenceAccelerator_.reset();
        last_ = false;
      }
    }
    else {
      if (convergenceCheck()) {
        if (verbose_)
          std::cout << "                      *** SCF converged ***" << std::endl;
        converged_ = true;
        lastEnergy_ = currentEnergy_;
        break;
      }
    }

    r_ = deltaE / (deltaEpredicted);

    // Reject step
    if (deltaE > 0 && !last_) {
      lastStepRejected = true;

      ++numRejectionsInARow;

      if (verbose_) {
        std::cout << "  .----------------------------------" << std::string(12, '-') << ".\n";
        std::cout << "  |  Steps rejected in a row      :  " << std::setw(2) << numRejectionsInARow
                  << std::string(10, ' ') << "|\n";
        std::cout << "  |  actual/predicted energy ratio:  " << std::setw(10) << std::setprecision(3) << r_ << "  |\n";
        std::cout << "  |  Current trust radius:           " << std::setw(10) << std::setprecision(3) << trustRadius_ << "  |\n";
      }
      ++numberOfRejections_;
      trustRadius_ *= 0.7;
      if (trustRadius_ < 1e-7) {
        throw std::runtime_error("Trust radius too small");
      }

      if (verbose_) {
        std::cout << "  |  New trust radius:               " << std::setw(10) << std::setprecision(3) << trustRadius_ << "  |\n";
        std::cout << "  '----------------------------------" << std::string(12, '-') << "'\n";
      }
      // Back up - begin
      data_->D = backUpData.D;
      data_->D_OAO = backUpData.D_OAO;
      data_->C = backUpData.C;
      data_->C_OAO = backUpData.C_OAO;
      data_->hartreeFockData->F_OAO = backupHFdata.F_OAO;
      data_->hartreeFockData->G = backupHFdata.G;
      data_->hartreeFockData->F = backupHFdata.F;
      data_->hartreeFockData->J = backupHFdata.J;
      data_->hartreeFockData->K = backupHFdata.K;
      data_->hartreeFockData->I = backupHFdata.I;
      data_->hartreeFockData->L = backupHFdata.L;
      // Back up - end
    }
    // Accept step
    else if (!localRegion) {
      numRejectionsInARow = 0;
      if (0 < r_ && r_ <= 0.25) {
        trustRadius_ *= 0.7;
      }
      else if (r_ > 0.75) {
        if ((trustRadius_ * 1.2) <= trahSettings_->maxTrustRadius) {
          trustRadius_ *= 1.2;
        }
        else {
          trustRadius_ = trahSettings_->maxTrustRadius;
        }
      }
      lastStepRejected = false;
      lastEnergy_ = currentEnergy_;
    }
    else if (localRegion) {
      lastStepRejected = false;
      lastEnergy_ = currentEnergy_;
    }
  }

  if (verbose_)
    std::cout << "\nTotal number of microiterations = " << _microIter << "\n\n";

  projector_.finalize();
}

template<SymmetryType Symmetry>
auto SecondOrder<Symmetry>::printHeader() -> void {
  using std::cout;
  using std::endl;
  using std::setw;

  int count = 0;
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      ++count; // [F,D]
    }
    else {
      if (elem.second.msVector[0] > 0) {
        ++count; // [F,D]
      }
      if (elem.second.msVector[1] > 0) {
        ++count; // [F,D]
      }
    }
  }
  cout << std::right << endl;
  cout << std::string(15 + 10 + 21 + 19 + 2 * 14 + 10 + count * 14 + 10 + 6, '-') << endl;
  cout << setw(15) << "Iteration";
  cout << setw(10) << "Microit";
  cout << setw(21) << "Energy/Ha" << setw(19) << "Delta E/Ha";
  cout << setw(14) << "Gradient";
  cout << setw(14) << "Error";
  cout << setw(10) << "TRadius";
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      cout << setw(12) << "[F,D]:" << setw(2) << Utils::ElementInfo::symbol(elem.first);
    }
    else {
      if (elem.second.msVector[0] > 0) {
        cout << setw(12) << "[F,D]-a:" << setw(2) << Utils::ElementInfo::symbol(elem.first);
      }
      if (elem.second.msVector[1] > 0) {
        cout << setw(12) << "[F,D]-b:" << setw(2) << Utils::ElementInfo::symbol(elem.first);
      }
    }
  }
  cout << setw(10) << "Method";
  cout << endl;
  cout << std::string(15 + 10 + 21 + 19 + 2 * 14 + 10 + count * 14 + 10 + 6, '-') << endl;
}

template<SymmetryType Symmetry>
auto SecondOrder<Symmetry>::printIteration() -> void {
  using std::cout;
  using std::endl;
  using std::setprecision;
  using std::setw;

  cout << std::fixed << std::right;
  cout << setw(15) << _macroIter;
  cout << setw(10) << davidsonMicroIterations_;
  cout << setw(21) << setprecision(12) << currentEnergy_;
  cout << setw(19) << setprecision(12) << currentEnergy_ - lastEnergy_;
  cout << std::scientific << std::right;
  cout << setw(14) << setprecision(3) << gradientNorm_;
  cout << setw(14) << setprecision(3) << davidsonError_;
  cout << std::fixed << std::right;
  if (!localRegion) {
    cout << setw(10) << setprecision(3) << trustRadius_;
  }
  else {
    cout << setw(10) << str_optimizer_;
  }
  cout << std::scientific << std::right;
  for (auto const& elem : *molecule_) {
    if (elem.second.isRestricted) {
      cout << setw(14) << setprecision(3) << errorMap_[elem.first].restricted;
    }
    else {
      if (elem.second.msVector[0] > 0) {
        cout << setw(14) << setprecision(3) << errorMap_[elem.first].alpha;
      }
      if (elem.second.msVector[1] > 0) {
        cout << setw(14) << setprecision(3) << errorMap_[elem.first].beta;
      }
    }
  }
  cout << setw(10) << str_method_;
  cout << endl;
  cout << std::fixed << std::right;
}

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_SECONDORDER_H
