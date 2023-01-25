/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_STABILITYANALYSIS_H
#define KIWI_STABILITYANALYSIS_H

#include <Kiwi/HartreeFock/SecondOrder/ArhDavidson.h>
#include <Kiwi/HartreeFock/SecondOrder/ArhInterface.h>
#include <Kiwi/HartreeFock/SecondOrder/SigmaVectorEvaluator.h>
#include <Kiwi/KiwiOpt/Davidson.h>
#include <iomanip>

namespace Scine {
namespace Kiwi {

auto buildHessian(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data) -> void {
  std::shared_ptr<ArhInterface> ptrArhInterface_;
  std::shared_ptr<Arh::SigmaVectorEvaluator> ptrSigmaVectorEvaluator_;

  Data backUp;

  // bool lastStepRejected = false;
  // bool secondOrderRegion = false;
  // int numRejectionsInARow = 0;

  ptrArhInterface_ = std::make_shared<ArhInterface>(molecule, data);

  ptrArhInterface_->setUseExactHessian(true);

  // int numEigPairs = 1;
  int maxSubspaceDimension = 50;

  //
  ptrArhInterface_->setMaxSize(1);
  ptrArhInterface_->addIteration();
  ptrArhInterface_->evaluateGradient();
  ptrSigmaVectorEvaluator_ =
      std::make_shared<Arh::SigmaVectorEvaluator>(ptrArhInterface_, maxSubspaceDimension, false, false, false, false);

  auto dim = ptrSigmaVectorEvaluator_->fullDim;
  auto dimVector = ptrArhInterface_->getDim();

  Eigen::MatrixXd Hessian(dim, dim);

  std::cout << "Starting explicit Hessian evaluation.\n   ... this could take a while.\n\n";
  std::cout << "Progress\n";

  // for (auto row = 0; row < dim; ++row) {
  //   for (auto col = 0; col <= row; ++col) {
  //     std::cout << row << "/" << col << std::endl;
  //     Eigen::VectorXd rowVector = Eigen::VectorXd::Zero(dim);
  //     Eigen::VectorXd colVector = Eigen::VectorXd::Zero(dim);
  //     rowVector(row) = 1;
  //     colVector(col) = 1;

  //    std::cout << "<H(" << row << ")," << col << "> = " << std::fixed << std::setprecision(4);
  //    std::cout << ptrSigmaVectorEvaluator_->evaluate(rowVector).transpose() * colVector << std::endl;
  //    std::cout << "<" << row << ",H(" << col << ")> = " << std::fixed << std::setprecision(4);
  //    std::cout << rowVector.transpose() * ptrSigmaVectorEvaluator_->evaluate(colVector) << std::endl;
  //  }
  //}
  for (auto col = 0; col < dim; ++col) {
    Eigen::VectorXd colVector = Eigen::VectorXd::Zero(dim);
    colVector(col) = 1;

    Hessian.col(col) = ptrSigmaVectorEvaluator_->evaluate(colVector);
  }
  std::cout << "Finished\n";
  // std::cout << "Hessian\n" << std::fixed << std::setprecision(4) << Hessian-Hessian.transpose() << std::endl;

  std::cout << "Symmetry test = " << std::fixed << std::setprecision(10)
            << (Hessian - Hessian.transpose()).cwiseAbs().sum() << std::endl;

  // for (auto i=0UL; i<dimVector.size(); ++i) {
  //   for (auto j=0UL; j<=i; ++j) {
  //     int iBlockStart=0;
  //     int jBlockStart=0;
  //     for (auto it=0UL; it<i; it++) {
  //       iBlockStart += dimVector[it];
  //     }
  //     for (auto jt=0UL; jt<j; jt++) {
  //       jBlockStart += dimVector[jt];
  //     }
  //     std::cout << i << "/" << j << std::endl;
  //     std::cout << "Norm of block divided by size of block = " << Hessian.block(iBlockStart, jBlockStart,
  //     dimVector[i], dimVector[j]).norm() / (dimVector[i]*dimVector[j]) << std::endl; std::cout << "Hessian block\n" <<
  //     std::scientific << std::setprecision(4) << Hessian.block(iBlockStart, jBlockStart, dimVector[i], dimVector[j])
  //     << std::endl;
  //   }
  // }
}

/**
 * @brief Cacluates the lowest eigenvalue of the Hessian. If the eigenvalue is negative, a saddle point is encountered.
 * If it is 0, in principle the result is mathematically inconclusive, but in practice it can be safely assumed that
 * a local minimum is reached. In case the eigenvalue is positive, a local minimum is encountered.
 * @param molecule
 * @param data
 */
auto runStabilityAnalysis(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Data>& data) -> void {
  std::shared_ptr<ArhInterface> ptrArhInterface_;
  std::shared_ptr<Arh::SigmaVectorEvaluator> ptrSigmaVectorEvaluator_;

  Data backUp;

  ptrArhInterface_ = std::make_shared<ArhInterface>(molecule, data);

  ptrArhInterface_->setUseExactHessian(true);

  int maxSubspaceDimension = 100;

  //
  ptrArhInterface_->setMaxSize(1);
  ptrArhInterface_->addIteration();
  ptrArhInterface_->evaluateGradient();
  ptrSigmaVectorEvaluator_ =
      std::make_shared<Arh::SigmaVectorEvaluator>(ptrArhInterface_, maxSubspaceDimension, true, true, true, true);
  ptrSigmaVectorEvaluator_->addRandomTrialVector();

  auto thresh = ptrSigmaVectorEvaluator_->getGradientNorm();
  if (thresh > 1e-4) {
    thresh = 1e-4;
  }

  std::cout << "\n\n";
  std::cout << "   **********************\n";
  std::cout << "     Stability analysis " << std::endl;
  std::cout << "   **********************\n";
  std::cout << "\n\n";

  Arh::Davidson davidson(ptrSigmaVectorEvaluator_);
  davidson.setMaxIterations(maxSubspaceDimension);
  davidson.setThresh(thresh);

  davidson.compute(true);

  // buildHessian(molecule, data);
}

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_STABILITYANALYSIS_H
