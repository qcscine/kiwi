/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_PRECONDITIONEDDIRECTCONJUGATEGRADIENT_H
#define KIWI_PRECONDITIONEDDIRECTCONJUGATEGRADIENT_H

#include <Eigen/Eigenvalues>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {
namespace Optimization {

template<typename T>
class PreconditionedDirectConjugateGradient {
  using VectorT = Eigen::Matrix<T, Eigen::Dynamic, 1>;

 private:
  T thresh_ = static_cast<T>(double(1e-16));

  VectorT r_;
  VectorT z_;
  VectorT p_;
  VectorT Ap_;
  VectorT x_;

  VectorT b_;
  double rDotz_previous_;
  double rDotz_;
  int maxIterations_ = 100;
  int iterations_;

  T error_;

  double alpha_ = 0;
  double beta_ = 0;
  bool successful_ = false;

  bool wasSuccessful() const {
    return successful_;
  }

 public:
  PreconditionedDirectConjugateGradient(VectorT&& x) : x_(std::move(x)) {
    p_.resize(x_.size());
    r_.resize(x_.size());
  }

  /**
   * Guess is the trial vector, i.e., linearly transformed operator A.
   * This approach has the advantage that we do not require the specific knowledge of A, but just
   * the linearly transformed A.
   * A can also be very large, so it has never to be build specifically.
   *
   *      x = A x_0;
   *
   * @param x
   */
  auto updateGuess(VectorT&& x) -> void {
    x_ = std::move(x);
  }

  /**
   * Image means im(A) = { b \in B : b = Ax }
   *    Ax = b
   * @param b
   */
  auto updateImage(VectorT&& b) -> void {
    b_ = std::move(b);
  }

  VectorT&& getSolution() {
    return std::move(x_);
  }

  void setMaxIterations(int maxIterations) {
    maxIterations_ = maxIterations;
  }

  int getIterations() const {
    return iterations_;
  }

  T getError() const {
    return error_;
  }

  void setThresh(T thresh) {
    thresh_ = thresh;
  }

  template<class SigmaVectorEvaluator, typename MatrixType>
  auto optimize(SigmaVectorEvaluator& evaluator, const MatrixType& precond) -> void {
    successful_ = false;

    r_ = b_ - evaluator.evaluate(x_);
    z_ = precond * r_;
    rDotz_previous_ = r_.template dot(z_);
    p_ = z_;

    iterations_ = 1;
    while (true) {
      if (iterations_ == maxIterations_) {
        std::cout << "PreconditionedDirectConjugateGradient did not convergege after " << iterations_ << " iterations.\n";
        break;
      }

      Ap_ = evaluator.evaluate(p_);

      alpha_ = rDotz_previous_ / p_.template dot(Ap_);

      x_ += alpha_ * p_;
      r_ -= alpha_ * Ap_;

      error_ = r_.norm();

      if (error_ < thresh_) {
        successful_ = true;
        break;
      }

      z_ = precond * r_;
      rDotz_ = r_.template dot(z_);
      beta_ = rDotz_ / rDotz_previous_;

      p_ = z_ + beta_ * p_;

      rDotz_previous_ = rDotz_;
      ++iterations_;
    }
  }
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_PRECONDITIONEDDIRECTCONJUGATEGRADIENT_H
