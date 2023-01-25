/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_GDIIS_H
#define KIWI_GDIIS_H

#include <Eigen/Dense>
#include <array>
#include <memory>

namespace Scine {
namespace Kiwi {
namespace Optimization {

/**
 * @brief An implementation of the GDIIS optimization acceleration algorithm.
 */
class GDIIS {
 public:
  GDIIS() = default;

  /**
   * @brief Construct a new GDIIS.
   * @param hessianInverse A reference to the Hessian inverse to be used for the
   *                       generation of an actual step. The content of the
   *                       reference can continuosly be updated using e.g. the
   *                       BFGS scheme or it can also just remain the identity.
   * @param maxm The maximum number of old steps to be stored and to be
   *             extrapolated from. Upon reaching the maximum number of points,
   *             the oldest one will be replaced with the newly given one.
   */
  GDIIS(std::shared_ptr<Eigen::MatrixXd> hessianInverse, unsigned int maxm)
    : _invH(std::move(hessianInverse)), _maxm(maxm), _nParam(_invH->cols()), _x(_invH->cols(), maxm), _g(_invH->cols(), maxm) {
    _g.setZero();
    _x.setZero();
  };
  /**
   * @brief Store data into the GDIIS without extrapolation.
   * @param parameters The current parameters.
   * @param gradients  The current gradients (for the current parameters).
   */
  void store(Eigen::VectorXd& parameters, Eigen::VectorXd& gradients) {
    unsigned int current = _cycle % _maxm;
    _x.col(current) = parameters;
    _g.col(current) = gradients;
    _cycle++;
  }
  /**
   * @brief Resets the storage of the GDIIS.
   */
  void flush() {
    _cycle = 0;
  }
  /**
   * @brief Stores the new data and extrapolates to optimum parameters using all
   *        stored data.
   * @param parameters The current parameters.
   * @param gradients  The current gradients (for the current parameters).
   */
  void update(Eigen::VectorXd& parameters, Eigen::VectorXd& gradients) {
    /*
     * If the GDIIS is supposed to extrapolate with a sufficient cache,
     *   exit and return the original parameters.
     */
    if (_maxm < 2) {
      return;
    }
    unsigned int current = _cycle % _maxm;
    _cycle++;
    _x.col(current) = parameters;
    _g.col(current) = gradients;

    unsigned int n = _maxm;
    if (_cycle < _maxm) {
      n = _cycle;
    }

    Eigen::VectorXd extrapolation;
    if (_cycle < 2) {
      return;
    }

    Eigen::MatrixXd B(n + 1, n + 1);
    B.col(n).setOnes();
    B.row(n).setOnes();
    B(n, n) = 0.0;
    Eigen::MatrixXd dx(_g.rows(), n);
    for (unsigned int i = 0; i < n; i++) {
      dx.col(i) = -*_invH * _g.col(i);
    }
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = 0; j < i + 1; j++) {
        const double tmp(dx.col(i).dot(dx.col(j)));
        B(i, j) = tmp;
        B(j, i) = tmp;
      }
    }

    // Generate Coefficients
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
    rhs[n] = 1.0;
    Eigen::VectorXd coefficients = B.fullPivHouseholderQr().solve(rhs);
    coefficients /= coefficients.segment(0, n).sum();

    // Extrapolate
    Eigen::VectorXd tmpGrad = Eigen::VectorXd::Zero(_nParam);
    Eigen::VectorXd tmpParam = Eigen::VectorXd::Zero(_nParam);
    for (unsigned int i = 0; i < n; i++) {
      tmpGrad.noalias() += coefficients[i] * _g.col(i);
      tmpParam.noalias() += coefficients[i] * _x.col(i);
    }

    // Estimate expected Change; return and flush saved gradients if positive
    Eigen::VectorXd dxDIIS = (tmpParam - parameters).eval();
    double expChange = tmpGrad.dot(dxDIIS) / (tmpGrad.norm() * dxDIIS.norm());
    if (expChange > 1e-4) {
      this->flush();
      return;
    }

    parameters.noalias() = tmpParam;
    gradients.noalias() = tmpGrad;
  }

 private:
  std::shared_ptr<Eigen::MatrixXd> _invH;
  unsigned int _maxm;
  unsigned int _nParam;
  unsigned int _cycle = 0;
  Eigen::MatrixXd _x;
  Eigen::MatrixXd _g;
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_GDIIS_H
