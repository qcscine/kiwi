/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_ENGINE_H
#define KIWI_ENGINE_H

#include <Utils/DataStructures/BasisSet.h>
#include <Eigen/Dense>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace Scine {
namespace Kiwi {
namespace ParticleDensity {

/**
 * @class Engine
 * @brief Handles the calculation of the delta-function integral at a point (x,y,z) between two AOs with spherical or
 * cartesian GTOs.
 */
class Engine {
 public:
  Engine(const Utils::Integrals::BasisSet& basis) : _basis(basis), _dimension(basis.nbf()) {
  }

  auto evaluateDensity(double x, double y, double z) -> Eigen::MatrixXd {
    auto dens = density(x, y, z);
    return dens * dens.transpose();
  }

  auto evaluateOrbitals(double x, double y, double z) -> Eigen::VectorXd {
    return density(x, y, z);
  }

 private:
  const Utils::Integrals::BasisSet& _basis;
  std::size_t _dimension;

  inline static auto gaussian(const double a, const double x, const double y, const double z) -> double {
    return std::exp(-a * x * x - a * y * y - a * z * z);
  }

  inline static auto cartesianPrefactor(const int& lx, const int& ly, const int& lz, const double& x, const double& y,
                                        const double& z) -> double {
    return std::pow(x, lx) * std::pow(y, ly) * std::pow(z, lz);
  }

  inline auto density(const double& x, const double& y, const double& z) const -> Eigen::VectorXd {
    Eigen::VectorXd dens(_dimension);

    int index = 0;
    for (const auto& shell : _basis) {
      const auto& O = shell.getShift();
      const auto& alpha = shell.getVecAlpha();
      const auto& contr = shell.getVecCoeffs();
      auto dx = x - O(0);
      auto dy = y - O(1);
      auto dz = z - O(2);
      //
      // phi(r) = ( sum_i C_i exp(-a_i r^2 ) ) *
      //
      // S = sum_i C_i exp(-a_i r^2 )
      double S = 0.0;
      for (auto a = 0; a < int(alpha.size()); ++a) {
        S += contr.at(a) * gaussian(alpha.at(a), dx, dy, dz);
      }
      if (shell.isPureSolid()) {
        // S_l^m(x,y,z)
        densitySpherical(dens, index, shell, dx, dy, dz, S);
      }
      else {
        // x^i y^j z^k
        densityCartesian(dens, index, shell, dx, dy, dz, S);
      }
    }

    return dens;
  }

  inline static auto densityCartesian(Eigen::VectorXd& dens, int& index, const Utils::Integrals::Shell& shell,
                                      const double& dx, const double& dy, const double& dz, const double& S) -> void {
    int l = int(shell.l());
    // official mapping from libint manual
    for (int i = 0; i <= l; ++i) {
      int lx = l - i; // exponent of x
      for (int j = 0; j <= i; ++j) {
        int ly = i - j; // exponent of y
        int lz = j;     // exponent of z
        // 2.
        // phi(r) = S * x^lx y^ly z^lz
        dens(index) = S * cartesianPrefactor(lx, ly, lz, dx, dy, dz);
        ++index;
      }
    }
  }

  inline static auto densitySpherical(Eigen::VectorXd& dens, int& index, const Utils::Integrals::Shell& shell,
                                      const double& dx, const double& dy, const double& dz, const double& S) -> void {
    namespace bg = boost::geometry;
    bg::model::point<double, 3, bg::cs::cartesian> point_cart(dx, dy, dz);
    bg::model::point<double, 3, bg::cs::spherical<bg::radian>> point_sph;
    bg::transform(point_cart, point_sph);
    const double& phi = bg::get<0>(point_sph);   // phi in [0,2 pi]
    const double& theta = bg::get<1>(point_sph); // theta in [0, pi]
    const double& r = bg::get<2>(point_sph);     // r = sqrt(x^2+y^2+z^2)

    //
    // phi(r) = ( sum_i C_i exp(-a_i r^2 ) ) * Y_l^m(x,y,z)
    //
    int l = int(shell.l());
    // Standard ordering
    // For real solid harmonics, see
    // Helgaker -- Molecular electronic structure theory, 6.4.1, Eqns. 6.4.19-21
    for (int m = -l; m <= l; ++m) {
      // 2. Y_l^m(x,y,z)
      // m=0
      if (m == 0) {
        dens(index) =
            S * std::sqrt((4 * M_PI) / (2 * l + 1)) * std::pow(r, l) * boost::math::spherical_harmonic_r(l, m, theta, phi);
      }
      else if (m < 0) {
        dens(index) = S * std::pow(-1, m) * std::sqrt((8 * M_PI) / (2 * l + 1)) * std::pow(r, l) *
                      boost::math::spherical_harmonic_i(l, m, theta, phi);
      }
      else if (m > 0) {
        dens(index) = S * std::pow(-1, m) * std::sqrt((8 * M_PI) / (2 * l + 1)) * std::pow(r, l) *
                      boost::math::spherical_harmonic_r(l, m, theta, phi);
      }
      ++index;
    }
  }
};

} // namespace ParticleDensity
} // namespace Kiwi
} // namespace Scine
#endif // KIWI_ENGINE_H
