/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_AO2MOHELPER_H
#define KIWI_AO2MOHELPER_H

#include <Eigen/Dense>
#include <iostream>

namespace Scine {
namespace Kiwi {

/**
 * @class AO2MOHelper
 * @brief Implements the AO to MO transformation with a scaling: 5*N^4 + 4*N^5
 *
 * In this implementation of the AO to MO transformation, we make use of the fast matrix multiplication of Eigen.
 * The price to pay here are 5 super-matrix transpositions, which however scale with N^4.
 *
 * TODO: access elements via pointer.
 *
 */
class AO2MOHelper {
 private:
  Eigen::MatrixXd& AO;
  Eigen::MatrixXd& C1;
  Eigen::MatrixXd& C2;
  long ao1 = 0L;
  long ao1sq = 0L;
  long ao2 = 0L;
  long mo1 = 0L;
  long mo2 = 0L;
  long ao1mo2 = 0L;
  long ao1ao2 = 0L;
  long mo1mo2 = 0L;
  long mo2sq = 0L;
  Eigen::MatrixXd MO;

  auto transform4() -> Eigen::MatrixXd;
  auto transform3(Eigen::MatrixXd&& AO3MO1) -> Eigen::MatrixXd;
  auto transform2(Eigen::MatrixXd&& AO2MO2) -> Eigen::MatrixXd;
  auto transform1(Eigen::MatrixXd&& AO1MO3) -> Eigen::MatrixXd;
  auto transpose4() -> Eigen::MatrixXd;
  auto transpose3(Eigen::MatrixXd& AO3MO1) const -> Eigen::MatrixXd;
  auto transpose2(Eigen::MatrixXd& AO2MO2) const -> Eigen::MatrixXd;
  auto transpose1(Eigen::MatrixXd& AO1MO3) const -> Eigen::MatrixXd;
  auto finalize(Eigen::MatrixXd&& MOTMP) -> void;

 public:
  AO2MOHelper(Eigen::MatrixXd& ao, Eigen::MatrixXd& c1, Eigen::MatrixXd& c2)
    : AO(ao),
      C1(c1),
      C2(c2),
      ao1(C1.rows()),
      ao1sq(ao1 * ao1),
      ao2(C2.rows()),
      mo1(C1.cols()),
      mo2(C2.cols()),
      ao1mo2(ao1 * mo2),
      mo1mo2(mo1 * mo2) {
    ao1ao2 = ao1 * ao2;
    mo2sq = mo2 * mo2;
  }

  auto perform() -> void;

  auto getResult() -> Eigen::MatrixXd;
};

auto AO2MOHelper::perform() -> void {
  finalize(transform1(transform2(transform3(transform4()))));
  // MO = transform1(transform2(transform3(transform4())));
}

auto AO2MOHelper::transform4() -> Eigen::MatrixXd {
  Eigen::MatrixXd AO3MO1 = Eigen::MatrixXd::Zero(ao1 * ao1 * ao2, mo2);

  Eigen::MatrixXd TMP1 = transpose4();

  // (mu nu lambda | sigma) C_sigma,l  --> (mu nu lambda | l)
  AO3MO1 = TMP1 * C2;

  return AO3MO1;
}

auto AO2MOHelper::transpose4() -> Eigen::MatrixXd {
  Eigen::MatrixXd TMP1 = Eigen::MatrixXd::Zero(ao1 * ao1 * ao2, ao2);

  // 1. ( mu nu | lambda sigma ) -> ( mu nu lambda | sigma )
  for (auto mu = 0L; mu < ao1; ++mu) {
    for (auto nu = 0L; nu < ao1; ++nu) {
      for (auto lambda = 0L; lambda < ao2; ++lambda) {
        for (auto sigma = 0L; sigma < ao2; ++sigma) {
          TMP1(ao1ao2 * mu + ao2 * nu + lambda, sigma) = AO(ao1 * mu + nu, ao2 * lambda + sigma);
        }
      }
    }
  }

  return TMP1;
}

auto AO2MOHelper::transform3(Eigen::MatrixXd&& AO3MO1) -> Eigen::MatrixXd {
  Eigen::MatrixXd AO2MO2 = Eigen::MatrixXd::Zero(ao1 * ao1 * mo2, mo2);

  Eigen::MatrixXd TMP2 = transpose3(AO3MO1);

  // (mu nu l | lambda ) C_lambda,k --> ( mu nu l | k )
  AO2MO2 = TMP2 * C2;
  return AO2MO2;
}

auto AO2MOHelper::transpose3(Eigen::MatrixXd& AO3MO1) const -> Eigen::MatrixXd {
  Eigen::MatrixXd TMP2 = Eigen::MatrixXd::Zero(ao1 * ao1 * mo2, ao2);
  // 2. ( mu nu lambda | l ) --> (mu nu l | lambda )
  for (auto mu = 0L; mu < ao1; ++mu) {
    for (auto nu = 0L; nu < ao1; ++nu) {
      for (auto l = 0L; l < mo2; ++l) {
        for (auto lambda = 0L; lambda < ao2; ++lambda) {
          TMP2(ao1mo2 * mu + mo2 * nu + l, lambda) = AO3MO1(ao1ao2 * mu + ao2 * nu + lambda, l);
        }
      }
    }
  }

  return TMP2;
}

auto AO2MOHelper::transform2(Eigen::MatrixXd&& AO2MO2) -> Eigen::MatrixXd {
  Eigen::MatrixXd AO1MO3 = Eigen::MatrixXd::Zero(ao1 * mo2 * mo2, mo1);

  Eigen::MatrixXd TMP3 = transpose2(AO2MO2);

  // (mu k l | nu ) C_nu,j --> (mu k l | j )
  AO1MO3 = TMP3 * C1;

  return AO1MO3;
}

auto AO2MOHelper::transpose2(Eigen::MatrixXd& AO2MO2) const -> Eigen::MatrixXd {
  Eigen::MatrixXd TMP3 = Eigen::MatrixXd::Zero(ao1 * mo2 * mo2, ao1);

  // 3. ( mu nu l | k ) --> (mu k l | nu )
  for (auto mu = 0L; mu < ao1; ++mu) {
    for (auto nu = 0L; nu < ao1; ++nu) {
      for (auto k = 0L; k < mo2; ++k) {
        for (auto l = 0L; l < mo2; ++l) {
          TMP3(mo2sq * mu + mo2 * k + l, nu) = AO2MO2(ao1mo2 * mu + mo2 * nu + l, k);
        }
      }
    }
  }

  return TMP3;
}

auto AO2MOHelper::transform1(Eigen::MatrixXd&& AO1MO3) -> Eigen::MatrixXd {
  Eigen::MatrixXd TMPMO = Eigen::MatrixXd::Zero(mo1 * mo2 * mo2, mo1);
  Eigen::MatrixXd TMP4 = transpose1(AO1MO3);

  // (j k l | mu ) C_mu,i --> (j k l | i )
  TMPMO = TMP4 * C1;

  return TMPMO;
}

auto AO2MOHelper::transpose1(Eigen::MatrixXd& AO1MO3) const -> Eigen::MatrixXd {
  Eigen::MatrixXd TMP4 = Eigen::MatrixXd::Zero(mo1 * mo2 * mo2, ao1);

  // 4. ( mu k l | j ) --> (j k l | mu )
  for (auto mu = 0L; mu < ao1; ++mu) {
    for (auto j = 0L; j < mo1; ++j) {
      for (auto k = 0L; k < mo2; ++k) {
        for (auto l = 0L; l < mo2; ++l) {
          TMP4(mo2sq * j + mo2 * k + l, mu) = AO1MO3(mo2sq * mu + mo2 * k + l, j);
        }
      }
    }
  }

  return TMP4;
}

auto AO2MOHelper::finalize(Eigen::MatrixXd&& MOTMP) -> void {
  MO.resize(mo1 * mo1, mo2 * mo2);

  // 5. (j k l | i) --> ( i j | k l)
  for (auto i = 0L; i < mo1; ++i) {
    for (auto j = 0L; j < mo1; ++j) {
      for (auto k = 0L; k < mo2; ++k) {
        for (auto l = 0L; l < mo2; ++l) {
          MO(mo1 * i + j, mo2 * k + l) = MOTMP(mo2sq * j + mo2 * k + l, i);
        }
      }
    }
  }
}

// auto AO2MOHelper::transform4() -> Eigen::MatrixXd {
//
//  Eigen::MatrixXd TMP1 = Eigen::MatrixXd::Zero(ao1 * ao1, ao2 * mo2);
//
//  for (auto mu = 0L; mu < ao1; ++mu) {
//    for (auto nu = 0L; nu < ao1; ++nu) {
//      for (auto lambda = 0L; lambda < ao2; ++lambda) {
//        for (auto l = 0L; l < mo2; ++l) {
//          for (auto sigma = 0L; sigma < ao2; ++sigma) {
//            TMP1(ao1 * mu + nu, lambda*ao2 + l) += AO(ao1 * mu + nu, ao2 * lambda + sigma) * C2(sigma, l);
//          }
//        }
//      }
//    }
//  }
//
//  return TMP1;
//}
//
// auto AO2MOHelper::transform3(Eigen::MatrixXd&& TMP1) -> Eigen::MatrixXd {
//
//  Eigen::MatrixXd TMP2 = Eigen::MatrixXd::Zero(ao1 * ao1, mo2 * mo2);
//
//  for (auto mu = 0L; mu < ao1; ++mu) {
//    for (auto nu = 0L; nu < ao1; ++nu) {
//      for (auto k = 0L; k < mo2; ++k) {
//        for (auto l = 0L; l < mo2; ++l) {
//          for (auto lambda = 0L; lambda < ao2; ++lambda) {
//            TMP2(ao1 * mu + nu, k*mo2 + l) += TMP1(ao1 * mu + nu, ao2 * lambda + l) * C2(lambda, k);
//          }
//        }
//      }
//    }
//  }
//
//  return TMP2;
//}
//
// auto AO2MOHelper::transform2(Eigen::MatrixXd&& TMP2) -> Eigen::MatrixXd {
//
//  Eigen::MatrixXd TMP3 = Eigen::MatrixXd::Zero(ao1 * mo1, mo2 * mo2);
//
//  for (auto mu = 0L; mu < ao1; ++mu) {
//    for (auto j = 0L; j < mo1; ++j) {
//      for (auto k = 0L; k < mo2; ++k) {
//        for (auto l = 0L; l < mo2; ++l) {
//          for (auto nu = 0L; nu < ao1; ++nu) {
//            TMP3(ao1 * mu + j, k*mo2 + l) += TMP2(ao1 * mu + nu, mo2 * k + l) * C1(nu, j);
//          }
//        }
//      }
//    }
//  }
//
//  return TMP3;
//}
//
// auto AO2MOHelper::transform1(Eigen::MatrixXd&& TMP3) -> Eigen::MatrixXd {
//
//  Eigen::MatrixXd RES = Eigen::MatrixXd::Zero(mo1 * mo1, mo2 * mo2);
//
//  for (auto i = 0L; i < mo1; ++i) {
//    for (auto j = 0L; j < mo1; ++j) {
//      for (auto k = 0L; k < mo2; ++k) {
//        for (auto l = 0L; l < mo2; ++l) {
//          for (auto mu = 0L; mu < ao1; ++mu) {
//            RES(mo1 * i + j, k*mo2 + l) += TMP3(ao1 * mu + j, mo2 * k + l) * C1(mu, i);
//          }
//        }
//      }
//    }
//  }
//
//  return RES;
//}

auto AO2MOHelper::getResult() -> Eigen::MatrixXd {
  return std::move(MO);
}

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_AO2MOHELPER_H
