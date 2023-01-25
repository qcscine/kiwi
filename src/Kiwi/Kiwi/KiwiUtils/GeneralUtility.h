/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_GENERALUTILITY_H
#define KIWI_GENERALUTILITY_H

#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace Scine {
namespace Kiwi {

class Clock {
 public:
  /**
   * Deleted functions for the singleton pattern.
   * @{
   */
  Clock(const Clock&) = delete;
  Clock& operator=(const Clock&) = delete;
  //! @}
 private:
  using clock_time = std::chrono::high_resolution_clock::time_point;
  using chrono_sec = std::chrono::duration<double>;
  using chrono_milli = std::chrono::duration<double, std::milli>;

  std::vector<std::chrono::high_resolution_clock::time_point> uClock;

  std::vector<std::string> uTimestring;

  /**
   * @brief private constructor for singleton pattern
   */
  Clock(){};

 public:
  /**
   * @brief Getter for an instance of the Clock singleton.
   * The private constructor must initialize the Clock environment.
   */
  static Clock& getInstance() {
    static Clock clock;
    return clock;
  }

  inline void time(const std::string& name, bool reset = false, bool unit = true) {
    bool first = true;
    size_t index = 0;

    for (size_t i = 0; i < uTimestring.size(); i++) {
      if (name == uTimestring.at(i)) {
        if (reset) {
          uTimestring.erase(uTimestring.begin() + i);
        }
        first = false;
        index = i;
      }
    }
    if (first) {
      std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
      uClock.push_back(now);
      uTimestring.push_back(name);
    }
    else {
      std::chrono::duration<double, std::milli> res = std::chrono::high_resolution_clock::now() - uClock.at(index);
      std::cout << std::fixed << std::setprecision(3) << std::setw(16) << res.count() * 1.0e-3;
      if (unit)
        std::cout << " sec";
      std::cout << std::endl;
      uClock.erase(uClock.begin() + index);
      uTimestring.erase(uTimestring.begin() + index);
    }
  }

  static inline clock_time clock_now() {
    return std::chrono::high_resolution_clock::now();
  }
};

inline auto printAllowedKeywords(const std::vector<std::string>& vec) -> void {
  std::cout << "Allowed keywords are:\n";
  for (auto const& elem : vec) {
    std::cout << "\t" << elem << "\n";
  }
}

[[maybe_unused]] static Eigen::MatrixXd csvFileReader(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.good())
    throw std::runtime_error("CSV file " + filename + " not found.");
  std::vector<std::vector<double>> stdmat;
  std::string line;

  // Iterate through each line and split the content using the delimiter ",".
  while (getline(file, line)) {
    std::vector<std::string> vec_line;
    std::vector<double> vec_row;
    boost::algorithm::split(vec_line, line, boost::is_any_of(","));
    for (std::string const& elem : vec_line) {
      try {
        vec_row.push_back(boost::lexical_cast<double>(elem));
      }
      catch (const boost::bad_lexical_cast&) {
        throw std::runtime_error("Bad lexical cast while reading CSV file " + filename + ".");
      }
    }
    stdmat.push_back(vec_row);
  }

  file.close();

  int const N = stdmat.size();
  int const M = stdmat.begin()->size();

  Eigen::MatrixXd eigmat = Eigen::MatrixXd::Constant(N, M, 0.0);

  for (auto i = 0; i < N; ++i)
    for (auto j = 0; j < M; ++j) {
      eigmat(i, j) = stdmat[i][j];
    }
  return eigmat;
}

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_GENERALUTILITY_H
