/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Ao2MoTask.h"
#include "Keywords.h"
#include <Kiwi/KiwiUtils/AO2MO.h>
#include <Utils/IO/Yaml.h>
#include <yaml-cpp/yaml.h>

using namespace Scine;
using namespace Kiwi;

Ao2MoTask::Ao2MoTask(YAML::Node& input) {
  settings(input);
}

auto Ao2MoTask::name() const -> const std::string {
  return "Ao2Mo";
}

auto Ao2MoTask::settings(YAML::Node& input) -> void {
  Scine::Utils::checkYamlKeyRecognition(input, Scine::Kiwi::Keywords::ao2moSettings);

  if (input["print keywords"]) {
    auto printKeywords = input["print keywords"].as<bool>();
    if (printKeywords) {
      printAllowedKeywords(Scine::Kiwi::Keywords::ao2moSettings);
    }
  }
  if (input["write"]) {
    _write = input["write"].as<bool>();
  }
  if (input["thresh"]) {
    _integralThresh = input["thresh"].as<double>();
  }
}

auto Ao2MoTask::run() -> void {
  if (this->data->integralDirect) {
    std::cout << "\n\nAO2MO was called with integral-direct mode.\nIntegrals will be calculated and stored in memory "
                 "now.\n\n";
    this->data->integralDirect = false;
    this->data->twoBodyIntegrals();
    this->data->makeExchange();
  }

  AO2MO ao2mo(this->data, true);

  ao2mo.perform();

  if (_write) {
    ao2mo.write(_integralThresh);
  }
}
