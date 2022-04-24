// Author: Yipeng Sun, Svede Braun
// License: BSD 2-clause
// Last Change: Sun Apr 24, 2022 at 05:29 PM -0400

#pragma once

#include <string>

using std::string;

// We copy the string when passing it to the function, so the original one
// remains unchanged
string capitalize(string str) {
  for (auto& s : str) {
    s = toupper(s);
    break;
  }
  return str;
}
