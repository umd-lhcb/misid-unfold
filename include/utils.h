// Author: Yipeng Sun, Svede Braun
// License: BSD 2-clause
// Last Change: Sun Apr 24, 2022 at 08:23 PM -0400

#pragma once

#include <filesystem>
#include <string>

#include <TTree.h>

using std::string;

/////////////
// General //
/////////////

// We copy the string when passing it to the function, so the original one
// remains unchanged
string capitalize(string str) {
  for (auto& s : str) {
    s = toupper(s);
    break;
  }
  return str;
}

string absDirPath(string pathRaw) {
  auto path    = std::filesystem::path(pathRaw);
  auto dirPath = path.parent_path();
  return std::filesystem::absolute(dirPath).string();
}

//////////////////
// ROOT-related //
//////////////////

bool branchExists(TTree* tree, string brName) {
  auto brPtr = tree->GetListOfBranches()->FindObject(brName.data());

  if (brPtr == nullptr) return false;
  return true;
}
