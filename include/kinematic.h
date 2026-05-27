// Author: Yipeng Sun, Svende Braun, Lucas Meyer Garcia
// License: BSD 2-clause
// Last Change: Tue Sep 20, 2022 at 03:26 AM -0400

#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

using std::string, std::vector;

using ROOT::Math::PxPyPzEVector;
using ROOT::Math::PxPyPzMVector;
using ROOT::Math::XYZVector;

//////////////////////
// Rebuild momentum //
//////////////////////

PxPyPzEVector rebuildMu4Mom(const PxPyPzEVector&  v4Mu,
                            const vector<double>& smrFac,
                            const string&         mode = "PThetaPhi") {
  if (mode == "PThetaPhi") {
    // Get variations
    const double& rp     = smrFac[3];
    const double& dtheta = smrFac[4];
    const double& dphi   = smrFac[5];
    // Compute smeared vector
    double p2_old = v4Mu.P2();
    double p      = std::sqrt(p2_old) * rp;
    double theta  = v4Mu.Theta() + dtheta;
    double phi    = v4Mu.Phi() + dphi;
    double pt     = p * std::sin(theta);
    double px     = pt * std::cos(phi);
    double py     = pt * std::sin(phi);
    double pz     = p * std::cos(theta);
    double e_old  = v4Mu.E();
    double e      = std::sqrt(e_old * e_old - p2_old + p * p);
    return PxPyPzEVector(px, py, pz, e);
  } else if (mode == "PxPyPz") {
    // Get variations
    double rx = smrFac[0];
    double ry = smrFac[1];
    double rz = smrFac[2];
    // Compute smeared vector
    auto vec = PxPyPzMVector{};
    vec.SetPx(v4Mu.Px() * rx);
    vec.SetPy(v4Mu.Py() * ry);
    vec.SetPz(v4Mu.Pz() * rz);
    vec.SetM(v4Mu.M());
    return PxPyPzEVector(vec);
  } else {
    throw std::runtime_error(
        "rebuildMu4Mom: Unexpected misid smearing mode. Options are PxPyPz and "
        "PThetaPhi (default).");
  }
}

XYZVector buildBFlightDir(double endVtxX, double ownPvX, double endVtxY,
                          double ownPvY, double endVtxZ, double ownPvZ) {
  return XYZVector(endVtxX - ownPvX, endVtxY - ownPvY, endVtxZ - ownPvZ);
}

//////////////////////////////
// Rest frame approximation //
//////////////////////////////

PxPyPzEVector estB4Mom(const PxPyPzEVector& v4BReco,
                       const XYZVector&     v3BFlight) {
  // Used in TupleToolTauMuDiscrVars for both B0 and B+
  constexpr double mBRef = 5279.61;

  auto mB  = v4BReco.M();
  auto pzB = v4BReco.Pz();

  auto   bFlightDir = v3BFlight.Unit();
  double cosX       = bFlightDir.X();
  double cosY       = bFlightDir.Y();
  double cosZ       = bFlightDir.Z();

  double pBMag = (mBRef * pzB) / (mB * cosZ);
  return PxPyPzEVector(pBMag * cosX, pBMag * cosY, pBMag * cosZ,
                       std::sqrt(pBMag * pBMag + mBRef * mBRef));
}

// all in GeV(^2)!
// also removed all template parameters because RDataFrame doesn't like them.
Double_t m2Miss(const PxPyPzEVector& v4BEst, const PxPyPzEVector& v4BReco) {
  return (v4BEst - v4BReco).M2() * 1e-6;
}

Double_t el(const PxPyPzEVector& v4BEst, const PxPyPzEVector& v4Mu) {
  auto boost    = v4BEst.BoostToCM();
  auto v4MuRest = ROOT::Math::VectorUtil::boost(v4Mu, boost);
  return v4MuRest.E() * 1e-3;
}

Double_t q2(const PxPyPzEVector& v4BEst, const PxPyPzEVector& v4D) {
  return (v4BEst - v4D).M2() * 1e-6;
}

// in MeV!
Double_t calcBM(const PxPyPzEVector& v4BReco) { return v4BReco.M(); }

Double_t eta(const PxPyPzEVector& v4) { return v4.Eta(); }

Double_t P(const PxPyPzEVector& v4) { return v4.P(); }
