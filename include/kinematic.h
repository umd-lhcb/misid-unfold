// Author: Yipeng Sun, Svende Braun
// License: BSD 2-clause
// Last Change: Tue Sep 20, 2022 at 03:26 AM -0400

#pragma once

#include <vector>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TROOT.h>

using std::vector;

using ROOT::Math::LorentzVector;
using ROOT::Math::PtEtaPhiMVector;
using ROOT::Math::PxPyPzEVector;
using ROOT::Math::PxPyPzMVector;

using ROOT::Math::DisplacementVector3D;
using ROOT::Math::XYZVector;

#define K_M 493.677
#define PI_M 139.570
#define B_M 5279.34
#define B0_M 5279.65

//////////////////////
// Rebuild momentum //
//////////////////////

PxPyPzEVector rebuildMu4Mom(PxPyPzEVector v4Mu, vector<double> smrFac, TString mode = "PThetaPhi") {
  if (mode == "PxPyPz") {
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
  }
  else if (mode == "PThetaPhi") {
    // Get variations
    const double& rp = smrFac[3];
    const double& dtheta = smrFac[4];
    const double& dphi = smrFac[5];
    // Compute smeared vector
    double p = v4Mu.P() * rp;
    double theta = v4Mu.Theta() + dtheta;
    double pt = p * sin(theta);
    double pz = p * cos(theta);
    double eta = 0.5 * log((p+pz)/(p-pz));
    double phi = v4Mu.Phi() + dphi;
    auto vec = PtEtaPhiMVector( pt, eta, phi, v4Mu.M() );
    return PxPyPzEVector(vec);
  }
  else {
    throw std::runtime_error("rebuildMu4Mom: Unexpected misid smearing mode. Options are PxPyPz and PThetaPhi (default).");
  }
}

XYZVector buildBFlightDir(double endVtxX, double ownPvX, double endVtxY,
                          double ownPvY, double endVtxZ, double ownPvZ) {
  return XYZVector(endVtxX - ownPvX, endVtxY - ownPvY, endVtxZ - ownPvZ);
}

//////////////////////////////
// Rest frame approximation //
//////////////////////////////

PxPyPzEVector estB4Mom(PxPyPzEVector v4BReco, XYZVector v3BFlight,
                       double mBRef = B_M) {
  auto mB  = v4BReco.M();
  auto pzB = v4BReco.Pz();

  auto cosX = v3BFlight.Unit().X();
  auto cosY = v3BFlight.Unit().Y();
  auto cosZ = v3BFlight.Unit().Z();

  Double_t pBMag = (mBRef / mB) * pzB / cosZ;
  return PxPyPzEVector(pBMag * cosX, pBMag * cosY, pBMag * cosZ,
                       TMath::Sqrt(pBMag * pBMag + mBRef * mBRef));
}

// all in GeV(^2)!
// also removed all template parameters because RDataFrame doesn't like them.
Double_t m2Miss(PxPyPzEVector& v4BEst, PxPyPzEVector& v4BReco) {
  return (v4BEst - v4BReco).M2() / 1000 / 1000;
}

Double_t el(PxPyPzEVector& v4BEst, PxPyPzEVector& v4Mu) {
  auto boost    = v4BEst.BoostToCM();
  auto v4MuRest = ROOT::Math::VectorUtil::boost(v4Mu, boost);
  return v4MuRest.E() / 1000;
}

Double_t q2(PxPyPzEVector& v4BEst, PxPyPzEVector& v4D) {
  return (v4BEst - v4D).M2() / 1000 / 1000;
}

// in MeV!
Double_t calcBM(PxPyPzEVector& v4BReco) { return v4BReco.M(); }
