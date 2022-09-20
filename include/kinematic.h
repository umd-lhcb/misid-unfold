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

PxPyPzEVector rebuildMu4Mom(PxPyPzEVector v4Mu, vector<double> smrFac) {
  auto vec = PxPyPzMVector{};

  vec.SetPx(v4Mu.Px() * smrFac[0]);
  vec.SetPy(v4Mu.Py() * smrFac[1]);
  vec.SetPz(v4Mu.Pz() * smrFac[2]);
  vec.SetM(v4Mu.M());

  return PxPyPzEVector(vec);
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
