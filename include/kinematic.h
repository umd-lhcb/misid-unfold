// Author: Yipeng Sun, Svede Braun
// License: BSD 2-clause
// Last Change: Sun Apr 24, 2022 at 08:51 PM -0400

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

template <typename T>
PxPyPzEVector rebuildMu4Mom(LorentzVector<T> v4Mu, vector<double> smrFac,
                            double m = PI_M) {
  auto vec = PxPyPzEVector{};

  vec.SetPx(v4Mu.Px() * smrFac[0]);
  vec.SetPy(v4Mu.Py() * smrFac[1]);
  vec.SetPz(v4Mu.Pz() * smrFac[2]);
  vec.SetM(m);

  return vec;
}

XYZVector buildBFlightDir(double endVtxX, double ownPvX, double endVtxY,
                          double ownPvY, double endVtxZ, double ownPvZ) {
  return XYZVector(endVtxX - ownPvX, endVtxY - ownPvY, endVtxZ - ownPvZ);
}

//////////////////////////////
// Rest frame approximation //
//////////////////////////////

template <typename T1, typename T2>
PxPyPzEVector estB4Mom(LorentzVector<T1>         v4BReco,
                       DisplacementVector3D<T2>& v3BFlight,
                       double                    mBRef = B_M) {
  auto mB  = v4BReco.M();
  auto pzB = v4BReco.Pz();

  auto cosX = v3BFlight.Unit().X();
  auto cosY = v3BFlight.Unit().Y();
  auto cosZ = v3BFlight.Unit().Z();

  Double_t pBMag = (mBRef / mB) * pzB / cosZ;
  return PxPyPzEVector(pBMag * cosX, pBMag * cosY, pBMag * cosZ,
                       TMath::Sqrt(pBMag * pBMag + mBRef * mBRef));
}

template <typename T>
Double_t m2Miss(LorentzVector<T> v4BEst, LorentzVector<T> v4BReco) {
  return (v4BEst - v4BReco).M2();
}

template <typename T>
Double_t el(LorentzVector<T> v4BEst, LorentzVector<T> v4Mu) {
  auto boost    = v4BEst.BoostToCM();
  auto v4MuRest = ROOT::Math::VectorUtil::boost(v4Mu, boost);
  return v4MuRest.E();
}

template <typename T>
Double_t q2(LorentzVector<T> v4BEst, LorentzVector<T> v4D) {
  return (v4BEst - v4D).M2();
}
