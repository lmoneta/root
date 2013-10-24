// @(#)root/mathcore:$Id: TRandom5.cxx $
// Author: Pierre Schweitzer 20/09/2012

//////////////////////////////////////////////////////////////////////////
//
// TRandom5
//
// Random number generator class based on the MIXMAX generator from
//  K. Savvidy
//
// The period of the generator is 10^4682 for N=256, and
//                                10^1597 for N=88
//
// This implementation is only a wrapper around the real implemention, see mixmax.cxx and mixmax.h
// The generator, in C, is available also at hepforge: http://mixmax.hepforge.org
//
//	Generator described in 
//  N.Z.Akopov, G.K.Savvidy and N.G.Ter-Arutyunian, Matrix Generator of Pseudorandom Numbers, 
//	J.Comput.Phys. 97, 573 (1991); 
//	Preprint EPI-867(18)-86, Yerevan Jun.1986;
//
//////////////////////////////////////////////////////////////////////////

#include "TRandom5.h"

#include "RandomMixMax.h"


ClassImp(TRandom5)

//______________________________________________________________________________
TRandom5::TRandom5(UInt_t seed) 
{
//*-*-*-*-*-*-*-*-*-*-*default constructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ===================

   fRng = new ROOT::Math::RandomMixMax();
   SetName("Random5");
   SetTitle("Random number generator: MixMax");
   SetSeed(seed);
}

//______________________________________________________________________________
TRandom5::~TRandom5()
{
//*-*-*-*-*-*-*-*-*-*-*default destructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ==================
   if (fRng) delete fRng; 
}

//______________________________________________________________________________
Double_t TRandom5::Rndm(Int_t)
{
   return (*fRng)();
}

//______________________________________________________________________________
void TRandom5::RndmArray(Int_t n, Double_t *array)
{
// Return an array of n random numbers uniformly distributed in ]0,1]
fRng->RndmArray(n,array);
//   for(Int_t i=0; i<n; i++) array[i]=(*fRng)();
}

//______________________________________________________________________________
void TRandom5::RndmArray(Int_t n, Float_t *array)
{
  // Return an array of n random numbers uniformly distributed in ]0,1]

  for(Int_t i=0; i<n; i++) array[i]=(Float_t)(fRng->GetSeed()*0x1p-32);
}

//______________________________________________________________________________
UInt_t TRandom5::GetSeed() const
{
   return fRng->GetSeed();
}

//______________________________________________________________________________
void TRandom5::SetSeed(UInt_t seed)
{
   fRng->SetSeed(seed);
}

void TRandom5::SetSeed64(ULong64_t seed)
{
   fRng->SetSeed64(seed);
}

//______________________________________________________________________________
void TRandom5::SeedUniqueStream(UInt_t clusterID, UInt_t machineID, UInt_t runID, UInt_t  streamID)
{
   fRng->SeedUniqueStream(clusterID,  machineID,  runID,   streamID);
}
