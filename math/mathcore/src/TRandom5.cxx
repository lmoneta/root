// @(#)root/mathcore:$Id: TRandom5.cxx $
// Author: Pierre Schweitzer 20/09/2012

//////////////////////////////////////////////////////////////////////////
//
// TRandom5
//
// Random number generator class based on the Threefry2x32 generator from
// John K. Salmon
//
// The period of the generator is 2**128.
//
// This implementation is only a wrapper around the real implemention.
//
// For more information see:
// Salmon, J., Moraes, M., Dror, R., and Shaw, D. (2011). Parallel random numbers:
// as easy as 1, 2, 3. In Proceedings of 2011 International Conference for High Per-
// formance Computing, Networking, Storage and Analysis, SC ’11, pages 16:1–16:12,
// New York, NY, USA. ACM.
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

   for(Int_t i=0; i<n; i++) array[i]=(*fRng)();
}

//______________________________________________________________________________
void TRandom5::RndmArray(Int_t n, Float_t *array)
{
  // Return an array of n random numbers uniformly distributed in ]0,1]

  for(Int_t i=0; i<n; i++) array[i]=(Float_t)(*fRng)();
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
