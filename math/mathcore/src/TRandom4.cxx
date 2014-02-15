// @(#)root/mathcore:$Id: TRandom4.cxx $
// Author: Pierre Schweitzer 20/09/2012

//////////////////////////////////////////////////////////////////////////
//
// TRandom4
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

#include "TRandom4.h"

#include "Random123.h"


ClassImp(TRandom4)

//______________________________________________________________________________
TRandom4::TRandom4(UInt_t seed) 
{
//*-*-*-*-*-*-*-*-*-*-*default constructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ===================

   fRandom123 = new ROOT::Math::Random123();
   SetName("Random4");
   SetTitle("Random number generator: Threefry");
   SetSeed(seed);
}

//______________________________________________________________________________
TRandom4::~TRandom4()
{
//*-*-*-*-*-*-*-*-*-*-*default destructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ==================
   if (fRandom123) delete fRandom123; 
}

//______________________________________________________________________________
Double_t TRandom4::Rndm(Int_t)
{
   return (*fRandom123)();
}

//______________________________________________________________________________
void TRandom4::RndmArray(Int_t n, Double_t *array)
{
  // Return an array of n random numbers uniformly distributed in ]0,1]

   //for(Int_t i=0; i<n; i++) array[i]=(*fRandom123)();
   fRandom123->GenerateArray(array, array+n);
   //for(Int_t i=0; i<n; i++) array[i]=(*fRandom123)();
}

//______________________________________________________________________________
void TRandom4::RndmArray(Int_t n, Float_t *array)
{
  // Return an array of n random numbers uniformly distributed in ]0,1]

   fRandom123->GenerateArray(array, array+n);
//  for(Int_t i=0; i<n; i++) array[i]=(Float_t)(*fRandom123)();
}

//______________________________________________________________________________
UInt_t TRandom4::GetSeed() const
{
   return fRandom123->GetSeed();
}

//______________________________________________________________________________
void TRandom4::SetSeed(UInt_t seed)
{
   fRandom123->SetSeed(seed);
}
