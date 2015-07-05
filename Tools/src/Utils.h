#ifndef Utils_h__
#define Utils_h__

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <Math/VectorUtil.h>

typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;
typedef ROOT::Math::PtEtaPhiEVector PtEtaPhiEVector;


inline double squareDouble( const double & val )
{
  return val * val;
}

template<class T>
double MT( const T & lep, const T & neut )
{
  return sqrt(2*lep.pt*neut.pt*(1-cos(lep.phi-neut.phi)));
}


// inline double calcEoverP( const KFInParticle & part )
// {
//   /// |P| = Pt * cosh(eta)
//   return ( part.energy / ( part.pt * cosh( part.eta ) ) );

// }


//
/// theta star calculation... template for different usages
//
template<class T>
std::pair<double, int>
computeThetaStar( const T & topP4, const T & lepP4,
		  const T & bP4, const T & wP4  )
{
  // 1) Boost all particles to top rest frame
  ROOT::Math::XYZVector topRF = topP4.BoostToCM();

  T muTopRF = ROOT::Math::VectorUtil::boost(lepP4, topRF);
  T bTopRF  = ROOT::Math::VectorUtil::boost(bP4, topRF);
  T wTopRF  = ROOT::Math::VectorUtil::boost(wP4, topRF);

  // 2) Boost lepton & b to CM of W
  ROOT::Math::XYZVector wRF = wTopRF.BoostToCM();

  T muWRF = ROOT::Math::VectorUtil::boost( muTopRF, wRF );
  T bWRF  = ROOT::Math::VectorUtil::boost( -bTopRF, wRF );

  // 3) Angle between b and boosted lepton in top rest frame
  double theta1 = ROOT::Math::VectorUtil::Angle( muWRF, bWRF ) ;
  double theta  = ROOT::Math::VectorUtil::Angle( muWRF, wTopRF );

  int bugg = 0 ;
  if ( fabs( theta1 - theta ) > 1.e-04 ) {
    printf("\n mmmmmm Cos theta -b: %f  W: %f \n", cos(theta1), cos(theta) );
    bugg = 1;
  }

  return std::make_pair(theta, bugg);

}



template
std::pair<double, int>
computeThetaStar< ROOT::Math::PxPyPzEVector >
( const ROOT::Math::PxPyPzEVector & topP4,
  const ROOT::Math::PxPyPzEVector & lepP4,
  const ROOT::Math::PxPyPzEVector & bP4,
  const ROOT::Math::PxPyPzEVector & wP4 );

template
std::pair<double, int>
computeThetaStar< ROOT::Math::PtEtaPhiEVector >
( const ROOT::Math::PtEtaPhiEVector & topP4,
  const ROOT::Math::PtEtaPhiEVector & lepP4,
  const ROOT::Math::PtEtaPhiEVector & bP4,
  const ROOT::Math::PtEtaPhiEVector & wP4 );

#endif


