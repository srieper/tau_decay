#ifndef DELTAR__STUFF
#define DELTAR__STUFF
#include "Math/Vector4D.h"
#include "TLorentzVector.h"

namespace deltaR {
  const double pi=3.14159265358979323846264338328;

  template <class Vector1, class Vector2>
  inline double DeltaPhi( const Vector1 & v1, const Vector2 & v2) {
    double dphi = v2.Phi() - v1.Phi();
    if ( dphi > pi ) {
      dphi -= 2.0*pi;
    } else if ( dphi <= -pi ) {
      dphi += 2.0*pi;
    }
    return dphi;
  }

  template <class Vector1, class Vector2>
  inline double DeltaR2( const Vector1 & v1, const Vector2 & v2) {
    double dphi = DeltaPhi(v1,v2);
    double deta = v2.Eta() - v1.Eta();
    return dphi*dphi + deta*deta;
  }

  template <class Vector1, class Vector2>
  inline double DeltaR( const Vector1 & v1, const Vector2 & v2) {
    return std::sqrt( DeltaR2(v1,v2) );
  }
}
#endif
