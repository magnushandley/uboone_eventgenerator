#ifndef UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADERNUMI_H
#define UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADERNUMI_H

#include "FluxReader.h"

#include "TVector3.h"
#include "TRotation.h"

namespace bsim {
  class Dk2Nu;
}
namespace simb {
  class MCFlux;
}

namespace hpsgen {
  class FluxReaderNuMI : public FluxReader {
    public:
      FluxReaderNuMI(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& fRNG);
      virtual ~FluxReaderNuMI();
      bsim::Dk2Nu* get_dk2nu() const {return dk2nu; }
      void get_MCFlux(simb::MCFlux& flux);
    private:
      bool get_kaon_from_flux(TLorentzVector& kmom, TLorentzVector& kpos, int& kpdg, int& pi_type, double& weight);
      double K2nu_branching_ratio(const int ndecay);
      bsim::Dk2Nu *dk2nu;
      const TVector3 beampos;
      const double beamtime;
      const TRotation rot;
  };
}

#endif // UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADERNUMI_H
