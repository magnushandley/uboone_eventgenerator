#ifndef UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADERBNB_H
#define UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADERBNB_H

#include "FluxReader.h"

#include "TVector3.h"
#include "TRotation.h"

namespace simb {
  class MCFlux;
}

namespace hpsgen {
  class FluxReaderBNB : public FluxReader {
    public:
      FluxReaderBNB(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& fRNG);
      virtual ~FluxReaderBNB();
      void get_MCFlux(simb::MCFlux& flux);
    private: 
      bool get_kaon_from_flux(TLorentzVector& kmom, TLorentzVector& kpos, int& kpdg, int& pi_type, double& weight);
      double K2nu_branching_ratio(const int id, const int ntp, const double parent_mass, const double Necm);
      const TVector3 beampos;
      const double beamtime;
      const TRotation rot;
      
      struct branch {
        float beamwgt;
        int ntp;
        int npart;
        int id[20];
        float ini_pos[20][3];
        float ini_mom[20][3];
        float ini_eng[20];
        float ini_t[20];
        float fin_mom[20][3];
        float fin_pol[20][3];
      } fBranch;
  };
}

#endif // UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADERBNB_H
