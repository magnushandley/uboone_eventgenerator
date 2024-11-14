#ifndef UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADER_H
#define UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADER_H

#include "fhiclcpp/ParameterSet.h"
#include "TLorentzVector.h"
#include <vector>

class TChain;

namespace simb {
  class MCFlux;
}

namespace bsim {
  class Dk2Nu;
}
namespace CLHEP {
  class HepRandomEngine;
}



namespace hpsgen {
  class FluxReader {
    public:
      FluxReader(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& fRNG, const char* tree_name,
          const char* pot_branch_name, const char* pot_tree);
      virtual ~FluxReader();
      virtual double get_kaon(TLorentzVector& kmom, TLorentzVector& kpos, int& kpdg, int& pi_type);
      double POTSeen(const double weight = 1.) const;
      virtual void get_MCFlux(simb::MCFlux& flux) = 0;
    protected:
      void set_branch(void* obj, const char* bname);
      void get_current_entry();
      virtual bool get_kaon_from_flux(TLorentzVector& kmom, TLorentzVector& kpos, int& kpdg, int& pi_type, double& weight) = 0;
    private:
      void get_next_entry();
      void count_pot(const char* bname, const char* tname);
      const std::string fFluxLocation;
      const std::string fFluxAccessSchema;
      bool fUseCache;
      TChain *fFluxTree;
      long fCurrEntry;
      std::vector<double> fPOTperFluxFile;
      double fTotalPOTinFluxFiles;
      std::vector<unsigned long> fEntriesPerFluxFile;
      unsigned int fFullChainsSeen;
      unsigned int fMaxFluxFileMB; // for IFDH copies

      struct cache {
        TLorentzVector kmom;
        TLorentzVector kpos;
        int kpdg;
        int pi_type;
        double weight;
      };
      std::vector<cache> fCache;
      std::vector<long> fSelectedEntries;
      std::vector<long>::iterator fCurrSelectedEntry;
  };
}

#endif // UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_FLUXREADER_H
