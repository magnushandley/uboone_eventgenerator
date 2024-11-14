#ifndef UBCORE_EVENTGENERATOR_HEAVYNEUTRALLEPTONGEN_GENKINEMATICS_H
#define UBCORE_EVENTGENERATOR_HEAVYNEUTRALLEPTONGEN_GENKINEMATICS_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "fhiclcpp/ParameterSet.h"
#include "TLorentzVector.h"
#include <map>

namespace CLHEP {
  class HepRandomEngine;
}

namespace hpsgen {
  struct Geo {
    explicit Geo(art::ServiceHandle<geo::Geometry const>& g);
    const TVector3& centre() const { return det_centre; }
    const TVector3& half_dims() const { return det_half_dims; }
    const double max_path_length() const { return max_path_len; }
    private:
    TVector3 det_centre;
    TVector3 det_half_dims;
    double max_path_len;
  };
  struct PhysicalConstants {
    explicit PhysicalConstants(fhicl::ParameterSet const& p);
    double mass_elec() const { return p_mass_elec; }
    double mass_muon() const { return p_mass_muon; }
    double mass_pion_pm() const { return p_mass_pion_pm; }
    double mass_pion_0() const { return p_mass_pion_0; }
    double speed_light() const { return p_speed_light; }
    double higgs_vev() const { return p_higgs_vev; }
    double hbar() const { return p_hbar; }
    double ckm_ts() const { return p_ckm_ts; }
    double ckm_td() const { return p_ckm_td; }
    double ckm_ud() const { return p_ckm_ud; }
    double decayconst_pion() const { return p_decayconst_pion; }
    double lifetime_kaon_0() const { return p_lifetime_kaon_0; }
    double lifetime_kaon_pm() const { return p_lifetime_kaon_pm; }
    double mass_top() const { return p_mass_top; }
    double gfermi() const { return p_gfermi; }
    private:
      const double p_mass_elec; // GeV
      const double p_mass_muon; // GeV
      const double p_mass_pion_pm; // GeV
      const double p_mass_pion_0; // GeV
      const double p_mass_top; // GeV
      const double p_lifetime_kaon_0; // ns
      const double p_lifetime_kaon_pm; // ns
      const double p_speed_light; // cm/ns
      const double p_higgs_vev; // GeV
      const double p_hbar; // GeV ns
      const double p_ckm_td;
      const double p_ckm_ts;
      const double p_ckm_ud;
      const double p_decayconst_pion;
      const double p_gfermi; 
  };
  class GenKinematics {
    typedef CLHEP::HepRandomEngine rng;
    public:
      explicit GenKinematics(fhicl::ParameterSet const& p);
      ~GenKinematics();



      // bool generate(const TLorentzVector& kaon_decay_pos, const TLorentzVector& kaon_4mom, const double HNL_mass,
      //     const double model_theta, const int pion_type, const double flux_weight, const double max_weight, rng& rand, std::multimap<int,TLorentzVector>& result) const;

      bool generate(const TLorentzVector& kaon_decay_pos, const TLorentzVector& kaon_4mom, const int kaon_pdg, const double HNL_mass,
           const int fProdLepType, const int fDecayLepType ,const double flux_weight, const double max_weight, rng& rand, std::multimap<int,TLorentzVector>& result) const;

      bool generate_uniform(const TLorentzVector& kaon_decay_pos, const TLorentzVector& kaon_4mom, const int kaon_pdg,const double HNL_mass,
          rng& rand, std::multimap<int,TLorentzVector>& result) const;
      const PhysicalConstants& get_constants() const { return consts; }
      void update_geometry(art::ServiceHandle<geo::Geometry const>& g);
    private:

      TLorentzVector gen_random_scalar_mom(const double HNL_mass, const double kaon_mass,
          const double lep_mass, const TLorentzVector& kaon4mom, rng& rand) const;

      bool intersects_ray(const TVector3& orig, const TVector3& dir, double* lambdas) const;

      TLorentzVector gen_random_scalar_decay_pos(const TLorentzVector& scalar4mom, const TLorentzVector& kaonpos,
          const double model_tau, rng& rand, const double* lambdas, double& weight) const;
      
      TLorentzVector gen_random_scalar_decay_pos_uniform(const TLorentzVector& scalar4mom, const TLorentzVector& kaonpos,
          rng& rand, const double* lambdas, double& weight) const;

      bool pos_inside_detector(const TLorentzVector& scalar_dk_pos) const;

      double twobodylep_decay_width(const double HNL_mass, const double lep_mass) const;

      std::multimap<int,TLorentzVector> gen_daughters(const int kaon_pdg,const TLorentzVector& parent_mom, const int decay_lep_type , rng& rand) const;

 
      double kfactor(const double meson_mass,const double lep_mass, const double HNL_mass) const;


      PhysicalConstants consts;
      Geo *geo;
  };
}

#endif // UBCORE_EVENTGENERATOR_HEAVYNEUTRALLEPTONGEN_GENKINEMATICS_H
