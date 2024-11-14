#ifndef UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_GENKINEMATICS_H
#define UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_GENKINEMATICS_H

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
    private:
    TVector3 det_centre;
    TVector3 det_half_dims;
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
    double lifetime_kaon_0() const { return p_lifetime_kaon_0; }
    double lifetime_kaon_pm() const { return p_lifetime_kaon_pm; }
    double mass_top() const { return p_mass_top; }
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
  };
  class GenKinematics {
    typedef CLHEP::HepRandomEngine rng;
    public:
      explicit GenKinematics(fhicl::ParameterSet const& p);
      ~GenKinematics();
      bool generate(const TLorentzVector& kaon_decay_pos, const TLorentzVector& kaon_4mom, const double scalar_mass,
          const double model_theta, const int pion_type, const double flux_weight, const double max_weight, rng& rand, std::multimap<int,TLorentzVector>& result) const;
      const PhysicalConstants& get_constants() const { return consts; }
      void update_geometry(art::ServiceHandle<geo::Geometry const>& g);
    private:

      TLorentzVector gen_random_scalar_mom(const double scalar_mass, const double kaon_mass,
          const double pion_mass, const TLorentzVector& kaon4mom, rng& rand) const;

      bool intersects_ray(const TVector3& orig, const TVector3& dir, double* lambdas) const;

      TLorentzVector gen_random_scalar_decay_pos(const TLorentzVector& scalar4mom, const TLorentzVector& kaonpos,
          const double model_tau, rng& rand, const double* lambdas, double& weight) const;

      bool pos_inside_detector(const TLorentzVector& scalar_dk_pos) const;

      std::multimap<int,TLorentzVector> gen_daughters(const TLorentzVector& parent_mom, const double theta, rng& rand) const;

      double branching_ratio(const double kmass, const double scalar_mass, const double model_theta, const int ktype) const;

      double scalar_tau(const double scalar_mass, const double model_theta) const;
      double partial_gamma_lep(const double scalar_mass, const double lep_mass, const double model_theta) const;
      double partial_gamma_pi(const double scalar_mass, const double pi_mass, const double model_theta) const;

      PhysicalConstants consts;
      Geo *geo;
  };
}

#endif // UBCORE_EVENTGENERATOR_HIGGSPORTALSCALARGEN_GENKINEMATICS_H
