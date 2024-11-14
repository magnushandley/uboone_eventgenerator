#ifndef UBCORE_EVENTGENERATOR_HNLGEN_GENKINEMATICS_H
#define UBCORE_EVENTGENERATOR_HNLGEN_GENKINEMATICS_H

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "fhiclcpp/ParameterSet.h"
#include "TLorentzVector.h"
#include <map>

#undef save_neg_lams

#ifdef save_neg_lams
class TTree;
#endif

namespace CLHEP {
  class HepRandomEngine;
}

namespace hnlgen {
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
    inline double mass_elec() const { return p_mass_elec; }
    inline double mass_muon() const { return p_mass_muon; }
    inline double mass_pion_pm() const { return p_mass_pion_pm; }
    inline double mass_pion_0() const { return p_mass_pion_0; }
    inline double mass_kaon_pm() const { return p_mass_kaon_pm; }
    inline double mass_kaon_0() const { return p_mass_kaon_0; }
    inline double speed_light() const { return p_speed_light; }
    inline double higgs_vev() const { return p_higgs_vev; }
    inline double hbar() const { return p_hbar; }
    inline double ckm_ud() const { return p_ckm_ud; }
    inline double ckm_ts() const { return p_ckm_ts; }
    inline double ckm_td() const { return p_ckm_td; }
    inline double lifetime_kaon_0() const { return p_lifetime_kaon_0; }
    inline double lifetime_kaon_pm() const { return p_lifetime_kaon_pm; }
    inline double mass_top() const { return p_mass_top; }
    inline double gFermi2() const { return p_gFermi2 ; }
    inline double pion_decay_constant() const { return p_pion_decay_constant; }
    inline double sin_thW() const { return p_sin_thW; }
    private:
      const double p_mass_elec; // GeV
      const double p_mass_muon; // GeV
      const double p_mass_pion_pm; // GeV
      const double p_mass_pion_0; // GeV
      const double p_mass_kaon_pm; // GeV
      const double p_mass_kaon_0; // GeV
      const double p_mass_top; // GeV
      const double p_lifetime_kaon_0; // nsmlm2
      const double p_lifetime_kaon_pm; // ns
      const double p_speed_light; // cm/ns
      const double p_higgs_vev; // GeV
      const double p_hbar; // GeV ns
      const double p_ckm_ud;
      const double p_ckm_td;
      const double p_ckm_ts;
      const double p_gFermi2;
      const double p_pion_decay_constant;
      const double p_sin_thW;
  };
  struct ModelParameters {
    ModelParameters(fhicl::ParameterSet const& p);
    double hnl_mass, Ue4, Umu4, Utau4;
    enum class FermionNature { Dirac, Majorana } fermion_type;
    std::vector<std::string> modes;
    bool print_ratio;
    bool is_Majorana() const { return fermion_type==FermionNature::Majorana; }
    bool is_Dirac() const { return fermion_type==FermionNature::Dirac; }
    double sum_U2() const { return std::pow(Ue4,2) + std::pow(Umu4,2) + std::pow(Utau4,2); }
    bool simulate(const std::string& mode) const {
      return std::find(modes.begin(), modes.end(), mode) != modes.end();
    }
    void print() const;
  };
  double ____tmp;
  class GenKinematics {
    typedef CLHEP::HepRandomEngine rng;
    public:
      explicit GenKinematics(fhicl::ParameterSet const& p, rng& rand
#ifdef save_neg_lams
          ,TTree* tt_ = 0, double& t_pfn_ = ____tmp, double& t_pfp_ = ____tmp, double& t_mp_ = ____tmp, double& t_mn_ = ____tmp, double& t_snum_ = ____tmp, double& t_snup_ = ____tmp, double& t_cthm_ = ____tmp, double& t_cthp_ = ____tmp, double& t_klm_ = ____tmp, double& t_klp_ = ____tmp, double& t_integrand_ = ____tmp
#endif
          );
      ~GenKinematics();
      bool generate(const TLorentzVector& parent_decay_pos, const TLorentzVector& parent_4mom, const int parent_pdg, const double flux_weight, const double max_weight, rng& rand, std::multimap<int,TLorentzVector>& result) const;
      const PhysicalConstants& get_constants() const { return consts; }
      void update_geometry(art::ServiceHandle<geo::Geometry const>& g);
    private:

      std::pair<int,TLorentzVector> gen_random_hnl_mom(const int parent_pdg, const TLorentzVector& parent4mom, rng& rand) const;

      bool intersects_ray(const TVector3& orig, const TVector3& dir, double* lambdas) const;

      TLorentzVector gen_random_hnl_decay_pos(const TLorentzVector& hnl4mom, const TLorentzVector& kaonpos,
          const double model_tau, rng& rand, const double* lambdas, double& weight) const;

      bool pos_inside_detector(const TLorentzVector& hnl_dk_pos) const;

      std::multimap<int,TLorentzVector> gen_daughters(const TLorentzVector& parent_mom,const int parent_pdg, rng& rand) const;

      double parent_branching_ratio(const int parent_pdg) const;

      double hnl_tau(const int pdg) const;

      using decay_output_fn_t = std::function<std::vector<TLorentzVector>(rng& rand, const TLorentzVector&)>;
                                        // pdg
      using decay_prod_vec_t = std::vector<int>;

      //                               decay products  ,BR or gamma,  integrands
      using branch_info_t = std::tuple<decay_prod_vec_t, double, decay_output_fn_t>;
                                // (BR or Gamma)
      using decay_info_t = std::pair<double, std::vector<branch_info_t>>;
            // pdg         ,     
      std::map<int, decay_info_t> decay_branches;

      void calculate_branching_ratios(rng& rand);

      PhysicalConstants consts;
      Geo *geo;
      ModelParameters params;
      const size_t num_evals_numerical_integration;
      const bool use_s23_weighting;
      const bool gen_nu_lep_lep;
      const bool gen_nu_e_e;
      const bool gen_nu_e_mu;
      const bool gen_nu_mu_mu;
      const bool gen_nu_pi0;
#ifdef save_neg_lams
      TTree *tt;
      double& t_pfp;
      double& t_pfn;
      double& t_mp;
      double& t_mn;
      double& t_snum;
      double& t_snup;
      double& t_cthm;
      double& t_cthp;
      double& t_klm;
      double& t_klp;
      double& t_integrand;
#endif
  };
}

#endif // UBCORE_EVENTGENERATOR_HNLGEN_GENKINEMATICS_H
