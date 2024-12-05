#include "GenKinematics.h"
#include "EvtTimeFNALBeam.h"

#include "CLHEP/Random/RandFlat.h"

#ifdef save_neg_lams
#include <TTree.h>
#endif
//#include "CLHEP/Random/RandomEngine.h"

namespace pdgcodes {
  const int k_elec = 11;
  const int k_muon = 13;
  const int k_nu_e = 12;
  const int k_nu_mu = 14;
  const int k_nu_tau = 16;
  const int k_pion_0 = 111;
  const int k_pion_pm = 211;
  const int k_kaon_0 = 130;
  const int k_kaon_pm = 321;
  const int k_HNL_poshel = 91;
  const int k_HNL_neghel = 89;
}

// units are [GeV] [ns] [cm]
hnlgen::PhysicalConstants::PhysicalConstants(fhicl::ParameterSet const& p) 
: p_mass_elec(p.get<double>("mass_elec",0.5109989461e-3)) // PDG 2019
  ,p_mass_muon(p.get<double>("mass_muon",0.1056583745)) // PDG 2019
  ,p_mass_pion_pm(p.get<double>("mass_pion_pm",0.13957061)) // PDG 2019
  ,p_mass_pion_0(p.get<double>("mass_pion_0",0.1349770)) // PDG 2019
  ,p_mass_kaon_pm(p.get<double>("mass_kaon_pm",0.493677))
  ,p_mass_kaon_0(p.get<double>("mass_kaon_0",0.497611))
  ,p_mass_top(p.get<double>("mass_top",172.9)) // PDG 2019
  ,p_lifetime_kaon_0(p.get<double>("lifetime_kaon_0",51.16)) // PDG 2019
  ,p_lifetime_kaon_pm(p.get<double>("lifetime_kaon_pm",12.38)) // PDG 2019
  ,p_speed_light(p.get<double>("speed_light",29.9792458)) // PDG, exact
  ,p_higgs_vev(p.get<double>("higgs_vev",246.22)) // PDG
  ,p_hbar(p.get<double>("hbar",6.582119569e-16)) // PDG, exact
  ,p_ckm_ud(p.get<double>("ckm_ud",0.97370)) // PDG 2019
  ,p_ckm_td(p.get<double>("ckm_td",8.1e-3)) // PDG 2019
  ,p_ckm_ts(p.get<double>("ckm_ts",39.4e-3)) // PDG 2019
  ,p_gFermi2(p.get<double>("gFermi2",std::pow(1.1663787e-5,2)))
  ,p_pion_decay_constant(p.get<double>("pion_decay_constant",0.130))
  ,p_sin_thW(p.get<double>("sin_thW",0.480843)) // sin^2 => 0.23121
{
}

hnlgen::ModelParameters::ModelParameters(fhicl::ParameterSet const& p)
: hnl_mass{p.get<double>("model_hnl_mass")}
  ,Ue4{p.get<double>("model_U_e_4_angle")}
  ,Umu4{p.get<double>("model_U_mu_4_angle")}
  ,Utau4{p.get<double>("model_U_tau_4_angle")}
  ,fermion_type{[](auto const& s) -> FermionNature {
    if(s == "dirac") return FermionNature::Dirac;
    if(s == "majorana") return FermionNature::Majorana;
    throw cet::exception("Configuration")
      << "HNL fermionic nature '"<<s<<"' should be 'dirac' or 'majorana'" ;
  }(p.get<std::string>("model_hnl_fermion_nature"))}
  ,modes{p.get<std::vector<std::string>>("enabled_decay_modes")}
  ,print_ratio{p.get<bool>("print_ratio",false)}
{
  print();
}

void hnlgen::ModelParameters::ModelParameters::print() const {
  std::cout << Form("MODEL %s hnl mass %g Ue4 %g Umu4 %g Utau4 %g \n",fermion_type==FermionNature::Dirac?"dirac":"majorana",
      hnl_mass,Ue4,Umu4,Utau4);
}

hnlgen::Geo::Geo(art::ServiceHandle<geo::Geometry const>& g) {
  // Find boundary of active volume
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<g->NTPC(); ++i)
  {
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = g->TPC(i);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-g->DetHalfWidth(i))
      minx = world[0]-g->DetHalfWidth(i);
    if (maxx<world[0]+g->DetHalfWidth(i))
      maxx = world[0]+g->DetHalfWidth(i);
    if (miny>world[1]-g->DetHalfHeight(i))
      miny = world[1]-g->DetHalfHeight(i);
    if (maxy<world[1]+g->DetHalfHeight(i))
      maxy = world[1]+g->DetHalfHeight(i);
    if (minz>world[2]-g->DetLength(i)/2.)
      minz = world[2]-g->DetLength(i)/2.;
    if (maxz<world[2]+g->DetLength(i)/2.)
      maxz = world[2]+g->DetLength(i)/2.;
  }
  det_centre.SetXYZ(0.5*(minx+maxx),0.5*(miny+maxy),0.5*(minz+maxz));
  det_half_dims.SetXYZ(0.5*(maxx-minx),0.5*(maxy-miny),0.5*(maxz-minz));
  //std::cout << "Constructed detector geometry. Centre: "; det_centre.Print();
  //std::cout << " half-widths: "; det_half_dims.Print();
}

hnlgen::GenKinematics::GenKinematics(fhicl::ParameterSet const& p, rng& rand
#ifdef save_neg_lams
          ,TTree* tt_, double& t_pfn_, double& t_pfp_, double& t_mp_, double& t_mn_
          , double& t_snum_, double& t_snup_, double& t_cthm_, double& t_cthp_, double& t_klm_, double& t_klp_
          , double& t_integrand_
#endif
    ) 
  : consts(p), geo(0), params(p)
    ,num_evals_numerical_integration(p.get<size_t>("num_evals_numerical_integration",10000000))
    ,use_s23_weighting(p.get<bool>("use_s23_weighting",true))
    ,gen_nu_lep_lep(p.get<bool>("model_gen_nu_lep_lep",true))
    ,gen_nu_e_e(p.get<bool>("model_gen_nu_e_e",true))
    ,gen_nu_e_mu(p.get<bool>("model_gen_nu_e_mu",true))
    ,gen_nu_mu_mu(p.get<bool>("model_gen_nu_mu_mu",true))
    ,gen_nu_pi0(p.get<bool>("model_gen_nu_pi0",false))
#ifdef save_neg_lams
          ,tt(tt_), t_pfp(t_pfp_), t_pfn(t_pfn_), t_mp(t_mp_), t_mn(t_mn_)
  , t_snum(t_snum_), t_snup(t_snup_), t_cthm(t_cthm_), t_cthp(t_cthp_), t_klm(t_klm_), t_klp(t_klp_)
  ,t_integrand(t_integrand_)
#endif
{
  calculate_branching_ratios(rand);
}

hnlgen::GenKinematics::~GenKinematics() 
{
  if(geo) delete geo;
}


bool hnlgen::GenKinematics::generate(const TLorentzVector& parent_decay_pos, const TLorentzVector& parent_4mom,
    const int parent_pdg, const double flux_weight, const double max_weight,
          rng& rand, std::multimap<int,TLorentzVector>& result) const {
  
  auto hnl_pdg_4mom = gen_random_hnl_mom(parent_pdg, parent_4mom, rand);
  auto const& hnl_4mom = hnl_pdg_4mom.second;
  const int hnl_pdg = hnl_pdg_4mom.first;
  double lambdas[2];
  //Magnus test - REMOVE
  if(!intersects_ray(parent_decay_pos.Vect(), hnl_4mom.Vect(), lambdas)) return false;
  //lambdas[0] = 9746.;
  //lambdas[1] = 10765.;
  //std::cout<< "trial" << std::endl;
  double dk_weight = 0.;
  TLorentzVector hnl_dk_pos = gen_random_hnl_decay_pos(hnl_4mom, parent_decay_pos, hnl_tau(hnl_pdg),
      rand, lambdas, dk_weight);
  if(!pos_inside_detector(hnl_dk_pos)) return false;
  
  const double br_weight = parent_branching_ratio(parent_pdg);

  const double weight = dk_weight * br_weight * flux_weight;

  if(max_weight > 0.) {
    if(weight > max_weight) {
      throw art::Exception(art::errors::LogicError) << "weight "<<weight<<" > max_weight "<<max_weight<<std::endl
        <<" Modify max_weight (suggested "<<1.1*weight<<") and re-run."<<std::endl;
    }
    if(CLHEP::RandFlat::shoot(&rand, max_weight) > weight) return false;
  }

  result = gen_daughters(hnl_4mom, hnl_pdg, rand);
  // empty result if unmodelled (eg 3nu) decay mode is chosen
  if(result.empty()) {
    return false;
  }
  result.emplace(0,hnl_dk_pos);
  result.emplace(hnl_pdg, hnl_4mom);

  if(max_weight <= 0.) {
    result.emplace(99,TLorentzVector{dk_weight,br_weight,flux_weight,weight});
  }

  return true;
}

std::pair<int,TLorentzVector> hnlgen::GenKinematics::gen_random_hnl_mom(const int parent_pdg, const TLorentzVector& par4mom, rng& rand) const {
    auto decinfo = decay_branches.find(parent_pdg);
  if(decinfo == decay_branches.end()) {
    auto err =  art::Exception(art::errors::LogicError) << "Unable to find pdg "<<parent_pdg
      <<" in precomputed BRs [ valid pdgs: ";
    for(auto const& d : decay_branches) {
      err << d.first << ", ";
    }
    err << "]";
    throw err;    
  }
  const double br = decinfo->second.first;
  auto const& branches = decinfo->second.second;
  
  auto is_hnl = [](const int pdg) {
    return std::abs(pdg) == pdgcodes::k_HNL_poshel || std::abs(pdg) == pdgcodes::k_HNL_neghel;
  };
  
  size_t which_b = branches.size()-1;
  
  if(branches.size() > 1) {
    std::vector<double> brs;
    std::vector<double> partial_sums;
    std::transform(branches.begin(), branches.end(), std::back_inserter(brs), [](auto const& b) { return std::get<1>(b); });
    std::partial_sum(brs.begin(), brs.end(), std::back_inserter(partial_sums));
    const double u = CLHEP::RandFlat::shoot(&rand, br);
    for(size_t i = 0; i < branches.size(); ++i) {
      if(u < partial_sums[i]) {
        which_b = i;
        break;
      }
    }
  }

  auto const& pdgs = std::get<0>(branches[which_b]);
  auto dprods = std::get<2>(branches[which_b])(rand, par4mom);
  for(size_t i = 0; i < pdgs.size(); ++i) {
    if(is_hnl(pdgs[i])) {
      return std::make_pair(pdgs[i],dprods[i]);
    }
  }
  
  // shouldn't be here
  throw art::Exception(art::errors::LogicError) << "Did not produce HNL!";
  return {};

}

TLorentzVector hnlgen::GenKinematics::gen_random_hnl_decay_pos(const TLorentzVector& hnl4mom, const TLorentzVector& parpos,
    const double model_tau, rng& rand, const double* lambdas, double& weight) const {
  const double gamma = hnl4mom.Gamma();
  const double speed = hnl4mom.Beta() * consts.speed_light();
  const double a = std::min(lambdas[0],lambdas[1])/speed;
  const double b = std::max(lambdas[0],lambdas[1])/speed;
  const double probA =  std::exp(-a/gamma/model_tau);
  const double probB =  std::exp(-b/gamma/model_tau);
  weight = probA - probB; // integral along exponential;
  const double p0 = CLHEP::RandFlat::shoot(&rand);
  const double length = (a + b - gamma*model_tau*std::log(std::exp(b/gamma/model_tau)*(1-p0) + std::exp(a/gamma/model_tau)*p0))*speed;
  const TVector3 traj = length * hnl4mom.Vect().Unit();
  const double time_lab = traj.Mag() / speed;
  const TLorentzVector traj4(traj,time_lab);
  return parpos+traj4;
}

void hnlgen::GenKinematics::update_geometry(art::ServiceHandle<geo::Geometry const>& g) {
  if(geo) {
    delete geo;
  }
  geo = new Geo(g);
}

// there may be better ways of doing this, using geo::Geometry directly?
bool hnlgen::GenKinematics::intersects_ray(const TVector3& orig, const TVector3& dir, double* lambdas) const {
  //find detector center and half-dimensions
  const TVector3& det_centre = geo->centre(); const TVector3& det_half_dims = geo->half_dims();
  //unit vector in direction of hnl
  const TVector3& unit_dir = dir.Unit();
  int n_intersects = 0;
  unsigned int ilam = 0;
  for(int coord = 0; coord < 3; ++coord) {
    if(std::abs(unit_dir[coord])>0.) {
      for(int side = -1; side < 2; side += 2) { // side = -1 or +1
        const double plane = det_centre[coord] + side * det_half_dims[coord];
        const double lambda = (plane - orig[coord])/unit_dir[coord];
        if(lambda < 0) continue; // no backwards-going hnls
        bool intersects_planes[2] = {false, false};
        unsigned int iplane = 0;
        for(int other_coord = 0; other_coord < 3; ++other_coord) {
          if(other_coord == coord) continue;
          const double oth_plane = lambda * unit_dir[other_coord] + orig[other_coord];
          if(std::abs(oth_plane - det_centre[other_coord]) < det_half_dims[other_coord]) {
            intersects_planes[iplane]=true;
          }
          iplane++;
        }
        if(intersects_planes[0] && intersects_planes[1]) {
          n_intersects++;
          if(ilam < 2) {
            lambdas[ilam++] = lambda;
          }
        }
      }
    }
  }
  return n_intersects >= 2;
}


// there may be better ways of doing this, using geo::Geometry directly?
bool hnlgen::GenKinematics::pos_inside_detector(const TLorentzVector& hnl_dk_pos) const {
  const TVector3& det_centre = geo->centre(); const TVector3& det_half_dims = geo->half_dims();
  bool intersects_vol = true;
  for(int coord = 0; coord < 3; ++coord) {
    const double pos = hnl_dk_pos.Vect()[coord];
    if(std::abs(pos - det_centre[coord]) > det_half_dims[coord]) {
      intersects_vol = false;
      break;
    }
  }
  return intersects_vol;
}

std::multimap<int,TLorentzVector> hnlgen::GenKinematics::gen_daughters(const TLorentzVector& parent_mom, const int parent_pdg, rng& rand) const {
  auto decinfo = decay_branches.find(parent_pdg);
  if(decinfo == decay_branches.end()) {
    auto err =  art::Exception(art::errors::LogicError) << "Unable to find pdg "<<parent_pdg
      <<" in precomputed BRs [ valid pdgs: ";
    for(auto const& d : decay_branches) {
      err << d.first << ", ";
    }
    err << "]";
    throw err;    
  }
  const double br = decinfo->second.first;
  auto const& branches = decinfo->second.second;
  
  size_t which_b = branches.size()-1;
  
  std::multimap<int,TLorentzVector> ret;

  if(branches.size() > 0) {
    std::vector<double> brs;
    std::vector<double> partial_sums;
    std::transform(branches.begin(), branches.end(), std::back_inserter(brs), [](auto const& b) { return std::get<1>(b); });
    std::partial_sum(brs.begin(), brs.end(), std::back_inserter(partial_sums));
    const double u = CLHEP::RandFlat::shoot(&rand, br);
    for(size_t i = 0; i < branches.size(); ++i) {
      if(u < partial_sums[i]) {
        which_b = i;
        break;
      }
      if(i == branches.size()-1) {
        // u is larger than active BRs in decay_branches, meaning one of the 
        // unmodelled decay modes (eg 3nu) was actually selected
        // return empty multimap
        return ret;
      }
    }
  }


  auto const& pdgs = std::get<0>(branches[which_b]);
  auto dprods = std::get<2>(branches[which_b])(rand, parent_mom);

  for(size_t i = 0; i < pdgs.size(); ++i) {
    ret.emplace(pdgs[i], dprods[i]);
  }
  
  return ret;
  
}

double hnlgen::GenKinematics::hnl_tau(const int pdg) const {
  auto const hnl_info = decay_branches.find(pdg);
  if(hnl_info == decay_branches.end()) {
    auto err =  art::Exception(art::errors::LogicError) << "Unable to find pdg "<<pdg
      <<" in precomputed BRs [ valid pdgs: ";
    for(auto const& d : decay_branches) {
      err << d.first << ", ";
    }
    err << "]";
    throw err;
  }
  const double gamma = std::get<0>(hnl_info->second);
  return consts.hbar() / gamma;
}

double hnlgen::GenKinematics::parent_branching_ratio(const int parent_pdg) const {
  auto const par_info = decay_branches.find(parent_pdg);
  if(par_info == decay_branches.end()) {
    auto err =  art::Exception(art::errors::LogicError) << "Unable to find pdg "<<parent_pdg
      <<" in precomputed BRs [ valid pdgs: ";
    for(auto const& d : decay_branches) {
      err << d.first << ", ";
    }
    err << "]";
    throw err;
  }
  return std::get<0>(par_info->second);
}


namespace detail {
    template<class F>
    struct neg_fn_t {
        F f;
        template<class... Args>
        constexpr auto operator()(Args&&... args) &
            noexcept(noexcept(-std::invoke(f, std::forward<Args>(args)...)))
            -> decltype(-std::invoke(f, std::forward<Args>(args)...))
        {
            return -std::invoke(f, std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) const&
            noexcept(noexcept(-std::invoke(f, std::forward<Args>(args)...)))
            -> decltype(-std::invoke(f, std::forward<Args>(args)...))
        {
            return -std::invoke(f, std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) &&
            noexcept(noexcept(-std::invoke(std::move(f), std::forward<Args>(args)...)))
            -> decltype(-std::invoke(std::move(f), std::forward<Args>(args)...))
        {
            return -std::invoke(std::move(f), std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) const&&
            noexcept(noexcept(-std::invoke(std::move(f), std::forward<Args>(args)...)))
            -> decltype(-std::invoke(std::move(f), std::forward<Args>(args)...))
        {
            return -std::invoke(std::move(f), std::forward<Args>(args)...);
        }
    };
}
 
template<class F>
constexpr detail::neg_fn_t<std::decay_t<F>> negate_func(F&& f)
{
    return { std::forward<F>(f) };
}
 

namespace detail {
    template<class F1, class F2>
    struct add_fn_t {
        F1 f1;
        F2 f2;
        template<class... Args>
        constexpr auto operator()(Args&&... args) &
            noexcept(noexcept(std::invoke(f1, std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...)))
            -> decltype(std::invoke(f1, std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...))
        {
            return std::invoke(f1, std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) const&
            noexcept(noexcept(std::invoke(f1, std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...)))
            -> decltype(std::invoke(f1, std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...))
        {
            return std::invoke(f1, std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) &&
            noexcept(noexcept(std::invoke(std::move(f1), std::forward<Args>(args)...)
                  +std::invoke(std::move(f2), std::forward<Args>(args)...)))
            -> decltype(std::invoke(std::move(f1), std::forward<Args>(args)...)
                +std::invoke(std::move(f2), std::forward<Args>(args)...))
        {
            return std::invoke(std::move(f1), std::forward<Args>(args)...)+std::invoke(std::move(f2), std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) const&&
            noexcept(noexcept(std::invoke(std::move(f1), std::forward<Args>(args)...)
                  +std::invoke(std::move(f2), std::forward<Args>(args)...)))
            -> decltype(std::invoke(std::move(f1), std::forward<Args>(args)...)
                +std::invoke(std::move(f2), std::forward<Args>(args)...))
        {
            return std::invoke(std::move(f1), std::forward<Args>(args)...)+std::invoke(std::move(f2), std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) &&
            noexcept(noexcept(std::invoke(f1, std::forward<Args>(args)...)
                  +std::invoke(std::move(f2), std::forward<Args>(args)...)))
            -> decltype(std::invoke(f1, std::forward<Args>(args)...)
                +std::invoke(std::move(f2), std::forward<Args>(args)...))
        {
            return std::invoke(f1, std::forward<Args>(args)...)+std::invoke(std::move(f2), std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) const&&
            noexcept(noexcept(std::invoke(f1, std::forward<Args>(args)...)
                  +std::invoke(std::move(f2), std::forward<Args>(args)...)))
            -> decltype(std::invoke(f1, std::forward<Args>(args)...)
                +std::invoke(std::move(f2), std::forward<Args>(args)...))
        {
            return std::invoke(f1, std::forward<Args>(args)...)+std::invoke(std::move(f2), std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) &&
            noexcept(noexcept(std::invoke(std::move(f1), std::forward<Args>(args)...)
                  +std::invoke(f2, std::forward<Args>(args)...)))
            -> decltype(std::invoke(std::move(f1), std::forward<Args>(args)...)
                +std::invoke(f2, std::forward<Args>(args)...))
        {
            return std::invoke(std::move(f1), std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...);
        }
 
        template<class... Args>
        constexpr auto operator()(Args&&... args) const&&
            noexcept(noexcept(std::invoke(std::move(f1), std::forward<Args>(args)...)
                  +std::invoke(f2, std::forward<Args>(args)...)))
            -> decltype(std::invoke(std::move(f1), std::forward<Args>(args)...)
                +std::invoke(f2, std::forward<Args>(args)...))
        {
            return std::invoke(std::move(f1), std::forward<Args>(args)...)+std::invoke(f2, std::forward<Args>(args)...);
        }
    };
}

template<class F1, class F2>
constexpr detail::add_fn_t<std::decay_t<F1>,std::decay_t<F2>> add_funcs(F1&& f1, F2&& f2)
{
    return { std::forward<F1>(f1), std::forward<F2>(f2) };
}


void hnlgen::GenKinematics::calculate_branching_ratios(rng& rand) {
  /*
  auto sqrt_kallen = [](const double a, const double b, const double c) {
    return std::sqrt(std::pow(a-b-c,2)-4.*b*c);
  };
  */
  auto sqrt_kallen_1 = [](const double b, const double c) {
    return std::sqrt(std::pow(1.-b-c,2)-4.*b*c);
  };
  auto pseudoscalar_kappa = [sqrt_kallen_1](const int hel, const double U2, const double mP2, const double mN2, const double ml2) {
    const double xi_N = mN2/mP2;
    const double xi_l = ml2/mP2;
    const double sqkl = (mN2 < ml2) ? sqrt_kallen_1(xi_N,xi_l) : sqrt_kallen_1(xi_l,xi_N);
    return U2 * sqkl * (xi_l + xi_N + hel * (xi_N - xi_l) * sqkl) / (2. * xi_l * std::pow(1. - xi_l,2));
  };
  auto pdk = [](double a, double b, double c) {
    return 0.5*std::sqrt((a+b+c)*(a-b+c)*(a+b-c)*(a-b-c))/a;
  };
  auto make_uniform_2body_decayer = [pdk](int pdg1, int pdg2, double mP, double m1, double m2){
    return [mP,m1,m2,pdk](rng& rand, const TLorentzVector& parent){
      const double p_dec = pdk(mP,m1,m2);
      const double costh = CLHEP::RandFlat::shoot(&rand, -1, 1);
      const double sinth = std::sqrt(1.-costh*costh);
      const double phi = CLHEP::RandFlat::shoot(&rand, 2*M_PI);
      TLorentzVector p1{p_dec*sinth*std::cos(phi), p_dec*sinth*std::sin(phi), p_dec*costh, std::sqrt(p_dec*p_dec+m1*m1)};
      TLorentzVector p2{-p_dec*sinth*std::cos(phi), -p_dec*sinth*std::sin(phi), -p_dec*costh, std::sqrt(p_dec*p_dec+m2*m2)};
      p1.Boost(parent.BoostVector());
      p2.Boost(parent.BoostVector());
      return std::vector<TLorentzVector>{p1,p2};
    };
  };
  auto pseudoscalar_decays = [&consts=consts,&params=params,&decay_branches=decay_branches,
       pseudoscalar_kappa,make_uniform_2body_decayer]
    (const int pdg, const double lifetime, const double mass) {
    const double parent_gamma = consts.hbar() / lifetime;
    const double mP = mass;
    const double mP2 = mP*mP;
    const double mN = params.hnl_mass;
    const double mN2 = mN * mN;
    const double me = consts.mass_elec();
    const double me2 = me * me;
    const double mm = consts.mass_muon();
    const double mm2 = mm * mm;

    double tot_br = 0.;
    std::vector<branch_info_t> info_pos, info_neg;

    if(mP > me + mN && params.simulate("K->e") && std::abs(params.Ue4)>0.) {
      const double kappa_e_pos = pseudoscalar_kappa(+1, std::pow(params.Ue4,2),mP2,mN2,me2);
      const double gamma_e_pos = kappa_e_pos * parent_gamma;
      tot_br += gamma_e_pos / parent_gamma;
      std::cout << "Adding K->eN+ m= "<<params.hnl_mass<<" BR="<<gamma_e_pos / parent_gamma<<std::endl;

      branch_info_t info = std::make_tuple(decay_prod_vec_t{pdgcodes::k_HNL_poshel,-pdgcodes::k_elec},
          gamma_e_pos / parent_gamma,
          make_uniform_2body_decayer(pdgcodes::k_HNL_poshel,-pdgcodes::k_elec,mP,mN,me));
      info_pos.push_back(info);

      info = std::make_tuple(decay_prod_vec_t{(params.is_Dirac()?-1:1)*pdgcodes::k_HNL_neghel,pdgcodes::k_elec},
          gamma_e_pos / parent_gamma,
          make_uniform_2body_decayer((params.is_Dirac()?-1:1)*pdgcodes::k_HNL_neghel,pdgcodes::k_elec,mP,mN,me));
      info_neg.push_back(info);

      const double kappa_e_neg = pseudoscalar_kappa(-1, std::pow(params.Ue4,2),mP2,mN2,me2);
      const double gamma_e_neg = kappa_e_neg * parent_gamma;
      tot_br += gamma_e_neg / parent_gamma;
      
      std::cout << "Adding K->eN- m= "<<params.hnl_mass<<" BR="<<gamma_e_neg / parent_gamma<<std::endl;

      info = std::make_tuple(decay_prod_vec_t{pdgcodes::k_HNL_neghel,-pdgcodes::k_elec},
        gamma_e_neg / parent_gamma,
          make_uniform_2body_decayer(pdgcodes::k_HNL_neghel,-pdgcodes::k_elec,mP,mN,me));
      info_pos.push_back(info);

      info = std::make_tuple(decay_prod_vec_t{(params.is_Dirac()?-1:1)*pdgcodes::k_HNL_poshel,pdgcodes::k_elec},
        gamma_e_neg / parent_gamma,
          make_uniform_2body_decayer((params.is_Dirac()?-1:1)*pdgcodes::k_HNL_poshel,pdgcodes::k_elec,mP,mN,me));
      info_neg.push_back(info);
    }

    if(mP > mm + mN && params.simulate("K->mu") && std::abs(params.Umu4)>0.) {
      const double kappa_mu_pos = pseudoscalar_kappa(+1, std::pow(params.Umu4,2),mP2,mN2,mm2);
      const double gamma_mu_pos = kappa_mu_pos * parent_gamma;
      tot_br += gamma_mu_pos / parent_gamma;
      std::cout << "Adding K->muN+ m= "<<params.hnl_mass<<" BR="<<gamma_mu_pos / parent_gamma<<std::endl;

      branch_info_t info = std::make_tuple(decay_prod_vec_t{pdgcodes::k_HNL_poshel,-pdgcodes::k_muon},
          gamma_mu_pos / parent_gamma,
          make_uniform_2body_decayer(pdgcodes::k_HNL_poshel,-pdgcodes::k_muon,mP,mN,mm));
      info_pos.push_back(info);

      info = std::make_tuple(decay_prod_vec_t{(params.is_Dirac()?-1:1)*pdgcodes::k_HNL_neghel,pdgcodes::k_muon},
          gamma_mu_pos / parent_gamma,
          make_uniform_2body_decayer((params.is_Dirac()?-1:1)*pdgcodes::k_HNL_neghel,pdgcodes::k_muon,mP,mN,mm));
      info_neg.push_back(info);

      const double kappa_mu_neg = pseudoscalar_kappa(-1, std::pow(params.Umu4,2),mP2,mN2,mm2);
      const double gamma_mu_neg = kappa_mu_neg * parent_gamma;
      tot_br += gamma_mu_neg / parent_gamma;
      std::cout << "Adding K->muN- m= "<<params.hnl_mass<<" BR="<<gamma_mu_pos / parent_gamma<<std::endl;

      info = std::make_tuple(decay_prod_vec_t{pdgcodes::k_HNL_neghel,-pdgcodes::k_muon},
          gamma_mu_neg / parent_gamma,
          make_uniform_2body_decayer(pdgcodes::k_HNL_neghel,-pdgcodes::k_muon,mP,mN,mm));
      info_pos.push_back(info);

      info = std::make_tuple(decay_prod_vec_t{(params.is_Dirac()?-1:1)*pdgcodes::k_HNL_poshel,pdgcodes::k_muon},
          gamma_mu_neg / parent_gamma,
          make_uniform_2body_decayer((params.is_Dirac()?-1:1)*pdgcodes::k_HNL_poshel,pdgcodes::k_muon,mP,mN,mm));
      info_neg.push_back(info);
    }

    decay_branches[pdg]  = std::make_pair(tot_br,info_pos);
    decay_branches[-pdg] = std::make_pair(tot_br,info_neg);    
  };
  // kaons
  pseudoscalar_decays(pdgcodes::k_kaon_pm, consts.lifetime_kaon_pm(), consts.mass_kaon_pm());
  // pions
  // k0Ls
  // muons

  {
    auto k = decay_branches.find(pdgcodes::k_kaon_pm);
    if(k == decay_branches.end()
        || k->second.second.empty()
        || std::accumulate(k->second.second.begin(), k->second.second.end(), 0.,
          [](double a, auto const& b) { return a+std::get<1>(b);} ) <= 0.) {
      throw cet::exception("Configuration") << " no kaon decays are possible.";
  }
  }
  // HNLs
  {
    auto I1 = [sqrt_kallen_1](double x, double y) {
      return sqrt_kallen_1(x,y)*(std::pow(1.-x,2) - y*(1.+x));
    };
    auto I1_x0 = [](/*double x = 0.,*/ double y) {
      return std::pow(1.-y,2);
    };
    const double mpi0 = consts.mass_pion_0();
    const double mpi02 = mpi0*mpi0;
    const double mN = params.hnl_mass;
    const double mN2 = mN*mN;
    const double me = consts.mass_elec();
    const double me2 = me*me;
    const double mm = consts.mass_muon();
    const double mm2 = mm*mm;
    const double mpi = consts.mass_pion_pm();
    const double mpi2 = mpi*mpi;
    const double gF2 = consts.gFermi2();

    const double gamma_3nu = (params.is_Dirac() ? 0.5 : 1.) * params.sum_U2() * gF2 * std::pow(mN,5) / 96. / std::pow(M_PI,3);
    
    const double gamma_pi0_nu = (mN < mpi0) ? 0.
      : (params.is_Dirac() ? 0.5 : 1.) * params.sum_U2() * gF2
      * std::pow(consts.pion_decay_constant(),2) * std::pow(mN,3) / 16. / M_PI * I1_x0(mpi02/mN2);
    
    auto gamma_mes_lep = [Vud2=std::pow(consts.ckm_ud(),2),mN,mN2,mN3=mN*mN2,mP2=mpi2,
         fP2=std::pow(consts.pion_decay_constant(),2),I1,gF2]
      (const double U2, const double ml2) {
       return Vud2 * U2 * gF2 * fP2 * mN3 / 16. / M_PI * I1(ml2/mN2,mP2/mN2);
    };
    
    const double gamma_pi_e = (mN < me + mpi) ? 0.
      : gamma_mes_lep(std::pow(params.Ue4,2), me2);
    
    const double gamma_pi_mu = (mN < mm + mpi) ? 0.
      : gamma_mes_lep(std::pow(params.Umu4,2), mm2);

    double tot_gamma_pos_nu = gamma_3nu + gamma_pi0_nu + gamma_pi_e + gamma_pi_mu;
          std::cout << Form("Adding (any)N -> nu nu nu ; m= %g Gamma=%g \n",
               mN,gamma_3nu);
          std::cout << Form("Adding (any)N -> pi e ; m= %g Gamma=%g \n",
               mN,gamma_pi_e);
          std::cout << Form("Adding (any)N -> pi mu ; m= %g Gamma=%g \n",
               mN,gamma_pi_mu);
    double tot_gamma_neg_nu = tot_gamma_pos_nu, tot_gamma_pos_anu = tot_gamma_pos_nu, tot_gamma_neg_anu = tot_gamma_pos_nu;

    const double sin2_thW = std::pow(consts.sin_thW(),2);
    const double gL = -0.5+sin2_thW;
    const double gR = sin2_thW;
    const double Ue2 = std::pow(params.Ue4,2);
    const double Um2 = std::pow(params.Umu4,2);
    const double Ut2 = std::pow(params.Utau4,2);

    auto C1_nu = [Ue2,Um2,Ut2,gL,gL2=gL*gL](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return Ue2*((neg_type==pos_type?gL2:0.) + (neg_type==1?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Um2*((neg_type==pos_type?gL2:0.) + (neg_type==2?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Ut2*((neg_type==pos_type?gL2:0.) + (neg_type==3?1.:0.)*(1.+(neg_type==pos_type?gL:0.))); 
    };
    auto C2_nu = [Ue2,Um2,Ut2,gR2=gR*gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return (neg_type != pos_type) ? 0. : gR2 * (Ue2+Um2+Ut2);
    };
    auto C3_nu = [Ue2,Um2,Ut2,gL,gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return (neg_type != pos_type) ? 0. : gR
        * (Ue2*((pos_type==1?1.:0.)+gL)+Um2*((pos_type==2?1.:0.)+gL)+Ut2*((pos_type==3?1.:0.)+gL));
    };
    auto C4_nu = C1_nu;
    auto C5_nu = C2_nu;
    auto C6_nu = C3_nu;

    auto C1_anu = [Ue2,Um2,Ut2,gR2=gR*gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return (neg_type != pos_type) ? 0. : gR2 * (Ue2+Um2+Ut2);
    };
    auto C2_anu = [Ue2,Um2,Ut2,gL,gL2=gL*gL](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return Ue2*((neg_type==pos_type?gL2:0.) + (pos_type==1?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Um2*((neg_type==pos_type?gL2:0.) + (pos_type==2?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Ut2*((neg_type==pos_type?gL2:0.) + (pos_type==3?1.:0.)*(1.+(neg_type==pos_type?gL:0.))); 
    };
    auto C3_anu = [Ue2,Um2,Ut2,gL,gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return (neg_type != pos_type) ? 0. : gR
        * (Ue2*((neg_type==1?1.:0.)+gL)+Um2*((neg_type==2?1.:0.)+gL)+Ut2*((neg_type==3?1.:0.)+gL));
    };
    auto C4_anu = negate_func(C1_anu);
    auto C5_anu = negate_func(C2_anu);
    auto C6_anu = negate_func(C3_anu);

    auto C1_maj = [Ue2,Um2,Ut2,gL,gR,gsum=gL*gL+gR*gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return Ue2*((neg_type==pos_type?gsum:0.) + (neg_type==1?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Um2*((neg_type==pos_type?gsum:0.) + (neg_type==2?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Ut2*((neg_type==pos_type?gsum:0.) + (neg_type==3?1.:0.)*(1.+(neg_type==pos_type?gL:0.))); 
    };
    auto C2_maj = [Ue2,Um2,Ut2,gL,gsum=gL*gL+gR*gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return Ue2*((neg_type==pos_type?gsum:0.) + (pos_type==1?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Um2*((neg_type==pos_type?gsum:0.) + (pos_type==2?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Ut2*((neg_type==pos_type?gsum:0.) + (pos_type==3?1.:0.)*(1.+(neg_type==pos_type?gL:0.))); 
    };
    auto C3_maj = [Ue2,Um2,Ut2,gL,gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return (neg_type != pos_type) ? 0. : 2. * gR
        * (Ue2*((neg_type==1?1.:0.)+gL)+Um2*((neg_type==2?1.:0.)+gL)+Ut2*((neg_type==3?1.:0.)+gL));
    };
    auto C4_maj = [Ue2,Um2,Ut2,gL,gdiff=gL*gL-gR*gR](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return Ue2*((neg_type==pos_type?gdiff:0.) + (neg_type==1?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Um2*((neg_type==pos_type?gdiff:0.) + (neg_type==2?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Ut2*((neg_type==pos_type?gdiff:0.) + (neg_type==3?1.:0.)*(1.+(neg_type==pos_type?gL:0.))); 
    };
    auto C5_maj = [Ue2,Um2,Ut2,gL,gdiff=gR*gR-gL*gL](int neg_type /* 1=elec, 2=mu, 3=tau */, int pos_type) {
      return Ue2*((neg_type==pos_type?gdiff:0.) + (pos_type==1?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Um2*((neg_type==pos_type?gdiff:0.) + (pos_type==2?1.:0.)*(1.+(neg_type==pos_type?gL:0.)))
        + Ut2*((neg_type==pos_type?gdiff:0.) + (pos_type==3?1.:0.)*(1.+(neg_type==pos_type?gL:0.))); 
    };
    auto C6_maj = [](int , int){ return 0.; };

    auto make_partial_gamma_nu_lepm_lepp = [gF2,mN,mN2,sqkl1=sqrt_kallen_1](int hel, const double mlm2, const double mlp2, const double C1, const double C2, const double C3, const double C4, const double C5, const double C6) {
      const double xi_p = mlp2/mN2;
      const double xi_m = mlm2/mN2;
      return [pfac=gF2*std::pow(mN,5)/128./std::pow(M_PI,5),sqkl1,hel,C1,C2,C3,C4,C5,C6,xi_m,xi_p,mN]
        (const double s_nu_m, const double s_nu_p, const double cos_theta_m, const double cos_theta_p) {
          /* These equations changed based on arXiv:2104.05719 */
        const double A0_2 = C1*(s_nu_m)*(1.+xi_p-s_nu_m) + C2*(s_nu_p)*(1+xi_m-s_nu_p) + 2.*C3*std::sqrt(xi_m*xi_p)*(s_nu_m+s_nu_p-xi_m-xi_p);
        const double pfac_p = (C4*(s_nu_m)-2.*C6*std::sqrt(xi_m*xi_p));
        const double pfac_m = (C5*(s_nu_p)-2.*C6*std::sqrt(xi_m*xi_p));
        
        
        const double A1_2 =
          (std::abs(pfac_p) > 0. ?   pfac_p*sqkl1(s_nu_m,xi_p)*cos_theta_p : 0.)
          + (std::abs(pfac_m) > 0. ? pfac_m*sqkl1(s_nu_p,xi_m)*cos_theta_m : 0.);

        const double integrand = pfac * (A0_2+hel*A1_2);

        if(!(integrand > 0)) {
          // take care of floating point rounding, integrand could be small negative when it should be 0
          const double diff = std::abs(std::abs(A0_2)-std::abs(A1_2));
          const double avg = std::abs(A0_2)>0&&std::abs(A1_2)>0 ? 0.5 * (std::abs(A0_2) + std::abs(A1_2)) : 1.;
          if(!(diff/avg < 1e-15)) {
            /*
               const double Enu_star = (s_nu_m-xi_m)/(2*std::sqrt(s_nu_m)*mN);
               const double Epos_star = (1.-s_nu_m-xi_p)/(2.*std::sqrt(s_nu_m)*mN);
               const double pnu_star = Enu_star;
               const double ppos_star = std::sqrt(Epos_star*Epos_star-xi_p*mN*mN);
               const double s_min = std::pow(Enu_star + Epos_star,2)-std::pow(pnu_star+ppos_star,2);
               const double s_max = std::pow(Enu_star + Epos_star,2)-std::pow(pnu_star-ppos_star,2);
               if(s_nu_p < s_min || s_nu_p > s_max) {
               std::cerr << Form("s < s_min or s > s_max : (%g < %g < %g)\n",s_min, s_nu_p, s_max);
               }
               */
            auto kallen_1 = [](const double b, const double c) {
              return std::pow(1.-b-c,2)-4.*b*c;
            };
            throw art::Exception(art::errors::LogicError) << "Bad integrand_fn "<<Form("(%g, %g, %g, %g, %g, %g, %g)",pfac,A0_2,A1_2,(pfac * (A0_2+hel*A1_2)),diff,avg,diff/avg)
              <<" at points (s_nu_m,s_nu_p,cos_theta_m,cos_theta_p) = "
              <<Form("(%g, %g, %g, %g))",s_nu_m,s_nu_p,cos_theta_m,cos_theta_p)
              <<" -- meta = (hel,C1,C2,C3,C4,C5,C6,xi_m,xi_p,    sqkl1(s_nu_m,xi_p),sqkl1(s_nu_p,xi_m),    log(pfac_p),log(pfac_m)) = "
              <<Form("(%d, %g, %g, %g, %g, %g, %g, %g, %g, \n %g, %g, \n %g,%g)",
                  hel,C1,C2,C3,C4,C5,C6,xi_m,xi_p,
                  kallen_1(s_nu_m,xi_p),kallen_1(s_nu_p,xi_m),
                  std::log10(std::abs(pfac_p)),std::log10(std::abs(pfac_m)) )
              <<std::endl; 
          }
          std::cerr <<Form("Negative integrand! A0: %g; A1: %g; hel: %d; diff/avg: %g \n returning 0",A0_2,A1_2,hel,diff/avg);
          return 0.;
        }

#ifdef save_neg_lams
        auto kallen_1 = [](const double b, const double c) {
          return std::pow(1.-b-c,2)-4.*b*c;
        };
        return std::make_tuple(pfac * (A0_2+hel*A1_2),
            pfac_m,pfac_p,
            s_nu_m,s_nu_p,
            cos_theta_m,cos_theta_p,
            kallen_1(s_nu_p,xi_m),kallen_1(s_nu_m,xi_p));
#else
        return integrand;
#endif
      };
    };

    auto make_partial_gamma_nu_pi0 = [xi_p=mpi02/mN2](const int hel){
      auto diff_decay_rate = [xi_p,hel](double cos_theta_l) {
        return std::pow(xi_p-1,2)*(1.-hel*cos_theta_l);
      };
      return diff_decay_rate;
    };


#ifdef save_neg_lams
    using integrand_output_gamma_t = std::tuple<double,double,double,double,double,double,double,double,double>;
#else
    using integrand_output_gamma_t = double;
#endif
    using integrand_output_t = std::pair<integrand_output_gamma_t,std::vector<TLorentzVector>>;
    using integrand_func_t = std::function<integrand_output_t(double,double,double,double,double,bool)>;
    using integrand_func2D_t = std::function<integrand_output_t(double,double,bool)>;
    using minmax_s_t = std::pair<double,double>;
    using tot_gamma_t = double;
    using max_integrand_t = double;
#ifdef save_neg_lams
    //using tot_gamma_t = std::tuple<double,double,double>;
    using integrand_t = std::tuple<tot_gamma_t,max_integrand_t,integrand_func_t,minmax_s_t,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>;
#else
    using integrand_t = std::tuple<tot_gamma_t,max_integrand_t,integrand_func_t,minmax_s_t>;
    using integrand2D_t = std::tuple<tot_gamma_t,max_integrand_t,integrand_func2D_t>;
#endif

    auto make_2body_decay_integrand = [mN,max_evals=num_evals_numerical_integration,pdk]
      (rng& rand, auto gamma_fn, double mlep, double mpscal) -> integrand2D_t {

      double sum_integrand = 0.;
      double max_integrand = 0.;

      const double pdecay = pdk(mN,mlep,mpscal);
      const double eneL = std::sqrt(pdecay*pdecay + mlep*mlep);
      const double eneP = std::sqrt(pdecay*pdecay + mpscal*mpscal);
    

      auto integrand_fn = [pdecay,eneL,eneP,gamma_fn]
        (double costh_p, double phi_p, bool save_output = true) {
        
        std::vector<TLorentzVector> outputs;
        if(save_output) {
          const double sinth_p = std::sqrt(1.-costh_p*costh_p);
          outputs.push_back(TLorentzVector{pdecay*sinth_p*std::cos(phi_p),
          pdecay*sinth_p*std::sin(phi_p), pdecay*costh_p, eneL});
          outputs.push_back(TLorentzVector{-pdecay*sinth_p*std::cos(phi_p),
          -pdecay*sinth_p*std::sin(phi_p), -pdecay*costh_p, eneP});
        }
        
        return std::pair(gamma_fn(costh_p),outputs);
      };

      for(size_t i = 0; i < max_evals; ++i) {
        const double costh_p = CLHEP::RandFlat::shoot(&rand, -1., 1.);
        const double phi_p = CLHEP::RandFlat::shoot(&rand, 2*M_PI);
        try {
          const double integrand = integrand_fn(costh_p, phi_p, false).first;
          sum_integrand += integrand;
          if(integrand > max_integrand) max_integrand = 1.1*integrand;
        } catch(art::Exception& e) {
          //std::cerr<< "found error "<<e.explain_self()<<std::endl;
          if(e.explain_self().find("Bad integrand_fn") != std::string::npos) {
            //std::cerr<< "modifying error "<<std::endl;
            std::ostringstream ostr;
            ostr << "Bad integrand "
              << (integrand_fn(costh_p, phi_p, false).first)
              <<" at points (costh_p,phi_p) = "
              <<Form("(%g, %g))",costh_p,phi_p)
              <<" -- meta = (mN,mlep,mpscal) = "
              <<Form("(%g, %g, %g)",mN,mlep,mpscal)
              <<std::endl;
            throw art::Exception(art::errors::LogicError,ostr.str(), e);
          }
          else throw e;
        }
      }
      const double tot_gamma = sum_integrand / max_evals *  4. * M_PI;
      return std::make_tuple(tot_gamma, max_integrand, integrand_fn);
    };
    
    auto make_3body_decay_integrand = [mN,mN2,max_evals=num_evals_numerical_integration,s23weighting=use_s23_weighting]
      (rng& rand, auto gamma_fn, double mlm2, double mlp2) -> integrand_t {
      const double xi_p = mlp2/mN2;
      const double xi_m = mlm2/mN2;
      const double min_s = xi_m;
      const double max_s = std::pow(1. - std::sqrt(xi_p),2);

      double sum_integrand = 0.;
      //double sum_uniform_integrand = 0.;
      double max_integrand = 0.;
    
      // xi1, xi2 are squared reduced masses
      // pdk1 for pdk(1, xi1=xip, xi2)
      auto pdk1 = [xi_p,xi_p2=xi_p*xi_p](const double xi2) {
        return 0.5*std::sqrt(xi_p2 + std::pow(xi2-1.,2) - 2.*xi_p*(1.+xi2));
      };
      // pdk0 for pdk(xi0, xi2=xi_m, 0)
      auto pdk0 = [xi_m](const double xi0) {
        return 0.5*(xi0-xi_m)/std::sqrt(xi0);
      };

      auto integrand_fn = [mN,mN2,pdk1,pdk0,xi_p,xi_m,gamma_fn,max_evals,s23weighting]
        (double s_nu_lep_neg, double costh_p, double phi_p, double costh_n_rest, double phi_n_rest, bool save_output = true) {
        const double p_pos = mN * pdk1(s_nu_lep_neg);
        const double sinth_p = std::sqrt(1.-costh_p*costh_p);
        const double p_neg_nu_rest = mN * pdk0(s_nu_lep_neg);
        const double sinth_n_rest = std::sqrt(1.-costh_n_rest*costh_n_rest);
        
        std::vector<TLorentzVector> outputs;
        const TLorentzVector mom_nu_lep_neg{-p_pos*sinth_p*std::cos(phi_p),
          -p_pos*sinth_p*std::sin(phi_p), -p_pos*costh_p, std::sqrt(s_nu_lep_neg*mN2 + p_pos*p_pos)};

        const double Eneg_rest = std::sqrt(std::pow(p_neg_nu_rest,2) + xi_m*mN2);
        TLorentzVector mom_neg{p_neg_nu_rest * sinth_n_rest * std::cos(phi_n_rest),
          p_neg_nu_rest * sinth_n_rest * std::sin(phi_n_rest), p_neg_nu_rest * costh_n_rest, Eneg_rest};
        mom_neg.Boost(mom_nu_lep_neg.BoostVector());
        const double costh_m = mom_neg.Vect().Unit().CosTheta();

        if(save_output) {
          TLorentzVector mom_pos{p_pos*sinth_p*std::cos(phi_p),
          p_pos*sinth_p*std::sin(phi_p), p_pos*costh_p, std::sqrt(xi_p*mN2 + p_pos*p_pos)};
          TLorentzVector mom_nu{-p_neg_nu_rest * sinth_n_rest * std::cos(phi_n_rest),
          -p_neg_nu_rest * sinth_n_rest * std::sin(phi_n_rest), -p_neg_nu_rest * costh_n_rest, p_neg_nu_rest}; 
          mom_nu.Boost(mom_nu_lep_neg.BoostVector());
          outputs.push_back(mom_nu);
          outputs.push_back(mom_neg);
          outputs.push_back(mom_pos);
        }
        
        const double Eneg = mom_neg.E();
        const double s_nu_lep_pos = 1 + xi_m - 2*Eneg/mN; // s == (pnu + plp)^2 == (pN - plm)^2

        const double s23weight = !s23weighting ? 1.
          : [=](){
          const double E2star = mN*(s_nu_lep_neg + xi_m)/(2.*std::sqrt(s_nu_lep_neg));
          const double E3star = mN*(1. - s_nu_lep_neg - xi_p)/(2.*std::sqrt(s_nu_lep_neg));
          const double p2star = std::sqrt(E2star*E2star - xi_m*xi_m*mN2);
          const double p3star = std::sqrt(E3star*E3star - xi_p*xi_p*mN2);
          const double low_lim_s23 = std::pow(E2star + E3star,2)-std::pow(p2star+p3star,2);
          const double up_lim_s23 = std::pow(E2star + E3star,2)-std::pow(p2star-p3star,2);
          return up_lim_s23 - low_lim_s23; }();

        return std::pair(gamma_fn(s_nu_lep_neg, s_nu_lep_pos, costh_m, costh_p)*s23weight,outputs);
      };
      /*
      auto uniform_integrand_fn = [mN,mN2,xi_p,xi_m,s23weighting]
        (double s_nu_lep_neg) {

        const double s23weight = !s23weighting ? 1.
          : [=](){
          const double E2star = mN*(s_nu_lep_neg + xi_m)/(2.*std::sqrt(s_nu_lep_neg));
          const double E3star = mN*(1. - s_nu_lep_neg - xi_p)/(2.*std::sqrt(s_nu_lep_neg));
          const double p2star = std::sqrt(E2star*E2star - xi_m*xi_m*mN2);
          const double p3star = std::sqrt(E3star*E3star - xi_p*xi_p*mN2);
          const double low_lim_s23 = std::pow(E2star + E3star,2)-std::pow(p2star+p3star,2);
          const double up_lim_s23 = std::pow(E2star + E3star,2)-std::pow(p2star-p3star,2);
          return up_lim_s23 - low_lim_s23; }();

        return s23weight;
      };
      */

#ifdef save_neg_lams
      std::vector<double> pfps, pfns, snums,snups,thms,thps,klms,klps,kint;
      pfps.reserve(max_evals);
      pfns.reserve(max_evals);
      snums.reserve(max_evals);
      snups.reserve(max_evals);
      thms.reserve(max_evals);
      thps.reserve(max_evals);
      klms.reserve(max_evals);
      klps.reserve(max_evals);
      kint.reserve(max_evals);
#endif
      for(size_t i = 0; i < max_evals; ++i) {
        const double s_nu_lep_neg = CLHEP::RandFlat::shoot(&rand, min_s, max_s);
        const double costh_p = CLHEP::RandFlat::shoot(&rand, -1., 1.);
        const double phi_p = CLHEP::RandFlat::shoot(&rand, 2*M_PI);
        const double costh_n_rest = CLHEP::RandFlat::shoot(&rand, -1., 1.);
        const double phi_n_rest = CLHEP::RandFlat::shoot(&rand, 2*M_PI);
        try {
#ifdef save_neg_lams
          auto const& iii = integrand_fn(s_nu_lep_neg, costh_p, phi_p, costh_n_rest, phi_n_rest, false).first;
          const double integrand = std::get<0>(iii);
          pfns.push_back(std::get<1>(iii));
          pfps.push_back(std::get<2>(iii));
          snums.push_back(std::get<3>(iii));
          snups.push_back(std::get<4>(iii));
          thms.push_back(std::get<5>(iii));
          thps.push_back(std::get<6>(iii));
          klms.push_back(std::get<7>(iii));
          klps.push_back(std::get<8>(iii));
          kint.push_back(integrand);
#else
          const double integrand = integrand_fn(s_nu_lep_neg, costh_p, phi_p, costh_n_rest, phi_n_rest, false).first;
#endif
          sum_integrand += integrand;
          //const double uniform_integrand = uniform_integrand_fn(s_nu_lep_neg);
          //sum_uniform_integrand += uniform_integrand;
          if(integrand > max_integrand) max_integrand = 1.1*integrand;
        } catch(art::Exception& e) {
          //std::cerr<< "found error "<<e.explain_self()<<std::endl;
          if(e.explain_self().find("Bad integrand_fn") != std::string::npos) {
            //std::cerr<< "modifying error "<<std::endl;
            std::ostringstream ostr;
            ostr << "Bad integrand "
#ifdef save_neg_lams
              << std::get<0>(integrand_fn(s_nu_lep_neg, costh_p, phi_p, costh_n_rest, phi_n_rest, false).first)
#else
              << (integrand_fn(s_nu_lep_neg, costh_p, phi_p, costh_n_rest, phi_n_rest, false).first)
#endif
              <<" at points (s_nu_lep_neg,costh_p,phi_p,costh_n_rest,phi_n_rest) = "
              <<Form("(%g, %g, %g, %g, %g))",s_nu_lep_neg,costh_p,phi_p,costh_n_rest,phi_n_rest)
              <<" -- meta = (mN,pdk1,pdk0,xi_p,xi_m) = "
              <<Form("(%g, %g, %g, %g, %g)",mN,pdk1(s_nu_lep_neg),pdk0(s_nu_lep_neg),xi_p,xi_m)
              <<std::endl;
            throw art::Exception(art::errors::LogicError,ostr.str(), e);
          }
          else throw e;
        }
      }
      const double tot_gamma = sum_integrand / max_evals * (max_s - min_s) * std::pow(4 * M_PI, 2);
      //const double tot_gamma_uniform = sum_uniform_integrand / max_evals * (max_s - min_s);
#ifdef save_neg_lams
      return std::make_tuple(tot_gamma, max_integrand, integrand_fn, std::make_pair(min_s, max_s),
          pfns,pfps, snums,snups,thms,thps,klms,klps,kint);
#else
      return std::make_tuple(tot_gamma /* / tot_gamma_uniform*/, max_integrand, integrand_fn, std::make_pair(min_s, max_s));
#endif
    };
    

    auto make_twobody_final_state_function = [](const integrand2D_t& integrand_info) {
      auto integrand_fn = std::get<2>(integrand_info);
      const double max_integrand = std::get<1>(integrand_info);
      return [integrand_fn,max_integrand](rng& rand, const TLorentzVector& parmom) {
        while(true) {
          const double costh_p = CLHEP::RandFlat::shoot(&rand, -1., 1.);
          const double phi_p = CLHEP::RandFlat::shoot(&rand, 2*M_PI);
          auto const& integrand_out = integrand_fn(costh_p, phi_p, true);
          const double u = CLHEP::RandFlat::shoot(&rand, max_integrand);
          if(u  < integrand_out.first) {
            std::vector<TLorentzVector> outputs;
            std::transform(integrand_out.second.begin(), integrand_out.second.end(), std::back_inserter(outputs),
                [&parmom](const TLorentzVector& v) {
                  if(!(parmom.Vect().Mag() > 0.)) {
                    return v;
                  }
                  TVector3 new3v = v.Vect();
                  new3v.RotateUz(parmom.Vect().Unit());
                  TLorentzVector newv{new3v,v.E()};
                  newv.Boost(parmom.BoostVector());
                  return newv;
                });
            return outputs;
          }
        }
      };
    };

    auto make_final_state_function = [](const integrand_t& integrand_info) {
      const double min_s = std::get<3>(integrand_info).first;
      const double max_s = std::get<3>(integrand_info).second;
      auto integrand_fn = std::get<2>(integrand_info);
      const double max_integrand = std::get<1>(integrand_info);
      return [min_s,max_s,integrand_fn,max_integrand](rng& rand, const TLorentzVector& parmom) {
        while(true) {
          const double s_nu_lep_neg = CLHEP::RandFlat::shoot(&rand, min_s, max_s);
          const double costh_p = CLHEP::RandFlat::shoot(&rand, -1., 1.);
          const double phi_p = CLHEP::RandFlat::shoot(&rand, 2*M_PI);
          const double costh_n_rest = CLHEP::RandFlat::shoot(&rand, -1., 1.);
          const double phi_n_rest = CLHEP::RandFlat::shoot(&rand, 2*M_PI);
          auto const& integrand_out = integrand_fn(s_nu_lep_neg, costh_p, phi_p, costh_n_rest, phi_n_rest, true);
          const double u = CLHEP::RandFlat::shoot(&rand, max_integrand);
#ifdef save_neg_lams
          if(u  < std::get<0>(integrand_out.first)) {
#else
          if(u  < integrand_out.first) {
#endif
            std::vector<TLorentzVector> outputs;
            std::transform(integrand_out.second.begin(), integrand_out.second.end(), std::back_inserter(outputs),
                [&parmom](const TLorentzVector& v) {
                  if(!(parmom.Vect().Mag() > 0.)) {
                    return v;
                  }
                  TVector3 new3v = v.Vect();
                  new3v.RotateUz(parmom.Vect().Unit());
                  TLorentzVector newv{new3v,v.E()};
                  newv.Boost(parmom.BoostVector());
                  return newv;
                });
            return outputs;
          }
        }
      };
    };
      
    auto append_decays = [make_partial_gamma_nu_lepm_lepp,make_3body_decay_integrand,make_final_state_function,mN,gF2,print_ratio=params.print_ratio,sw2=sin2_thW,Ue2,Um2,Ut2,mj=params.is_Majorana()]
      (auto& decayinfo, int hel, int pdgnu, int pdgneg, int pdgpos,
       double mlm2, double mlp2, int flavm, int flavp,
       auto C1, auto C2, auto C3, auto C4, auto C5, auto C6,
       double& tot_gamma, rng& rand) {
        const double c1 = C1(flavm,flavp);
        const double c2 = C2(flavm,flavp);
        const double c3 = C3(flavm,flavp);
        const double c4 = C4(flavm,flavp);
        const double c5 = C5(flavm,flavp);
        const double c6 = C6(flavm,flavp);
        if(std::abs(c1)>0 || std::abs(c2) > 0 || std::abs(c3) > 0 || std::abs(c4) > 0 || std::abs(c5) >0 || std::abs(c6) > 0) {
          auto fn = make_partial_gamma_nu_lepm_lepp(hel,mlm2,mlp2,c1,c2,c3,c4,c5,c6);
          auto info = make_3body_decay_integrand(rand, fn, mlm2, mlp2);
          //const double gamma = std::get<0>(info);
          auto acoth = [](const double x) {
            return 0.5 * (std::log(1+1/x) - std::log(1-1./x));
          };
          const double xi1 = mlm2/mN/mN;
          const double xi2 = mlp2/mN/mN;
          const double i1_34_fullform =  flavm == flavp
            ? std::sqrt(1 - 4*xi1)*(1 - 2*xi1*(7 + xi1 + 6*std::pow(xi1,2))) + 12*std::pow(xi1,2)*(-1 + std::pow(xi1,2))*std::atanh((std::sqrt(1 - 4*xi1)*(-1 + 2*xi1))/(1 + 2*(-2 + xi1)*xi1))
            : std::sqrt(std::pow(-1 + xi1,2) - 2*(1 + xi1)*xi2 + std::pow(xi2,2))*(1 + xi1*(-7 + (-7 + xi1)*xi1) - 7*xi2 + (12 - 7*xi1)*xi1*xi2 - 7*(1 + xi1)*std::pow(xi2,2) + std::pow(xi2,3)) + 12*(std::pow(xi1,2) + (-2 + std::pow(xi1,2))*std::pow(xi2,2))*acoth((-1 + xi1 - xi2)/std::sqrt(std::pow(xi1,2) + std::pow(-1 + xi2,2) - 2*xi1*(1 + xi2))) + 12*std::pow(xi1,2)*(-1 + std::pow(xi2,2))*std::atanh(((-1 + xi2)*std::sqrt(std::pow(xi1,2) + std::pow(-1 + xi2,2) - 2*xi1*(1 + xi2)))/(std::pow(-1 + xi2,2) - xi1*(1 + xi2)))
            ;
          const double i1_34 = flavm == flavp 
            ? std::max(
              1 - 16*xi1 - 32*std::pow(xi1,3) + 32*std::pow(xi1,5) - 24*std::pow(xi1,2)*(-1 + std::log(xi1)) + 2*std::pow(xi1,4)*(1 + 12*std::log(xi1)),
              (32*std::pow(1 - 4*xi1,3.5))/35.
            )
            : std::max(xi1,xi2) < 0.96 ? 1 + 8*xi1*std::pow(-1 + xi2,3) - xi2*(8 + (-8 + xi2)*std::pow(xi2,2)) + 8*std::pow(xi1,3)*(1 + xi2)*(1 + 2*(std::pow(xi2,2) + std::pow(xi2,4))) + (16*std::pow(xi1,5)*std::pow(xi2,2)*(1 + 5*xi2*(1 + xi2*(3 + 7*xi2))))/5. + std::pow(xi1,4)*(-1 + 2*std::pow(xi2,2)*(3 + xi2*(8 + 3*xi2*(5 + 8*xi2)))) + (2*std::pow(xi1,2)*xi2*(-60 + xi2*(-30 + xi2*(40 + xi2*(15 + 8*xi2)))))/5. + 12*std::pow(xi1,2)*(-1 + std::pow(xi2,2))*std::log(xi1) + 12*(-1 + std::pow(xi1,2))*std::pow(xi2,2)*std::log(xi2) : i1_34_fullform
            ;
          /*
          const double i1_34_a = std::sqrt(1 - 4*xi1);
          const double i1_34_ba = (std::sqrt(1 - 4*xi1)*(-1 + 2*xi1));
          const double i1_34_bb = (1 + 2*(-2 + xi1)*xi1);
          const double i1_34_c = ((std::sqrt(1 - 4*xi1)*(-1 + 2*xi1))/(1 + 2*(-2 + xi1)*xi1));
          const double i1_34_d = std::atanh((std::sqrt(1 - 4*xi1)*(-1 + 2*xi1))/(1 + 2*(-2 + xi1)*xi1));
          */
          const double i1_43_fullform = flavm == flavp ? i1_34_fullform
            : std::sqrt(std::pow(-1 + xi2,2) - 2*(1 + xi2)*xi1 + std::pow(xi1,2))*(1 + xi2*(-7 + (-7 + xi2)*xi2) - 7*xi1 + (12 - 7*xi2)*xi2*xi1 - 7*(1 + xi2)*std::pow(xi1,2) + std::pow(xi1,3)) + 12*(std::pow(xi2,2) + (-2 + std::pow(xi2,2))*std::pow(xi1,2))*acoth((-1 + xi2 - xi1)/std::sqrt(std::pow(xi2,2) + std::pow(-1 + xi1,2) - 2*xi2*(1 + xi1))) + 12*std::pow(xi2,2)*(-1 + std::pow(xi1,2))*std::atanh(((-1 + xi1)*std::sqrt(std::pow(xi2,2) + std::pow(-1 + xi1,2) - 2*xi2*(1 + xi1)))/(std::pow(-1 + xi1,2) - xi2*(1 + xi1)))
            ;
          const double i1_43= flavm == flavp ? i1_34
            : std::max(xi1,xi2) < 0.96 ? 1 + 8*xi2*std::pow(-1 + xi1,3) - xi1*(8 + (-8 + xi1)*std::pow(xi1,2)) + 8*std::pow(xi2,3)*(1 + xi1)*(1 + 2*(std::pow(xi1,2) + std::pow(xi1,4))) + (16*std::pow(xi2,5)*std::pow(xi1,2)*(1 + 5*xi1*(1 + xi1*(3 + 7*xi1))))/5. + std::pow(xi2,4)*(-1 + 2*std::pow(xi1,2)*(3 + xi1*(8 + 3*xi1*(5 + 8*xi1)))) + (2*std::pow(xi2,2)*xi1*(-60 + xi1*(-30 + xi1*(40 + xi1*(15 + 8*xi1)))))/5. + 12*std::pow(xi2,2)*(-1 + std::pow(xi1,2))*std::log(xi2) + 12*(-1 + std::pow(xi2,2))*std::pow(xi1,2)*std::log(xi1) : i1_43_fullform
            ;
          const double i2_34_fullform = flavm == flavp
            ? 8*xi1*(std::sqrt(1 - 4*xi1)*(1 + (5 - 6*xi1)*xi1) + 6*xi1*(1 + 2*(-1 + xi1)*xi1)*(-2*std::log(1 + std::sqrt(1 - 4*xi1)) + std::log(4*xi1)))
            : 12*std::sqrt(xi1*xi2)*(-1./3.*((-2 + (-5 + xi1)*xi1 - 5*xi2 + 10*xi1*xi2 + std::pow(xi2,2))*std::sqrt(std::pow(-1 + xi1,2) - 2*(1 + xi1)*xi2 + std::pow(xi2,2))) + 4*xi1*xi2*(-1 + xi1 + xi2)*acoth((-1 + xi1 + xi2)/std::sqrt(std::pow(xi1,2) + std::pow(-1 + xi2,2) - 2*xi1*(1 + xi2))) + 2*(-xi1 + xi2)*std::atanh(((xi1 - xi2)*std::sqrt(std::pow(xi1,2) + std::pow(-1 + xi2,2) - 2*xi1*(1 + xi2)))/(std::pow(xi1,2) + (-1 + xi2)*xi2 - xi1*(1 + 2*xi2))) + (-xi1 - xi2 + 2*xi1*xi2)*std::log((-1 + xi1 + xi2 - std::sqrt(std::pow(-1 + xi1,2) - 2*(1 + xi1)*xi2 + std::pow(xi2,2)))/(-1 + xi1 + xi2 + std::sqrt(std::pow(-1 + xi1,2) - 2*(1 + xi1)*xi2 + std::pow(xi2,2)))));
          const double i2_34 = flavm == flavp 
            ? (xi1<0.15 ? 8*xi1*((-1 + xi1)*(-1 + 2*xi1*(-2 + xi1 + 5*std::pow(xi1,2))) + 6*xi1*(1 + 2*(-1 + xi1)*xi1)*std::log(xi1))
            : (32*std::pow(1 - 4*xi1,3.5)*(1 + 8*xi1))/105.)
            :4*std::sqrt(xi1*xi2)*(2 + xi1*(3 + (-6 + xi1)*xi1) + 3*xi2 + std::pow(xi1,2)*(-9 + xi1*(4 + xi1))*xi2 + (-6 + xi1*(-9 + 2*xi1*(6 + xi1*(3 + 2*xi1))))*std::pow(xi2,2) + (1 + 2*xi1*(2 + xi1*(3 + xi1*(4 + 5*xi1))))*std::pow(xi2,3) + xi1*(1 + 2*xi1*(2 + 5*xi1*(1 + 2*xi1)))*std::pow(xi2,4) + 6*xi1*(1 + xi2*(-2 + xi1 + xi2))*std::log(xi1) + 6*xi2*(1 + xi1*(-2 + xi1 + xi2))*std::log(xi2))
            ;
          const double calc_gamma_durham = gF2 * std::pow(mN,5) / (192. * M_PI * M_PI * M_PI) * (c1 * i1_34 + c2 * i1_43 + c3 * i2_34); /* FIXME!!!! gamma from equations, not MC integration */
          [&](){return i1_43_fullform*i2_34_fullform*calc_gamma_durham;};
          //const double calc_gamma_durham_fullform = gF2 * std::pow(mN,5) / (192. * M_PI * M_PI * M_PI) * (c1 * i1_34_fullform + c2 * i1_43_fullform + c3 * i2_34_fullform); /* FIXME!!!! gamma from equations, not MC integration */

          /*
          if(calc_gamma_durham < 0) {
            std::cerr << Form("ERROR! -ve calc_gamma_durham (c1,c2,c3) (%g,%g,%g) (i1_34,i1_43,i2_34) (%g,%g,%g) (full-) (%g,%g,%g) -- (%g,%g,%g,%g)\n",
              c1,c2,c3,
              i1_34,i1_43,i2_34,
              i1_34_fullform,i1_43_fullform,i2_34_fullform,
              (-4096*std::sqrt(1 - 4*xi1)*std::pow(-0.25 + xi1,3.5))/(35.*std::sqrt(-1 + 4*xi1)),
              std::sqrt(1 - 4*xi1),
              std::pow(-0.25 + xi1,3.5),
              std::sqrt(-1 + 4*xi1)
            );
            //throw 0;
          }
          */
          /*
          if(std::isnan(calc_gamma_durham_fullform) || std::isinf(calc_gamma_durham_fullform)) {
            std::cerr << Form("NAN (%g,%g,%g) (%g,%g,%g,%g,%g,%g,%g,%g) \n",i1_34_fullform,i1_43_fullform,i2_34_fullform,
              std::sqrt(1 - 4*xi1),
              (1 - 2*xi1*(7 + xi1 + 6*std::pow(xi1,2))),
               12*std::pow(xi1,2),
               (-1 + std::pow(xi1,2)),
               std::atanh((std::sqrt(1 - 4*xi1)*(-1 + 2*xi1))/(1 + 2*(-2 + xi1)*xi1)),
               std::sqrt(1 - 4*xi1)*(-1 + 2*xi1),
               std::sqrt(1 - 4*xi1),
               (1 + 2*(-2 + xi1)*xi1)
            );
            throw 0;
          }
          */

#ifdef USE_GENIE_EQUATIONS
          /* USING GENIE EQUATIONS (2007.03701), MARCH 2023 */
          //const double sw_bit = flavm == flavp ? sw2 : 0.;
          const double Lx = std::log((1.-3*xi1 - (1.-xi1)*std::sqrt(1-4*xi1))/(xi1*(1+std::sqrt(1-4*xi1))));
          const double c1_bit = .25*(1-4.*sw2+8.*sw2*sw2);
          const double c2_bit = .5*(2.*sw2*sw2-sw2);
          const double f1_bit = (1.-14.*xi1 - 2*xi1*xi1 -12.*xi1*xi1*xi1)*std::sqrt(1-4.*xi1)+12.*xi1*xi1*(xi1*xi1-1.)*Lx;
          const double f2_bit = 4.*(xi1*(2.+10.*xi1-12.*xi1*xi1)*std::sqrt(1-4*xi1)+6.*xi1*xi1*(1.-2*xi1+2*xi1*xi1)*Lx);
          const double xM = std::max(xi1,xi2);
          const double xM_bit = (1.-8*xM+8*xM*xM*xM - xM*xM*xM*xM - 12*xM*xM*std::log(xM));
          const double e_bit = Ue2 * ((flavm==flavp)?((c1_bit+(flavm==1?2:0)*sw2)*f1_bit + (c2_bit+(flavm==1?1:0)*sw2)*f2_bit):(flavm==1?xM_bit:0.));
          const double mu_bit = Um2 * ((flavm==flavp)?((c1_bit+(flavm==2?2:0)*sw2)*f1_bit + (c2_bit+(flavm==2?1:0)*sw2)*f2_bit):(flavm==2?xM_bit:0.));
          const double tau_bit = Ut2 * ((flavm==flavp)?((c1_bit+(flavm==3?2:0)*sw2)*f1_bit + (c2_bit+(flavm==3?1:0)*sw2)*f2_bit):(flavm==3?xM_bit:0.));
          const double calc_gamma =  gF2 * std::pow(mN,5) / (192. * M_PI * M_PI * M_PI) * (e_bit+mu_bit+tau_bit) * (mj&&(flavm==flavp)?2.:1.);
#else
          const double calc_gamma = calc_gamma_durham;
#endif

          decayinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{pdgnu,pdgneg,pdgpos},
                calc_gamma,
                make_final_state_function(info)
                ));
          tot_gamma += calc_gamma;
          
          const char* flavs[] = {"", "e","mu","tau"};
          //std::cout << Form("Adding %sN%s -> nu %s- %s+ ; Gamma=%g (%g,%g,%g) m= %g calc_gamma=%g\n",
          //    pdgnu<0?"~":"", hel>0?"+":"-", flavs[flavm], flavs[flavp], gamma,i1_34,i1_43,i2_34,  /*i1_34_a,i1_34_ba,i1_34_bb,i1_34_c,i1_34_d,*/  mN,calc_gamma);
          if(print_ratio) {
          std::cout << Form("Adding %sN%s -> nu %s- %s+ ; m= %g Gamma=%g OldGamma=%g ratio=%g \n",
              pdgnu<0?"~":"", hel>0?"+":"-", flavs[flavm], flavs[flavp], mN, calc_gamma, std::get<0>(info), calc_gamma/std::get<0>(info));
          }
          else {
          std::cout << Form("Adding %sN%s -> nu %s- %s+ ; m= %g Gamma=%g \n",
              pdgnu<0?"~":"", hel>0?"+":"-", flavs[flavm], flavs[flavp], mN, calc_gamma);
          }
#ifdef save_neg_lams
          return std::make_tuple(std::get<4>(info),std::get<5>(info),std::get<6>(info),std::get<7>(info),std::get<8>(info),std::get<9>(info),std::get<10>(info),std::get<11>(info),std::get<12>(info));
        }
        else {
          std::cerr << "All Cs are 0: flavm:"<<flavm<<" flavp:"<<flavp<<std::endl;
          std::vector<double> tmp;
          return std::make_tuple(tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp);
#endif
        }
      };

    std::map<int, decay_info_t> dummymap{};
    auto& dummybranch = dummymap[0];

    // hnl - poshel
    if(params.is_Dirac()) {
      const int hel = +1;
      // particle
      auto& posnuinfo = decay_branches[pdgcodes::k_HNL_poshel];

      // nu pi0
      if(mN > mpi0 && gen_nu_pi0) {
        posnuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(hel), 0., mpi0))
                ));
        const int pdgnu = pdgcodes::k_nu_e;
          std::cout << Form("Adding %sN%s -> nu pi0 ; m= %g Gamma=%g \n",
              pdgnu<0?"~":"", hel>0?"+":"-", mN,gamma_pi0_nu);
      }
      
      // e-e+
      if(mN > 2. * me) {
#ifdef save_neg_lams
        auto const& ret = 
#endif
        append_decays((gen_nu_lep_lep|gen_nu_e_e)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_e,pdgcodes::k_elec,-pdgcodes::k_elec,
            me2, me2, 1, 1, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_pos_nu, rand);
#ifdef save_neg_lams
        if(tt) {
        t_mn = me;
        t_mp = me;
        auto const& v1 = std::get<0>(ret);
        auto const& v2 = std::get<1>(ret);
        auto const& v3 = std::get<2>(ret);
        auto const& v4 = std::get<3>(ret);
        auto const& v5 = std::get<4>(ret);
        auto const& v6 = std::get<5>(ret);
        auto const& v7 = std::get<6>(ret);
        auto const& v8 = std::get<7>(ret);
        auto const& v9 = std::get<8>(ret);
        for(size_t i = 0; i < std::min(v1.size(),v2.size()); ++i) {
          t_pfn = v1[i];
          t_pfp = v2[i];
          t_snum = v3[i];
          t_snup = v4[i];
          t_cthm = v5[i];
          t_cthp = v6[i];
          t_klm = v7[i];
          t_klp = v8[i];
          t_integrand = v9[i];
          tt->Fill();
        }
        }
#endif
      }

      if(mN > me + mm) {
        // e- mu+
#ifdef save_neg_lams
        {
        auto const& ret = 
#endif
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_mu,pdgcodes::k_elec,-pdgcodes::k_muon,
            me2, mm2, 1, 2, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_pos_nu, rand);
#ifdef save_neg_lams
        if(tt) {
        t_mn = me;
        t_mp = mm;
        auto const& v1 = std::get<0>(ret);
        auto const& v2 = std::get<1>(ret);
        auto const& v3 = std::get<2>(ret);
        auto const& v4 = std::get<3>(ret);
        auto const& v5 = std::get<4>(ret);
        auto const& v6 = std::get<5>(ret);
        auto const& v7 = std::get<6>(ret);
        auto const& v8 = std::get<7>(ret);
        auto const& v9 = std::get<8>(ret);
        for(size_t i = 0; i < std::min(v1.size(),v2.size()); ++i) {
          t_pfn = v1[i];
          t_pfp = v2[i];
          t_snum = v3[i];
          t_snup = v4[i];
          t_cthm = v5[i];
          t_cthp = v6[i];
          t_klm = v7[i];
          t_klp = v8[i];
          t_integrand = v9[i];
          tt->Fill();
        }
        }
        }
#endif

        // e+ mu-
#ifdef save_neg_lams
        {
        auto const& ret = 
#endif
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_e,pdgcodes::k_muon,-pdgcodes::k_elec,
            mm2, me2, 2, 1, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_pos_nu, rand);
#ifdef save_neg_lams
        if(tt) {
        t_mn = mm;
        t_mp = me;
        auto const& v1 = std::get<0>(ret);
        auto const& v2 = std::get<1>(ret);
        auto const& v3 = std::get<2>(ret);
        auto const& v4 = std::get<3>(ret);
        auto const& v5 = std::get<4>(ret);
        auto const& v6 = std::get<5>(ret);
        auto const& v7 = std::get<6>(ret);
        auto const& v8 = std::get<7>(ret);
        auto const& v9 = std::get<8>(ret);
        for(size_t i = 0; i < std::min(v1.size(),v2.size()); ++i) {
          t_pfn = v1[i];
          t_pfp = v2[i];
          t_snum = v3[i];
          t_snup = v4[i];
          t_cthm = v5[i];
          t_cthp = v6[i];
          t_klm = v7[i];
          t_klp = v8[i];
          t_integrand = v9[i];
          tt->Fill();
        }
        }
        }
#endif
      }
      
#ifdef save_neg_lams
        auto bla = [](auto...){};
        bla(C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj);
#else
      // mu+ mu-
      if(mN > 2. * mm) {
        append_decays((gen_nu_lep_lep|gen_nu_mu_mu)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_mu,pdgcodes::k_muon,-pdgcodes::k_muon,
            mm2, mm2, 2, 2, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_pos_nu, rand);
      }

      // anti-particle
      auto& posanuinfo = decay_branches[-pdgcodes::k_HNL_poshel];
      
      // nu pi0
      if(mN > mpi0 && gen_nu_pi0) {
        posanuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{-pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(-1), 0., mpi0))
                ));
        const int pdgnu = -pdgcodes::k_nu_e;
          std::cout << Form("Adding %sN%s -> nu pi0 ; m= %g Gamma=%g \n",
              pdgnu<0?"~":"", "-", mN,gamma_pi0_nu);
      }
      // e-e+
      if(mN > 2. * me) {
        append_decays((gen_nu_lep_lep|gen_nu_e_e)?posanuinfo:dummybranch, +1, -pdgcodes::k_nu_e,pdgcodes::k_elec,-pdgcodes::k_elec,
            me2, me2, 1, 1, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_pos_anu, rand);
      }

      if(mN > me + mm) {
        // e- mu+
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?posanuinfo:dummybranch, +1, -pdgcodes::k_nu_mu,pdgcodes::k_elec,-pdgcodes::k_muon,
            me2, mm2, 1, 2, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_pos_anu, rand);

        // e+ mu-
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?posanuinfo:dummybranch, +1, -pdgcodes::k_nu_e,pdgcodes::k_muon,-pdgcodes::k_elec,
            mm2, me2, 2, 1, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_pos_anu, rand);
      }
      
      // mu+ mu-
      if(mN > 2. * mm) {
        append_decays((gen_nu_lep_lep|gen_nu_mu_mu)?posanuinfo:dummybranch, +1, -pdgcodes::k_nu_mu,pdgcodes::k_muon,-pdgcodes::k_muon,
            mm2, mm2, 2, 2, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_pos_anu, rand);
      }
    }
    else {
      auto& posnuinfo = decay_branches[pdgcodes::k_HNL_poshel];
      
      // nu pi0
      if(mN > mpi0 && gen_nu_pi0) {
        posnuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                0.5*gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(+1), 0., mpi0))
                ));
        posnuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{-pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                0.5*gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(-1), 0., mpi0))
                ));
          std::cout << Form("Adding N+ -> nu pi0 ; m= %g Gamma=%g \n",
              mN,gamma_pi0_nu);
      }

      // e-e+
      if(mN > 2. * me) {
        append_decays((gen_nu_lep_lep|gen_nu_e_e)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_e,pdgcodes::k_elec,-pdgcodes::k_elec,
            me2, me2, 1, 1, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_pos_nu, rand);
      }

      if(mN > me + mm) {
        // e- mu+
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_mu,pdgcodes::k_elec,-pdgcodes::k_muon,
            me2, mm2, 1, 2, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_pos_nu, rand);

        // e+ mu-
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_e,pdgcodes::k_muon,-pdgcodes::k_elec,
            mm2, me2, 2, 1, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_pos_nu, rand);
      }
      
      // mu+ mu-
      if(mN > 2. * mm) {
        append_decays((gen_nu_lep_lep|gen_nu_mu_mu)?posnuinfo:dummybranch, +1, pdgcodes::k_nu_mu,pdgcodes::k_muon,-pdgcodes::k_muon,
            mm2, mm2, 2, 2, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_pos_nu, rand);
      }
      
    }

    // hnl - neghel
    if(params.is_Dirac()) {
      // particle
      auto& negnuinfo = decay_branches[pdgcodes::k_HNL_neghel];

      const int hel = -1;

      // nu pi0
      if(mN > mpi0 && gen_nu_pi0) {
        negnuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(hel), 0., mpi0))
                ));
        const int pdgnu = pdgcodes::k_nu_e;
          std::cout << Form("Adding %sN%s -> nu pi0 ; m= %g Gamma=%g \n",
              pdgnu<0?"~":"", hel>0?"+":"-", mN,gamma_pi0_nu);
      }
      // e-e+
      if(mN > 2. * me) {
        append_decays((gen_nu_lep_lep|gen_nu_e_e)?negnuinfo:dummybranch, -1, pdgcodes::k_nu_e,pdgcodes::k_elec,-pdgcodes::k_elec,
            me2, me2, 1, 1, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_neg_nu, rand);
      }

      if(mN > me + mm) {
        // e- mu+
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?negnuinfo:dummybranch, -1, pdgcodes::k_nu_mu,pdgcodes::k_elec,-pdgcodes::k_muon,
            me2, mm2, 1, 2, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_neg_nu, rand);

        // e+ mu-
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?negnuinfo:dummybranch, -1, pdgcodes::k_nu_e,pdgcodes::k_muon,-pdgcodes::k_elec,
            mm2, me2, 2, 1, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_neg_nu, rand);
      }
      
      // mu+ mu-
      if(mN > 2. * mm) {
        append_decays((gen_nu_lep_lep|gen_nu_mu_mu)?negnuinfo:dummybranch, -1, pdgcodes::k_nu_mu,pdgcodes::k_muon,-pdgcodes::k_muon,
            mm2, mm2, 2, 2, C1_nu, C2_nu, C3_nu, C4_nu, C5_nu, C6_nu,
            tot_gamma_neg_nu, rand);
      }

      // anti-particle
      auto& neganuinfo = decay_branches[-pdgcodes::k_HNL_neghel];
      
      // nu pi0
      if(mN > mpi0 && gen_nu_pi0) {
        neganuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{-pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(-hel), 0., mpi0))
                ));
        const int pdgnu = -pdgcodes::k_nu_e;
          std::cout << Form("Adding %sN%s -> nu pi0 ; m= %g Gamma=%g \n",
              pdgnu<0?"~":"", -hel>0?"+":"-", mN,gamma_pi0_nu);
      }

      // e-e+
      if(mN > 2. * me) {
        append_decays((gen_nu_lep_lep|gen_nu_e_e)?neganuinfo:dummybranch, -1, -pdgcodes::k_nu_e,pdgcodes::k_elec,-pdgcodes::k_elec,
            me2, me2, 1, 1, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_neg_anu, rand);
      }

      if(mN > me + mm) {
        // e- mu+
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?neganuinfo:dummybranch, -1, -pdgcodes::k_nu_mu,pdgcodes::k_elec,-pdgcodes::k_muon,
            me2, mm2, 1, 2, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_neg_anu, rand);

        // e+ mu-
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?neganuinfo:dummybranch, -1, -pdgcodes::k_nu_e,pdgcodes::k_muon,-pdgcodes::k_elec,
            mm2, me2, 2, 1, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_neg_anu, rand);
      }
      
      // mu+ mu-
      if(mN > 2. * mm) {
        append_decays((gen_nu_lep_lep|gen_nu_mu_mu)?neganuinfo:dummybranch, -1, -pdgcodes::k_nu_mu,pdgcodes::k_muon,-pdgcodes::k_muon,
            mm2, mm2, 2, 2, C1_anu, C2_anu, C3_anu, C4_anu, C5_anu, C6_anu,
            tot_gamma_neg_anu, rand);
      }
    }
    else {
      auto& negnuinfo = decay_branches[pdgcodes::k_HNL_neghel];
      

      // nu pi0
      if(mN > mpi0 && gen_nu_pi0) {
        negnuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                0.5*gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(-1), 0., mpi0))
                ));
        negnuinfo.second.push_back(std::make_tuple(
                decay_prod_vec_t{-pdgcodes::k_nu_e,pdgcodes::k_pion_0},
                0.5*gamma_pi0_nu,
                make_twobody_final_state_function(make_2body_decay_integrand(rand, make_partial_gamma_nu_pi0(+1), 0., mpi0))
                ));
          std::cout << Form("Adding ~N- -> nu pi0 ; m= %g Gamma=%g \n",
               mN,gamma_pi0_nu);
      }

      // e-e+
      if(mN > 2. * me) {
        append_decays((gen_nu_lep_lep|gen_nu_e_e)?negnuinfo:dummybranch, -1, -pdgcodes::k_nu_e,pdgcodes::k_elec,-pdgcodes::k_elec,
            me2, me2, 1, 1, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_neg_nu, rand);
      }

      if(mN > me + mm) {
        // e- mu+
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?negnuinfo:dummybranch, -1, -pdgcodes::k_nu_mu,pdgcodes::k_elec,-pdgcodes::k_muon,
            me2, mm2, 1, 2, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_neg_nu, rand);

        // e+ mu-
        append_decays((gen_nu_lep_lep|gen_nu_e_mu)?negnuinfo:dummybranch, -1, -pdgcodes::k_nu_e,pdgcodes::k_muon,-pdgcodes::k_elec,
            mm2, me2, 2, 1, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_neg_nu, rand);
      }
      
      // mu+ mu-
      if(mN > 2. * mm) {
        append_decays((gen_nu_lep_lep|gen_nu_mu_mu)?negnuinfo:dummybranch, -1, -pdgcodes::k_nu_mu,pdgcodes::k_muon,-pdgcodes::k_muon,
            mm2, mm2, 2, 2, C1_maj, C2_maj, C3_maj, C4_maj, C5_maj, C6_maj,
            tot_gamma_neg_nu, rand);
      }
      
#endif
    }

    if(params.is_Dirac()) {
      std::get<0>(decay_branches[pdgcodes::k_HNL_poshel]) = tot_gamma_pos_nu;
      std::cout << Form("Adding N+ -> X ; m= %g TotGamma=%g lifetime=%g \n",mN,tot_gamma_pos_nu,hnl_tau(pdgcodes::k_HNL_poshel));
      std::get<0>(decay_branches[-pdgcodes::k_HNL_poshel]) = tot_gamma_pos_anu;
      std::cout << Form("Adding ~N+ -> X ; m= %g TotGamma=%g lifetime=%g \n",mN,tot_gamma_pos_anu,hnl_tau(-pdgcodes::k_HNL_poshel));
      std::get<0>(decay_branches[pdgcodes::k_HNL_neghel]) = tot_gamma_neg_nu;
      std::cout << Form("Adding N- -> X ; m= %g TotGamma=%g lifetime=%g \n",mN,tot_gamma_neg_nu,hnl_tau(pdgcodes::k_HNL_neghel));
      std::get<0>(decay_branches[-pdgcodes::k_HNL_neghel]) = tot_gamma_neg_anu;
      std::cout << Form("Adding ~N- -> X ; m= %g TotGamma=%g lifetime=%g \n",mN,tot_gamma_neg_anu,hnl_tau(-pdgcodes::k_HNL_neghel));
      std::sort(std::get<1>(decay_branches[pdgcodes::k_HNL_poshel]).begin(), 
          std::get<1>(decay_branches[pdgcodes::k_HNL_poshel]).end(), [](auto const& a, auto const& b){
          return std::get<1>(a) > std::get<1>(b);
          });
      std::sort(std::get<1>(decay_branches[pdgcodes::k_HNL_neghel]).begin(), 
          std::get<1>(decay_branches[pdgcodes::k_HNL_neghel]).end(), [](auto const& a, auto const& b){
          return std::get<1>(a) > std::get<1>(b);
          });
      std::sort(std::get<1>(decay_branches[-pdgcodes::k_HNL_poshel]).begin(), 
          std::get<1>(decay_branches[-pdgcodes::k_HNL_poshel]).end(), [](auto const& a, auto const& b){
          return std::get<1>(a) > std::get<1>(b);
          });
      std::sort(std::get<1>(decay_branches[-pdgcodes::k_HNL_neghel]).begin(), 
          std::get<1>(decay_branches[-pdgcodes::k_HNL_neghel]).end(), [](auto const& a, auto const& b){
          return std::get<1>(a) > std::get<1>(b);
          });
    }
    else {
      std::get<0>(decay_branches[pdgcodes::k_HNL_poshel]) = tot_gamma_pos_nu;
      std::cout << Form("Adding (~)N+ -> X ; m= %g TotGamma=%g lifetime=%g \n",mN,tot_gamma_pos_nu,hnl_tau(pdgcodes::k_HNL_poshel));
      std::get<0>(decay_branches[pdgcodes::k_HNL_neghel]) = tot_gamma_neg_nu;
      std::cout << Form("Adding (~)N- -> X ; m= %g TotGamma=%g lifetime=%g \n",mN,tot_gamma_neg_nu,hnl_tau(pdgcodes::k_HNL_neghel));
      std::sort(std::get<1>(decay_branches[pdgcodes::k_HNL_poshel]).begin(), 
          std::get<1>(decay_branches[pdgcodes::k_HNL_poshel]).end(), [](auto const& a, auto const& b){
          return std::get<1>(a) > std::get<1>(b);
          });
      std::sort(std::get<1>(decay_branches[pdgcodes::k_HNL_neghel]).begin(), 
          std::get<1>(decay_branches[pdgcodes::k_HNL_neghel]).end(), [](auto const& a, auto const& b){
          return std::get<1>(a) > std::get<1>(b);
          });
    }
  }
}
