#include "GenKinematics.h"

#include "CLHEP/Random/RandFlat.h"

//#include "CLHEP/Random/RandomEngine.h"

namespace pdg {
  const int k_elec = 11;
  const int k_muon = 13;
  const int k_pion_0 = 111;
  const int k_pion_pm = 211;
  const int k_scalar = 54;
}

// units are [GeV] [ns] [cm]
hpsgen::PhysicalConstants::PhysicalConstants(fhicl::ParameterSet const& p) 
: p_mass_elec(p.get<double>("mass_elec",0.5109989461e-3)), // PDG 2019
  p_mass_muon(p.get<double>("mass_muon",0.1056583745)), // PDG 2019
  p_mass_pion_pm(p.get<double>("mass_pion_pm",0.13957061)), // PDG 2019
  p_mass_pion_0(p.get<double>("mass_pion_0",0.1349770)), // PDG 2019
  p_mass_top(p.get<double>("mass_top",172.9)), // PDG 2019
  p_lifetime_kaon_0(p.get<double>("lifetime_kaon_0",51.16)), // PDG 2019
  p_lifetime_kaon_pm(p.get<double>("lifetime_kaon_pm",12.38)), // PDG 2019
  p_speed_light(p.get<double>("speed_light",29.9792458)), // PDG, exact
  p_higgs_vev(p.get<double>("higgs_vev",246.22)), // PDG
  p_hbar(p.get<double>("hbar",6.582119569e-16)), // PDG, exact
  p_ckm_td(p.get<double>("ckm_td",8.1e-3)), // PDG 2019
  p_ckm_ts(p.get<double>("ckm_ts",39.4e-3)) // PDG 2019
{
}

hpsgen::Geo::Geo(art::ServiceHandle<geo::Geometry const>& g) {
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

hpsgen::GenKinematics::GenKinematics(fhicl::ParameterSet const& p) 
  : consts(p), geo(0)
{
}

hpsgen::GenKinematics::~GenKinematics() 
{
  if(geo) delete geo;
}


bool hpsgen::GenKinematics::generate(const TLorentzVector& kaon_decay_pos, const TLorentzVector& kaon_4mom, const double scalar_mass,
          const double model_theta, const int pion_type, const double flux_weight, const double max_weight,
          rng& rand, std::multimap<int,TLorentzVector>& result) const {
  const double pimass = (pion_type == 0 ? consts.mass_pion_0() : consts.mass_pion_pm());
  TLorentzVector scalar4mom = gen_random_scalar_mom(scalar_mass, kaon_4mom.M(), pimass, kaon_4mom, rand);
  double lambdas[2];
  if(!intersects_ray(kaon_decay_pos.Vect(), scalar4mom.Vect(), lambdas)) return false;
  double dk_weight = 0.;
  TLorentzVector scalar_dk_pos = gen_random_scalar_decay_pos(scalar4mom, kaon_decay_pos, scalar_tau(scalar_mass, model_theta),
      rand, lambdas, dk_weight);
  if(!pos_inside_detector(scalar_dk_pos)) return false;
  
  const double br_weight = branching_ratio(kaon_4mom.M(), scalar_mass, model_theta, pion_type);

  const double weight = dk_weight * br_weight * flux_weight;

  if(max_weight > 0.) {
    if(weight > max_weight) {
      throw art::Exception(art::errors::LogicError) << "weight "<<weight<<" > max_weight "<<max_weight<<std::endl
        <<" Modify max_weight (suggested "<<1.1*weight<<") and re-run."<<std::endl;
    }
    if(CLHEP::RandFlat::shoot(&rand, max_weight) > weight) return false;
  }

  result = gen_daughters(scalar4mom, model_theta, rand);
  result.emplace(0,scalar_dk_pos);
  result.emplace(pdg::k_scalar, scalar4mom);

  if(max_weight < 0.) {
    result.emplace(99,TLorentzVector{dk_weight,br_weight,flux_weight,weight});
  }

  return true;
}

TLorentzVector hpsgen::GenKinematics::gen_random_scalar_mom(const double scalar_mass, const double kaon_mass, const double pion_mass,
    const TLorentzVector& kaon4mom, rng& rand) const {
  if(kaon_mass - pion_mass - scalar_mass <= 0.) {
    throw art::Exception(art::errors::LogicError) << "scalar mass too heavy K " << kaon_mass << " pi "<<pion_mass <<" S " << scalar_mass<<std::endl;
  }
  const double costh = CLHEP::RandFlat::shoot(&rand, -1, 1);
  const double sinth = std::sqrt(1.-costh*costh);
  double phi = CLHEP::RandFlat::shoot(&rand, 2*M_PI);

  const double totE = (kaon_mass*kaon_mass-pion_mass*pion_mass+scalar_mass*scalar_mass)/(2.*kaon_mass);
  const double mom = std::sqrt(totE*totE-scalar_mass*scalar_mass);

  TLorentzVector s4mom(mom*sinth*std::sin(phi),mom*sinth*std::cos(phi),mom*costh,totE);

  s4mom.Boost(kaon4mom.BoostVector());

  return s4mom;
}

TLorentzVector hpsgen::GenKinematics::gen_random_scalar_decay_pos(const TLorentzVector& scalar4mom, const TLorentzVector& kaonpos,
    const double model_tau, rng& rand, const double* lambdas, double& weight) const {
  const double gamma = scalar4mom.Gamma();
  const double speed = scalar4mom.Beta() * consts.speed_light();
  const double a = std::min(lambdas[0],lambdas[1])/speed;
  const double b = std::max(lambdas[0],lambdas[1])/speed;
  const double probA =  std::exp(-a/gamma/model_tau);
  const double probB =  std::exp(-b/gamma/model_tau);
  weight = probA - probB; // integral along exponential;
  const double p0 = CLHEP::RandFlat::shoot(&rand);
  const double length = (a + b - gamma*model_tau*std::log(std::exp(b/gamma/model_tau)*(1-p0) + std::exp(a/gamma/model_tau)*p0))*speed;
  const TVector3 traj = length * scalar4mom.Vect().Unit();
  const double time_lab = traj.Mag() / speed;
  const TLorentzVector traj4(traj,time_lab);
  return kaonpos+traj4;
}

void hpsgen::GenKinematics::update_geometry(art::ServiceHandle<geo::Geometry const>& g) {
  if(geo) {
    delete geo;
  }
  geo = new Geo(g);
}

// there may be better ways of doing this, using geo::Geometry directly?
bool hpsgen::GenKinematics::intersects_ray(const TVector3& orig, const TVector3& dir, double* lambdas) const {
  const TVector3& det_centre = geo->centre(); const TVector3& det_half_dims = geo->half_dims();
  const TVector3& unit_dir = dir.Unit();
  int n_intersects = 0;
  unsigned int ilam = 0;
  for(int coord = 0; coord < 3; ++coord) {
    if(std::abs(unit_dir[coord])>0.) {
      for(int side = -1; side < 2; side += 2) { // side = -1 or +1
        const double plane = det_centre[coord] + side * det_half_dims[coord];
        const double lambda = (plane - orig[coord])/unit_dir[coord];
        if(lambda < 0) continue; // no backwards-going scalars
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
bool hpsgen::GenKinematics::pos_inside_detector(const TLorentzVector& scalar_dk_pos) const {
  const TVector3& det_centre = geo->centre(); const TVector3& det_half_dims = geo->half_dims();
  bool intersects_vol = true;
  for(int coord = 0; coord < 3; ++coord) {
    const double pos = scalar_dk_pos.Vect()[coord];
    if(std::abs(pos - det_centre[coord]) > det_half_dims[coord]) {
      intersects_vol = false;
      break;
    }
  }
  return intersects_vol;
}

std::multimap<int,TLorentzVector> hpsgen::GenKinematics::gen_daughters(const TLorentzVector& parent_mom, const double model_theta, rng& rand) const {
  const double scalar_mass = parent_mom.M();
  const double gammas[4] = {partial_gamma_lep(scalar_mass,consts.mass_elec(),model_theta),
    partial_gamma_lep(scalar_mass,consts.mass_muon(),model_theta),
    partial_gamma_pi(scalar_mass,consts.mass_pion_0(),model_theta),
    2*partial_gamma_pi(scalar_mass,consts.mass_pion_pm(),model_theta)};
  const double totgamma = gammas[0]+gammas[1]+gammas[2]+gammas[3];
  int pdg, anti_pdg;
  double md;
  if(totgamma > gammas[0]) {
    const double r = CLHEP::RandFlat::shoot(&rand);
    if(r < gammas[0]/totgamma) {
      pdg = pdg::k_elec;
      anti_pdg = -pdg;
      md = consts.mass_elec();
    }
    else if(r < (gammas[0]+gammas[1])/totgamma) {
      pdg = pdg::k_muon;
      anti_pdg = -pdg;
      md = consts.mass_muon();
    }
    else if(r < (gammas[0]+gammas[1]+gammas[2])/totgamma) {
      pdg = pdg::k_pion_0;
      anti_pdg = pdg; // pi0 particle==antiparticle
      md = consts.mass_pion_0();
    }
    else {
      pdg = pdg::k_pion_pm;
      anti_pdg = -pdg;
      md = consts.mass_pion_pm();
    }
  }
  else {
    pdg = pdg::k_elec;
    anti_pdg = -pdg;
    md = consts.mass_elec();
  }
  const double phi = CLHEP::RandFlat::shoot(&rand, 2.*M_PI);
  const double costh = CLHEP::RandFlat::shoot(&rand, -1., 1.);
  const double th = std::acos(costh);
  const double ene = scalar_mass / 2.;
  const double mom = std::sqrt(ene*ene - md*md);
  TLorentzVector d1(mom*cos(phi)*sin(th), mom*sin(phi)*sin(th), mom*cos(th), ene);
  TLorentzVector d2(-mom*cos(phi)*sin(th), -mom*sin(phi)*sin(th), -mom*cos(th), ene);
  d1.Boost(parent_mom.BoostVector());
  d2.Boost(parent_mom.BoostVector());
  return {{pdg,d1},{anti_pdg,d2}};
}

double hpsgen::GenKinematics::partial_gamma_pi(const double scalar_mass, const double pi_mass, const double model_theta) const {
  if(scalar_mass < 2.*pi_mass) return 0.;
  // eq 4 of arXiv:1909.11670v1 without factor of three in orginal document.
  const double g = 2.*scalar_mass*scalar_mass/9. + 11.*pi_mass*pi_mass/9.;
  return model_theta*model_theta*g*g/(32.*M_PI*consts.higgs_vev()*consts.higgs_vev()*scalar_mass)*std::sqrt(1.-4.*pi_mass*pi_mass/scalar_mass/scalar_mass);
}

double hpsgen::GenKinematics::partial_gamma_lep(const double scalar_mass, const double lep_mass, const double model_theta) const {
  if(scalar_mass < 2.*lep_mass) return 0.;
  // eq 3 of arXiv:1909.11670v1
  return model_theta*model_theta * lep_mass*lep_mass*scalar_mass/(8.*M_PI*consts.higgs_vev()*consts.higgs_vev())*std::pow(1.-4.*lep_mass*lep_mass/scalar_mass/scalar_mass,1.5);
}

double hpsgen::GenKinematics::scalar_tau(const double scalar_mass, const double model_theta) const {
  const double gamma = partial_gamma_lep(scalar_mass,consts.mass_elec(),model_theta)
    + partial_gamma_lep(scalar_mass,consts.mass_muon(),model_theta)
    + partial_gamma_pi(scalar_mass,consts.mass_pion_0(),model_theta)
    + 2*partial_gamma_pi(scalar_mass,consts.mass_pion_pm(),model_theta);
  return consts.hbar() / gamma;
}

double hpsgen::GenKinematics::branching_ratio(const double kmass, const double scalar_mass, const double model_theta, const int ktype) const {
  auto kallen_lambda = [](double a, double b, double c)->double {
    return a*a + b*b + c*c -2*a*b -2*a*c -2*b*c;
  };
  const double vev = consts.higgs_vev();
  const double Vtd = consts.ckm_td();
  const double Vts = consts.ckm_ts();
  const double tmass = consts.mass_top();
  const double k_lifetime = (ktype == 0 ? consts.lifetime_kaon_0() : consts.lifetime_kaon_pm());
  const double pimass = (ktype==0 ? consts.mass_pion_0() : consts.mass_pion_pm());
  // equation (5) of arXiv:1909.11670v1
  const double branch_ratio = model_theta * model_theta / (16*M_PI*kmass)
    * std::pow(3 * Vtd * Vts * tmass*tmass * kmass*kmass / (32*M_PI*M_PI * vev*vev*vev), 2) 
    * std::sqrt(kallen_lambda(1., scalar_mass*scalar_mass/kmass/kmass, pimass*pimass/kmass/kmass))
    / (consts.hbar() / k_lifetime);
  return branch_ratio;
}
