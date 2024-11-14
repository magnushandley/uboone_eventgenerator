#include "GenKinematics.h"

#include "CLHEP/Random/RandFlat.h"

//#include "CLHEP/Random/RandomEngine.h"

namespace pdg {
  const int k_elec = 11;
  const int k_muon = 13;
  //const int k_pion_0 = 111;
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
  p_ckm_ts(p.get<double>("ckm_ts",39.4e-3)), // PDG 2019
  p_ckm_ud(p.get<double>("ckm_ud",0.97420)), // PDG 2019
  p_decayconst_pion(p.get<double>("decayconst_pion",0.13041)), //PDG [GeV]
  p_gfermi(p.get<double>("gfermi",1.16637e-5)) //PDG [GeV]
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
  max_path_len = 2.*det_half_dims.Mag();
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



bool hpsgen::GenKinematics::generate(const TLorentzVector& kaon_decay_pos, const TLorentzVector& kaon_4mom, const int kaon_pdg, const double HNL_mass,
 const int fProdLepType, const int fDecayLepType ,const double flux_weight, const double max_weight,
          rng& rand, std::multimap<int,TLorentzVector>& result) const {


  int prod_lep_type = 0;
  int decay_lep_type = 0;
  double split_weight = 1; //weight which accounts for splitting in both prod and decay

  if(fProdLepType==2){
    prod_lep_type = (CLHEP::RandFlat::shoot(&rand) > 0.5  ? 0 : 1); //0 mu 1 e, can change to diff prop later
    split_weight *= 2.; //we are going to split prod in half so account for that  
  } 
  else{
    prod_lep_type = fProdLepType;
  }


  if(fDecayLepType==2){
    decay_lep_type = (CLHEP::RandFlat::shoot(&rand) > 0.5  ? 0 : 1); //0 mu 1 e, can change to diff prop later
    split_weight *= 2.; //we are going to split decay in half so account for that
  } 
  else{
    decay_lep_type = fDecayLepType;
  }


  const double prod_lepmass = (prod_lep_type == 0 ? consts.mass_muon() : consts.mass_elec());
  const double decay_lepmass = (decay_lep_type == 0 ? consts.mass_muon() : consts.mass_elec());
  
  // const double pimass = (pion_type == 0 ? consts.mass_pion_0() : consts.mass_pion_pm()); //delete eventually

  
  if(kaon_4mom.M() - prod_lepmass - HNL_mass <= 0.) return false; //is this allowd
 
  TLorentzVector scalar4mom = gen_random_scalar_mom(HNL_mass, kaon_4mom.M(), prod_lepmass, kaon_4mom, rand);
  double lambdas[2];
  if(!intersects_ray(kaon_decay_pos.Vect(), scalar4mom.Vect(), lambdas)) return false;
  double dk_weight = 0.;
 
  TLorentzVector scalar_dk_pos = gen_random_scalar_decay_pos_uniform(scalar4mom, kaon_decay_pos, rand, lambdas, dk_weight);

  //     rand, lambdas, dk_weight); //can go back to non uniform but wont matter
  if(!pos_inside_detector(scalar_dk_pos)) return false;

  double br_weight = 1;//if e like need to scale down from mu like to e like br  
  if(prod_lep_type==1) br_weight=br_weight*2.488e-5; //https://pdg.lbl.gov/2018/listings/rpp2018-list-K-plus-minus.pdf

  const double kfactor_weight = kfactor(kaon_4mom.M(), prod_lepmass , HNL_mass);


    // const double weight = dk_weight * br_weight * flux_weight;

  // std::cout<<"kfactor_weight*br_weight: "<<kfactor_weight*br_weight<<std::endl; //(1/scalar4mom.Gamma()) is decay prob (shape)

  br_weight=br_weight*kfactor_weight;


  const double decay_len=consts.hbar()/twobodylep_decay_width(HNL_mass, decay_lepmass)
        *consts.speed_light()*scalar4mom.Beta()*scalar4mom.Gamma(); //[cm]

 
  dk_weight=dk_weight/decay_len; //previously decay weight is just the length of path through detector
  const double weight = dk_weight * flux_weight * br_weight * split_weight; 

  if(max_weight > 0.) {
    if(weight > max_weight) {
      throw art::Exception(art::errors::LogicError) << "weight "<<weight<<" > max_weight "<<max_weight<<std::endl
        <<" Modify max_weight (suggested "<<1.1*weight<<") and re-run."<<std::endl;
    }
    if(CLHEP::RandFlat::shoot(&rand, max_weight) > weight) return false;
  }


  result = gen_daughters(kaon_pdg,scalar4mom,decay_lep_type,rand);//needs kaon parent info to know what kind of HNL we have * consts.speed_light()*scalar4mom.Beta()



  result.emplace(0,scalar_dk_pos); 
  result.emplace(pdg::k_scalar, scalar4mom);

  if(max_weight < 0.) {
    result.emplace(99,TLorentzVector{dk_weight,br_weight,flux_weight,weight});
  }

  return true;
}


bool hpsgen::GenKinematics::generate_uniform(const TLorentzVector& kaon_decay_pos, const TLorentzVector& kaon_4mom,const int kaon_pdg, const double HNL_mass,
           rng& rand, std::multimap<int,TLorentzVector>& result) const {
  const double max_weight = geo->max_path_length();


  const int prod_lep_type = (CLHEP::RandFlat::shoot(&rand) > 0.5  ? 0 : 1); //0 mu 1 e, can change to diff prop later

  
  // const double pimass = (pion_type == 0 ? consts.mass_pion_0() : consts.mass_pion_pm()); //delete eventually

  const double lepmass = (prod_lep_type == 0 ? consts.mass_muon() : consts.mass_elec());
  if(kaon_4mom.M() - lepmass - HNL_mass <= 0.) return false; //is this allowd

  TLorentzVector scalar4mom = gen_random_scalar_mom(HNL_mass, kaon_4mom.M(), lepmass, kaon_4mom, rand);
  double lambdas[2];
  if(!intersects_ray(kaon_decay_pos.Vect(), scalar4mom.Vect(), lambdas)) return false;
  double dk_weight = 0.;
  TLorentzVector scalar_dk_pos = gen_random_scalar_decay_pos_uniform(scalar4mom, kaon_decay_pos, rand, lambdas, dk_weight);
  if(!pos_inside_detector(scalar_dk_pos)) return false;
  
  const double weight = dk_weight; // in this case, dk_weight is just the path length

  if(weight > max_weight) {
    throw art::Exception(art::errors::LogicError) << "weight "<<weight<<" > max_weight "<<max_weight<<std::endl
      <<" Modify max_weight (suggested "<<1.1*weight<<") and re-run."<<std::endl;
  }
  if(CLHEP::RandFlat::shoot(&rand, max_weight) > weight) return false;
  const int decay_lep_type = (CLHEP::RandFlat::shoot(&rand) > 0.5  ? 0 : 1); //0 mu 1 e, can change to diff prop later
  result = gen_daughters(kaon_pdg,scalar4mom,decay_lep_type, rand);
  result.emplace(0,scalar_dk_pos);
  result.emplace(pdg::k_scalar, scalar4mom);

  return true;
}

TLorentzVector hpsgen::GenKinematics::gen_random_scalar_mom(const double HNL_mass, const double kaon_mass, const double lep_mass,
    const TLorentzVector& kaon4mom, rng& rand) const {
  if(kaon_mass - lep_mass - HNL_mass <= 0.) {
    throw art::Exception(art::errors::LogicError) << "scalar mass too heavy K " << kaon_mass << " lep "<<lep_mass <<" S " << HNL_mass<<std::endl;
  }
  const double costh = CLHEP::RandFlat::shoot(&rand, -1, 1);
  const double sinth = std::sqrt(1.-costh*costh);
  double phi = CLHEP::RandFlat::shoot(&rand, 2*M_PI);

  const double totE = (kaon_mass*kaon_mass-lep_mass*lep_mass+HNL_mass*HNL_mass)/(2.*kaon_mass);
  const double mom = std::sqrt(totE*totE-HNL_mass*HNL_mass);

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

TLorentzVector hpsgen::GenKinematics::gen_random_scalar_decay_pos_uniform(const TLorentzVector& scalar4mom, const TLorentzVector& kaonpos,
    rng& rand, const double* lambdas, double& weight) const {
  const double speed = scalar4mom.Beta() * consts.speed_light();
  const double a = std::min(lambdas[0],lambdas[1]);
  const double b = std::max(lambdas[0],lambdas[1]);
  weight = b-a;
  const double p0 = CLHEP::RandFlat::shoot(&rand);
  const double length = a + weight*p0;
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



std::multimap<int,TLorentzVector> hpsgen::GenKinematics::gen_daughters(const int kaon_pdg, const TLorentzVector& parent_mom, const int decay_lep_type, rng& rand) const {
  const double HNL_mass = parent_mom.M();

  int pdg, anti_pdg;
  double lep_mass;
  
  // const int d_prod_lep_type = (CLHEP::RandFlat::shoot(&rand) > 0.5  ? 0 : 1); //0 mu 1 e, can change to diff prop later
  int kaon_sign = (kaon_pdg > 0 ? 1 : -1); //positive kaons produce lep- pi+ (postive pdg)
  if(decay_lep_type==1) {
      pdg = pdg::k_elec*kaon_sign;
      anti_pdg = pdg::k_pion_pm*kaon_sign;
      lep_mass = consts.mass_elec();
    }
  else{
      pdg = pdg::k_muon*kaon_sign;
      anti_pdg = pdg::k_pion_pm*kaon_sign;
      lep_mass = consts.mass_muon();
    }


  const double maj = CLHEP::RandFlat::shoot(&rand);
  if(maj < 0.5) { //flip in half of cases, should make this fhicl config
      pdg = pdg*-1;
      anti_pdg = anti_pdg*-1;
    }

  const double phi = CLHEP::RandFlat::shoot(&rand, 2.*M_PI);
  const double costh = CLHEP::RandFlat::shoot(&rand, -1., 1.);
  const double th = std::acos(costh);


  const double en_lep = ((HNL_mass*HNL_mass) + (lep_mass*lep_mass) - (consts.mass_pion_pm() * consts.mass_pion_pm())) / (2.*HNL_mass);
  const double en_pi = (HNL_mass-en_lep);
  const double mom = std::sqrt(en_lep*en_lep - lep_mass*lep_mass); //equal opposite in rf


  TLorentzVector d1(mom*cos(phi)*sin(th), mom*sin(phi)*sin(th), mom*cos(th), en_lep);
  TLorentzVector d2(-mom*cos(phi)*sin(th), -mom*sin(phi)*sin(th), -mom*cos(th), en_pi);
  d1.Boost(parent_mom.BoostVector());
  d2.Boost(parent_mom.BoostVector());
  return {{pdg,d1},{anti_pdg,d2}};
}


double hpsgen::GenKinematics::twobodylep_decay_width(const double HNL_mass, const double lep_mass) const {
  auto kallen_lambda = [](double a, double b, double c)->double {
    return a*a + b*b + c*c -2*a*b -2*a*c -2*b*c;
  };
  auto I1_factor = [kallen_lambda](double x,double y)->double{
    return ((1+x-y)*(1+x)-(4*x))*std::sqrt(kallen_lambda(1,x,y));
  };


  const double Vud = consts.ckm_ud();
  const double decayconst = consts.decayconst_pion();
  const double gFermi = consts.gfermi();
  const double pimass = consts.mass_pion_pm();

  if(HNL_mass - lep_mass - pimass <= 0.) return false;
  // std::cout<<"I1_factor: "<<I1_factor(lep_mass*lep_mass/HNL_mass/HNL_mass,pimass*pimass/HNL_mass/HNL_mass)<<std::endl;
  double width = 2*(gFermi*gFermi)/(16*M_PI)
    *(decayconst*decayconst)
    *(Vud*Vud)*std::pow(HNL_mass,3)
    *I1_factor(lep_mass*lep_mass/HNL_mass/HNL_mass,pimass*pimass/HNL_mass/HNL_mass); //times two for majoran
  return width;
}




double hpsgen::GenKinematics::kfactor(const double meson_mass,const double lep_mass, const double HNL_mass) const {
  auto kallen_lambda = [](double a, double b, double c)->double {
    return a*a + b*b + c*c -2*a*b -2*a*c -2*b*c;
  };
  auto Rho = [kallen_lambda](double a, double b)->double {
    double fl = kallen_lambda(1.,a,b);
    if(fl<0) return 0.;
    double result = (a+b-(a-b)*(a-b))*std::sqrt(fl);
    if(result>0) return result;
    else return 0.;
  };

  double num = Rho((lep_mass*lep_mass)/(meson_mass*meson_mass),(HNL_mass*HNL_mass)/(meson_mass*meson_mass));
  double den = (lep_mass*lep_mass)/(meson_mass*meson_mass) * (1.-(lep_mass*lep_mass)/(meson_mass*meson_mass)) * (1.-(lep_mass*lep_mass)/(meson_mass*meson_mass));

  return num/den;
}

