#include "FluxReaderBNB.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "GENIE/Framework/ParticleData/PDGUtils.h"

hpsgen::FluxReaderBNB::FluxReaderBNB(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& fRNG)
  : FluxReader(p, fRNG, "h101", "10000", "h120"),
  beampos(1.24325*100, -0.0093*100, -463.363525*100), // /* from /uboone/app/users/zarko/windowRotation/correctWindow.c */
  beamtime(-1561.8610), // time of flight of neutrinos from t0=POT to detector
  rot([](){
      /* beam rotation, from /uboone/app/users/zarko/windowRotation/correctWindow.c */
      TRotation r;
      r.RotateX(0.016/10.712);
      r.RotateY(0.036/10.712);
      return r;
      }())
{
  set_branch(&fBranch.beamwgt, "beamwgt");
  set_branch(&fBranch.ntp,     "ntp");
  set_branch(&fBranch.npart,   "npart");
  set_branch( fBranch.id,      "id");
  set_branch( fBranch.ini_pos, "ini_pos");
  set_branch( fBranch.ini_mom, "ini_mom");
  set_branch( fBranch.ini_eng, "ini_eng");
  set_branch( fBranch.ini_t,   "ini_t");
  set_branch( fBranch.fin_mom, "fin_mom");
  set_branch( fBranch.fin_pol, "fin_pol");
}

hpsgen::FluxReaderBNB::~FluxReaderBNB() {
}

bool hpsgen::FluxReaderBNB::get_kaon_from_flux(TLorentzVector& kmom, TLorentzVector& kpos, int& kpdg, int& pi_type, double& weight) {
  if(!(fBranch.npart > 1 && fBranch.id[1]>=10 && fBranch.id[1] <= 12)) return false; // only select kaon decays
  // id 10 == K0L; id 11 = K+, id 12 = K-
  kpdg = (fBranch.id[1] == 10 ? 130 : (fBranch.id[1] == 12 ? -321 : 321));
  pi_type = (fBranch.id[1] == 10 ? 0 : (fBranch.id[1] == 12 ? -1 : 1));
  kpos.SetVect(rot * TVector3(fBranch.ini_pos[0][0],fBranch.ini_pos[0][1],fBranch.ini_pos[0][2]) + beampos);
  kpos.SetT(fBranch.ini_t[0] + beamtime);
  // there is no fin_eng, so we need a roundabout method for the final kaon momentum
  const double kmass = TLorentzVector(fBranch.ini_mom[1][0],fBranch.ini_mom[1][1],fBranch.ini_mom[1][2],fBranch.ini_eng[1]).M();
  kmom.SetVectM(rot * TVector3(fBranch.fin_mom[1][0],fBranch.fin_mom[1][1],fBranch.fin_mom[1][2]), kmass);
  if(kmom.M() < 0.4) {
    kmom.Print();
    throw cet::exception("LogicError")<<"Found a massless kaon! mass="<<kmass<<std::endl;
  }

  const double Necm = [&]() {
    double Nenergy = fBranch.ini_eng[0] ;
    double Ndxdz   = fBranch.ini_mom[0][0]/fBranch.ini_mom[0][2] ;
    double Ndydz   = fBranch.ini_mom[0][1]/fBranch.ini_mom[0][2] ;

    double ppenergy = fBranch.ini_eng[1];
    double pdPx     = fBranch.fin_mom[1][0] ;
    double pdPy     = fBranch.fin_mom[1][1] ;
    double pdPz     = fBranch.fin_mom[1][2] ;

    double ppdxdz   = fBranch.ini_mom[1][0]/fBranch.ini_mom[1][2] ;
    double ppdydz   = fBranch.ini_mom[1][1]/fBranch.ini_mom[1][2] ;
    double pppz     = fBranch.ini_mom[1][2] ;
    double parent_mass=sqrt(ppenergy*ppenergy-
        pppz*pppz*(ppdxdz*ppdxdz +
          ppdydz*ppdydz +
          1.));

    double parent_energy = sqrt(pdPx*pdPx +
        pdPy*pdPy +
        pdPz*pdPz + 
        parent_mass*parent_mass);
    double gamma         = parent_energy / parent_mass;
    double beta[3];
    beta[0] = pdPx/parent_energy;
    beta[1] = pdPy/parent_energy;
    beta[2] = pdPz/parent_energy;

    double partial = fBranch.ini_mom[0][2] * gamma * ( beta[0] * Ndxdz + 
        beta[1] * Ndydz + 
        beta[2] );

    double Necm = gamma * Nenergy - partial;
    return Necm;
  }();

  weight = fBranch.beamwgt  / K2nu_branching_ratio(fBranch.id[1],fBranch.ntp,kmass,Necm);

  return true;
}


// branching ratios taken from the GEANT4 version used by bnb beamline sim in MCC9,
// located in g4bnb/src/BooNEDecayPhysics.cc
double hpsgen::FluxReaderBNB::K2nu_branching_ratio(const int id, const int ntp, const double parent_mass, const double Necm) {
  /* ndecay:
   *
   1  KL0 → νe  + π− + e+
   2  KL0 → ν ̄e + π+ + e−
   3  KL0 → νμ  + π− + μ+
   4  KL0 → ν ̄μ + π+ + μ−
   5  K+  → νμ  + μ+
   6  K+  → νe  + π0 + e+
   7  K+  → νμ  + π0 + μ+
   8  K−  → ν ̄μ + μ−
   9  K−  → ν ̄e + π0 + e−
   10 K−  → ν ̄μ + π0 + μ−
   11 μ+  → ν ̄μ + νe + e+
   12 μ−  → ν  + ν ̄e + e−
   13 π+  → νμ  + μ+
   14 π−  → ν ̄μ + μ−
   */
  const double muon_mass = 0.105658389;
  if (id == 10 && ntp == 1) {
    //ndecay = 1;
    return 0.20275; // over 2 because BR=0.4055 is to both charge states
  }
  else if (id == 10 && ntp == 2) {
    //ndecay = 2;
    return 0.20275;
  }
  else if (id == 10 && ntp == 3) {
    //ndecay = 3;
    return 0.1352;
  }
  else if (id == 10 && ntp == 4) {
    //ndecay = 4;
    return 0.1352;
  }
  else if (id == 11 && ntp == 3) {
    //check if it is a two or three body decay
    if (fabs((parent_mass*parent_mass-muon_mass*muon_mass)/(2.*parent_mass)-Necm)/Necm <= 0.001){
      //two body decay (numu + mu+)
      //ndecay = 5;
      return 0.6356;
    }
    else{
      //three body decay (numu + pi0 + mu+)
      //ndecay = 7;
      return 0.03352;
    }
  } else if (id == 11 && ntp == 1) {
    //ndecay = 6;
    return 0.0507;
  }
  else if (id == 12 && ntp == 4) {
    if (fabs((parent_mass*parent_mass-muon_mass*muon_mass)/(2.*parent_mass)-Necm)/Necm <= 0.001){
      //two body decay (numu + mu+)
      //ndecay = 8;
      return 0.6356;
    }
    else{
      //three body decay (numu + pi0 + mu+)
      //ndecay = 10;
      return 0.03352;
    }
  } else if (id == 12 && ntp == 2) {
    //ndecay = 9;
    return 0.0507;
  }
  else if (id == 5 ) {
    //ndecay = 11;
    return 1.0;
  }
  else if (id == 6 ) {
    //ndecay = 12;
    return 1.0;
  }
  else if (id == 8 ) {
    //ndecay = 13;
    return 1. - 1.23e-4;
  }
  else if (id == 9 ) {
    //ndecay = 14;
    return 1. - 1.23e-4;
  }
  return 1e100;   // shouldn't get here, but if so, weight will be highly suppressed
}

// mostly copied from BooNEtoGSimple.cxx
void hpsgen::FluxReaderBNB::get_MCFlux(simb::MCFlux& flux) {
  get_current_entry();
  flux.Reset();
  flux.fFluxType = simb::kSimple_Flux;

  // maintained variable names from gnumi ntuples
  // see http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/[/v19/output_gnumi.html]


  if ( fBranch.ntp == 1 ){
    flux.fntype = 12; //nue
    //fentry->ntype = 12; //nue
  }
  else if ( fBranch.ntp == 2 ){
    flux.fntype = -12; //nuebar
    //	fentry->ntype = -12; //nue
  }
  else if ( fBranch.ntp == 3 ){
    flux.fntype = 14; //numu
    //	fentry->ntype = 14; //nue
  }
  else if ( fBranch.ntp == 4 ){
    flux.fntype = -14; //numubar
    //	fentry->ntype = -14; //nue
  }
  else{
    std::cerr<<"Neutrino type not recognized!!! ntp = "<< fBranch.ntp <<std::endl;
  }
  flux.fdk2gen = 0.; // L of neutrino. We don't care.
  flux.fnenergyn = flux.fnenergyf = fBranch.ini_eng[0];

  flux.frun      = 0;
  flux.fevtno    = 0; // not needed
  flux.ftpx      = fBranch.ini_mom[fBranch.npart-2][0];
  flux.ftpy      = fBranch.ini_mom[fBranch.npart-2][1];
  flux.ftpz      = fBranch.ini_mom[fBranch.npart-2][2];
  flux.ftptype   = fBranch.id[fBranch.npart-2] != 0 ? genie::pdg::GeantToPdg(fBranch.id[fBranch.npart-2]) : 0;   // converted to PDG
  flux.fvx       = fBranch.ini_pos[0][0];
  flux.fvy       = fBranch.ini_pos[0][1];
  flux.fvz       = fBranch.ini_pos[0][2];

  const double Nenergy = fBranch.ini_eng[0] ;
  const double Ndxdz   = fBranch.ini_mom[0][0]/fBranch.ini_mom[0][2] ;
  const double Ndydz   = fBranch.ini_mom[0][1]/fBranch.ini_mom[0][2] ;
  //double Npz     = fBranch.ini_mom[0][2] ;

  const double ppenergy = fBranch.ini_eng[1];
  const double pdPx     = fBranch.fin_mom[1][0] ;
  const double pdPy     = fBranch.fin_mom[1][1] ;
  const double pdPz     = fBranch.fin_mom[1][2] ;

  //double ppvx     = fBranch.ini_pos[1][0] ;
  //double ppvy     = fBranch.ini_pos[1][1] ;
  //double ppvz     = fBranch.ini_pos[1][2] ;    

  const double ppdxdz   = fBranch.ini_mom[1][0]/fBranch.ini_mom[1][2] ;
  const double ppdydz   = fBranch.ini_mom[1][1]/fBranch.ini_mom[1][2] ;
  const double pppz     = fBranch.ini_mom[1][2] ;

  //Get the neutrino energy in the parent decay cm
  const double parent_mass=sqrt(ppenergy*ppenergy-
      pppz*pppz*(ppdxdz*ppdxdz +
        ppdydz*ppdydz +
        1.));

  const double parent_energy = sqrt(pdPx*pdPx +
      pdPy*pdPy +
      pdPz*pdPz + 
      parent_mass*parent_mass);

  const double gamma         = parent_energy / parent_mass;
  double beta[3];
  beta[0] = pdPx/parent_energy;
  beta[1] = pdPy/parent_energy;
  beta[2] = pdPz/parent_energy;

  const double partial = fBranch.ini_mom[0][2] * gamma * ( beta[0] * Ndxdz + 
      beta[1] * Ndydz + 
      beta[2] );

  const double Necm = gamma * Nenergy - partial;
  if (fBranch.id[1] == 10 && fBranch.ntp == 1) 
    flux.fndecay = 1;
  else if (fBranch.id[1] == 10 && fBranch.ntp == 2) 
    flux.fndecay = 2;
  else if (fBranch.id[1] == 10 && fBranch.ntp == 3) 
    flux.fndecay = 3;
  else if (fBranch.id[1] == 10 && fBranch.ntp == 4) 
    flux.fndecay = 4;
  else if (fBranch.id[1] == 11 && fBranch.ntp == 3) {
    //check if it is a two or three body decay
    if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-Necm)/Necm <= 0.001)
      //two body decay (numu + mu+)
      flux.fndecay = 5;
    else
      //three body decay (numu + pi0 + mu+)
      flux.fndecay = 7;
  } else if (fBranch.id[1] == 11 && fBranch.ntp == 1) 
    flux.fndecay = 6;
  else if (fBranch.id[1] == 12 && fBranch.ntp == 4) {
    if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-Necm)/Necm <= 0.001)
      //two body decay (numu + mu+)
      flux.fndecay = 8;
    else
      //three body decay (numu + pi0 + mu+)
      flux.fndecay = 10;
  } else if (fBranch.id[1] == 12 && fBranch.ntp == 2) 
    flux.fndecay = 9;
  else if (fBranch.id[1] == 5 ) 
    flux.fndecay = 11;
  else if (fBranch.id[1] == 6 ) 
    flux.fndecay = 12;
  else if (fBranch.id[1] == 8 ) 
    flux.fndecay = 13;
  else if (fBranch.id[1] == 9 ) 
    flux.fndecay = 14;
  flux.fppmedium = 0;

  flux.fpdpx     = fBranch.fin_mom[1][0];
  flux.fpdpy     = fBranch.fin_mom[1][1];
  flux.fpdpz     = fBranch.fin_mom[1][2];

  double apppz = fBranch.ini_mom[1][2];
  if ( TMath::Abs(fBranch.ini_mom[1][2]) < 1.0e-30 ) apppz = 1.0e-30;
  flux.fppdxdz   = fBranch.ini_mom[1][0] / apppz;
  flux.fppdydz   = fBranch.ini_mom[1][1] / apppz;
  flux.fpppz     = fBranch.ini_mom[1][2];

  flux.fptype    = fBranch.id[1] != 0 ? genie::pdg::GeantToPdg(fBranch.id[1]) : 0;;


  // anything useful stuffed into vdbl or vint?
  // need to check the metadata  auxintname, auxdblname

  flux.ftgen     = fBranch.npart;
  flux.ftgptype  = -9999.;


  if ( fBranch.id[1] == 5 ||
      fBranch.id[1] == 6) {
    flux.fmupare  = fBranch.ini_eng[2];
    flux.fmuparpx = fBranch.fin_mom[2][0];
    flux.fmuparpy = fBranch.fin_mom[2][1];
    flux.fmuparpz = fBranch.fin_mom[2][2];
  } else {
    flux.fmupare  = -9999.;
    flux.fmuparpx = -9999.;
    flux.fmuparpy = -9999.;
    flux.fmuparpz = -9999.;
  }
  flux.fnecm     = Necm;
  flux.fnimpwt   = fBranch.beamwgt;
  flux.fnwtnear = flux.fnwtfar = -9999.;




}
