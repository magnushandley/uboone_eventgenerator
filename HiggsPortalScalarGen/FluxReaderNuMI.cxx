#include "FluxReaderNuMI.h"

#include "dk2nu/tree/dk2nu.h"
#include "nusimdata/SimulationBase/MCFlux.h"

hpsgen::FluxReaderNuMI::FluxReaderNuMI(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& fRNG)
  : FluxReader(p, fRNG, "dk2nuTree", "pots", "dkmetaTree"), dk2nu(0),
  beampos(-31387.58422, -3316.402543, -60100.2414), // NuMI flux note
  beamtime(-2276.3), // time of flight of neutrinos from t0=POT to detector
  rot([](){
      TRotation R;
      // Rotation matrix using the 0,0,0 position for MicroBooNE (beam to det input) // From NuMI flux note
      const TVector3 Rot_row_x = {  0.92103853804025681562,    0.022713504803924120662,  0.38880857519374290021  };
      const TVector3 Rot_row_y = {  4.6254001262154668408e-05, 0.99829162468141474651,  -0.058427989452906302359 };
      const TVector3 Rot_row_z = { -0.38947144863934973769,    0.053832413938664107345,  0.91946400794392302291  };

      R.RotateAxes(Rot_row_x, Rot_row_y, Rot_row_z); // Also inverts so now to det to beam
      R.Invert(); // go back to beam to det
      return R;
      }())
{
  set_branch(&dk2nu, "dk2nu");
}

hpsgen::FluxReaderNuMI::~FluxReaderNuMI() {
  if(dk2nu) delete dk2nu;
}

bool hpsgen::FluxReaderNuMI::get_kaon_from_flux(TLorentzVector& kmom, TLorentzVector& kpos, int& kpdg, int& pi_type, double& weight) {
  if(dk2nu->decay.ndecay <= 0 || dk2nu->decay.ndecay > 10) return false; // only select kaon decays
  //std::cerr << "DEBUG "<<__FILE__<<" "<<__LINE__<<" " <<dk2nu->decay.ndecay<<" "<<dk2nu->decay.ptype<<std::endl;
  kpdg = dk2nu->decay.ptype;
  pi_type = (dk2nu->decay.ndecay < 5 ? 0 : (kpdg < 0 ? -1 : 1));
  kpos.SetVect(rot * TVector3(dk2nu->decay.vx,dk2nu->decay.vy,dk2nu->decay.vz) + beampos);
  kpos.SetT(dk2nu->ancestor.back().startt + beamtime);
  const double kmass = TLorentzVector(
      std::abs(dk2nu->decay.pppz) > 0 ? dk2nu->decay.pppz*dk2nu->decay.ppdxdz : 0.,
      std::abs(dk2nu->decay.pppz) > 0 ? dk2nu->decay.pppz*dk2nu->decay.ppdydz : 0.,
      dk2nu->decay.pppz,
      dk2nu->decay.ppenergy
      ).M();
  kmom.SetVectM(rot * TVector3(dk2nu->decay.pdpx,dk2nu->decay.pdpy,dk2nu->decay.pdpz),kmass);
  if(kmom.M() < 0.4) {
    dk2nu->Print();
    kmom.Print();
    throw cet::exception("LogicError")<<"Found a massless kaon! mass="<<kmass<<std::endl;
  }

    const double tot_K2nu_branching_ratio = (dk2nu->decay.ndecay < 5)
    ?
    K2nu_branching_ratio(1) + K2nu_branching_ratio(2) + K2nu_branching_ratio(3) + K2nu_branching_ratio(4)
    :
    K2nu_branching_ratio(5) + K2nu_branching_ratio(6) + K2nu_branching_ratio(7)
    ;

  weight = dk2nu->decay.nimpwt  / tot_K2nu_branching_ratio;//K2nu_branching_ratio(dk2nu->decay.ndecay);

  return true; }





// branching ratios taken from the GEANT4 version used by g4numi in MCC9,
// located in /cvmfs/minerva.opensciencegrid.org/minerva/beamsim/x86_64/geant4
double hpsgen::FluxReaderNuMI::K2nu_branching_ratio(const int ndecay) {
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
  if (ndecay == 1) {
    return 0.2020; 
  }
  else if (ndecay == 2) {
    return 0.2020;
  }
  else if (ndecay == 3) {
    return 0.1348;
  }
  else if (ndecay == 4) {
    return 0.1348;
  }
  else if (ndecay == 5) {
    return 0.6339;
  }
  else if (ndecay == 6) {
    return 0.0493;
  }
  else if (ndecay == 7) {
    return 0.0330;
  }
  else if (ndecay == 8) {
    return 0.6339;
  }
  else if (ndecay == 9) {
    return 0.0493;
  }
  else if (ndecay == 10) {
    return 0.0330;
  }
  else if (ndecay == 11) {
    return 1.0;
  }
  else if (ndecay == 12) {
    return 1.0;
  }
  else if (ndecay == 13) {
    return 1.0;
  }
  else if (ndecay == 14) {
    return 1.0;
  }
  return 1e100; // shouldn't get here, but if so, weight will be highly suppressed
}

void hpsgen::FluxReaderNuMI::get_MCFlux(simb::MCFlux& flux) {
  get_current_entry();
  flux.Reset();
  flux.fFluxType = simb::kDk2Nu;

  flux.frun      = dk2nu->job;
  flux.fevtno    = dk2nu->potnum;

  // ignore vector<bsim::NuRay> (see nuchoice above)

  // bsim::Decay object
  flux.fnorig    = dk2nu->decay.norig;
  flux.fndecay   = dk2nu->decay.ndecay;
  flux.fntype    = dk2nu->decay.ntype;
  flux.fppmedium = dk2nu->decay.ppmedium;
  flux.fptype    = dk2nu->decay.ptype;

  flux.fvx       = dk2nu->decay.vx;
  flux.fvy       = dk2nu->decay.vy;
  flux.fvz       = dk2nu->decay.vz;
  flux.fpdpx     = dk2nu->decay.pdpx;
  flux.fpdpy     = dk2nu->decay.pdpy;
  flux.fpdpz     = dk2nu->decay.pdpz;

  flux.fppdxdz   = dk2nu->decay.ppdxdz;
  flux.fppdydz   = dk2nu->decay.ppdydz;
  flux.fpppz     = dk2nu->decay.pppz;
  flux.fppenergy = dk2nu->decay.ppenergy;

  flux.fmuparpx  = dk2nu->decay.muparpx;
  flux.fmuparpy  = dk2nu->decay.muparpy;
  flux.fmuparpz  = dk2nu->decay.muparpz;
  flux.fmupare   = dk2nu->decay.mupare;

  flux.fnecm     = dk2nu->decay.necm;
  flux.fnimpwt   = dk2nu->decay.nimpwt;

  // no place for:  vector<bsim::Ancestor>

  // production vertex of nu parent
  flux.fppvx      = dk2nu->ppvx;
  flux.fppvy      = dk2nu->ppvy;
  flux.fppvz      = dk2nu->ppvz;

  // bsim::TgtExit object
  flux.ftvx      = dk2nu->tgtexit.tvx;
  flux.ftvy      = dk2nu->tgtexit.tvy;
  flux.ftvz      = dk2nu->tgtexit.tvz;
  flux.ftpx      = dk2nu->tgtexit.tpx;
  flux.ftpy      = dk2nu->tgtexit.tpy;
  flux.ftpz      = dk2nu->tgtexit.tpz;
  flux.ftptype   = dk2nu->tgtexit.tptype;   // converted to PDG
  flux.ftgen     = dk2nu->tgtexit.tgen;

  // ignore vector<bsim::Traj>

}
