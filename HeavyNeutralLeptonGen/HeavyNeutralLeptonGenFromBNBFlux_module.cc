////////////////////////////////////////////////////////////////////////
// Class:       HeavyNeutralLeptonGenFromBNBFlux
// Plugin Type: producer (art v3_01_02)
// File:        HeavyNeutralLeptonGenFromBNBFlux_module.cc
//
// Generated at Fri Mar 20 11:25:20 2020 by Pawel Guzowski using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <memory>

#include "TTree.h"

#include "GenKinematics.h"
#include "FluxReaderBNB.h"

namespace hpsgen {
  class HeavyNeutralLeptonGenFromBNBFlux;
}


class hpsgen::HeavyNeutralLeptonGenFromBNBFlux : public art::EDProducer {
public:
  explicit HeavyNeutralLeptonGenFromBNBFlux(fhicl::ParameterSet const& p);
  ~HeavyNeutralLeptonGenFromBNBFlux();
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HeavyNeutralLeptonGenFromBNBFlux(HeavyNeutralLeptonGenFromBNBFlux const&) = delete;
  HeavyNeutralLeptonGenFromBNBFlux(HeavyNeutralLeptonGenFromBNBFlux&&) = delete;
  HeavyNeutralLeptonGenFromBNBFlux& operator=(HeavyNeutralLeptonGenFromBNBFlux const&) = delete;
  HeavyNeutralLeptonGenFromBNBFlux& operator=(HeavyNeutralLeptonGenFromBNBFlux&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run& r) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

private:

  // Declare member data here.
  CLHEP::HepRandomEngine& fRNG;
  GenKinematics *fKinHelper;
  FluxReaderBNB *fFluxHelper;

  const std::string fScalarParams;
  const std::vector<double> fScalarMass;
  const double fModelTheta;
  const double fMaxWeight;

  const int fProdLepType;
  const int fDecayLepType;


  double fPrevTotPOT;
  double fPrevTotGoodPOT;
  
  const std::vector<int> fSelectKaonPDGs;
  const std::string fSelectKaons;
  const double fCutKaonMom;

  const double fGlobalTimeOffset;
  const double fBeamWindowDuration;

  const std::vector<double> fBeamWindowCut;
  const std::vector<int> fFinalStateCut;

  const bool fMajoranaDecay;

  bool fTreeOnlyMode;

  TTree *fEventTree;
  double fEventTree_kaon_mom_x;
  double fEventTree_kaon_mom_y;
  double fEventTree_kaon_mom_z;
  double fEventTree_kaon_energy;
  double fEventTree_kaon_decay_x;
  double fEventTree_kaon_decay_y;
  double fEventTree_kaon_decay_z;
  double fEventTree_kaon_decay_t;
  double fEventTree_scalar_mom_x;
  double fEventTree_scalar_mom_y;
  double fEventTree_scalar_mom_z;
  double fEventTree_scalar_energy;
  double fEventTree_scalar_decay_x;
  double fEventTree_scalar_decay_y;
  double fEventTree_scalar_decay_z;
  double fEventTree_scalar_decay_t;
  double fEventTree_daughter1_mom_x;
  double fEventTree_daughter1_mom_y;
  double fEventTree_daughter1_mom_z;
  double fEventTree_daughter1_energy;
  double fEventTree_daughter2_mom_x;
  double fEventTree_daughter2_mom_y;
  double fEventTree_daughter2_mom_z;
  double fEventTree_daughter2_energy;
  int    fEventTree_daughter1_pdg;
  int    fEventTree_daughter2_pdg;
//testcomment
  double fEventTree_weight;
  double fEventTree_flux_weight;
  double fEventTree_decay_weight;
  double fEventTree_branching_ratio_weight;
  int    fEventTree_kaon_pdg;
  bool   fEventTree_selected;

  //added by Magnus
  double fEventTree_time_shift;
  double fEventTree_unshifted_time;

  TTree *fSubRunTree;
  double fSubRunTree_totpot;
  ULong64_t   fSubRunTree_n_kaons_read;
  ULong64_t   fSubRunTree_n_scalars_gen;
  int    fSubRunTree_n_scalar_decays_in_detector;

};


hpsgen::HeavyNeutralLeptonGenFromBNBFlux::HeavyNeutralLeptonGenFromBNBFlux(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fRNG(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "hpsgen", p, "RNGSeed")),
  fKinHelper(new GenKinematics(p)), fFluxHelper(new FluxReaderBNB(p, fRNG)),
  fScalarParams(p.get<std::string>("scalar_params","fixed")),
  fScalarMass(
      [&p](){
        try {
        return p.get<std::vector<double>>("HNL_mass",{0.125});
        }
        catch(...){
        return std::vector<double>{p.get<double>("HNL_mass",0.125)};
        }
      }()
  ),
  fModelTheta(p.get<double>("model_theta",1e-5)),
  fMaxWeight(p.get<double>("max_weight",0.)),
  fProdLepType(p.get<int>("prod_lep_type",2)),
  fDecayLepType(p.get<int>("decay_lep_type",2)),
  fSelectKaonPDGs(p.get<std::vector<int>>("select_kaon_pdgs",{})),
  fSelectKaons(p.get<std::string>("select_kaon_decay_type","")),
  fCutKaonMom(p.get<double>("cut_kaon_mom",0.)),
  // note these are not interaction times or trigger time, but neutrino times, in line with e.g GENIE:
  fGlobalTimeOffset(p.get<double>("global_time_offset",3125.)), // time of bnb window start
  fBeamWindowDuration(p.get<double>("beam_window_duration",1600)), // bnb window duration
  // this is the interaction time:
  fBeamWindowCut(p.get<std::vector<double>>("beam_window_cut",{})), // cut on final interaction time, eg for HNL window [4800,5400]
  fFinalStateCut(p.get<std::vector<int>>("final_state_cut",{})), // cut on final state, 
  fMajoranaDecay(p.get<bool>("majorana_decay",true)), //both sign decays or only right signed
  fTreeOnlyMode(p.get<bool>("tree_only_mode",false)) // produce only the TFile tree, not the artroot output to process faster
{
  produces< sumdata::RunData, art::InRun >();
  produces< std::map<std::string,double> , art::InRun >("generatorConfig");
  produces< sumdata::POTSummary, art::InSubRun >();
  if(!fTreeOnlyMode) {
    produces< std::vector<simb::MCTruth> >();
    produces< std::vector<simb::MCFlux>  >();
    produces< art::Assns<simb::MCTruth, simb::MCFlux> >();
  }
  
  if(fScalarMass.empty()) {
    throw cet::exception("Configuration") << "Need to supply a model scalar mass";
  }
  if(fScalarParams != "fixed" && fScalarParams != "random") {
    throw cet::exception("Configuration") << "scalar_params should be 'fixed' or 'random' only";
  }
  if(fScalarParams == "random" && fScalarMass.size() != 2 && fMaxWeight != 0) {
    throw cet::exception("Configuration") << "Need to supply a model scalar mass range [low,high] and max_weight should be 0";
  }
  if(!fSelectKaons.empty() && fSelectKaons != "kdar" && fSelectKaons != "kdif") {
    throw cet::exception("Configuration") << "select_kaon_decay_type should be 'kdar' or 'kdif', or not defined";
  }
  if(!fBeamWindowCut.empty() && fBeamWindowCut.size() != 2) {
     throw cet::exception("Configuration") << "beam_window_cut should be [ t1, t2 ]";
  }
  if(!fFinalStateCut.empty()) {
    const int first = fFinalStateCut.front();
    for(auto const& v : fFinalStateCut) {
      if(v * first <= 0) { 
        throw cet::exception("Configuration") << " all final state numbers should have the same sign";
      }
    }
  }
  if(!fSelectKaonPDGs.empty()) {
    for(auto const& v : fSelectKaonPDGs) {
      if(v != 130 && std::abs(v) != 321) {
        throw cet::exception("Configuration") << " select_kaon_pdgs should be a choice of 130, 321, -321";
      }
    }
  }

  art::ServiceHandle<art::TFileService> tfs; 
  fEventTree = tfs->make<TTree>("event_tree","tree of events");
  fEventTree->Branch("kaon_mom_x",&fEventTree_kaon_mom_x);
  fEventTree->Branch("kaon_mom_y",&fEventTree_kaon_mom_y);
  fEventTree->Branch("kaon_mom_z",&fEventTree_kaon_mom_z);
  fEventTree->Branch("kaon_energy",&fEventTree_kaon_energy);
  fEventTree->Branch("kaon_decay_x",&fEventTree_kaon_decay_x);
  fEventTree->Branch("kaon_decay_y",&fEventTree_kaon_decay_y);
  fEventTree->Branch("kaon_decay_z",&fEventTree_kaon_decay_z);
  fEventTree->Branch("kaon_decay_t",&fEventTree_kaon_decay_t);
  fEventTree->Branch("scalar_mom_x",&fEventTree_scalar_mom_x);
  fEventTree->Branch("scalar_mom_y",&fEventTree_scalar_mom_y);
  fEventTree->Branch("scalar_mom_z",&fEventTree_scalar_mom_z);
  fEventTree->Branch("scalar_energy",&fEventTree_scalar_energy);
  fEventTree->Branch("scalar_decay_x",&fEventTree_scalar_decay_x);
  fEventTree->Branch("scalar_decay_y",&fEventTree_scalar_decay_y);
  fEventTree->Branch("scalar_decay_z",&fEventTree_scalar_decay_z);
  fEventTree->Branch("scalar_decay_t",&fEventTree_scalar_decay_t);
  fEventTree->Branch("daughter1_mom_x",&fEventTree_daughter1_mom_x);
  fEventTree->Branch("daughter1_mom_y",&fEventTree_daughter1_mom_y);
  fEventTree->Branch("daughter1_mom_z",&fEventTree_daughter1_mom_z);
  fEventTree->Branch("daughter1_energy",&fEventTree_daughter1_energy);
  fEventTree->Branch("daughter2_mom_x",&fEventTree_daughter2_mom_x);
  fEventTree->Branch("daughter2_mom_y",&fEventTree_daughter2_mom_y);
  fEventTree->Branch("daughter2_mom_z",&fEventTree_daughter2_mom_z);
  fEventTree->Branch("daughter2_energy",&fEventTree_daughter2_energy);
  fEventTree->Branch("weight",&fEventTree_weight);
  fEventTree->Branch("flux_weight",&fEventTree_flux_weight);
  fEventTree->Branch("decay_weight",&fEventTree_decay_weight);
  fEventTree->Branch("branching_ratio_weight",&fEventTree_branching_ratio_weight);
  fEventTree->Branch("selected",&fEventTree_selected);
  
  fEventTree->Branch("kaon_pdg",&fEventTree_kaon_pdg);
  fEventTree->Branch("daughter1_pdg",&fEventTree_daughter1_pdg);
  fEventTree->Branch("daughter2_pdg",&fEventTree_daughter2_pdg);

  fSubRunTree = tfs->make<TTree>("subrun_tree","");
  fSubRunTree->Branch("tot_pot",&fSubRunTree_totpot);
  fSubRunTree->Branch("n_kaons_read",&fSubRunTree_n_kaons_read);
  fSubRunTree->Branch("n_scalars_gen",&fSubRunTree_n_scalars_gen);
  fSubRunTree->Branch("n_scalar_decays_in_detector",&fSubRunTree_n_scalar_decays_in_detector);
  fSubRunTree_n_kaons_read = 0;
  fSubRunTree_n_scalars_gen = 0;
  fSubRunTree_n_scalar_decays_in_detector = 0;

  //added by Magnus
  fEventTree->Branch("time_shift",&fEventTree_time_shift);
  fEventTree->Branch("unshifted_time",&fEventTree_unshifted_time);
}

hpsgen::HeavyNeutralLeptonGenFromBNBFlux::~HeavyNeutralLeptonGenFromBNBFlux()
{
  delete fKinHelper;
  delete fFluxHelper;
}

void hpsgen::HeavyNeutralLeptonGenFromBNBFlux::produce(art::Event& e)
{
  double HNL_mass = 0.;
  if(fScalarParams == "fixed") {
    HNL_mass = fScalarMass.front();
    }
  else if(fScalarParams == "random") {
    HNL_mass = CLHEP::RandFlat::shoot(&fRNG,fScalarMass.front(),fScalarMass.back());
    std::cout << "Choosing scalar mass "<<HNL_mass<< " from range ["<<fScalarMass.front()<<","<<fScalarMass.back() <<"] "<<std::endl;
  }
  TLorentzVector kaon_4mom, kaon_pos;
  int pion_type;
  int kaon_pdg;
  while(true) {
    const double flux_weight = fFluxHelper->get_kaon(kaon_4mom,kaon_pos,kaon_pdg,pion_type);

    fSubRunTree_n_kaons_read++;

    if(!fSelectKaonPDGs.empty()) {
      if(std::find(fSelectKaonPDGs.begin(), fSelectKaonPDGs.end(), kaon_pdg) == fSelectKaonPDGs.end()) {
        continue;
      }
    }
    if(fSelectKaons == "kdif") {
      if(kaon_4mom.Vect().Mag() < fCutKaonMom) continue;
    }
    if(fSelectKaons == "kdar") {
      if(kaon_4mom.Vect().Mag() > fCutKaonMom) continue;
    }
    std::multimap<int,TLorentzVector> res;
    fSubRunTree_n_scalars_gen++;
    const bool passes = (fScalarParams == "random") ?
      fKinHelper->generate_uniform(kaon_pos,kaon_4mom,kaon_pdg, HNL_mass,  fRNG, res) :
      fKinHelper->generate(kaon_pos,kaon_4mom, kaon_pdg, HNL_mass, fProdLepType, fDecayLepType,flux_weight, fMaxWeight, fRNG, res);
    if(passes) {
      
      const TLorentzVector& dk_pos = res.find(0)->second;
      const TLorentzVector& scalar_mom = res.find(54)->second;
      auto d1ptr = [&res]() {
      for(auto i = res.begin(); i != res.end(); ++i) {
          auto const& v = *i;
          if((abs(v.first) == 13) | (abs(v.first) == 11)) return i; //want the first entry to be lepton. if i do it old way get the first neg pdg
          // if(v.first != 54 && v.first != 99 && v.first != 0) return i;
        }
        return res.end();
      }();
      auto const& d1 = *d1ptr;
      auto const& d2ptr = [&res,&d1ptr]() {
        for(auto i = res.begin(); i != res.end(); ++i) {
          auto const& v = *i;
          if(v.first != 54 && v.first != 99 && v.first != 0 && i != d1ptr) return i;
        };
        return res.end();
      }();
      auto const& d2 = *d2ptr;

      bool selected = true;

      if(!fFinalStateCut.empty()) {
        if(fFinalStateCut.front() > 0) {
          if(std::find(fFinalStateCut.begin(), fFinalStateCut.end(), std::abs(d1.first)) == fFinalStateCut.end()) {
            selected = false;
          }
        }
        else if(fFinalStateCut.front() < 0) {
          if(std::find(fFinalStateCut.begin(), fFinalStateCut.end(), -std::abs(d1.first)) != fFinalStateCut.end()) {
            selected = false;
          }
        }
      }

      //Added by Magnus
      auto const seed = 123;
      auto urbg = std::mt19937 {seed};  
      double const mu = 0.0; 
      double const sigma = 1.308; //ns
      auto norm = std::normal_distribution<double>{mu,sigma};
      auto gauss = norm(urbg);
      double bunch_spacing = 18.936; //ns

      //std::cout << "fRNG: "<< fRNG << std::endl;

      double time_shift_crude = fGlobalTimeOffset+CLHEP::RandFlat::shoot(&fRNG,fBeamWindowDuration);
      double time_shift = ((int) (time_shift_crude / bunch_spacing)) * bunch_spacing;
      time_shift += gauss;
      
      TLorentzVector shift_to_detector_time(0.,0.,0.,time_shift);

      const double dk_t = (dk_pos+shift_to_detector_time).T();
      if(!fBeamWindowCut.empty() && ( dk_t < fBeamWindowCut.front() || dk_t > fBeamWindowCut.back() )) {
        selected=false;
      }

        if(!fMajoranaDecay) {
          if(d1.first/kaon_pdg < 0) { //if daughter and kaon pdg codes are oppesite signs it is a non allowed decay
            selected = false;
          }
        }

      if(fMaxWeight < 0.) {
        auto const& r99 = res.find(99);
        if(r99 == res.end()) {
          throw cet::exception("LogicError") << "there should be a weight lorentz vector" <<  std::endl;
        }
        fEventTree_weight = r99->second.T();
        fEventTree_decay_weight = r99->second.X();
        fEventTree_branching_ratio_weight = r99->second.Y();
        fEventTree_flux_weight = r99->second.Z();
        //const double qq = (fScalarParams == "random" ? model_theta : 0.);
      }

      
      fEventTree_kaon_mom_x = kaon_4mom.X();
      fEventTree_kaon_mom_y = kaon_4mom.Y();
      fEventTree_kaon_mom_z = kaon_4mom.Z();
      fEventTree_kaon_energy = kaon_4mom.T();
      fEventTree_kaon_decay_x = kaon_pos.X();
      fEventTree_kaon_decay_y = kaon_pos.Y();
      fEventTree_kaon_decay_z = kaon_pos.Z();
      fEventTree_kaon_decay_t = (kaon_pos+shift_to_detector_time).T();
      fEventTree_scalar_mom_x = scalar_mom.X();
      fEventTree_scalar_mom_y = scalar_mom.Y();
      fEventTree_scalar_mom_z = scalar_mom.Z();
      fEventTree_scalar_energy = scalar_mom.E();
      fEventTree_scalar_decay_x = dk_pos.X();
      fEventTree_scalar_decay_y = dk_pos.Y();
      fEventTree_scalar_decay_z = dk_pos.Z();
      fEventTree_scalar_decay_t = (dk_pos+shift_to_detector_time).T();
      fEventTree_daughter1_mom_x = d1.second.X();
      fEventTree_daughter1_mom_y = d1.second.Y();
      fEventTree_daughter1_mom_z = d1.second.Z();
      fEventTree_daughter1_energy = d1.second.E();
      fEventTree_daughter2_mom_x = d2.second.X();
      fEventTree_daughter2_mom_y = d2.second.Y();
      fEventTree_daughter2_mom_z = d2.second.Z();
      fEventTree_daughter2_energy = d2.second.E();
      fEventTree_kaon_pdg = kaon_pdg;
      fEventTree_daughter1_pdg = d1.first;
      fEventTree_daughter2_pdg = d2.first;
      fEventTree_selected = selected;

      //Added by Magnus
      fEventTree_time_shift = time_shift;
      fEventTree_unshifted_time = dk_pos.T();

      fEventTree->Fill();

      fSubRunTree_n_scalar_decays_in_detector++;

      if(!selected) continue;

      if(!fTreeOnlyMode) {

        std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>(1));
        simb::MCTruth& truth = truthcol->back();

        simb::MCParticle kaon(1,kaon_pdg,"beamline",-1,kaon_4mom.M(),0);
        kaon.AddTrajectoryPoint(kaon_pos+shift_to_detector_time,kaon_4mom);

        simb::MCParticle scalar(2,54,"decay",1,HNL_mass,2);
        scalar.AddTrajectoryPoint(kaon_pos,scalar_mom);
        scalar.AddTrajectoryPoint(dk_pos+shift_to_detector_time,scalar_mom);

        simb::MCParticle dgt1(3,d1.first,"decay",2,d1.second.M(),1);
        dgt1.AddTrajectoryPoint(dk_pos+shift_to_detector_time,d1.second);

        simb::MCParticle dgt2(4,d2.first,"decay",2,d2.second.M(),1);
        dgt2.AddTrajectoryPoint(dk_pos+shift_to_detector_time,d2.second);

        truth.Add(kaon);
        truth.Add(scalar);
        truth.Add(dgt1);
        truth.Add(dgt2);

        if(fMaxWeight < 0.) {
          truth.SetNeutrino(0,0,0,0,0,0,fEventTree_weight,fEventTree_decay_weight,fEventTree_branching_ratio_weight,fEventTree_flux_weight);
        }
        truth.SetOrigin(simb::kUnknown);

        std::unique_ptr< std::vector<simb::MCFlux> > mcfluxcol(new std::vector<simb::MCFlux>(1));
        fFluxHelper->get_MCFlux(mcfluxcol->back());

        std::unique_ptr< art::Assns<simb::MCTruth, simb::MCFlux> > tfassn(new art::Assns<simb::MCTruth, simb::MCFlux>);
        util::CreateAssn(*this, e, *truthcol, *mcfluxcol, *tfassn, mcfluxcol->size()-1, mcfluxcol->size());

        e.put(std::move(truthcol));
        e.put(std::move(mcfluxcol));
        e.put(std::move(tfassn));
      }

      return;
    }
  }
}

void hpsgen::HeavyNeutralLeptonGenFromBNBFlux::beginJob()
{
  fPrevTotPOT = 0.;
  fPrevTotGoodPOT = 0.;
}

void hpsgen::HeavyNeutralLeptonGenFromBNBFlux::beginRun(art::Run& r)
{
  art::ServiceHandle<geo::Geometry const> geo;
  r.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
  fKinHelper->update_geometry(geo);
  std::unique_ptr<std::map< std::string, double >> mc_config(new std::map< std::string, double >);
  (*mc_config)["mass_elec"] = fKinHelper->get_constants().mass_elec();
  (*mc_config)["mass_muon"] = fKinHelper->get_constants().mass_muon();
  (*mc_config)["mass_pion_pm"] = fKinHelper->get_constants().mass_pion_pm();
  (*mc_config)["mass_pion_0"] = fKinHelper->get_constants().mass_pion_0();
  (*mc_config)["speed_light"] = fKinHelper->get_constants().speed_light();
  (*mc_config)["higgs_vev"] = fKinHelper->get_constants().higgs_vev();
  (*mc_config)["hbar"] = fKinHelper->get_constants().hbar();
  (*mc_config)["ckm_ts"] = fKinHelper->get_constants().ckm_ts();
  (*mc_config)["ckm_td"] = fKinHelper->get_constants().ckm_td();
  (*mc_config)["lifetime_kaon_0"] = fKinHelper->get_constants().lifetime_kaon_0();
  (*mc_config)["lifetime_kaon_pm"] = fKinHelper->get_constants().lifetime_kaon_pm();
  (*mc_config)["mass_top"] = fKinHelper->get_constants().mass_top();
  if(fScalarParams == "fixed") {
    (*mc_config)["model_theta"] = fModelTheta;
  }
  r.put(std::move(mc_config),"generatorConfig");
}

void hpsgen::HeavyNeutralLeptonGenFromBNBFlux::beginSubRun(art::SubRun& sr)
{
  fPrevTotPOT = fFluxHelper->POTSeen(fMaxWeight);
  fPrevTotGoodPOT = fFluxHelper->POTSeen(fMaxWeight);
}

void hpsgen::HeavyNeutralLeptonGenFromBNBFlux::endSubRun(art::SubRun& sr)
{
  auto p = std::make_unique<sumdata::POTSummary>();
  p->totpot = fFluxHelper->POTSeen(fMaxWeight) - fPrevTotPOT;
  p->totgoodpot = fFluxHelper->POTSeen(fMaxWeight) - fPrevTotGoodPOT;
  fSubRunTree_totpot = p->totpot;
  fSubRunTree->Fill();
  sr.put(std::move(p));
}

DEFINE_ART_MODULE(hpsgen::HeavyNeutralLeptonGenFromBNBFlux)
