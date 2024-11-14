////////////////////////////////////////////////////////////////////////
// Class:       HNLGenFromBNBFlux
// Plugin Type: producer (art v3_01_02)
// File:        HNLGenFromBNBFlux_module.cc
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

namespace hnlgen {
  class HNLGenFromBNBFlux;
}


class hnlgen::HNLGenFromBNBFlux : public art::EDProducer {
public:
  explicit HNLGenFromBNBFlux(fhicl::ParameterSet const& p);
  ~HNLGenFromBNBFlux();
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HNLGenFromBNBFlux(HNLGenFromBNBFlux const&) = delete;
  HNLGenFromBNBFlux(HNLGenFromBNBFlux&&) = delete;
  HNLGenFromBNBFlux& operator=(HNLGenFromBNBFlux const&) = delete;
  HNLGenFromBNBFlux& operator=(HNLGenFromBNBFlux&&) = delete;

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

  const double fModelUe4;
  const double fModelUmu4;
  const double fModelUtau4;
  const bool fIsMajorana;
  const double fMaxWeight;

  double fPrevTotPOT;
  double fPrevTotGoodPOT;
  
  const std::vector<int> fSelectKaonPDGs;
  const std::string fSelectKaons;
  const double fCutKaonMom;

  const double fGlobalTimeOffset;
  const double fBeamWindowDuration;

  const std::vector<double> fBeamWindowCut;
  const std::vector<int> fFinalStateCut;

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
  double fEventTree_hnl_mom_x;
  double fEventTree_hnl_mom_y;
  double fEventTree_hnl_mom_z;
  double fEventTree_hnl_energy;
  double fEventTree_hnl_decay_x;
  double fEventTree_hnl_decay_y;
  double fEventTree_hnl_decay_z;
  double fEventTree_hnl_decay_t;
  double fEventTree_daughter1_mom_x;
  double fEventTree_daughter1_mom_y;
  double fEventTree_daughter1_mom_z;
  double fEventTree_daughter1_energy;
  double fEventTree_daughter2_mom_x;
  double fEventTree_daughter2_mom_y;
  double fEventTree_daughter2_mom_z;
  double fEventTree_daughter2_energy;
  double fEventTree_weight;
  double fEventTree_flux_weight;
  double fEventTree_decay_weight;
  double fEventTree_branching_ratio_weight;
  int    fEventTree_daughter1_pdg;
  int    fEventTree_daughter2_pdg;
  int    fEventTree_kaon_pdg;
  bool   fEventTree_selected;

  TTree *fSubRunTree;
  double fSubRunTree_totpot;
  ULong64_t   fSubRunTree_n_kaons_read;
  ULong64_t   fSubRunTree_n_hnl_gen;
  int    fSubRunTree_n_hnl_decays_in_detector;

};


hnlgen::HNLGenFromBNBFlux::HNLGenFromBNBFlux(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fRNG(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "hnlgen", p, "RNGSeed")),
  fKinHelper(new GenKinematics(p,fRNG)), fFluxHelper(new FluxReaderBNB(p, fRNG)),
  fModelUe4(p.get<double>("model_U_e_4_angle")),
  fModelUmu4(p.get<double>("model_U_mu_4_angle")),
  fModelUtau4(p.get<double>("model_U_tau_4_angle")),
  fIsMajorana{[](auto const& s) {
    if(s == "dirac") return false;
    if(s == "majorana") return true;
    throw cet::exception("Configuration")
      << "HNL fermionic nature '"<<s<<"' should be 'dirac' or 'majorana'" ;
  }(p.get<std::string>("model_hnl_fermion_nature"))},
  fMaxWeight(p.get<double>("max_weight",0.)),
  fSelectKaonPDGs(p.get<std::vector<int>>("select_kaon_pdgs",{})),
  fSelectKaons(p.get<std::string>("select_kaon_decay_type","")),
  fCutKaonMom(p.get<double>("cut_kaon_mom",0.)),
  // note these are not interaction times or trigger time, but neutrino times, in line with e.g GENIE:
  fGlobalTimeOffset(p.get<double>("global_time_offset",3125.)), // time of bnb window start
  fBeamWindowDuration(p.get<double>("beam_window_duration",1600)), // bnb window duration
  // this is the interaction time:
  fBeamWindowCut(p.get<std::vector<double>>("beam_window_cut",{})), // cut on final interaction time, eg for HNL window [4800,5400]
  fFinalStateCut(p.get<std::vector<int>>("final_state_cut",{})), // cut on final state, 
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
  fEventTree->Branch("hnl_mom_x",&fEventTree_hnl_mom_x);
  fEventTree->Branch("hnl_mom_y",&fEventTree_hnl_mom_y);
  fEventTree->Branch("hnl_mom_z",&fEventTree_hnl_mom_z);
  fEventTree->Branch("hnl_energy",&fEventTree_hnl_energy);
  fEventTree->Branch("hnl_decay_x",&fEventTree_hnl_decay_x);
  fEventTree->Branch("hnl_decay_y",&fEventTree_hnl_decay_y);
  fEventTree->Branch("hnl_decay_z",&fEventTree_hnl_decay_z);
  fEventTree->Branch("hnl_decay_t",&fEventTree_hnl_decay_t);
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
  fSubRunTree->Branch("n_hnl_gen",&fSubRunTree_n_hnl_gen);
  fSubRunTree->Branch("n_hnl_decays_in_detector",&fSubRunTree_n_hnl_decays_in_detector);
  fSubRunTree_n_kaons_read = 0;
  fSubRunTree_n_hnl_gen = 0;
  fSubRunTree_n_hnl_decays_in_detector = 0;
}

hnlgen::HNLGenFromBNBFlux::~HNLGenFromBNBFlux()
{
  delete fKinHelper;
  delete fFluxHelper;
}

void hnlgen::HNLGenFromBNBFlux::produce(art::Event& e)
{
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
    fSubRunTree_n_hnl_gen++;
    if(fKinHelper->generate(kaon_pos, kaon_4mom,kaon_pdg, flux_weight, fMaxWeight, fRNG, res)) {
      
      const TLorentzVector& dk_pos = res.find(0)->second;
      auto hnl_it = std::find_if(res.begin(), res.end(), [](auto& r) {
          return r.first == 89 || r.first == 91 || r.first == -89 || r.first==-91;
          });
      if(hnl_it == res.end()) continue;
      const TLorentzVector& hnl_mom = hnl_it->second;
      auto d1ptr = [&res]() {
        for(auto i = res.begin(); i != res.end(); ++i) {
          auto const& v = *i;
          if(v.first != 89 && v.first != 91 && v.first != -89 && v.first != -91
              && v.first != 12 && v.first != -12 && v.first != 99 && v.first != 0) return i;
        }
        return res.end();
      }();
      auto const& d1 = *d1ptr;
      auto const& d2ptr = [&res,&d1ptr]() {
        for(auto i = res.begin(); i != res.end(); ++i) {
          auto const& v = *i;
          if(v.first != 89 && v.first != 91 && v.first != -89 && v.first != -91
              && ((v.first != 12 && v.first != -12) || d1ptr->first == 111) && v.first != 99 && v.first != 0 && i != d1ptr) return i;
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
      
      if(fMaxWeight <= 0.) {
        auto const& r99 = res.find(99);
        if(r99 == res.end()) {
          throw cet::exception("LogicError") << "there should be a weight lorentz vector" <<  std::endl;
        }
        fEventTree_weight = r99->second.T();
        fEventTree_decay_weight = r99->second.X();
        fEventTree_branching_ratio_weight = r99->second.Y();
        fEventTree_flux_weight = r99->second.Z();
      }
      
      //Magnus edit for ns timing

      //Details from https://github.com/NuSoftHEP/nutools/blob/v2_18_01/nutools/EventGeneratorBase/GENIE/EvtTimeFNALBeam.cxx
      
      auto const seed = 123;
      auto urbg = std::mt19937 {seed};  
      double const mu = 0.0; 
      double const sigma = 0.750; //ns
      auto norm = std::normal_distribution<double>{mu,sigma};
      auto gauss = norm(urbg);
      double bunch_spacing = 18.831; //ns

      //std::cout << "fRNG: "<< fRNG << std::endl;

      double time_shift_crude = fGlobalTimeOffset+CLHEP::RandFlat::shoot(&fRNG,fBeamWindowDuration);
      double time_shift = ((int) (time_shift_crude / bunch_spacing)) * bunch_spacing;
      time_shift += gauss;

      TLorentzVector shift_to_detector_time(0.,0.,0.,time_shift);

      fEventTree_kaon_mom_x = kaon_4mom.X();
      fEventTree_kaon_mom_y = kaon_4mom.Y();
      fEventTree_kaon_mom_z = kaon_4mom.Z();
      fEventTree_kaon_energy = kaon_4mom.T();
      fEventTree_kaon_decay_x = kaon_pos.X();
      fEventTree_kaon_decay_y = kaon_pos.Y();
      fEventTree_kaon_decay_z = kaon_pos.Z();
      fEventTree_kaon_decay_t = (kaon_pos+shift_to_detector_time).T();
      fEventTree_hnl_mom_x = hnl_mom.X();
      fEventTree_hnl_mom_y = hnl_mom.Y();
      fEventTree_hnl_mom_z = hnl_mom.Z();
      fEventTree_hnl_energy = hnl_mom.E();
      fEventTree_hnl_decay_x = dk_pos.X();
      fEventTree_hnl_decay_y = dk_pos.Y();
      fEventTree_hnl_decay_z = dk_pos.Z();
      fEventTree_hnl_decay_t = (dk_pos+shift_to_detector_time).T();
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
      fEventTree->Fill();

      fSubRunTree_n_hnl_decays_in_detector++;
      
      if(!selected) continue;

      if(!fTreeOnlyMode) {

        std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>(1));
        simb::MCTruth& truth = truthcol->back();

        simb::MCParticle kaon(1,kaon_pdg,"beamline",-1,kaon_4mom.M(),0);
        kaon.AddTrajectoryPoint(kaon_pos+shift_to_detector_time,kaon_4mom);

        simb::MCParticle hnl(2,hnl_it->first,"decay",1,hnl_mom.M(),2);
        hnl.AddTrajectoryPoint(kaon_pos+shift_to_detector_time,hnl_mom);
        hnl.AddTrajectoryPoint(dk_pos+shift_to_detector_time,hnl_mom);

        simb::MCParticle dgt1(3,d1.first,"decay",2,d1.second.M(),1);
        dgt1.AddTrajectoryPoint(dk_pos+shift_to_detector_time,d1.second);

        simb::MCParticle dgt2(4,d2.first,"decay",2,d2.second.M(),1);
        dgt2.AddTrajectoryPoint(dk_pos+shift_to_detector_time,d2.second);

        truth.Add(kaon);
        truth.Add(hnl);
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

void hnlgen::HNLGenFromBNBFlux::beginJob()
{
  fPrevTotPOT = 0.;
  fPrevTotGoodPOT = 0.;
}

void hnlgen::HNLGenFromBNBFlux::beginRun(art::Run& r)
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
    (*mc_config)["model_Ue4"] = fModelUe4;
    (*mc_config)["model_Umu4"] = fModelUmu4;
    (*mc_config)["model_Utau4"] = fModelUtau4;
    (*mc_config)["model_majorana"] = fIsMajorana ? 1.:0.;
  r.put(std::move(mc_config),"generatorConfig");
}

void hnlgen::HNLGenFromBNBFlux::beginSubRun(art::SubRun& sr)
{
  fPrevTotPOT = fFluxHelper->POTSeen(fMaxWeight);
  fPrevTotGoodPOT = fFluxHelper->POTSeen(fMaxWeight);
}

void hnlgen::HNLGenFromBNBFlux::endSubRun(art::SubRun& sr)
{
  auto p = std::make_unique<sumdata::POTSummary>();
  p->totpot = fFluxHelper->POTSeen(fMaxWeight) - fPrevTotPOT;
  p->totgoodpot = fFluxHelper->POTSeen(fMaxWeight) - fPrevTotGoodPOT;
  fSubRunTree_totpot = p->totpot;
  fSubRunTree->Fill();
  sr.put(std::move(p));
}

DEFINE_ART_MODULE(hnlgen::HNLGenFromBNBFlux)
