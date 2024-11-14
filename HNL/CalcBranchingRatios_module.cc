////////////////////////////////////////////////////////////////////////
// Class:       CalcBranchingRatios
// Plugin Type: analyzer (art v3_01_02)
// File:        CalcBranchingRatios_module.cc
//
// Generated at Tue Dec  7 16:05:09 2021 by Pawel Guzowski using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nutools/RandomUtils/NuRandomService.h"

#include "GenKinematics.h"

#ifdef save_neg_lams
#include <TFile.h>
#include <TTree.h>
#endif

namespace hnlgen {
  class CalcBranchingRatios;
}


class hnlgen::CalcBranchingRatios : public art::EDAnalyzer {
public:
  explicit CalcBranchingRatios(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CalcBranchingRatios(CalcBranchingRatios const&) = delete;
  CalcBranchingRatios(CalcBranchingRatios&&) = delete;
  CalcBranchingRatios& operator=(CalcBranchingRatios const&) = delete;
  CalcBranchingRatios& operator=(CalcBranchingRatios&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  CLHEP::HepRandomEngine& fRNG;

};


hnlgen::CalcBranchingRatios::CalcBranchingRatios(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  ,fRNG(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "hnlgen", p, "RNGSeed"))
  // More initializers here.
{
  fhicl::ParameterSet p_backup = p;
  //p_backup.put_or_replace<double>("model_U_mu_4_angle",1e-4);
  //p_backup.put_or_replace<double>("model_U_e_4_angle",0.);
  //p_backup.put_or_replace<double>("model_U_tau_4_angle",0.);
  //p_backup.put_or_replace<std::string>("model_hnl_fermion_nature","dirac");
  p_backup.put_or_replace<std::vector<std::string>>("enabled_decay_modes",{"K->e","K->mu"});
  p_backup.put_or_replace<bool>("print_ratio",true);
  
#ifdef save_neg_lams
  TFile *fout = new TFile("fout.root","recreate");
  TTree *fTree = new TTree("t","t");
  double t_mN, t_mp, t_mn, t_pfp, t_pfn, t_snum, t_snup, t_cthm, t_cthp, t_klm, t_klp, t_int;
  fTree->Branch("mN",&t_mN);
  fTree->Branch("mp",&t_mp);
  fTree->Branch("mn",&t_mn);
  fTree->Branch("pfp",&t_pfp);
  fTree->Branch("pfn",&t_pfn);
  fTree->Branch("snum",&t_snum);
  fTree->Branch("snup",&t_snup);
  fTree->Branch("cthm",&t_cthm);
  fTree->Branch("cthp",&t_cthp);
  fTree->Branch("klm",&t_klm);
  fTree->Branch("klp",&t_klp);
  fTree->Branch("int",&t_int);
#endif
  
  for(int m = 1; m < 500; ++m) {
    fhicl::ParameterSet new_p = p_backup;
    new_p.put_or_replace<double>("model_hnl_mass",m*1e-3);
    try {
#ifdef save_neg_lams
      t_mN = m*1e-3;
#endif
    new GenKinematics(new_p,fRNG
#ifdef save_neg_lams
        ,fTree,t_pfn,t_pfp,t_mp,t_mn, t_snum, t_snup, t_cthm, t_cthp, t_klm, t_klp,t_int
#endif
        );
    } catch(cet::exception e) {
      if(e.category() == "Configuration" && e.explain_self().find("no kaon decays are possible") != std::string::npos) {
        // too heavy for kaon decays
      }
      else throw e;
    }
  }
#ifdef save_neg_lams
  fout->cd();
  fTree->Write();
  delete fout;
#endif
  //throw 1;
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void hnlgen::CalcBranchingRatios::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(hnlgen::CalcBranchingRatios)
