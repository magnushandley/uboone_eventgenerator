#include "FluxReader.h"

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"

#include <numeric>
#include <glob.h>

#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "IFDH_service.h"

hpsgen::FluxReader::FluxReader(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& fRNG,
    const char* tree_name, const char* pot_branch_name, const char* pot_tree) : 
  fFluxLocation(p.get<std::string>("flux_location","")),
  fFluxAccessSchema(p.get<std::string>("flux_access_schema","direct")),
  fUseCache(p.get<bool>("use_flux_cache",true)),
  fFluxTree(new TChain(tree_name)),
  fCurrEntry(-1), fPOTperFluxFile{}, fTotalPOTinFluxFiles(0.), fEntriesPerFluxFile{}, fFullChainsSeen(0),
  fMaxFluxFileMB(p.get<unsigned int>("ifdh_max_size",1000)),
  fCache{}, fSelectedEntries{}
{
  if(fFluxLocation.empty()) {
    throw cet::exception("Configuration") << "must supply a flux file location";
  }
  if("direct" != fFluxAccessSchema && "ifdh" != fFluxAccessSchema && "xrootd" != fFluxAccessSchema) {
    throw cet::exception("Configuration") << "flux_access_schema should be 'direct', 'ifdh' or 'xrootd'";
  }
  if(fFluxAccessSchema == "xrootd" && fFluxLocation.compare(0,6,"/pnfs/") != 0) {
    throw cet::exception("Configuration") << "Can only use xrootd schema with files on /pnfs/";
  }
  
  if(fFluxAccessSchema == "direct" || fFluxAccessSchema == "xrootd") {
    glob_t glob_result;
    std::map<double,std::string> filenames;
    int ret = glob(fFluxLocation.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(ret != 0) {
      globfree(&glob_result);
      throw cet::exception("Configuration") << " unable to access '"<<fFluxLocation<<"'. Maybe you should use IFDH access?"<<std::endl;
      /*
      if(fFluxAccessSchema == "direct" || fFluxLocation.compare(0,6,"/pnfs/") != 0) {
        throw cet::exception("Configuration")<<"flux file globbing returned an error "<<ret<<std::endl;
      }
      const size_t slash_loc = fFluxLocation.find_last_of('/');
      const std::string spath = fFluxLocation.substr(0,slash_loc);
      const std::string fpath = fFluxLocation.substr(slash_loc+1,std::string::npos);
      art::ServiceHandle<IFDH> ifdh;
      auto const& results = ifdh->findMatchingFiles(spath,fpath);
      for(auto const& r : results) {
        std::string fn = r.first;
        //std::cout << "ifdh found filename: " << fn << std::endl;
        if(fn.compare(0,6,"/pnfs/") == 0) {
          fn.replace(0,6,"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/");
        }
        else {
          throw cet::exception("Configuration") << "Unable to find xrootd-able path for input flux file "<<fn<<std::endl;
        }
        filenames[CLHEP::RandFlat::shoot(&fRNG)] = fn;
      }
      */
    }
    else {
      for(size_t i = 0; i < glob_result.gl_pathc; ++i) {
        std::string fn = glob_result.gl_pathv[i];
        if(fFluxAccessSchema == "xrootd") {
          if(fn.compare(0,6,"/pnfs/") == 0) {
            fn.replace(0,6,"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/");
          }
          else {
            globfree(&glob_result);
            throw cet::exception("Configuration") << "Unable to find xrootd-able path for input flux file "<<fn<<std::endl;
          }
        }
        // make a random ordering of the filenames
        filenames[CLHEP::RandFlat::shoot(&fRNG)] = fn;
      }
      globfree(&glob_result);
    }
    for(auto const& f : filenames) {
      fFluxTree->Add(f.second.c_str());
    }
  }
  else if(fFluxAccessSchema == "ifdh") {
    //throw cet::exception("Configuration") << "ifdh not implemented yet." << std::endl;
    const size_t slash_loc = fFluxLocation.find_last_of('/');
    const std::string spath = fFluxLocation.substr(0,slash_loc);
    const std::string fpath = fFluxLocation.substr(slash_loc+1,std::string::npos);
    art::ServiceHandle<IFDH> ifdh;
    auto const& results = ifdh->findMatchingFiles(spath,fpath);
    std::map<double, std::pair<std::string,long>> filenames_and_sizes;
    for(auto const& r : results) {
      filenames_and_sizes[CLHEP::RandFlat::shoot(&fRNG)] = r;
    }
    long long int tot_size = 0;
    const long long int max_size = (long long int)fMaxFluxFileMB * 1024ll * 1024ll;
    std::vector<std::pair<std::string,long>> selected_list;
    for(auto const& fs: filenames_and_sizes) {
      selected_list.push_back(fs.second);
      const long size = fs.second.second;
      tot_size += size;
      if(tot_size > max_size) break;
    }
    auto const& locals = ifdh->fetchSharedFiles(selected_list);
    for(auto const& l : locals) {
      std::cout << "fetched "<<l.first<<std::endl;
      fFluxTree->Add(l.first.c_str());
    }
  }
  std::cout << "hpsgen::FluxReader::FluxReader : please wait, initializing flux (can take some time) ... "<<std::endl;
  const unsigned long nentries = fFluxTree->GetEntries();
  if(nentries < 1) {
    throw cet::exception("Configuration") << "There are no flux entries."<<std::endl;
  }
  count_pot(pot_branch_name, pot_tree);
}

hpsgen::FluxReader::~FluxReader() {
  delete fFluxTree;
}

double hpsgen::FluxReader::POTSeen(const double weight /* = 1. */) const {
  if(fCurrEntry < 0) return 0.;
  double partial_pot = 0.;
  const size_t curr_file = fFluxTree->GetTreeNumber();
  for(size_t i = 0; i < curr_file; ++i) partial_pot += fPOTperFluxFile[i];
  partial_pot += (double)(1 + fCurrEntry - fFluxTree->GetTreeOffset()[curr_file]) / fEntriesPerFluxFile[curr_file] * fPOTperFluxFile[curr_file];
  if(weight > 0) return (fFullChainsSeen * fTotalPOTinFluxFiles + partial_pot) / weight;
  return fFullChainsSeen * fTotalPOTinFluxFiles + partial_pot;
}

void hpsgen::FluxReader::set_branch(void* obj, const char* bname) {
  fFluxTree->SetBranchAddress(bname, obj);
}

void hpsgen::FluxReader::get_next_entry() {
  if(fCache.empty()) {
    if(fUseCache) {
      fCache.resize(fFluxTree->GetEntries());
      fSelectedEntries.reserve(fFluxTree->GetEntries());
    }
    else fCache.resize(1);
  }
  if(fUseCache && fFullChainsSeen > 0) {
    fCurrSelectedEntry++;
    if(fCurrSelectedEntry >= fSelectedEntries.end()) {
      fCurrSelectedEntry = fSelectedEntries.begin();
      fFullChainsSeen++;
    }
    fCurrEntry = *fCurrSelectedEntry;
    return;
  }
  fCurrEntry++;
  if(fCurrEntry >= fFluxTree->GetEntries()) {
    fCurrEntry = 0;
    fFullChainsSeen++;
    if(fUseCache) {
      fCurrSelectedEntry = fSelectedEntries.begin();
      fCurrSelectedEntry--; // this will be ++ed again on next entry
      return get_next_entry();
    }
  }
  fFluxTree->GetEntry(fCurrEntry);
  auto& cache = fUseCache ? fCache.at(fCurrEntry) : fCache.front();
  if(!get_kaon_from_flux(cache.kmom, cache.kpos, cache.kpdg, cache.pi_type, cache.weight)){
    return get_next_entry();
  }
  if(fUseCache) fSelectedEntries.push_back(fCurrEntry);
}

void hpsgen::FluxReader::get_current_entry() {
  if(fUseCache && fFullChainsSeen > 0) fFluxTree->GetEntry(fCurrEntry);
}

void hpsgen::FluxReader::count_pot(const char* bname, const char* tname) {
  const size_t n_trees = fFluxTree->GetNtrees();
  fTotalPOTinFluxFiles = 0.;
  fPOTperFluxFile.resize(n_trees);
  fEntriesPerFluxFile.resize(n_trees);
  for(size_t i = 0; i < n_trees; ++i) {
    TFile *f = TFile::Open(fFluxTree->GetListOfFiles()->At(i)->GetTitle());
    TTree *t = (TTree*)f->Get(tname);
    if(!t) {
      throw cet::exception("Configuration") << "cannot read pot counting tree in flux files";
    }
    t->Draw("1",bname,"goff");
    const double pot = std::accumulate(t->GetW(),t->GetW()+t->GetSelectedRows(), 0.);
    fPOTperFluxFile[i] = pot;
    fTotalPOTinFluxFiles += pot;
    fEntriesPerFluxFile[i] = (i < n_trees-1 ? fFluxTree->GetTreeOffset()[i+1]-fFluxTree->GetTreeOffset()[i] : fFluxTree->GetEntries()-fFluxTree->GetTreeOffset()[i] );
    delete f;
    //std::cerr <<i<< " pot: "<<pot << std::endl;
  }
  std::cout << "hpsgen::FluxReader::count_pot : read "<<n_trees<<" files, with total pot = "<<fTotalPOTinFluxFiles<<std::endl;
}

double hpsgen::FluxReader::get_kaon(TLorentzVector& kmom, TLorentzVector& kpos, int& kpdg, int& pi_type) {
  get_next_entry();
  auto const& centry = fUseCache ? fCache.at(fCurrEntry) : fCache.front();
  kmom = centry.kmom;
  kpos = centry.kpos;
  kpdg = centry.kpdg;
  pi_type = centry.pi_type;
  return centry.weight;
}
