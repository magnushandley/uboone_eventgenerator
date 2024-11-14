/**
 * @file HepMCFileGen_module.cc
 * @brief Producer generating Monte Carlo truth record in LArSoft format from a text file in HepMC format
 * @date Dec 2019
 * @author Marco Del Tutto
 * @comment adapted from TxtFileGen by Brian Rebel
 */
/**
 * @class evgen::HepMCFileGen
 *
 *  This module assumes that the input file has the hepevt format for
 *  each event to be simulated. In brief each event contains at least two
 *  lines. The first line contains two entries, the event number (which is
 *  ignored in ART/LArSoft) and the number of particles in the event. Each
 *  following line containes 15 entries to describe each particle. The entries
 *  are:
 *
 *  1.  status code (should be set to 1 for any particle to be tracked, others
 *      won't be tracked)
 *  2.  the pdg code for the particle
 *  3.  the entry of the first mother for this particle in the event,
 *      0 means no mother
 *  4.  the entry of the second mother for this particle in the event,
 *      0 means no mother
 *  5. the entry of the first daughter for this particle in the event,
 *      0 means no daughter
 *  6. the entry of the second daughter for this particle in the event,
 *      0 means no daughter
 *  7. x component of the particle momentum
 *  8. y component of the particle momentum
 *  9. z component of the particle momentum
 *  10. energy of the particle
 *  11. mass of the particle
 *  12. x position of the particle initial position
 *  13. y position of the particle initial position
 *  14. z position of the particle initial position
 *  15. time of the particle production
 *
 *  For example, if you want to simulate a single muon with a 5 GeV energy
 *  moving only in the z direction, the entry would be (see onemuon.hepmc):
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  0 1
 *  1 13 0 0 0 0 0. 0. 1.0 5.0011 0.105 1.0 1.0 1.0 0.0
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  Or if you want to simulate a muon neutrino event (see oneneutrino.hepmc): 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  0 3
 *  0 14 0 0 0 0 0.00350383 0.002469 0.589751 0.589766 0 208.939 63.9671 10.9272 4026.32
 *  1 13 1 0 0 0 -0.168856 -0.0498011 0.44465 0.489765 105.658 208.939 63.9671 10.9272 4026.32
 *  1 2212 1 0 0 0 0.151902 -0.124578 0.0497377 0.959907 938.272 208.939 63.9671 10.9272 4026.32
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  here the first particle is the initial state neutrino (status code 0, meaning initial state, 
 *  not to be prooagated by GEANT). The other particles are the initial state particles.
 *  For the neutrino, write ccnc in place of 1st daugther, and mode (qe, res, ...) 
 *  in place of 2nd daugther.
 *
 *  There are some assumptions that go into using this format that may not
 *  be obvious.  The first is that only particles with status code = 1
 *  are tracked in the LArSoft/Geant4 combination making the mother daughter
 *  relations somewhat irrelevant. That also means that you should let
 *  Geant4 handle any decays.
 *
 *  The units in LArSoft are cm for distances and ns for time.
 *  The use of `TLorentzVector` below does not imply space and time have the same units
 *   (do not use `TLorentzVector::Boost()`).
 */
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <glob.h>
#include <cstdlib>  // for unsetenv()

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/search_path.h"
#include "cetlib/getenv.h"
#include "cetlib/split_path.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TLorentzVector.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "ifdh.h"       // use bare
#undef USE_IFDH_SERVICE // ifdh for now

// #ifndef NO_IFDH_LIB
//   #define USE_IFDH_SERVICE 1
//   // IFDHC
//   #ifdef USE_IFDH_SERVICE
//     #include "IFDH_service.h"
//   #else
//     // bare IFDHC
//     #include "ifdh.h"
//   #endif
// #else
//   #undef USE_IFDH_SERVICE
//   // nothing doing ... use ifdef to hide any reference that might need header
//   #include <cassert>
// #endif

namespace evgen {
  class HepMCFileGen;
}

class evgen::HepMCFileGen : public art::EDProducer {
public:
  explicit HepMCFileGen(fhicl::ParameterSet const & p);

  void produce(art::Event & e)        override;
  void beginJob()               		  override;
  void beginRun(art::Run & run) 		  override;
  void beginSubRun(art::SubRun& sr)     override;
  void endSubRun(art::SubRun& sr)     override;

private:
  void ExpandInputFilePatternsDirect();
  void ExpandInputFilePatternsIFDH();
  void open_file();


  //kinematic variable calculators

  // neutron mass
  const double M_n = 0.939565420;

  std::vector<simb::MCParticle> fParticle_vect;
  TLorentzVector q(std::vector<simb::MCParticle> P_vect);
  double Theorist_Q2(std::vector<simb::MCParticle> P_vect);
  double Theorist_W(std::vector<simb::MCParticle> P_vect);
  double Theorist_X(std::vector<simb::MCParticle> P_vect);
  double Theorist_Y(std::vector<simb::MCParticle> P_vect);

  double Experimentalist_Q2(std::vector<simb::MCParticle> P_vect);
  double Experimentalist_W(std::vector<simb::MCParticle> P_vect);
  double Experimentalist_X(std::vector<simb::MCParticle> P_vect);
  double Experimentalist_Y(std::vector<simb::MCParticle> P_vect);

  std::ifstream* fInputFile;

  //added by C Thorpe - metadata file
  std::ifstream* fMetadataFile;
  std::string fMetadataFileName; //File Containing metadata  
  bool UseMetadataFile;

  std::string              fFileSearchPaths; ///< colon separated set of path stems (to be set)
  std::vector<std::string> fFilePatterns;    ///< wildcard patterns files containing histograms or ntuples, or txt (to be set)
  std::vector<std::string> fSelectedFiles;   ///< flux files selected after wildcard expansion and subset selection
  std::string              fFileCopyMethod;  ///< "DIRECT" = old direct access method, otherwise = ifdh approach schema ("" okay)
  
  double         fEventsPerPOT;     ///< Number of events per POT (to be set)
  int            fEventsPerSubRun;  ///< Keeps track of the number of processed events per subrun

  ifdh_ns::ifdh* fIFDH;             ///< (optional) flux file handling

  int fMaxEvents;
  int fTotalEvents;

  //misc
  bool isNu(int pdg);
  bool isLepton(int pdg);

};

//------------------------------------------------------------------------------
evgen::HepMCFileGen::HepMCFileGen(fhicl::ParameterSet const & p)
  : EDProducer{p}
  , fInputFile(0)
  , fMetadataFileName{p.get<std::string>("MetadataFileName","")} //if no metadata file name given, revert to fhicl EventsPerPOT 
  , fFileSearchPaths{p.get<std::string>("FileSearchPaths")}
  , fFilePatterns{p.get<std::vector<std::string>>("FilePatterns")}
  , fFileCopyMethod{p.get<std::string>("FluxCopyMethod","DIRECT")}
  , fEventsPerPOT{p.get<double>("EventsPerPOT", -1.)}
  , fEventsPerSubRun(0)
  , fMaxEvents{p.get<int>("MaxEvents",100)}

{
  srand (time(0));

  if(fMetadataFileName != "") UseMetadataFile = true;
  else UseMetadataFile = false;

  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();
  produces< sumdata::POTSummary, art::InSubRun >();
  //std::cout << "Metadata file name is: " << fMetadataFileName << std::endl;
}

//------------------------------------------------------------------------------
void evgen::HepMCFileGen::open_file()
{

  int random_file_index = rand() / double(RAND_MAX) * fSelectedFiles.size(); 

  mf::LogInfo("HepMCFileGen")
      << "Opening file " << fSelectedFiles.at(random_file_index) << std::endl;;

  fInputFile = new std::ifstream(fSelectedFiles.at(random_file_index).c_str(), std::fstream::in);

  // check that the file is a good one
  if( !fInputFile->good() )
    throw cet::exception("HepMCFileGen") << "input text file "
          << fSelectedFiles.at(random_file_index)
          << " cannot be read.\n";
}


//------------------------------------------------------------------------------
void evgen::HepMCFileGen::beginJob()
{

  if (fFileCopyMethod == "DIRECT")       ExpandInputFilePatternsDirect();
  else if (fFileCopyMethod == "IFDH")    ExpandInputFilePatternsIFDH();
  else {
    throw cet::exception("HepMCFileGen") << "FluxCopyMethod "
          << fFileCopyMethod
          << " not supported.\n";
  }

  open_file();
  
  //open metadata file
  if( UseMetadataFile ){

	  //TODO: Add more robust way of locating this file
	  std::cout << "Metadata file is: " << fFileSearchPaths + "/" + fMetadataFileName << std::endl;
	  std::string MetadataFileFullName = fFileSearchPaths + "/" + fMetadataFileName;
	  fMetadataFile = new std::ifstream(MetadataFileFullName, std::fstream::in); 

	  //check metadata file is ok
	  if( !fMetadataFile->good() ) throw cet::exception("HepMCFileGen") << "Metadata file cannot be read in produce" << std::endl;

	  //set sentinel value
	  double POT_Per_Event=-1;

	  std::string line;
	  while(getline(*fMetadataFile,line)){

		if(line.find("POT_Per_Event=") != std::string::npos){
                          std::size_t equalspos = line.find("=");
                          POT_Per_Event = atof(line.substr(equalspos+1).c_str());
                  }
	  }

	if(POT_Per_Event == -1) throw cet::exception("HepMCFileGen") << "Metadata file " << MetadataFileFullName << " does not contain any line POT_Per_Event=" << std::endl;

  fEventsPerPOT = 1.0/POT_Per_Event;

  }

  //reset the total event counter
  fTotalEvents = 0;

}

//------------------------------------------------------------------------------
void evgen::HepMCFileGen::beginRun(art::Run& run)
{
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(geo->DetectorName()));
  }

//------------------------------------------------------------------------------
void evgen::HepMCFileGen::beginSubRun(art::SubRun& sr)
  {

	//reset event counter for this subrun
	fEventsPerSubRun = 0;  

    return;
  }

//------------------------------------------------------------------------------
void evgen::HepMCFileGen::endSubRun(art::SubRun& sr)
  {
 
   //generate the POT summary

   //std::cout << "Using EventsPerPOT = " << fEventsPerPOT << std::endl;
   //std::cout << "Events in this subrun = " <<  fEventsPerSubRun << std::endl;
   std::cout << "Subrun POT = " << fEventsPerSubRun / fEventsPerPOT << std::endl;

    auto p = std::make_unique<sumdata::POTSummary>();
    p->totpot     = fEventsPerSubRun / fEventsPerPOT;
    p->totgoodpot = fEventsPerSubRun / fEventsPerPOT;
    sr.put(std::move(p));

    return;
  }

//------------------------------------------------------------------------------
void evgen::HepMCFileGen::produce(art::Event & e)
{

  //throw exception if trying to process more events than limit
  if(fTotalEvents >= fMaxEvents) throw cet::exception("HepMCFileGen") << "Limit of " << fMaxEvents << " events/job exceeded"  << std::endl;  


  // check that the file is still good
  if( !fInputFile->good() ) {
    open_file();
  }

  if( !fInputFile->good() ) {
    throw cet::exception("HepMCFileGen") << "input text file "
					<< " cannot be read in produce().\n";
  }


  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth truth;

  // declare the variables for reading in the event record
  int            event          = 0;
  unsigned short nParticles 	  = 0;
  int            status         = 0;
  int 	 	 pdg            = 0;
  int 	 	 firstMother    = 0;
  int 	 	 secondMother   = 0;
  int 	 	 firstDaughter  = 0;
  int 	 	 secondDaughter = 0;
  double 	 xMomentum      = 0.;
  double 	 yMomentum   	= 0.;
  double 	 zMomentum   	= 0.;
  double 	 energy      	= 0.;
  double 	 mass        	= 0.;
  double 	 xPosition   	= 0.;
  double 	 yPosition   	= 0.;
  double 	 zPosition   	= 0.;
  double 	 time        	= 0.;

  bool set_neutrino = false;

  // neutrino
  int ccnc = -1, mode = -1, itype = -1, target = -1, nucleon = -1, quark = -1;
  double w = -1, x = -1, y = -1, qsqr = -1;

  // read in line to get event number and number of particles
  std::string oneLine;
  std::getline(*fInputFile, oneLine);
  std::istringstream inputLine;
  inputLine.str(oneLine);

 inputLine >> event >> nParticles;

 //clear the particle vector
 
 fParticle_vect.clear();

  // now read in all the lines for the particles
  // in this interaction. only particles with
  // status = 1 get tracked in Geant4. see GENIE GHepStatus
  for(unsigned short i = 0; i < nParticles; ++i){
    std::getline(*fInputFile, oneLine);
    inputLine.clear();
    inputLine.str(oneLine);

    inputLine >> status >> pdg
	      >> firstMother >> secondMother >> firstDaughter >> secondDaughter
	      >> xMomentum   >> yMomentum    >> zMomentum     >> energy >> mass
	      >> xPosition   >> yPosition    >> zPosition     >> time;

    TLorentzVector pos(xPosition, yPosition, zPosition, time);
    TLorentzVector mom(xMomentum, yMomentum, zMomentum, energy);

    simb::MCParticle part(i, pdg, "primary", firstMother, mass, status);
    part.AddTrajectoryPoint(pos, mom);

    //initial neutrino will always have status code 0
    if ( (abs(pdg) == 14 || abs(pdg) == 12) && status == 0 ) {
      set_neutrino = true;
      ccnc = firstDaughter; // for the neutrino we write ccnc in place of 1st daugther
      mode = secondDaughter; // for the neutrino we write mode in place of 2nd daugther
      itype = -1;
      target = nucleon = quark = w = x = y = qsqr = -1;
    } 

    truth.Add(part);
    fParticle_vect.push_back(part);
  }


/*
  //calculate kinematic variables
  double th_qsqr = Theorist_Q2(fParticle_vect);
  double th_w = Theorist_W(fParticle_vect);
  double th_x = Theorist_X(fParticle_vect);
  double th_y = Theorist_Y(fParticle_vect);

//  std::cout << "Theorist Q2 = " << th_qsqr << std::endl;
//  std::cout << "Theorist W = " << th_w << std::endl;
//  std::cout << "Theorist x = " << th_x << std::endl;
//  std::cout << "Theorist y = " << th_y << std::endl;
*/


  double  ex_qsqr = Experimentalist_Q2(fParticle_vect);
  double  ex_w = Experimentalist_W(fParticle_vect);
  double  ex_x = Experimentalist_X(fParticle_vect);
  double  ex_y = Experimentalist_Y(fParticle_vect);

//  std::cout << "Experimentalist Q2 = " << ex_qsqr << std::endl;
//  std::cout << "Experimentalist W = " << ex_w << std::endl;
//  std::cout << "Experimentalist x = " << ex_x << std::endl;
//  std::cout << "Experimentalist y = " << ex_y << std::endl;

// std::cout << "ex_y/th_y = " << ex_y/th_y << std::endl;

  qsqr = ex_qsqr;
  w = ex_w;
  x = ex_x;
  y = ex_y;


  if (set_neutrino) {

    truth.SetNeutrino(ccnc,
                      mode,
                      itype,
                      target,
                      nucleon,
                      quark,
                      w,
                      x,
                      y,
                      qsqr);

    // set the neutrino information in MCTruth
    truth.SetOrigin(simb::kBeamNeutrino);
    // truth.SetGeneratorInfo(simb::Generator_t::kGENIE,
    //                          __GENIE_RELEASE__,
    //                          {{"tune", fTuneName}});
  }

  truthcol->push_back(truth);


  e.put(std::move(truthcol));

  fEventsPerSubRun++;

  fTotalEvents++;

  return;
}


void evgen::HepMCFileGen::ExpandInputFilePatternsDirect() {

  std::vector<std::string> dirs;
  cet::split_path(fFileSearchPaths,dirs);

  glob_t g;
  int flags = GLOB_TILDE;   // expand ~ home directories

  std::ostringstream patterntext;  // for info/error messages
  std::ostringstream dirstext;     // for info/error messages

  std::vector<std::string>::const_iterator uitr = fFilePatterns.begin();
  int ipatt = 0;

  for ( ; uitr != fFilePatterns.end(); ++uitr, ++ipatt ) {
    std::string userpattern = *uitr;
    patterntext << "\n\t" << userpattern;
    std::vector<std::string>::const_iterator ditr = dirs.begin();
    for ( ; ditr != dirs.end(); ++ditr ) {
      std::string dalt = *ditr;
      // if non-null, does it end with a "/"?  if not add one
      size_t len = dalt.size();
      if ( len > 0 && dalt.rfind('/') != len-1 ) dalt.append("/");
      if ( uitr == fFilePatterns.begin() ) dirstext << "\n\t" << dalt;
      std::string filepatt = dalt + userpattern;
      glob(filepatt.c_str(),flags,NULL,&g);
      if ( g.gl_pathc > 0 ) flags |= GLOB_APPEND; // next glob() will append to list
    }  // loop over FluxSearchPaths dirs
  }  // loop over user patterns

  std::ostringstream paretext;
  std::ostringstream flisttext;
  int nfiles = g.gl_pathc;
  if ( nfiles == 0 ) {
    paretext << "\n  expansion resulted in a null list for flux files";
  } else { 
    // some sets of files should be left in order
    // and no size limitations imposed ... just copy the list
    paretext << "\n  list of files will be processed in order";
    for (int i=0; i<nfiles; ++i) {
      std::string afile(g.gl_pathv[i]);
      fSelectedFiles.push_back(afile);
      flisttext << "[" << std::setw(3) << i << "] "
                << afile << "\n";
    }
    mf::LogInfo("HepMCFileGen")
      << "ExpandFilePatternsDirect initially found " << nfiles
      << " files for user patterns:"
      << patterntext.str() << "\n  using FileSearchPaths of: "
      << dirstext.str() <<  "\n" << paretext.str();
      //<< "\"" << cet::getenv("FW_SEARCH_PATH") << "\"";
    mf::LogDebug("HepMCFileGen") << "\n" << flisttext.str();
    // done with glob list
    globfree(&g);
  }
}


void evgen::HepMCFileGen::ExpandInputFilePatternsIFDH() {

#ifdef NO_IFDH_LIB
    std::ostringstream fmesg;
    std::string marker = "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    fmesg << marker
          << __FILE__ << ":" << __LINE__
          << "\nno IFDH implemented on this platform\n"
          << marker;
    // make sure the message goes everywhere
    std::cout << fmesg.str() << std::flush;
    std::cerr << fmesg.str();
    throw cet::exception("Attempt to use ifdh class") << fmesg.str();
    assert(0);
#else
    // if "method" just an identifier and not a scheme then clear it
    if ( fFileCopyMethod.find("IFDH") == 0 ) fFileCopyMethod = "";
  #ifdef USE_IFDH_SERVICE
    art::ServiceHandle<IFDH> ifdhp;
  #else
    if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
  #endif
    std::string spaths = fFileSearchPaths;
    const char* ifdh_debug_env = std::getenv("IFDH_DEBUG_LEVEL");
    if ( ifdh_debug_env ) {
      mf::LogInfo("HepMCFileGen") << "IFDH_DEBUG_LEVEL: " << ifdh_debug_env;
      fIFDH->set_debug(ifdh_debug_env);
    }
    // filenames + size
    std::vector<std::pair<std::string,long>> partiallist, fulllist, selectedlist, locallist;
    std::ostringstream patterntext;  // stringification of vector of patterns
    std::ostringstream fulltext;     // for info on full list of matches
    std::ostringstream selectedtext; // for info on selected files
    std::ostringstream localtext;    // for info on local files
    fulltext << "search paths: " << spaths;
    //std::vector<std::string>::const_iterator uitr = fFilePatterns.begin();
    // loop over possible patterns
    // IFDH handles various stems but done a list of globs
    size_t ipatt=0;
    auto uitr = fFilePatterns.begin();
    for ( ; uitr != fFilePatterns.end(); ++uitr, ++ipatt ) {
      std::string userpattern = *uitr;
      patterntext << "\npattern [" << std::setw(3) << ipatt << "] " << userpattern;
      fulltext    << "\npattern [" << std::setw(3) << ipatt << "] " << userpattern;
  #ifdef USE_IFDH_SERVICE
      partiallist = ifdhp->findMatchingFiles(spaths,userpattern);
  #else
      partiallist = fIFDH->findMatchingFiles(spaths,userpattern);
  #endif
      fulllist.insert(fulllist.end(),partiallist.begin(),partiallist.end());
      // make a complete list ...
      fulltext << " found " << partiallist.size() << " files";
      for (auto pitr = partiallist.begin(); pitr != partiallist.end(); ++pitr) {
        fulltext << "\n  " << std::setw(10) << pitr->second << " " << pitr->first;
      }
      partiallist.clear();
    }  // loop over user patterns
    size_t nfiles = fulllist.size();
    mf::LogInfo("HepMCFileGen")
      << "ExpandFilePatternsIFDH initially found " << nfiles << " files";
    mf::LogDebug("HepMCFileGen")
      << fulltext.str();
    if ( nfiles == 0 ) {
      selectedtext << "\n  expansion resulted in a null list for flux files";
    } else {
      // some sets of files should be left in order
      // and no size limitations imposed ... just copy the list
      selectedtext << "\n  list of files will be processed in order";
      selectedlist.insert(selectedlist.end(),fulllist.begin(),fulllist.end());
    } 
    mf::LogInfo("HepMCFileGen") << "Fetching files." << std::endl;
    // have a selected list of remote files
    // get paths to local copies
  // #ifdef USE_IFDH_SERVICE
  //   locallist = ifdhp->fetchSharedFiles(selectedlist,fFileCopyMethod);
  // #else
  //   locallist = fIFDH->fetchSharedFiles(selectedlist,fFileCopyMethod);
  // #endif
    localtext << "final list of files:";
    size_t i=0;
    for (auto litr = selectedlist.begin(); litr != selectedlist.end(); ++litr, ++i) {
        fSelectedFiles.push_back(litr->first);
        localtext << "\n\t[" << std::setw(3) << i << "]\t" << litr->first;
      }
    mf::LogInfo("HepMCFileGen") << localtext.str();
    // no null path allowed for at least these
#endif  // 'else' code only if NO_IFDH_LIB not defined
  }


//------------------------------------------------------------------------------
// Kinematic Variable Calculators 
//------------------------------------------------------------------------------

TLorentzVector evgen::HepMCFileGen::q(std::vector<simb::MCParticle> P_vect){

   std::vector<TLorentzVector> P4_nu; //neutrino 4 momentum
   std::vector<TLorentzVector> P4_lepton; //outgoing lepton

   for(size_t i_part=0;i_part<P_vect.size();i_part++){

      //grab a particle from the vector
      simb::MCParticle Part = P_vect.at(i_part);
      TLorentzVector this_P4(Part.Px(),Part.Py(),Part.Pz(),Part.E());	  

      //if an initial state neutrino
      if(Part.StatusCode() == 0 && isNu(Part.PdgCode()))  P4_nu.push_back(this_P4);

      //if a final state nu/lepton
      if(Part.StatusCode() == 1 && (isNu(Part.PdgCode()) || isLepton(Part.PdgCode()))) P4_lepton.push_back(this_P4);
   }

   //check you have the right number of particles
   if(P4_nu.size() == 1 && P4_lepton.size() == 1){
      return P4_nu.at(0) - P4_lepton.at(0); 
   }
   else { std::cout << "Initial state neutrinos/final state leptons != 1, returning 0,0,0,0 for q" << std::endl; return TLorentzVector(0,0,0,0); }

}


// Squared 4 Momentum Transfer

double evgen::HepMCFileGen::Theorist_Q2(std::vector<simb::MCParticle> P_vect){

return (-1)*q(P_vect)*q(P_vect);

}




// Hadronic invariant mass

double evgen::HepMCFileGen::Theorist_W(std::vector<simb::MCParticle> P_vect){

   std::vector<TLorentzVector> P4_nu; //neutrino 4 momentum
   std::vector<TLorentzVector> P4_lepton; //outgoing lepton
   std::vector<TLorentzVector> P4_nucleon; //struck nucleon

   for(size_t i_part=0;i_part<P_vect.size();i_part++){

      //grab a particle from the vector
      simb::MCParticle Part = P_vect.at(i_part);
      TLorentzVector this_P4(Part.Px(),Part.Py(),Part.Pz(),Part.E());	  

      //if an initial state neutrino
      if(Part.StatusCode() == 0 && isNu(Part.PdgCode()))  P4_nu.push_back(this_P4);

      //if a final state nu/lepton
      if(Part.StatusCode() == 1 && (isNu(Part.PdgCode()) || isLepton(Part.PdgCode()))) P4_lepton.push_back(this_P4);

      //if struck nucleon (or other struck particle)
      if(Part.StatusCode() == 11) P4_nucleon.push_back(this_P4);
   }


   //starting 4 momentum
   TLorentzVector P4_initial = {0.0,0.0,0.0,0.0};
   for(size_t i_nu=0;i_nu<P4_nu.size();i_nu++) P4_initial += P4_nu.at(i_nu);
   for(size_t i_nuc=0;i_nuc<P4_nucleon.size();i_nuc++) P4_initial += P4_nucleon.at(i_nuc);

   //std::cout << "Initial 4 momentum: " << P4_initial.E() << "  " << P4_initial.Px() << "  " << P4_initial.Py() << "  " << P4_initial.Pz() << std::endl;
   //std::cout << "W = " << sqrt(P4_initial*P4_initial) << std::endl;
   return sqrt(P4_initial*P4_initial);
}




// Bjorken x

double evgen::HepMCFileGen::Theorist_X(std::vector<simb::MCParticle> P_vect){

   std::vector<TLorentzVector> P4_nucleon; //struck nucleon

   for(size_t i_part=0;i_part<P_vect.size();i_part++){

      //grab a particle from the vector
      simb::MCParticle Part = P_vect.at(i_part);
      TLorentzVector this_P4(Part.Px(),Part.Py(),Part.Pz(),Part.E());	  

      //if struck nucleon (or other struck particle)
      if(Part.StatusCode() == 11) P4_nucleon.push_back(this_P4);
   }



   TLorentzVector P4_initial_nuc = {0.0,0.0,0.0,0.0};
   for(size_t i_nuc=0;i_nuc<P4_nucleon.size();i_nuc++) P4_initial_nuc += P4_nucleon.at(i_nuc);

   //calculate Bjorken x
   return Theorist_Q2(P_vect)/(2*q(P_vect)*P4_initial_nuc);
}




// Bjorken y

double evgen::HepMCFileGen::Theorist_Y(std::vector<simb::MCParticle> P_vect){

   std::vector<TLorentzVector> P4_nu; //neutrino 4 momentum
   std::vector<TLorentzVector> P4_lepton; //outgoing lepton
   std::vector<TLorentzVector> P4_nucleon; //struck nucleon

   for(size_t i_part=0;i_part<P_vect.size();i_part++){

      //grab a particle from the vector
      simb::MCParticle Part = P_vect.at(i_part);
      TLorentzVector this_P4(Part.Px(),Part.Py(),Part.Pz(),Part.E());	  

      //if an initial state neutrino
      if(Part.StatusCode() == 0 && isNu(Part.PdgCode()))  P4_nu.push_back(this_P4);

      //if a final state nu/lepton
      if(Part.StatusCode() == 1 && (isNu(Part.PdgCode()) || isLepton(Part.PdgCode()))) P4_lepton.push_back(this_P4);

      //if struck nucleon (or other struck particle)
      if(Part.StatusCode() == 11) P4_nucleon.push_back(this_P4);
   }

   TLorentzVector P4_initial_nuc = {0.0,0.0,0.0,0.0};
   for(size_t i_nuc=0;i_nuc<P4_nucleon.size();i_nuc++) P4_initial_nuc += P4_nucleon.at(i_nuc);

   //calculate Bjorken y
   return q(P_vect)*P4_initial_nuc/(P4_nu.at(0)*P4_initial_nuc);

}



// Squared 4 Momentum Transfer

double evgen::HepMCFileGen::Experimentalist_Q2(std::vector<simb::MCParticle> P_vect){

return (-1)*q(P_vect)*q(P_vect);

}




// Hadronic invariant mass

double evgen::HepMCFileGen::Experimentalist_W(std::vector<simb::MCParticle> P_vect){

return sqrt(M_n*M_n+2*M_n*q(P_vect).E()-Experimentalist_Q2(P_vect));

}




// Bjorken x

double evgen::HepMCFileGen::Experimentalist_X(std::vector<simb::MCParticle> P_vect){


return Experimentalist_Q2(P_vect)/(2*M_n*q(P_vect).E());

}




// Bjorken y

double evgen::HepMCFileGen::Experimentalist_Y(std::vector<simb::MCParticle> P_vect){

   std::vector<TLorentzVector> P4_nu; //neutrino 4 momentum
   std::vector<TLorentzVector> P4_lepton; //outgoing lepton

   for(size_t i_part=0;i_part<P_vect.size();i_part++){

      //grab a particle from the vector
      simb::MCParticle Part = P_vect.at(i_part);
      TLorentzVector this_P4(Part.Px(),Part.Py(),Part.Pz(),Part.E());	  

      //if an initial state neutrino
      if(Part.StatusCode() == 0 && isNu(Part.PdgCode()))  P4_nu.push_back(this_P4);

      //if a final state nu/lepton
      if(Part.StatusCode() == 1 && (isNu(Part.PdgCode()) || isLepton(Part.PdgCode()))) P4_lepton.push_back(this_P4);
   }

     
   //check you have the right number of particles
   if(P4_nu.size() == 1 && P4_lepton.size() == 1){
      return 1-(P4_lepton.at(0).E()/P4_nu.at(0).E()); 
   }
   else return -1;

}


//------------------------------------------------------------------------------
// Miscellaneous 
//------------------------------------------------------------------------------

bool evgen::HepMCFileGen::isNu(int pdg){

return abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16;

}

bool evgen::HepMCFileGen::isLepton(int pdg){

return abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 15;

}

DEFINE_ART_MODULE(evgen::HepMCFileGen)

















