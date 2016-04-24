#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "../src/MuonPogTree.h"
#include "../src/Utils.h"
#include "tdrstyle.C"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <vector>
#include <regex>
#include <map>


#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/lexical_cast.hpp>

// Helper classes defintion *****
// 1. SampleConfig : configuration class containing sample parameters
// 2. AlgoConfig   : configuration class containing macro logic parameters
// ******************************

namespace muon_pog {
 
  class SampleConfig {

  public :

    // config parameters (public for direct access)

    TString fileName;  
    TString sampleName;  
    Float_t cSection;
    Float_t nEvents;
    Bool_t applyReweighting;
        
    SampleConfig() {};
    
#ifndef __MAKECINT__ // avoid CINT vs boost problems   
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};

  private:

  };

  class AlgoConfig {

  public :
    
    // config parameters (public for direct access)
    
    Float_t muon_minPt;      
    TString muon_ID;
   
    AlgoConfig() {};
    
#ifndef __MAKECINT__ // avoid CINT vs boost problems 
    AlgoConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~AlgoConfig() {};
    
  private:

  };
  
}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file (fill sample configs and algo one)
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, AlgoConfig & algoConfig,
		   std::vector<SampleConfig> & sampleConfigs);  
}



// The main program******** *****
// 1. Get configuration file and produces configurations
// 2. Create Plotters and loop on the event to fill them
// 3. Writes results in configurable outuput file
// ******************************

int main(int argc, char* argv[]){
  using namespace muon_pog;


  if (argc != 3) 
    {
      std::cout << "Usage : "
		<< argv[0] << " PATH_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR\n";
      exit(100);
    }

  std::string configFile(argv[1]);
  
  std::cout << "[" << argv[0] << "] Using config file " << configFile << std::endl;

  // Output directory
  TString dirName = argv[2];
  system("mkdir -p " + dirName);
  TFile* outputFile = TFile::Open(dirName + "/results.root","RECREATE"); // CB find a better name for output file  

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();
 
  AlgoConfig algoConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,algoConfig,sampleConfigs);

  std::map<TString, std::map<TString, TH1F *> > histos;
  for (auto sampleConfig : sampleConfigs)
    {
      TString histoName = "hPfIso_" + sampleConfig.sampleName;
      histos[sampleConfig.sampleName]["hPfIso"] = new TH1F(histoName,histoName+";rel PF isolation (dBeta); entries",100,0,5);
    }      

  for (auto sampleConfig : sampleConfigs)
    {

      TString fileName = sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileName.Data() << std::endl;   
  
      // Initialize pointers to summary and full event structure

      muon_pog::Event*   ev   = new muon_pog::Event();

      TTree* tree;
      TBranch* evBranch;

      // Open file, get tree, set branches

      TFile* inputFile = TFile::Open(fileName,"READONLY");
      tree = (TTree*)inputFile->Get("MUONPOGTREE");
      if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);

      evBranch = tree->GetBranch("event");
      evBranch->SetAddress(&ev);

      // Watch number of entries
      int nEntries = tree->GetEntriesFast();
      if (sampleConfig.nEvents > 0 && sampleConfig.nEvents < nEntries)
	nEntries = sampleConfig.nEvents;
      std::cout << "[" << argv[0] << "] Number of entries being processed = " << nEntries << std::endl;

      int nFilteredEvents = 0;

      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;

	  if (iEvent % 25000 == 0 )
	    std::cout << "[" << argv[0] << "] processing event : " << iEvent << "\r" << std::flush;
	    
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;

	  for (auto & muon : ev->muons)
	    {
	      if (muon.pt > algoConfig.muon_minPt &&
		  hasGoodId(muon,algoConfig.muon_ID))
		histos[sampleConfig.sampleName]["hPfIso"]->Fill(muon.isoPflow04,weight);
	    }
	}
      
      delete ev;
      inputFile->Close();
      std::cout << std::endl;
	   
    }
  
  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;

}


muon_pog::AlgoConfig::AlgoConfig(boost::property_tree::ptree::value_type & vt)
{

  try
    {
      muon_minPt = vt.second.get<Float_t>("muon_minPt");
      muon_ID    = vt.second.get<std::string>("muon_ID");
    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[AlgoConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }

}

muon_pog::SampleConfig::SampleConfig(boost::property_tree::ptree::value_type & vt)
{
  
  try
    {
      fileName     = TString(vt.second.get<std::string>("fileName").c_str());
      sampleName   = TString(vt.first.c_str());
      cSection = vt.second.get<Float_t>("cSection");
      nEvents = vt.second.get<Float_t>("nEvents"); //CB do we really need this? can't we take nEvents from the file itself?
      applyReweighting = vt.second.get<Bool_t>("applyReweighting");
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}

void muon_pog::parseConfig(const std::string configFile, muon_pog::AlgoConfig & algoConfig,
			   std::vector<muon_pog::SampleConfig> & sampleConfigs)
{

  boost::property_tree::ptree pt;
  
  try
    {
      boost::property_tree::ini_parser::read_ini(configFile, pt);
    }
  catch (boost::property_tree::ini_parser::ini_parser_error iniParseErr)
    {
      std::cout << "[parseConfig] Can't open : " << iniParseErr.filename()
		<< "\n\tin line : " << iniParseErr.line()
		<< "\n\thas error :" << iniParseErr.message()
		<< std::endl;
      throw std::runtime_error("Bad INI parsing");
    }
  
  for( auto vt : pt )
    {
      if (vt.first.find("Algo") != std::string::npos)
	algoConfig = muon_pog::AlgoConfig(vt);
      else
	sampleConfigs.push_back(muon_pog::SampleConfig(vt));
    }
}
