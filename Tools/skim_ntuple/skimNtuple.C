#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1I.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "../src/MuonPogTree.h"
#include "../src/Utils.h"

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
// 1. SampleConfig : configuration class containing sample information
// 2. TagAndProbeConfig : configuration class containing TnP cuts information
// 3. Plotter : class containing the plot definition and defining the plot filling 
//              for a given sample <= CB modify this to add new variables
// ******************************

namespace muon_pog {
  
  class SampleConfig {

  public :

    // config parameters (public for direct access)

    TString fileName;  
    TString outFileName;  
 
    SampleConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};

  };

  class SkimConfig {
    
  public :
    
    // config parameters (public for direct access)
    
    Float_t muon_minPt;
    Int_t muon_nMuonsAboveCut;
    std::string hlt_path; 
    std::vector<int> runs;

    bool clear_muonHitsAndMatches;
    bool clear_genParticles;
    bool clear_hlt;
   
    SkimConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    SkimConfig(boost::property_tree::ptree::value_type & vt); 
#endif
    
    ~SkimConfig() {};
    
  private:
    std::vector<Float_t> toArrayF(const std::string & entries); 
    std::vector<int> toArrayI(const std::string & entries);
    
  };
  
}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file
// 1. comparisonPlot : make a plot overlayng data and MC for a given plot
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, SkimConfig & skimConfig,
		   std::vector<SampleConfig> & sampleConfigs);
  
}



// The main program******** *****
// 1. Get configuration file and produces configurations
// 2. Create Plotters and loop on the event to fill them
// 3. Writes results in cnfigurable outuput file
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

  TString dirName = argv[2];
  system("mkdir -p " + dirName);
  
  std::cout << "[" << argv[0] << "] Using config file " << configFile << std::endl;

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  SkimConfig skimConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,skimConfig,sampleConfigs);

  for (auto sampleConfig : sampleConfigs)
    {

      TString fileName = sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileName.Data() << std::endl;  
      
      TFile* outputFile = TFile::Open(dirName + "/" + sampleConfig.outFileName,"RECREATE");
      // Initialize pointers to summary and full event structure
      
      outputFile->mkdir("MuonPogTree")->cd();

      int splitBranches = 2;
      TTree* outputTree = new TTree("MUONPOGTREE","Muon POG skimmed Tree");      

      muon_pog::Event* ev = new muon_pog::Event();
      muon_pog::Event* outputEv = new muon_pog::Event();

      TChain* tree;

      outputTree->Branch("event",&outputEv,64000,splitBranches);

      TH1I * eventHisto = new TH1I("eventHisto","n events before and after filter",2,0.5,2.5);

      // Open file, get tree, set branches
      //TFile* inputFile = TFile::Open(fileName,"READONLY");
      //tree = (TTree*)inputFile->Get("MUONPOGTREE");
      //if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);
      tree = openFileOrDir(fileName.Data());

      tree->SetBranchAddress("event", &ev);

      std::cout << "[" << argv[0] << "] Number of entries = " 
		<< tree->GetEntries() << std::endl;

      int nFilteredEvents = 0;
      for (Long64_t iEvent=0; iEvent<tree->GetEntries(); ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;
          if (iEvent % 25000 == 0 )
            std::cout << "[" << argv[0] << "] processing event : " 
		      << iEvent << "\r" << std::flush;
          tree->GetEvent(iEvent); 
	  
	  bool isGoodRun = false;
	  if (skimConfig.runs.size() == 1 && 
	      skimConfig.runs.at(0) <= 0)
	    isGoodRun = true;

	  if (!isGoodRun)
	    {
	      for(auto & run : skimConfig.runs)
		{
		  if (run == ev->runNumber)
		    {
		      isGoodRun = true;
		      break;
		    }
		}
	    }

	  if (!isGoodRun)
	    continue;

	  Int_t nGoodMuons = 0;
	  
	  for(auto & muon : ev->muons)
	    {
	      if (muon.pt         >= skimConfig.muon_minPt ||
		  muon.pt_tuneP   >= skimConfig.muon_minPt ||
		  muon.pt_tracker >= skimConfig.muon_minPt)
		{
		  nGoodMuons++;
		}
	    }

	  if (nGoodMuons <skimConfig.muon_nMuonsAboveCut)
	    continue;

	  nFilteredEvents++;
	  (*outputEv) = (*ev);

	  if (skimConfig.clear_hlt)
	    {
	      outputEv->hlt.triggers.clear();
	      outputEv->hlt.objects.clear();
	    }

	  if (skimConfig.clear_genParticles)
	    {
	      outputEv->genParticles.clear();
	    }
	  
	  if (skimConfig.clear_muonHitsAndMatches)
	    {
	      for(auto & muon : outputEv->muons)
		{
		  muon.hits.clear();
		  muon.matches.clear();
		}
	    }
	  
	  outputTree->Fill();
	
	
	}

      std::cout << "\n[" << argv[0] << "] FilterEfficiency = " << (Float_t(nFilteredEvents)/tree->GetEntries()) << std::endl;
      
      eventHisto->SetBinContent(1,tree->GetEntries());
      eventHisto->SetBinContent(2,nFilteredEvents);
      eventHisto->Write();

      outputTree->Write();

      delete ev;
      delete outputEv;
      //inputFile->Close();
      
      outputFile->Write();
      outputFile->Close();

    }
  
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;

}


muon_pog::SkimConfig::SkimConfig(boost::property_tree::ptree::value_type & vt)
{

  try
    {

      hlt_path = vt.second.get<std::string>("hlt_path");
      muon_nMuonsAboveCut = vt.second.get<Int_t>("muon_nMuonsAboveCut");
      muon_minPt = vt.second.get<Float_t>("muon_minPt");
      runs = toArrayI(vt.second.get<std::string>("runs"));
      clear_muonHitsAndMatches = vt.second.get<bool>("clear_muonHitsAndMatches");
      clear_genParticles = vt.second.get<bool>("clear_genParticles");
      clear_hlt = vt.second.get<bool>("clear_hlt");

    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[SkimConfig] Can' t get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }

}

muon_pog::SampleConfig::SampleConfig(boost::property_tree::ptree::value_type & vt)
{
 
 try
    {
      fileName     = TString(vt.second.get<std::string>("fileName").c_str());
      outFileName  = TString(vt.second.get<std::string>("outFileName").c_str());
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[SampleConfig] Cant'g get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}



std::vector<int> muon_pog::SkimConfig::toArrayI(const std::string& entries)
{

  std::vector<int> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(atoi(item.c_str()));
  return result;

}


std::vector<Float_t> muon_pog::SkimConfig::toArrayF(const std::string& entries)
{
  
  std::vector<Float_t> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(atof(item.c_str()));
  return result;

}


//Functions

void muon_pog::parseConfig(const std::string configFile, muon_pog::SkimConfig & skimConfig,
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
      if (vt.first.find("Cut") != std::string::npos)
	skimConfig = muon_pog::SkimConfig(vt);
      else
	sampleConfigs.push_back(muon_pog::SampleConfig(vt));
    }
}
