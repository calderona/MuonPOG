#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
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
    TString sampleName;  
    Float_t cSection;
    Int_t nEvents;
    Bool_t applyReweighting;
        
    SampleConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};
    
  };

  class TagAndProbeConfig {

  public :
    
    // config parameters (public for direct access)
    
    Float_t pair_minInvMass;
    Float_t pair_maxInvMass;    
    
    Float_t tag_minPt;      

    TString     tag_ID;
    Float_t     tag_isoCut;
    Float_t     tag_hltDrCut;
    std::string tag_hltFilter;
    
    std::string muon_trackType; // applies to both, tag and probe     
  
    Float_t probe_minPt;      
    std::vector<TString> probe_IDs;
    std::vector<std::pair<TString,TString> > probe_etaBins;

    std::vector<Float_t> pu_weights;
     
    std::string hlt_path; 
   
    TagAndProbeConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    TagAndProbeConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~TagAndProbeConfig() {};
    
  private:
    std::vector<TString> toArray(const std::string & entries); 
    std::vector<Float_t> toArrayF(const std::string & entries); 
    std::vector<std::pair<TString,TString> > toPairArray(const std::vector<TString> &,
							 const std::vector<TString> &); 
  
  };

  class Observable {

  public :
    
    Observable() {};

    Observable(TString hName, TString sampleTag, TString xTitle, TString yTitle,
	       Int_t nBins, Float_t min, Float_t max, bool kinPlots);
    
    ~Observable() { m_plots.clear(); };

    // void fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Float_t mcScale, Int_t PV);
    void fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Int_t PV);
    
    std::vector<TH1 *> & plots() { return m_plots; }
    
  private :
    
    std::vector<TH1 *> m_plots;

  };

  class Plotter {

  public :

    enum HistoType { KIN=0, MOM, ID, ISO, TIME, CONT, EFF};
      
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, int nVtx, float weight, int run);

    std::map<Plotter::HistoType, std::map<TString, muon_pog::Observable> > m_plots;
    TagAndProbeConfig m_tnpConfig;
    SampleConfig m_sampleConfig;
        
  };
  
}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file
// 1. comparisonPlot : make a plot overlayng data and MC for a given plot
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, TagAndProbeConfig & tpConfig,
		   std::vector<SampleConfig> & sampleConfigs);
  
  void comparisonPlots(std::vector<Plotter> & plotters,
		       TFile *outFile, TString &  outputDir);

  void copyPhp(const TString &  outputDir);

  // void setTProfY(TProfile &prof1, TProfile &prof2);
  void addUnderFlow(TH1 &hist);
  void addOverFlow(TH1 &hist);
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
 
  TagAndProbeConfig tnpConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,tnpConfig,sampleConfigs);

  std::vector<Plotter> plotters;

  for (auto sampleConfig : sampleConfigs)
    {

      Plotter plotter(tnpConfig, sampleConfig);
      plotter.book(outputFile);
      
      plotters.push_back(plotter);
    }
 
  for (auto & plotter : plotters)
    {

      TString fileName = plotter.m_sampleConfig.fileName;
      std::cout << "[" << argv[0] << "] Processing file "
		<< fileName.Data() << std::endl;  
  
      // Initialize pointers to summary and full event structure

      muon_pog::Event* ev = new muon_pog::Event();
      TTree* tree;
      TBranch* evBranch;

      // Open file, get tree, set branches

      TFile* inputFile = TFile::Open(fileName,"READONLY");
      tree = (TTree*)inputFile->Get("MUONPOGTREE");
      if (!tree) inputFile->GetObject("MuonPogTree/MUONPOGTREE",tree);

      evBranch = tree->GetBranch("event");
      evBranch->SetAddress(&ev);

      // Watch number of entries
      int nEntries = plotter.m_sampleConfig.nEvents > 0 ? plotter.m_sampleConfig.nEvents : tree->GetEntriesFast();
      std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

      if (nEntries < plotter.m_sampleConfig.nEvents || plotter.m_sampleConfig.nEvents < 0)
	plotter.m_sampleConfig.nEvents = nEntries;
      std::cout << "[" << argv[0] << "] Number of entries that will be processed = " << plotter.m_sampleConfig.nEvents << std::endl;

      int nFilteredEvents = 0;

      for (Long64_t iEvent=0; iEvent<plotter.m_sampleConfig.nEvents; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;

	  if (iEvent % 25000 == 0 )
	    std::cout << "[" << argv[0] << "] processing event : " << iEvent << "\r" << std::flush;	  

	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;

	  if(plotter.m_sampleConfig.applyReweighting==true)
	    weight *= ev->nVtx < 60 ? tnpConfig.pu_weights[ev->nVtx] : 0;
	    
	  plotter.fill(ev->muons, ev->hlt, ev->nVtx, weight, ev->runNumber);
	}
      
      delete ev;
      inputFile->Close();
    }
  
  muon_pog::comparisonPlots(plotters,outputFile,dirName);
  muon_pog::copyPhp(dirName);
  
  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;

}


muon_pog::TagAndProbeConfig::TagAndProbeConfig(boost::property_tree::ptree::value_type & vt)
{

  try
    {

      hlt_path = vt.second.get<std::string>("hlt_path");

      pair_minInvMass = vt.second.get<Float_t>("pair_minInvMass");
      pair_maxInvMass = vt.second.get<Float_t>("pair_maxInvMass");

      tag_minPt  = vt.second.get<Float_t>("tag_minPt");
      tag_ID     = vt.second.get<std::string>("tag_muonID");
      tag_isoCut = vt.second.get<Float_t>("tag_isoCut"); // CB for now just comb reliso dBeta R04

      tag_hltFilter = vt.second.get<std::string>("tag_hltFilter");
      tag_hltDrCut  = vt.second.get<Float_t>("tag_hltDrCut");
      
      muon_trackType = vt.second.get<std::string>("muon_trackType");

      probe_minPt  = vt.second.get<Float_t>("probe_minPt");
      probe_IDs    = toArray(vt.second.get<std::string>("probe_muonIDs"));
      probe_etaBins = toPairArray(toArray(vt.second.get<std::string>("probe_etaMin")),
				  toArray(vt.second.get<std::string>("probe_etaMax")));

      pu_weights = toArrayF(vt.second.get<std::string>("pu_weights"));
  
    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can' t get data : has error : "
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
      nEvents  = vt.second.get<Int_t>("nEvents");
      applyReweighting = vt.second.get<Bool_t>("applyReweighting");
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Cant'g get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}


std::vector<TString> muon_pog::TagAndProbeConfig::toArray(const std::string& entries)
{
  
  std::vector<TString> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(TString(item));
  return result;

}

std::vector<Float_t> muon_pog::TagAndProbeConfig::toArrayF(const std::string& entries)
{
  
  std::vector<Float_t> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(atof(item.c_str()));
  return result;

}


std::vector<std::pair<TString,TString> > muon_pog::TagAndProbeConfig::toPairArray(const std::vector<TString> & etaMin,
										  const std::vector<TString> & etaMax)
{

  std::vector<std::pair<TString,TString> > result;

  std::vector<TString>::const_iterator etaMinIt  = etaMin.begin();
  std::vector<TString>::const_iterator etaMinEnd = etaMin.end();

  std::vector<TString>::const_iterator etaMaxIt  = etaMax.begin();
  std::vector<TString>::const_iterator etaMaxEnd = etaMax.end();
  
  for (; etaMinIt != etaMinEnd || etaMaxIt != etaMaxEnd; ++etaMinIt, ++etaMaxIt)
    {
      result.push_back(std::make_pair((*etaMinIt),
				      (*etaMaxIt)));
    }

  return result;

}

muon_pog::Observable::Observable(TString hName, TString sampleTag, TString xTitle, TString yTitle,
				 Int_t nBins, Float_t min, Float_t max, bool kinPlots)
{
  m_plots.push_back(new TH1F("h" + hName + "_" + sampleTag, hName + " ;" + xTitle + ";" + yTitle, nBins, min, max));
  if (kinPlots)
    { // CB book here for other plots vs kin variables
      m_plots.push_back(new TProfile("h" + hName + "VsEta_"      + sampleTag, hName + " vs #eta;    #eta;"          + xTitle, 24, -2.4, 2.4, min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPhi_"      + sampleTag, hName + " vs #phi;    #phi;"          + xTitle, 24, -TMath::Pi(),TMath::Pi(), min, max));     
      m_plots.push_back(new TProfile("h" + hName + "VsPt_"       + sampleTag, hName + " vs p_{T};   p_{T} (GeV);"   + xTitle, 50,  0., 150., min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPV_"       + sampleTag, hName + " vs PV;      # of PV;"       + xTitle, 60,  0., 60., min, max));
    }
}

void muon_pog::Observable::fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Int_t PV)
{
  
  m_plots.at(0)->Fill(value, weight);
  // CB fill here other plots vs kin variables
  if (m_plots.size() > 1)
    {   
      ((TProfile*)m_plots.at(1))->Fill(muonTk.Eta(), value, weight);
      ((TProfile*)m_plots.at(2))->Fill(muonTk.Phi(), value, weight);
      ((TProfile*)m_plots.at(3))->Fill(muonTk.Pt(),  value, weight);
      ((TProfile*)m_plots.at(4))->Fill(PV,           value, weight);
    }
}

void muon_pog::Plotter::book(TFile *outFile)
{

  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);

  outFile->mkdir(sampleTag+"/efficiencies");      
  outFile->mkdir(sampleTag+"/timing");
  outFile->mkdir(sampleTag+"/id_variables");
  outFile->mkdir(sampleTag+"/kinematical_variables");
  outFile->mkdir(sampleTag+"/momentum_variables");
  outFile->mkdir(sampleTag+"/isolation");
  outFile->mkdir(sampleTag+"/control");
 
  outFile->cd(sampleTag+"/efficiencies");
 
  //Medium ID N-1 plots
  m_plots[EFF]["Medium_Numerator_eta"]       = muon_pog::Observable("Medium_Numerator_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_Numerator_pt"]        = muon_pog::Observable("Medium_Numerator_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_Numerator_phi"]       = muon_pog::Observable("Medium_Numerator_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);

  m_plots[EFF]["Medium_Step0_eta"]           = muon_pog::Observable("Medium_Step0_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_Step0_pt"]            = muon_pog::Observable("Medium_Step0_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_Step0_phi"]           = muon_pog::Observable("Medium_Step0_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);

  m_plots[EFF]["Medium_Step1_eta"]           = muon_pog::Observable("Medium_Step1_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_Step1_pt"]            = muon_pog::Observable("Medium_Step1_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_Step1_phi"]           = muon_pog::Observable("Medium_Step1_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Medium_Step2_eta"]           = muon_pog::Observable("Medium_Step2_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_Step2_pt"]            = muon_pog::Observable("Medium_Step2_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_Step2_phi"]           = muon_pog::Observable("Medium_Step2_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Medium_Step1not2_eta"]       = muon_pog::Observable("Medium_Step1not2_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_Step1not2_pt"]        = muon_pog::Observable("Medium_Step1not2_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_Step1not2_phi"]       = muon_pog::Observable("Medium_Step1not2_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Medium_Step2not1_eta"]       = muon_pog::Observable("Medium_Step2not1_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_Step2not1_pt"]        = muon_pog::Observable("Medium_Step2not1_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_Step2not1_phi"]       = muon_pog::Observable("Medium_Step2not1_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Medium_isGlobal_eta"]        = muon_pog::Observable("Medium_isGlobal_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_isGlobal_pt"]         = muon_pog::Observable("Medium_isGlobal_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_isGlobal_phi"]        = muon_pog::Observable("Medium_isGlobal_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
     
  m_plots[EFF]["Medium_glbNormChi2_eta"]     = muon_pog::Observable("Medium_glbNormChi2_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_glbNormChi2_pt"]      = muon_pog::Observable("Medium_glbNormChi2_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_glbNormChi2_phi"]     = muon_pog::Observable("Medium_glbNormChi2_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
   
  m_plots[EFF]["Medium_trkStaChi2_eta"]      = muon_pog::Observable("Medium_trkStaChi2_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_trkStaChi2_pt"]       = muon_pog::Observable("Medium_trkStaChi2_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_trkStaChi2_phi"]      = muon_pog::Observable("Medium_trkStaChi2_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
   
  m_plots[EFF]["Medium_trkKink_eta"]         = muon_pog::Observable("Medium_trkKink_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_trkKink_pt"]          = muon_pog::Observable("Medium_trkKink_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_trkKink_phi"]         = muon_pog::Observable("Medium_trkKink_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
   
  m_plots[EFF]["Medium_muSegmCompL_eta"]     = muon_pog::Observable("Medium_muSegmCompL_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Medium_muSegmCompL_pt"]      = muon_pog::Observable("Medium_muSegmCompL_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Medium_muSegmCompL_phi"]     = muon_pog::Observable("Medium_muSegmCompL_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
 
  //Tight ID N-1 plots
  m_plots[EFF]["Tight_Numerator_eta"]        = muon_pog::Observable("Tight_Numerator_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_Numerator_pt"]         = muon_pog::Observable("Tight_Numerator_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_Numerator_phi"]        = muon_pog::Observable("Tight_Numerator_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);

  m_plots[EFF]["Tight_isGlobal_eta"]         = muon_pog::Observable("Tight_isGlobal_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_isGlobal_pt"]          = muon_pog::Observable("Tight_isGlobal_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_isGlobal_phi"]         = muon_pog::Observable("Tight_isGlobal_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);

  m_plots[EFF]["Tight_isPF_eta"]             = muon_pog::Observable("Tight_isPF_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_isPF_pt"]              = muon_pog::Observable("Tight_isPF_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_isPF_phi"]             = muon_pog::Observable("Tight_isPF_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);

  m_plots[EFF]["Tight_glbNormChi2_eta"]      = muon_pog::Observable("Tight_glbNormChi2_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_glbNormChi2_pt"]       = muon_pog::Observable("Tight_glbNormChi2_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_glbNormChi2_phi"]      = muon_pog::Observable("Tight_glbNormChi2_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Tight_glbMuonValidHits_eta"] = muon_pog::Observable("Tight_glbMuonValidHits_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_glbMuonValidHits_pt"]  = muon_pog::Observable("Tight_glbMuonValidHits_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_glbMuonValidHits_phi"] = muon_pog::Observable("Tight_glbMuonValidHits_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Tight_trkMuonMatchedStations_eta"]   = muon_pog::Observable("Tight_trkMuonMatchedStations_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_trkMuonMatchedStations_pt"]    = muon_pog::Observable("Tight_trkMuonMatchedStations_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_trkMuonMatchedStations_phi"]   = muon_pog::Observable("Tight_trkMuonMatchedStations_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
   
  m_plots[EFF]["Tight_dxy_eta"]              = muon_pog::Observable("Tight_dxy_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_dxy_pt"]               = muon_pog::Observable("Tight_dxy_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_dxy_phi"]              = muon_pog::Observable("Tight_dxy_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
   
  m_plots[EFF]["Tight_dz_eta"]               = muon_pog::Observable("Tight_dz_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_dz_pt"]                = muon_pog::Observable("Tight_dz_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_dz_phi"]               = muon_pog::Observable("Tight_dz_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Tight_trkPixelValidHits_eta"]= muon_pog::Observable("Tight_trkPixelValidHits_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_trkPixelValidHits_pt"] = muon_pog::Observable("Tight_trkPixelValidHits_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_trkPixelValidHits_phi"]= muon_pog::Observable("Tight_trkPixelValidHits_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);
  
  m_plots[EFF]["Tight_trkTrackerLayerWithMeas_eta"]  = muon_pog::Observable("Tight_trkTrackerLayerWithMeas_eta", sampleTag, "#eta", "# entries", 48, -2.4, 2.4, false);
  m_plots[EFF]["Tight_trkTrackerLayerWithMeas_pt"]   = muon_pog::Observable("Tight_trkTrackerLayerWithMeas_pt", sampleTag, "p_{T} (GeV)", "# entries", 75, 0., 150., false);
  m_plots[EFF]["Tight_trkTrackerLayerWithMeas_phi"]  = muon_pog::Observable("Tight_trkTrackerLayerWithMeas_phi", sampleTag, "#phi", "# entries", 48, -TMath::Pi(), TMath::Pi(), false);

  outFile->cd(sampleTag+"/timing");
  
  //std::cout << sampleTag << "  Begin: " << m_plots.size()<< std::endl;

  m_plots[TIME]["STAmuonTime"]       = muon_pog::Observable("STAmuonTime", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["STAmuonTimeBarrel"] = muon_pog::Observable("STAmuonTimeBarrel", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["STAmuonTimeEndcap"] = muon_pog::Observable("STAmuonTimeEndcap", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);

  m_plots[TIME]["UnbSTAmuonTime"] = muon_pog::Observable("UnbSTAmuonTime", sampleTag, "time (ns)", "# entries", 400, -200., 200., true);
  m_plots[TIME]["UnbSTAmuonTimeBarrel"] = muon_pog::Observable("UnbSTAmuonTimeBarrel", sampleTag, "time (ns)", "# entries", 400, -200., 200., true);
  m_plots[TIME]["UnbSTAmuonTimeEndcap"] = muon_pog::Observable("UnbSTAmuonTimeEndcap", sampleTag, "time (ns)", "# entries", 400, -200., 200., true);

  for (auto etaBin : m_tnpConfig.probe_etaBins)
    {
      
      TString etaTag = "_etaMin" + etaBin.first + "_etaMax" + etaBin.second;

      outFile->cd(sampleTag+"/id_variables");

      //Id var
      
      m_plots[ID]["NHitsGLB" + etaTag] = muon_pog::Observable("NHitsGLB" + etaTag, sampleTag, "# hits", "# entries", 80, 0., 80., true);
      m_plots[ID]["NHitsTRK" + etaTag] = muon_pog::Observable("NHitsTRK" + etaTag, sampleTag, "# hits", "# entries", 40, 0., 40., true);
      m_plots[ID]["NHitsSTA" + etaTag] = muon_pog::Observable("NHitsSTA" + etaTag, sampleTag, "# hits", "# entries", 60, 0., 60., true);
      
      m_plots[ID]["Chi2GLB" + etaTag]  = muon_pog::Observable("Chi2GLB_" + etaTag, sampleTag, "chi2/ndof", "# entries", 50, 0., 100., true);
      m_plots[ID]["Chi2TRK" + etaTag]  = muon_pog::Observable("Chi2TRK_" + etaTag, sampleTag, "chi2/ndof", "# entries", 100, 0., 50., true);
      
      m_plots[ID]["NMatchedStation" + etaTag]  = muon_pog::Observable("NMatchedStation_" + etaTag, sampleTag,"# stations", "# entries", 10, 0., 10., true);
      m_plots[ID]["NMuonValidHitsGLB" + etaTag]= muon_pog::Observable("NMuonValidHitsGLB_" + etaTag, sampleTag,"# hits", "# entries", 60, 0., 60., true);
      m_plots[ID]["PixelHitsTRK" + etaTag]    = muon_pog::Observable("PixelHitsTRK" + etaTag, sampleTag,"# hits", "# entries", 10, 0., 10., true);
      m_plots[ID]["PixelLayersTRK" + etaTag]  = muon_pog::Observable("PixelLayersTRK" + etaTag, sampleTag,"# layers", "# entries", 10, 0., 10., true);
      m_plots[ID]["TrackerLayersTRK" + etaTag]= muon_pog::Observable("TrackerLayersTRK" + etaTag, sampleTag,"# layers"," # entries", 30, 0., 30., true);
      
      m_plots[ID]["HitFractionTRK" + etaTag]  = muon_pog::Observable("HitFractionTRK" + etaTag, sampleTag, "fraction", " # entries", 20, 0., 1., true);
      m_plots[ID]["TrkStaChi2" + etaTag]      = muon_pog::Observable("TrkStaChi2" + etaTag, sampleTag, "chi2", "# entries", 50, 0., 100., true);
      m_plots[ID]["TrkKink" + etaTag]         = muon_pog::Observable("TrkKink" + etaTag, sampleTag, "prob.", "# entries", 100, 0., 250., true);
      m_plots[ID]["SegmentComp" + etaTag]     = muon_pog::Observable("SegmentComp" + etaTag, sampleTag,"prob.", "segmentComp", 100, 0., 1., true);
      
      m_plots[ID]["Dxy" + etaTag]             = muon_pog::Observable("Dxy" + etaTag, sampleTag, "d_{xy} (cm)", "# entries", 50,-0.5,0.5, true);
      m_plots[ID]["Dz" + etaTag]              = muon_pog::Observable("Dz" + etaTag, sampleTag, "d_{z} (cm)", "# entries", 60,-1.5,1.5, true);    

      m_plots[ID]["qOverPtTrkSta" + etaTag]      = muon_pog::Observable("qOverPtTrkSta" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.05,0.05, true); 
      // m_plots[ID]["qOverPtTrkStaOverPt" + etaTag]      = muon_pog::Observable("qOverPtTrkStaOverPt" + etaTag, sampleTag, "(q/p_{T}^{sta} - q/p_{T}^{trk})/q/p_{T}^{trk}", "# entries", 50,-5.,5., true); 
      // m_plots[ID]["qOverPtTrkStaPlus" + etaTag]      = muon_pog::Observable("qOverPtTrkStaPlus" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.03,0.03, true); 
      m_plots[ID]["qOverPtTrkSta200" + etaTag]      = muon_pog::Observable("qOverPtTrkSta200" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 
      // m_plots[ID]["qOverPtTrkSta200Plus" + etaTag]      = muon_pog::Observable("qOverPtTrkSta200Plus" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 

      m_plots[ID]["qOverPtTrkGlb" + etaTag]      = muon_pog::Observable("qOverPtTrkGlb" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.03,0.03, true); 
      // m_plots[ID]["qOverPtTrkGlbPlus" + etaTag]      = muon_pog::Observable("qOverPtTrkGlbPlus" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.03,0.03, true); 
      m_plots[ID]["qOverPtTrkGlb200" + etaTag]      = muon_pog::Observable("qOverPtTrkGlb200" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 
      // m_plots[ID]["qOverPtTrkGlb200Plus" + etaTag]      = muon_pog::Observable("qOverPtTrkGlb200Plus" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 
      m_plots[TIME]["TightMuonTime" + etaTag] = muon_pog::Observable("TightMuonTime", sampleTag + etaTag, "time (ns)", "# entries", 200, -100., 100., true);

      for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  outFile->cd(sampleTag+"/kinematical_variables");

	  m_plots[KIN]["ProbePt" + etaTag + IDTag]  = muon_pog::Observable("ProbePt" + etaTag + IDTag, sampleTag, "p_{T} (GeV)", "# entries", 75,0.,150., true);
	  m_plots[KIN]["ProbeEta" + etaTag + IDTag] = muon_pog::Observable("hProbeEta_" + etaTag + IDTag, sampleTag, "#eta", "# entries", 48,-2.4, 2.4, true);
	  m_plots[KIN]["ProbePhi" + etaTag + IDTag] = muon_pog::Observable("hProbePhi_" + etaTag + IDTag, sampleTag, "#phi", "# entries", 48,-TMath::Pi(),TMath::Pi(), true); 

	  outFile->cd(sampleTag+"/momentum_variables");

	  m_plots[MOM]["goodMuMass" + etaTag + IDTag]  = muon_pog::Observable("goodMuMass" + etaTag + IDTag, sampleTag, "invariant mass (GeV)", "# entries", 30,85.,115., false);
	  m_plots[MOM]["goodMuMassPlus" + etaTag + IDTag]  = muon_pog::Observable("goodMuMassPlus" + etaTag + IDTag, sampleTag, "invariant mass (GeV)", "# entries", 20 ,86.5,96.5, true);
	  m_plots[MOM]["goodMuMassMinus" + etaTag + IDTag]  = muon_pog::Observable("goodMuMassMinus" + etaTag + IDTag, sampleTag, "invariant mass (GeV)", "# entries", 20,86.5,96.5, true);   
	}
  
      outFile->cd(sampleTag+"/isolation");

      for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  m_plots[ISO]["ChHadIso" + etaTag + IDTag]    = muon_pog::Observable("ChHadIso_"    + etaTag + IDTag, sampleTag, "charged had. iso 0.4", "# entries", 50,0.,10., true);
	  m_plots[ISO]["ChHadIsoPU" + etaTag + IDTag]  = muon_pog::Observable("ChHadIsoPU_"  + etaTag + IDTag, sampleTag, "PU Charged had. iso 0.4", "# entries", 50,0.,10., true);
	  m_plots[ISO]["PhotonIso" + etaTag + IDTag]   = muon_pog::Observable("PhotonIso_"   + etaTag + IDTag, sampleTag, "photon iso 0.4", " # entries", 50,0.,10., true);
	  m_plots[ISO]["NeutralIso" + etaTag + IDTag]  = muon_pog::Observable("NeutralIso_"  + etaTag + IDTag, sampleTag, "neutral had. iso 0.4", "# entries", 50,0.,10., true);
	  m_plots[ISO]["DBetaRelIso" + etaTag + IDTag] = muon_pog::Observable("DBetaRelIso_" + etaTag + IDTag, sampleTag, "PFIso 0.4 (dBeta)", "# entries", 50,0.,2., true);
	}
    }

  outFile->cd(sampleTag+"/control");
  
  m_plots[CONT]["01-invMass"] = muon_pog::Observable("invMass", sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
  m_plots[CONT]["02-dilepPt"] = muon_pog::Observable("dilepPt", sampleTag ,"p_{T} (GeV)", "# entries", 100,0.,200., false);
  m_plots[CONT]["03-nVertices"] = muon_pog::Observable("nVertices", sampleTag ,"# vertices", "# entries", 60,0.,60., false);
  m_plots[CONT]["04-runNumber"] = muon_pog::Observable("runNumber", sampleTag ,"run number", "# entries", 9000,253000.,262000., false);

  m_plots[CONT]["99-invMassInRange"] = muon_pog::Observable("invMassInRange", sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);

  
  //  std::cout << sampleTag << "  End: " << m_plots.size()<< std::endl;
  

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, int nVtx, float weight, Int_t run)
{
  
  TLorentzVector emptyTk;
    
  //muon timing only
  for (auto & muon : muons)
    {
      if (muon.isStandAlone && ((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) ||
				(fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)))
	m_plots[TIME]["STAmuonTime"].fill(muon.muonTime,emptyTk,weight,nVtx);
      
      if (muon.isStandAlone && fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20)
	m_plots[TIME]["STAmuonTimeBarrel"].fill(muon.muonTime,emptyTk,weight,nVtx);
      
      if (muon.isStandAlone && fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11)
	m_plots[TIME]["STAmuonTimeEndcap"].fill(muon.muonTime,emptyTk,weight,nVtx);
    }
  
  if (!muon_pog::pathHasFired(hlt,m_tnpConfig.hlt_path)) return;

  std::vector<const muon_pog::Muon *> tagMuons;
  
  for (auto & muon : muons)
    {
      if (muon_pog::hasGoodId(muon,m_tnpConfig.tag_ID) && 
	  muon_pog::muonTk(muon,m_tnpConfig.muon_trackType).Pt() > 
	  muon_pog::hasFilterMatch(muon,hlt,
				   m_tnpConfig.tag_hltFilter,
				   m_tnpConfig.tag_hltDrCut) &&
	  muon.isoPflow04 < m_tnpConfig.tag_isoCut)
	tagMuons.push_back(&muon);
    }
  
  std::vector<const muon_pog::Muon *> probeMuons;
  

  for (auto & muon : muons)
    {
      for (auto tagMuonPointer : tagMuons)
	{
	  const muon_pog::Muon & tagMuon = *tagMuonPointer;
	  
	  if ( tagMuonPointer != &muon && 
	       muon_pog::chargeFromTrk(tagMuon,m_tnpConfig.muon_trackType) *
	       muon_pog::chargeFromTrk(muon,m_tnpConfig.muon_trackType) == -1)
	    {

	      // General Probe Muons	      
	      if((muon.isGlobal || muon.isTrackerArb) &&  // CB minimal cuts on potental probe 
		 muon.pt > m_tnpConfig.probe_minPt)
		{
		  TLorentzVector tagMuTk(muon_pog::muonTk(tagMuon,m_tnpConfig.muon_trackType));
		  TLorentzVector muTk(muon_pog::muonTk(muon,m_tnpConfig.muon_trackType));

		  if((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) || (fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)) 
		    m_plots[TIME]["UnbSTAmuonTime"].fill(muon.muonTime,muTk,weight,nVtx);
		  
		  if(fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20) 
		    m_plots[TIME]["UnbSTAmuonTimeBarrel"].fill(muon.muonTime,muTk,weight,nVtx);
		
		  if(fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11) 
		    m_plots[TIME]["UnbSTAmuonTimeEndcap"].fill(muon.muonTime,muTk,weight,nVtx);
		  
		  Float_t mass = (tagMuTk+muTk).M();
		  	  
		  // CB Fill control plots
		  m_plots[CONT]["01-invMass"].fill(mass,emptyTk,weight,nVtx);
		  if ( mass > m_tnpConfig.pair_minInvMass &&
		       mass < m_tnpConfig.pair_maxInvMass )
		    {
		      m_plots[CONT]["99-invMassInRange"].fill(mass,emptyTk,weight,nVtx);
		      
		      Float_t dilepPt = (tagMuTk+muTk).Pt();
		      m_plots[CONT]["02-dilepPt"].fill(dilepPt,emptyTk,weight,nVtx);
		      m_plots[CONT]["03-nVertices"].fill(nVtx,emptyTk,weight,nVtx);
		      
		      for (auto etaBin : m_tnpConfig.probe_etaBins)
			{
			  if (muTk.Eta() > etaBin.first.Atof() &&
			      muTk.Eta() < etaBin.second.Atof() )
			    {

			      TString etaTag = "_etaMin" + etaBin.first + "_etaMax" + etaBin.second;
	
			      for (auto & probe_ID : m_tnpConfig.probe_IDs)
				{
				  TString IDTag = "_" + probe_ID;

				  if(hasGoodId(muon,probe_ID) && muon.isoPflow04 < 0.25 && mass > 60 && mass < 115)
				    {
				    m_plots[MOM]["goodMuMass" + etaTag + IDTag].fill(mass, muTk, weight, nVtx);
				    if (mass > 86.5 && mass < 96.5)
				      {
					if (chargeFromTrk(muon,m_tnpConfig.muon_trackType) > 0)
					  m_plots[MOM]["goodMuMassPlus" + etaTag + IDTag].fill(mass, muTk, weight, nVtx);
					else
					  m_plots[MOM]["goodMuMassMinus" + etaTag + IDTag].fill(mass, muTk, weight, nVtx);
				      }
				    }
				}
			    }
			}
		      
		      probeMuons.push_back(&muon);
		      
		      //FP: muons to compute the MediumID N-1 Eff.
		      if(muon.isMedium){
			m_plots[EFF]["Medium_Numerator_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Medium_Numerator_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Medium_Numerator_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }

		      if(muon.trkValidHitFrac > 0.8 && muon.isLoose)
			{
			  m_plots[EFF]["Medium_Step0_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			  m_plots[EFF]["Medium_Step0_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			  m_plots[EFF]["Medium_Step0_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  
			  if(muon.isGlobal && muon.glbNormChi2 < 3. && muon.trkStaChi2 < 12. && muon.trkKink < 20. && muon.muSegmComp > 0.303){
			    m_plots[EFF]["Medium_Step1_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_Step1_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_Step1_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			    
			    if(muon.muSegmComp <= 0.451){ //FP here the sign is inverted
			      m_plots[EFF]["Medium_Step1not2_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			      m_plots[EFF]["Medium_Step1not2_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			      m_plots[EFF]["Medium_Step1not2_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			    }
			  }//end step1
			  else if(muon.muSegmComp > 0.451){
			    m_plots[EFF]["Medium_Step2not1_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_Step2not1_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_Step2not1_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  }

			
			  //N-1 of step1
			  if(muon.glbNormChi2 < 3. && muon.trkStaChi2 < 12. && muon.trkKink < 20. && muon.muSegmComp > 0.303){
			    m_plots[EFF]["Medium_isGlobal_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_isGlobal_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_isGlobal_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  }
			  
			  if(muon.isGlobal && muon.trkStaChi2 < 12. && muon.trkKink < 20. && muon.muSegmComp > 0.303 ){
			    m_plots[EFF]["Medium_glbNormChi2_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_glbNormChi2_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_glbNormChi2_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  }
			  
			  if(muon.isGlobal && muon.glbNormChi2 < 3. && muon.trkKink < 20. && muon.muSegmComp > 0.303 ){
			    m_plots[EFF]["Medium_trkStaChi2_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_trkStaChi2_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_trkStaChi2_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  }
			  
			  if(muon.isGlobal && muon.glbNormChi2 < 3. && muon.trkStaChi2 < 12. && muon.muSegmComp > 0.303 ){
			    m_plots[EFF]["Medium_trkKink_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_trkKink_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_trkKink_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  }
			  
			  if(muon.isGlobal && muon.glbNormChi2 < 3. && muon.trkStaChi2 < 12. && muon.trkKink < 20.) {
			    m_plots[EFF]["Medium_muSegmCompL_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_muSegmCompL_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_muSegmCompL_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  }
					  
			  //step2  
			  if(muon.muSegmComp > 0.451){
			    m_plots[EFF]["Medium_Step2_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_Step2_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			    m_plots[EFF]["Medium_Step2_phi"].fill(muon.phi,emptyTk,weight,nVtx);
			  }
			}
		      
		     
		      //FP: muons to compute the TightID N-1 efficiency
		      if(muon.isTight){
			m_plots[EFF]["Tight_Numerator_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_Numerator_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_Numerator_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      
		      if(muon.isPF && muon.glbNormChi2 < 10. && muon.glbMuonValidHits > 0 && muon.trkMuonMatchedStations > 1 && fabs(muon.dxyBest) < 0.2 &&
		      	 fabs(muon.dzBest) < 0.5 && muon.trkPixelValidHits > 0 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_isGlobal_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_isGlobal_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_isGlobal_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      		      
		      if(muon.isGlobal && muon.glbNormChi2 < 10. && muon.glbMuonValidHits > 0 && muon.trkMuonMatchedStations > 1 && fabs(muon.dxyBest) < 0.2 &&
		      	 fabs(muon.dzBest) < 0.5 && muon.trkPixelValidHits > 0 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_isPF_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_isPF_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_isPF_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		     
  		      if(muon.isGlobal && muon.isPF && muon.glbMuonValidHits > 0 && muon.trkMuonMatchedStations > 1 && fabs(muon.dxyBest) < 0.2 &&
		      	 fabs(muon.dzBest) < 0.5 && muon.trkPixelValidHits > 0 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_glbNormChi2_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_glbNormChi2_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_glbNormChi2_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      	
		      if(muon.isGlobal && muon.isPF && muon.glbNormChi2 < 10. && muon.trkMuonMatchedStations > 1 && fabs(muon.dxyBest) < 0.2 &&
		      	 fabs(muon.dzBest) < 0.5 && muon.trkPixelValidHits > 0 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_glbMuonValidHits_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_glbMuonValidHits_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_glbMuonValidHits_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      			      
		      if(muon.isGlobal && muon.isPF && muon.glbNormChi2 < 10. && muon.glbMuonValidHits > 0  && fabs(muon.dxyBest) < 0.2 &&
		      	 fabs(muon.dzBest) < 0.5 && muon.trkPixelValidHits > 0 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_trkMuonMatchedStations_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_trkMuonMatchedStations_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_trkMuonMatchedStations_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      			      
		      if(muon.isGlobal && muon.isPF && muon.glbNormChi2 < 10. && muon.glbMuonValidHits > 0 && muon.trkMuonMatchedStations > 1 &&
		      	 fabs(muon.dzBest) < 0.5 && muon.trkPixelValidHits > 0 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_dxy_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_dxy_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_dxy_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      
		      if(muon.isGlobal && muon.isPF && muon.glbNormChi2 < 10. && muon.glbMuonValidHits > 0 && muon.trkMuonMatchedStations > 1 && 
			 fabs(muon.dxyBest) < 0.2 && muon.trkPixelValidHits > 0 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_dz_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_dz_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_dz_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      			      
		      if(muon.isGlobal && muon.isPF && muon.glbNormChi2 < 10. && muon.glbMuonValidHits > 0 && muon.trkMuonMatchedStations > 1 && 
			 fabs(muon.dxyBest) < 0.2 && fabs(muon.dzBest) < 0.5 && muon.trkTrackerLayersWithMeas > 5){
			m_plots[EFF]["Tight_trkPixelValidHits_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_trkPixelValidHits_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_trkPixelValidHits_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      			      
		      if(muon.isGlobal && muon.isPF && muon.glbNormChi2 < 10. && muon.glbMuonValidHits > 0 && muon.trkMuonMatchedStations > 1 && 
			 fabs(muon.dxyBest) < 0.2 && fabs(muon.dzBest) < 0.5 && muon.trkPixelValidHits > 0 ){
			m_plots[EFF]["Tight_trkTrackerLayerWithMeas_eta"].fill(muon.eta,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_trkTrackerLayerWithMeas_pt"].fill(muon.pt,emptyTk,weight,nVtx);
			m_plots[EFF]["Tight_trkTrackerLayerWithMeas_phi"].fill(muon.phi,emptyTk,weight,nVtx);
		      }
		      		      
		      continue; // CB If a muon is already a probe don't loo on other tags
		    }
		}
	    }
	}
    }
  
  // m_plots[CONT]["04_nProbesVsnTags"]->Fill(tagMuons.size(),probeMuons.size());

  m_plots[CONT]["04-runNumber"].fill(run,emptyTk,weight,nVtx);
		     
  
  for (auto probeMuonPointer : probeMuons)
    {
      const muon_pog::Muon & probeMuon = *probeMuonPointer;
      			  
      for (auto etaBin : m_tnpConfig.probe_etaBins)
	{

	  TLorentzVector probeMuTk(muon_pog::muonTk(probeMuon,m_tnpConfig.muon_trackType));
	  
	  if (probeMuTk.Eta() > etaBin.first.Atof() &&
	      probeMuTk.Eta() < etaBin.second.Atof() )
	    {
	      
	      TString etaTag = "_etaMin" + etaBin.first + "_etaMax" + etaBin.second;

	      //id var
	      if (probeMuon.isGlobal)
		{
		  m_plots[ID]["NHitsGLB"  + etaTag].fill(probeMuon.nHitsGlobal, probeMuTk ,weight,nVtx);
		  m_plots[ID]["Chi2GLB"  + etaTag].fill(probeMuon.glbNormChi2, probeMuTk, weight, nVtx);
		  m_plots[ID]["NMuonValidHitsGLB"  + etaTag].fill(probeMuon.glbMuonValidHits, probeMuTk, weight, nVtx);
		  m_plots[ID]["TrkStaChi2"  + etaTag].fill(probeMuon.trkStaChi2, probeMuTk, weight, nVtx);       
		  m_plots[ID]["TrkKink"  + etaTag].fill(probeMuon.trkKink, probeMuTk, weight, nVtx);          
		}

	      if (probeMuon.isTrackerArb && probeMuon.isGlobal)    
		{
		  m_plots[ID]["SegmentComp"  + etaTag].fill(probeMuon.muSegmComp, probeMuTk, weight, nVtx);      
		  m_plots[ID]["NHitsTRK"  + etaTag].fill(probeMuon.nHitsTracker, probeMuTk, weight, nVtx);
		  m_plots[ID]["Chi2TRK"  + etaTag].fill(probeMuon.trkNormChi2, probeMuTk, weight, nVtx);
		  m_plots[ID]["PixelHitsTRK"  + etaTag].fill(probeMuon.trkPixelValidHits, probeMuTk, weight, nVtx);
		  m_plots[ID]["PixelLayersTRK"  + etaTag].fill(probeMuon.trkPixelLayersWithMeas, probeMuTk, weight, nVtx); 
		  m_plots[ID]["TrackerLayersTRK"  + etaTag].fill(probeMuon.trkTrackerLayersWithMeas, probeMuTk, weight, nVtx); 
		  m_plots[ID]["HitFractionTRK"  + etaTag].fill(probeMuon.trkValidHitFrac, probeMuTk, weight, nVtx);   
		  m_plots[ID]["Dxy"  + etaTag].fill(probeMuon.dxy, probeMuTk, weight, nVtx);
		  m_plots[ID]["Dz"  + etaTag].fill(probeMuon.dz, probeMuTk, weight, nVtx);			
		}

	      if (probeMuon.isStandAlone || probeMuon.isGlobal) 
		{
		  m_plots[ID]["NHitsSTA"  + etaTag].fill(probeMuon.nHitsStandAlone, probeMuTk, weight, nVtx);
		}

	      if (probeMuon.isTrackerArb)
		{
		  m_plots[ID]["NMatchedStation"  + etaTag].fill(probeMuon.trkMuonMatchedStations, probeMuTk, weight, nVtx);
		}


	      if (probeMuon.isGlobal)
		{
		  Float_t qOverPtTrk = probeMuon.charge_tracker / probeMuon.pt_tracker;
		  Float_t qOverPtSta = probeMuon.charge_standalone / probeMuon.pt_standalone;
		  Float_t qOverPtGlb = probeMuon.charge_global / probeMuon.pt_global;

		  m_plots[ID]["qOverPtTrkSta" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
		  // m_plots[ID]["qOverPtTrkStaOverPt" + etaTag].fill((qOverPtSta - qOverPtTrk)/qOverPtTrk, probeMuTk, weight, nVtx);  
		  m_plots[ID]["qOverPtTrkGlb" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  

		  // if (probeMuon.charge_tracker > 0 )
		  //   {
		  //     m_plots[ID]["qOverPtTrkGlbPlus" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  
		  //     m_plots[ID]["qOverPtTrkStaPlus" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
		  //   }
		  if (probeMuon.pt_tracker > 200)
		    {
		      m_plots[ID]["qOverPtTrkSta200" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
		      m_plots[ID]["qOverPtTrkGlb200" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  
		      // if (probeMuon.charge_tracker > 0 )
		      // 	{
		      // 	  m_plots[ID]["qOverPtTrkSta200Plus" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
		      // 	  m_plots[ID]["qOverPtTrkGlb200Plus" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  
		      // 	}
		    }
		}

	      if( hasGoodId(probeMuon,"TIGHT") && probeMuon.isoPflow04 < 0.25 &&
		  ( (fabs(probeMuon.eta) < 1.2  && probeMuon.nHitsStandAlone > 20) || 
		    (fabs(probeMuon.eta) >= 1.2 && probeMuon.nHitsStandAlone > 11)
		    ))
		{
		  m_plots[TIME]["TightMuonTime" + etaTag].fill(probeMuon.muonTime, probeMuTk, weight, nVtx);
		}

	      for (auto & probe_ID : m_tnpConfig.probe_IDs)
	      	{
	      	  TString IDTag = "_" + probe_ID;
		  
	      	  if(hasGoodId(probeMuon,probe_ID)) 
	      	    {	 
	      	      m_plots[KIN]["ProbePt" + etaTag + IDTag].fill(probeMuTk.Pt(), probeMuTk, weight, nVtx);
	      	      m_plots[KIN]["ProbeEta" + etaTag + IDTag].fill(probeMuTk.Eta(), probeMuTk, weight, nVtx);
	      	      m_plots[KIN]["ProbePhi" + etaTag + IDTag].fill(probeMuTk.Phi(), probeMuTk, weight, nVtx);
		      
	      	      // Fill isolation plots for muons passign a given identification (programmable from cfg)
	      	      m_plots[ISO]["PhotonIso" + etaTag + IDTag].fill(probeMuon.photonIso, probeMuTk, weight, nVtx);
	      	      m_plots[ISO]["ChHadIso" + etaTag + IDTag].fill(probeMuon.chargedHadronIso, probeMuTk, weight, nVtx);
	      	      m_plots[ISO]["ChHadIsoPU" + etaTag + IDTag].fill(probeMuon.chargedHadronIsoPU, probeMuTk, weight, nVtx);
	      	      m_plots[ISO]["NeutralIso" + etaTag + IDTag].fill(probeMuon.neutralHadronIso, probeMuTk, weight, nVtx);
	      	      m_plots[ISO]["DBetaRelIso" + etaTag + IDTag].fill(probeMuon.isoPflow04, probeMuTk, weight, nVtx);
	      	    }
		  }
	    }
	}
    }
}


//Functions

void muon_pog::parseConfig(const std::string configFile, muon_pog::TagAndProbeConfig & tpConfig,
			   std::vector<muon_pog::SampleConfig> & sampleConfigs)
{

  boost::property_tree::ptree pt;
  
  try
    {
      boost::property_tree::ini_parser::read_ini(configFile, pt);
    }
  catch (boost::property_tree::ini_parser::ini_parser_error iniParseErr)
    {
      std::cout << "[TagAndProbeConfig] Can't open : " << iniParseErr.filename()
		<< "\n\tin line : " << iniParseErr.line()
		<< "\n\thas error :" << iniParseErr.message()
		<< std::endl;
      throw std::runtime_error("Bad INI parsing");
    }

  for( auto vt : pt )
    {
      if (vt.first.find("TagAndProbe") != std::string::npos)
	tpConfig = muon_pog::TagAndProbeConfig(vt);
      else
	sampleConfigs.push_back(muon_pog::SampleConfig(vt));
    }
}

void muon_pog::comparisonPlots(std::vector<muon_pog::Plotter> & plotters,
			       TFile *outFile, TString &  outputDir)
{

  outFile->cd("/");
  outFile->mkdir("comparison");
  
  outFile->mkdir("comparison/control");
  outFile->mkdir("comparison/timing");
  outFile->mkdir("comparison/isolation");
  outFile->mkdir("comparison/id_variables");
  outFile->mkdir("comparison/kinematical_variables");
  outFile->mkdir("comparison/momentum_variables");
  outFile->mkdir("comparison/efficiencies");
  
  outFile->cd("comparison");
  std::ofstream integrals(outputDir + "/Result.txt");
 
  std::vector<std::pair<Plotter::HistoType,TString> > plotTypesAndNames;
  for (auto & plotPairType : plotters.at(0).m_plots)
    {
      for (auto & plotPairName : plotPairType.second)
	{
	  plotTypesAndNames.push_back(std::make_pair(plotPairType.first, plotPairName.first));
	}
    }
 
  for (auto & plotTypeAndName : plotTypesAndNames)
    {

      TString outputDirMap[7] {"/comparison/kinematical_variables", "/comparison/momentum_variables", "/comparison/id_variables", 
	  "/comparison/isolation", "/comparison/timing", "/comparison/control", "/comparison/efficiencies"};
      
      Plotter::HistoType plotType = plotTypeAndName.first;
      TString plotName = plotTypeAndName.second;
      outFile->cd(outputDirMap[plotType]);

      std::stringstream sPlotName(plotName.Data());
      std::string plotTitle;
      std::getline(sPlotName, plotTitle, '_');

      std::vector<TH1 *>::size_type nPlots = plotters.at(0).m_plots[plotType][plotName].plots().size();

      for (std::vector<TH1 *>::size_type iPlot=0; iPlot < nPlots; ++ iPlot)
	{
	  TLegend *leg = new TLegend(0.65,0.7,0.95,0.95);
	  leg->SetBorderSize(0);
	  leg->SetLineWidth(0);
	  leg->SetFillColor(0);
	  leg->SetFillStyle(0);
	  
	  THStack hMc(plotName,"");
	  TProfile * pMc = 0;
	  TH1 *proj = 0;
	  Float_t errors[1000] = {0.};
	  TH1 * hData = 0;

	  int colorMap[5] {kGreen+1, kAzure+7, kGray+1, kRed, kOrange};
	  TString folderMap[5] {"", "/VsEta", "/VsPhi", "/VsPt", "/VsPV"};
	 
	  int iColor = 0;
	  
	  float integralData = 0;
	  float integralMC   = 0;
	  
	  for (auto & plotter : plotters)
	    {
	      if(plotter.m_sampleConfig.sampleName.Contains("Data"))
		{
		  if(plotName.Contains("muonTime"))
		    {
		      TH1* plot = plotter.m_plots[plotType][plotName].plots().at(0);
		      integralData = plot->Integral(plot->FindBin(-2.5),plot->FindBin(2.5)); 
		    }
		  else
		    integralData = plotter.m_plots[Plotter::CONT]["99-invMassInRange"].plots().at(0)->Integral();
		}
	      else
		{
		  if(plotName.Contains("muonTime"))
		    {
		      TH1* plot = plotter.m_plots[plotType][plotName].plots().at(0);
		      integralMC += (plot->Integral(plot->FindBin(-2.5),plot->FindBin(2.5)) *
				     plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents);
		    }
		  else
		    {
		      integralMC += (plotter.m_plots[Plotter::CONT]["99-invMassInRange"].plots().at(0)->Integral() *
				     plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents);
		      
		    }
		}
	    }
	  
	  for (auto & plotter : plotters)
	    {        
	      TH1* plot = plotter.m_plots[plotType][plotName].plots().at(iPlot);
	      if(plotter.m_sampleConfig.sampleName.Contains("Data"))
		{
		  hData = plot;
		  //FP Add OF and UF for iso variables only
		  if (!plot->IsA()->InheritsFrom("TProfile"))
		    {
		      hData->Sumw2();
		      if(!plotName.Contains("mass")){
			addOverFlow(*hData);		     
			addUnderFlow(*hData);
		      }
		    }
		  
		  Double_t stat = hData->GetEntries();	
		  leg->AddEntry(hData,Form(plotter.m_sampleConfig.sampleName+" [%9.0f]",stat),"LP");  
		}
	      else
		{
		  plot->SetFillColor(colorMap[iColor]);
		  plot->SetMarkerStyle(0);
		 		 
		  float scale = plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents;
		  
		  if (!plot->IsA()->InheritsFrom("TProfile"))
		    {
		      //FP Add OF and UF for iso variables only
		      if(!plotName.Contains("mass")){
			addOverFlow(*plot);
			//if(plotName.Contains("Iso") || plotName.Contains("Dxy") || plotName.Contains("Dz"))
			addUnderFlow(*plot);
		      }
		      
		      float scaleToData = integralData/integralMC;

		      if(!plotName.Contains("runNumber"))		      
			plot->Scale(scale*scaleToData);
	     
		      hMc.Add(plot);
		      leg->AddEntry(plot, plotter.m_sampleConfig.sampleName, "LF"); 
		      iColor++;
		    }
		  else
		    {
		      if (!pMc)
		       	{
			  pMc = (TProfile*)(plot->Clone("_pClone"));
			  pMc->Clear();
			  pMc->SetFillColor(0);
			  pMc->SetMarkerStyle(26);
			  pMc->SetMarkerColor(colorMap[iColor]);
			  
			  pMc->Add(plot,plotName.Contains("runNumber") ? 1 : scale);
			  // leg->AddEntry(pMc, "Weighted sum of MCs", "LP"); 
			  leg->AddEntry(pMc, plotter.m_sampleConfig.sampleName, "LP");
			}
		      else
		  	{
			  pMc->Add(plot,plotName.Contains("runNumber") ? 1 :scale);
			}
		    
		      iColor++;
		    }
		}
	    }
	  
	  TString histoName = hData->GetName();
	  
	  TCanvas *canvas = new TCanvas("c"+histoName, "c"+histoName, 500, 500);
	  canvas->cd();
	  
	  if(pMc)
	    {
	      if(plotName.Contains("goodMuMass"))
		{
		  if(histoName.Contains("VsPhi"))
		    hData->GetYaxis()->SetRangeUser(89.,93.);
		  else
		    hData->GetYaxis()->SetRangeUser(90.5,91.5);
		}
	      else
		{
		  Double_t max = hData->GetMaximum() > pMc->GetMaximum() ? hData->GetMaximum() : pMc->GetMaximum();
		  Double_t min = hData->GetMinimum() < pMc->GetMinimum() ? hData->GetMinimum() : pMc->GetMinimum();

		  Float_t range = max - min;
		  hData->GetYaxis()->SetRangeUser(min-range*0.2, max+range*0.2);
		}
	      
	      hData->Draw();
	      pMc->Draw("same");
	    }
	  else
	    {
	      canvas->SetLogy(1);
	      hData->Draw();
	      hMc.Draw("samehist");
	      hData->Draw("same");
	    }  
	  
	  leg->Draw();
	  
	  canvas->Update();
	  canvas->Write();
	  system("mkdir -p " + outputDir + outputDirMap[plotType] + "/" + plotTitle + folderMap[iPlot] +  "/no_ratio/");
	  canvas->SaveAs(outputDir + outputDirMap[plotType] + "/" + plotTitle + folderMap[iPlot] + "/no_ratio/c" + histoName + ".png");
	  
	  //Canvas with ratio plots
	  TCanvas *ratioCanvas = new TCanvas("rc"+histoName, "rc"+histoName, 500, 700);
	  ratioCanvas->Divide(1,2);
	  ratioCanvas->cd(1);
	  TPad *plotPad = (TPad*)ratioCanvas->GetPad(1);
	  plotPad->SetPad(0.,0.2,1.,1.);
	  
	  hData->GetXaxis()->SetLabelSize(0.);
	  hData->GetXaxis()->SetTitleSize(0.);
	  
	  if(pMc)
	    {
	      if(plotName.Contains("goodMuMass"))
		{
		  if(histoName.Contains("VsPhi"))
		    hData->GetYaxis()->SetRangeUser(89.,93.);
		  else
		    hData->GetYaxis()->SetRangeUser(90.5,91.5);
		}
	      else
		{
		  Double_t max = hData->GetMaximum() > pMc->GetMaximum() ? hData->GetMaximum() : pMc->GetMaximum();
		  Double_t min = hData->GetMinimum() < pMc->GetMinimum() ? hData->GetMinimum() : pMc->GetMinimum();

		  Double_t range = max-min;
		  hData->GetYaxis()->SetRangeUser(min-range*0.2, max+range*0.2);
		}
	      
	      hData->Draw();
	      pMc->Draw("same");
	    }
	  else
	    {
	      plotPad->SetLogy(1);
	      hData->Draw();
	      hMc.Draw("samehist");
	      hData->Draw("same");
	    }
	  
	  leg->Draw();
	  
	  // //Get Integral before/after ID cuts
	  // std::vector<TString> etaBinMin = {"0.0","2.1"};
	  // std::vector<TString> etaBinMax = {"2.4","2.4"};
	  
	  // for (int i = 0; i < 2; ++i){
	  //   TString etaTag = "_fEtaMin" + fEtaBinMin[i] + "_fEtaMax" + fEtaBinMax[i];	    

	  //   if(plotName == "HitFractionTRK" + etaTag && !pMc)
	  //     {
	  // 	//the cut is > 0.8
	  // 	double x1 = hData->FindBin(0.8);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_HitFractionTRK_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
	  // 	double Data_HitFractionTRK_b = hData->Integral(1,x1)/hData->Integral(1,x2);
	  // 	double MC_HitFractionTRK_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_HitFractionTRK_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: HitFractionTRK_a = " << Data_HitFractionTRK_a << "  HitFractionTRK_b = " << Data_HitFractionTRK_b << "  --> " << Data_HitFractionTRK_a+Data_HitFractionTRK_b << std::endl;
	  // 	integrals << "MC  : HitFractionTRK_a = " << MC_HitFractionTRK_a << "  HitFractionTRK_b = " << MC_HitFractionTRK_b << "  --> " << MC_HitFractionTRK_a+MC_HitFractionTRK_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	  //   if(plotName == "Chi2GLB" + etaTag && !pMc) // chi2<3
	  //     {
	  // 	double x1 = hData->FindBin(3.);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_Chi2GLB_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
	  // 	double Data_Chi2GLB_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_Chi2GLB_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_Chi2GLB_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: Chi2GLB_a = " << Data_Chi2GLB_a << "  Chi2GLB_b = " << Data_Chi2GLB_b << "  --> " << Data_Chi2GLB_a+Data_Chi2GLB_b << std::endl;
	  // 	integrals << "MC  : Chi2GLB_a = " << MC_Chi2GLB_a << "  Chi2GLB_b = " << MC_Chi2GLB_b << "  --> " << MC_Chi2GLB_a+MC_Chi2GLB_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	    
	  //   if(plotName == "TrkStaChi2" + etaTag && !pMc)
	  //     {
	  // 	double x1 = hData->FindBin(12.);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_TrkStaChi2_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
	  // 	double Data_TrkStaChi2_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_TrkStaChi2_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_TrkStaChi2_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	integrals << std::endl;
	  // 	integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: TrkStaChi2_a = " << Data_TrkStaChi2_a << "  TrkStaChi2_b = " << Data_TrkStaChi2_b << "  --> " << Data_TrkStaChi2_a+Data_TrkStaChi2_b << std::endl;
	  // 	integrals << "MC  : TrkStaChi2_a = " << MC_TrkStaChi2_a << "  TrkStaChi2_b = " << MC_TrkStaChi2_b << "  --> " << MC_TrkStaChi2_a+MC_TrkStaChi2_b << std::endl;
	  // 	integrals << std::endl;
		
	  //     }
	    
	  //   if(plotName == "TrkKink" + etaTag && !pMc)
	  //     {
	  // 	double x1 = hData->FindBin(20.);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_TrkKink_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
	  // 	double Data_TrkKink_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_TrkKink_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_TrkKink_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: TrkKink_a = " << Data_TrkKink_a << "  TrkKink_b = " << Data_TrkKink_b << "  --> " << Data_TrkKink_a+Data_TrkKink_b << std::endl;
	  // 	integrals << "MC  : TrkKink_a = " << MC_TrkKink_a << "  TrkKink_b = " << MC_TrkKink_b << "  --> " << MC_TrkKink_a+MC_TrkKink_b << std::endl;
	  // 	integrals << std::endl;
		
	  //     }
	    
	  //   if(plotName == "SegmentComp" + etaTag && !pMc)
	  //     {
	  // 	double x0  = hData->FindBin(.303);
	  // 	double x1  = hData->FindBin(.451);
	  // 	double x2  = hData->GetNbinsX();
		
	  // 	double Data_SegmentCompLoose_a = hData->Integral(x0+1,x2)/hData->Integral(1,x2);
	  // 	double Data_SegmentCompTight_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
	  // 	double Data_SegmentCompLoose_b = hData->Integral(1,x0)/hData->Integral(1,x2);
	  // 	double Data_SegmentCompTight_b = hData->Integral(1,x1)/hData->Integral(1,x2);
	  // 	double MC_SegmentCompLoose_a = ((TH1*)hMc.GetStack()->Last())->Integral(x0+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_SegmentCompTight_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_SegmentCompLoose_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x0)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_SegmentCompTight_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: SegmentCompLoose_a = " << Data_SegmentCompLoose_a << "  SegmentCompLoose_b = " << Data_SegmentCompLoose_b << "  --> " << Data_SegmentCompLoose_a+Data_SegmentCompLoose_b << std::endl;
	  // 	integrals << "MC  : SegmentCompLoose_a = " << MC_SegmentCompLoose_a << "  SegmentCompLoose_b = " << MC_SegmentCompLoose_b << "  --> " << MC_SegmentCompLoose_a+MC_SegmentCompLoose_b << std::endl;
	  // 	integrals << "Data: SegmentCompTight_a = " << Data_SegmentCompTight_a << "  SegmentCompTight_b = " << Data_SegmentCompTight_b << "  --> " << Data_SegmentCompLoose_a+Data_SegmentCompLoose_b << std::endl;
	  // 	integrals << "MC  : SegmentCompTight_a = " << MC_SegmentCompTight_a << "  SegmentCompTight_b = " << MC_SegmentCompTight_b << "  --> " << MC_SegmentCompLoose_a+MC_SegmentCompLoose_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	  //   //Tight ID
	  //   if(plotName == "Chi2GLB" + etaTag && !pMc)
	  //     {
	  // 	double x1 = hData->FindBin(10.);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_Chi2GLB_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
	  // 	double Data_Chi2GLB_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_Chi2GLB_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_Chi2GLB_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: Chi2GLB_a = " << Data_Chi2GLB_a << "  Chi2GLB_b = " << Data_Chi2GLB_b << "  --> " << Data_Chi2GLB_a+Data_Chi2GLB_b << std::endl;
	  // 	integrals << "MC  : Chi2GLB_a = " << MC_Chi2GLB_a << "  Chi2GLB_b = " << MC_Chi2GLB_b << "  --> " << MC_Chi2GLB_a+MC_Chi2GLB_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	  //   if(plotName == "NMuonValidHitsGLB" + etaTag && !pMc)
	  //     {
	  // 	double x1 = hData->FindBin(0.);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_NMuonValidHitsGLB_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
	  // 	double Data_NMuonValidHitsGLB_b = hData->Integral(1,x1)/hData->Integral(1,x2);
	  // 	double MC_NMuonValidHitsGLB_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_NMuonValidHitsGLB_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: NMuonValidHitsGLB_a = " << Data_NMuonValidHitsGLB_a << "  NMuonValidHitsGLB_b = " << Data_NMuonValidHitsGLB_b << "  --> " << Data_NMuonValidHitsGLB_a+Data_NMuonValidHitsGLB_b << std::endl;
	  // 	integrals << "MC  : NMuonValidHitsGLB_a = " << MC_NMuonValidHitsGLB_a << "  NMuonValidHitsGLB_b = " << MC_NMuonValidHitsGLB_b << "  --> " << MC_NMuonValidHitsGLB_a+MC_NMuonValidHitsGLB_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	  //   if(plotName == "NMatchedStation" + etaTag && !pMc)
	  //     {
	  // 	double x1 = hData->FindBin(1.);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_NMatchedStation_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
	  // 	double Data_NMatchedStation_b = hData->Integral(1,x1)/hData->Integral(1,x2);
	  // 	double MC_NMatchedStation_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_NMatchedStation_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: NMatchedStation_a = " << Data_NMatchedStation_a << "  NMatchedStation_b = " << Data_NMatchedStation_b << "  --> " << Data_NMatchedStation_a+Data_NMatchedStation_b << std::endl;
	  // 	integrals << "MC  : NMatchedStation_a = " << MC_NMatchedStation_a << "  NMatchedStation_b = " << MC_NMatchedStation_b << "  --> " << MC_NMatchedStation_a+MC_NMatchedStation_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	  //   if(plotName == "Dxy" + etaTag && !pMc)
	  //     {
	  // 	double mx1 = hData->FindBin(-0.2);
	  // 	double x1  = hData->FindBin(.2);
	  // 	double x2  = hData->GetNbinsX();
	  // 	double Data_Dxy_a = (hData->Integral(x1,x2) + hData->Integral(1,mx1))/hData->Integral(1,x2);  
	  // 	double Data_Dxy_b = hData->Integral(mx1+1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_Dxy_a = (((TH1*)hMc.GetStack()->Last())->Integral(x1,x2) + ((TH1*)hMc.GetStack()->Last())->Integral(1,mx1))/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_Dxy_b = ((TH1*)hMc.GetStack()->Last())->Integral(mx1+1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: Dxy_a = " << Data_Dxy_a << "  Dxy_b = " << Data_Dxy_b << "  --> " << Data_Dxy_a+Data_Dxy_b << std::endl;
	  // 	integrals << "MC  : Dxy_a = " << MC_Dxy_a << "  Dxy_b = " << MC_Dxy_b << "  --> " << MC_Dxy_a+MC_Dxy_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	  //   if(plotName == "Dz" + etaTag && !pMc)
	  //     {
	  // 	double mx1 = hData->FindBin(-0.5);
	  // 	double x1 = hData->FindBin(.5);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_Dz_a = (hData->Integral(x1,x2) + hData->Integral(1,mx1))/hData->Integral(1,x2); //the cut is above th 
	  // 	double Data_Dz_b = hData->Integral(mx1+1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_Dz_a = (((TH1*)hMc.GetStack()->Last())->Integral(x1,x2) + ((TH1*)hMc.GetStack()->Last())->Integral(1,mx1))/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_Dz_b = ((TH1*)hMc.GetStack()->Last())->Integral(mx1+1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: Dz_a = " << Data_Dz_a << "  Dz_b = " << Data_Dz_b << "  --> " << Data_Dz_a+Data_Dz_b << std::endl;
	  // 	integrals << "MC  : Dz_a = " << MC_Dz_a << "  Dz_b = " << MC_Dz_b << "  --> " << MC_Dz_a+MC_Dz_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	    
	  //   if(plotName == "PixelHitsTRK" + etaTag && !pMc)
	  //     {
	  // 	double x1 = hData->FindBin(1);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_PixelHitsTRK_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
	  // 	double Data_PixelHitsTRK_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_PixelHitsTRK_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_PixelHitsTRK_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: PixelHitsTRK_a = " << Data_PixelHitsTRK_a << "  PixelHitsTRK_b = " << Data_PixelHitsTRK_b << "  --> " << Data_PixelHitsTRK_a+Data_PixelHitsTRK_b << std::endl;
	  // 	integrals << "MC  : PixelHitsTRK_a = " << MC_PixelHitsTRK_a << "  PixelHitsTRK_b = " << MC_PixelHitsTRK_b << "  --> " << MC_PixelHitsTRK_a+MC_PixelHitsTRK_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	    
	  //   if(plotName == "TrackerLayersTRK" + etaTag && !pMc)
	  //     {
	  // 	double x1 = hData->FindBin(6);
	  // 	double x2 = hData->GetNbinsX();
	  // 	double Data_TrackerLayersTRK_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
	  // 	double Data_TrackerLayersTRK_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
	  // 	double MC_TrackerLayersTRK_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
	  // 	double MC_TrackerLayersTRK_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
	  // 	integrals << std::endl;
	  // 	integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
	  // 	integrals << "Data: TrackerLayersTRK_a = " << Data_TrackerLayersTRK_a << "  TrackerLayersTRK_b = " << Data_TrackerLayersTRK_b << "  --> " << Data_TrackerLayersTRK_a+Data_TrackerLayersTRK_b << std::endl;
	  // 	integrals << "MC  : TrackerLayersTRK_a = " << MC_TrackerLayersTRK_a << "  TrackerLayersTRK_b = " << MC_TrackerLayersTRK_b << "  --> " << MC_TrackerLayersTRK_a+MC_TrackerLayersTRK_b << std::endl;
	  // 	integrals << std::endl;
	  //     }
	  // }
	  
	  
	  ratioCanvas->cd(2);
	  TPad *ratioPad = (TPad*)ratioCanvas->GetPad(2);
	  ratioPad->SetPad(0.,0.,1.,0.31);
	  ratioPad->SetFillStyle(4000);
	  ratioPad->SetBottomMargin(0.2);	  

	  TH1 *hRatio = (TH1*)hData->Clone("_dataClone");
	  hRatio->SetTitle(" ");

	  hRatio->GetXaxis()->SetLabelSize(0.1);
	  hRatio->GetXaxis()->SetTitleSize(0.1);
	  hRatio->GetXaxis()->SetTitleOffset(.85);

	  hRatio->GetYaxis()->SetLabelSize(0.07);
	  hRatio->GetYaxis()->SetTitleSize(0.1);
	  hRatio->GetYaxis()->SetTitleOffset(.6);
	  hRatio->GetYaxis()->SetTitle("Data/MC");
	  hRatio->GetYaxis()->SetRangeUser(0.,2.);
	  
	  if(pMc)
	    hRatio->Divide((TH1*)pMc);
	  else
	    hRatio->Divide((TH1*)hMc.GetStack()->Last());
	  
	  hRatio->Draw();

	  Double_t Xmax = hRatio->GetXaxis()->GetXmax();
	  Double_t Xmin = hRatio->GetXaxis()->GetXmin();
	 
	  TLine *l = new TLine(Xmin,1,Xmax,1);
	  l->SetLineColor(1); 
	  l->Draw("same"); 
      
	  if(plotName == "03-nVertices" && !pMc)
	    {
	      Int_t nbins = hRatio->GetNbinsX();
	 
	      std::cout << "[comparisonPlots] : reweight w.r.t. data sample: " << std::endl;
	      std::cout << "pu_weights =";
	      for (Int_t i=1;i<=nbins;i++)
		{
		  std::cout << hRatio->GetBinContent(i) << ",";
		}
	      std::cout << std::endl;
	    }
	        
	  ratioCanvas->Write();
	  ratioCanvas->SaveAs(outputDir+ outputDirMap[plotType] + "/" + plotTitle + folderMap[iPlot] + "/rc" + histoName + ".png");
      
	  delete canvas;
	  delete ratioCanvas;
	}
    }
}

void muon_pog::copyPhp(const TString &  outputDir)
{

  system("cp index.php " + outputDir);
  
  boost::filesystem::directory_iterator dirIt(outputDir.Data());
  boost::filesystem::directory_iterator dirEnd;
  for (;dirIt != dirEnd; ++ dirIt)
  {
    if (boost::filesystem::is_directory(dirIt->status()))
      copyPhp(TString(dirIt->path().string()));
  }
  
}
