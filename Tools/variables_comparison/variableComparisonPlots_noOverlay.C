#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "../src/MuonPogTree.h"
#include "tdrstyle.C"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <map>

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
    Float_t eventi;
    Int_t applyReweighting;

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

    std::string tag_ID;  
    Float_t     tag_isoCut;
    Float_t     tag_hltDrCut;
    std::string tag_hltFilter;
    
    std::string muon_trackType; // applies to both, tag and probe     
  
    std::string probe_ID;  
    std::vector<TString> probe_fEtaMin;
    std::vector<TString> probe_fEtaMax;
    
    std::string hlt_path; 
   
    TagAndProbeConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    TagAndProbeConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~TagAndProbeConfig() {};
    
  private:
    std::vector<TString> toArray(const std::string& entries); 
  
  };

  class Plotter {

  public :
    
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, int nVtx, float weight);

    std::map<TString,TH1 *> m_plots;
    TagAndProbeConfig m_tnpConfig;
    SampleConfig m_sampleConfig;

  private :

    bool hasGoodId(const muon_pog::Muon & muon,
		   std::string leg);
    bool hasFilterMatch(const muon_pog::Muon & muon,
			const muon_pog::HLT  & hlt);
    Int_t chargeFromTrk(const muon_pog::Muon & muon);
    TLorentzVector muonTk(const muon_pog::Muon & muon);    
    
  };

}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file
// 1. comparisonPlot : make a plot overlayng data and MC for a given plot
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, TagAndProbeConfig & tpConfig,
		   std::vector<SampleConfig> & sampleConfigs);
  
  void comparisonPlot(TFile *outFile, TString plotName,
		      std::vector<Plotter> & plotters);

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
		<< argv[0] << " PAT_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR\n";
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
 
  for (auto plotter : plotters)
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
      int nEntries = tree->GetEntriesFast();
      std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

      int nFilteredEvents = 0;

      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
	//   for (Long64_t iEvent=0; iEvent<10000; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;
	  
	  if(iEvent%10000 == 0) printf("[%s] Processing event %8d/%8d [%4.1f%]\n", argv[0], iEvent, nEntries, float(iEvent)/float(nEntries)*100); 

	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;
	  
	  //weight fo Data to Data or DY comparison 
	  //float PUweight[60] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

	  // Weights for Data to DY comparison 
	  float PUweight_DY[60] = {0., 4.20351, 2.32244, 0.818425, 1.09792, 1.15235, 1.25735, 1.17727, 1.34312, 1.40308, 1.37595, 1.34705, 1.24496, 1.2293, 1.04211, 0.984, 0.839735, 0.69296, 0.65199, 0.508954, 0.496283, 0.409652, 0.362298, 0.314055, 0.290668, 0.17501, 0.172355, 0.124814, 0.167055, 0, 0.147691, 0.162636, 0, 0, 0, 0, 0, 0, 0, 0., 0., 0, 0., 0., 0., 0., 0., 0., 0, 0., 0, 0., 0., 0., 0., 0., 0., 0., 0., 0.}; 
	  
	  // Weights for Data to Data comparison 
	  float PUweight_Data50[60] = {0., 1.07875, 3.33432, 5.6095, 8.52213, 4.48899, 4.53211, 3.44043, 2.71373, 2.28532, 1.92331, 1.6078, 1.26329, 1.11807, 0.883009, 0.808497, 0.660469, 0.517963, 0.484396, 0.408176, 0.391339, 0.312365, 0.281628, 0.256951, 0.212307, 0.164809, 0.137559, 0.128981, 0.179792, 0, 0.134844, 0.154107, 0, 0, 0, 0, 0, 0., 0., 0., 0., 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

	  float *PUweight = std::string(plotter.m_sampleConfig.sampleName.Data()).find("Data") == std::string::npos ? PUweight_DY : PUweight_Data50; 

	  //if(applyReweighting || std::string(plotter.m_sampleConfig.sampleName.Data()).find("Data") == std::string::npos){
	  if(plotter.m_sampleConfig.applyReweighting==1){
	    if(ev->nVtx < 60)
	      weight *= PUweight[ev->nVtx];
	    else 
	      weight *=0;
	  }
   
	  plotter.fill(ev->muons, ev->hlt, ev->nVtx, weight);
	  
	}
      
      delete ev;
      delete evBranch;
      
      inputFile->Close();
      
    }

  outputFile->cd("/");
  /*
  outputFile->mkdir("comparison");
  outputFile->cd("comparison");

  muon_pog::comparisonPlot(outputFile,"nVertices",plotters);

  //timing
  muon_pog::comparisonPlot(outputFile,"STAmuonTime",plotters);
  muon_pog::comparisonPlot(outputFile,"STAmuonTimeBarrel",plotters);
  muon_pog::comparisonPlot(outputFile,"STAmuonTimeEndcap",plotters);

  muon_pog::comparisonPlot(outputFile,"UnbSTAmuonTime",plotters);
  muon_pog::comparisonPlot(outputFile,"UnbSTAmuonTimeBarrel",plotters);
  muon_pog::comparisonPlot(outputFile,"UnbSTAmuonTimeEndcap",plotters);

  // muon_pog::comparisonPlot(outputFile,"GLBmuonTime",plotters);
  // muon_pog::comparisonPlot(outputFile,"GLBmuonTimeBarrel",plotters);
  // muon_pog::comparisonPlot(outputFile,"GLBmuonTimeEndcap",plotters);

  /////
  muon_pog::comparisonPlot(outputFile,"invMass",plotters);
  muon_pog::comparisonPlot(outputFile,"dilepPt",plotters);

  std::vector<TString>::const_iterator fEtaMinIt  = tnpConfig.probe_fEtaMin.begin();
  std::vector<TString>::const_iterator fEtaMinEnd = tnpConfig.probe_fEtaMin.end();

  std::vector<TString>::const_iterator fEtaMaxIt  = tnpConfig.probe_fEtaMax.begin();
  std::vector<TString>::const_iterator fEtaMaxEnd = tnpConfig.probe_fEtaMax.end();
  
  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
    {
      TString etaTag = "_fEtaMin" + (*fEtaMinIt) + "_fEtaMax" + (*fEtaMaxIt);

      //Id variables
      muon_pog::comparisonPlot(outputFile,"nHitsGLB" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"nHitsTRK" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"nHitsSTA" + etaTag,plotters);

      muon_pog::comparisonPlot(outputFile,"nChi2GLB" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"nChi2TRK" + etaTag,plotters);

      muon_pog::comparisonPlot(outputFile,"matchedStation" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"muonValidHitsGLB" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"PixelHitsTRK" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"PixelLayersTRK" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"TrackerLayersTRK" + etaTag,plotters);

      muon_pog::comparisonPlot(outputFile,"HitFractionTRK" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"TrkStaChi2" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"TrkKink" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"segmentComp" + etaTag,plotters);

      muon_pog::comparisonPlot(outputFile,"probeDxy" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"probeDz" + etaTag,plotters);

      muon_pog::comparisonPlot(outputFile,"probePt" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"probeEta" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"probePhi" + etaTag,plotters);

      muon_pog::comparisonPlot(outputFile,"chHadIso" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"photonIso" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"neutralIso" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"dBetaRelIso" + etaTag,plotters);
    }
  */
  
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

      probe_ID     = vt.second.get<std::string>("probe_muonID");
      probe_fEtaMin = toArray(vt.second.get<std::string>("probe_fEtaMin"));
      probe_fEtaMax = toArray(vt.second.get<std::string>("probe_fEtaMax"));
      

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
      eventi = vt.second.get<Float_t>("eventi");
      applyReweighting = vt.second.get<Int_t>("applyReweighting");
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

void muon_pog::Plotter::book(TFile *outFile)
{

  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);

  std::vector<TString>::const_iterator fEtaMinIt  = m_tnpConfig.probe_fEtaMin.begin();
  std::vector<TString>::const_iterator fEtaMinEnd = m_tnpConfig.probe_fEtaMin.end();

  std::vector<TString>::const_iterator fEtaMaxIt  = m_tnpConfig.probe_fEtaMax.begin();
  std::vector<TString>::const_iterator fEtaMaxEnd = m_tnpConfig.probe_fEtaMax.end();

  m_plots["STAmuonTime" + sampleTag]        = new TH1F("STAmuonTime_" + sampleTag ,"STAmuonTime; time (ns) ; # entries", 400, -200., 200.);
  m_plots["STAmuonTimeBarrel" + sampleTag]  = new TH1F("STAmuonTimeBarrel_" + sampleTag ,"STAmuonTimeBarrel; time (ns) ; # entries", 400, -200., 200.);
  m_plots["STAmuonTimeEndcap" + sampleTag]  = new TH1F("STAmuonTimeEndcap_" + sampleTag ,"STAmuonTimeEndcap; time (ns) ; # entries", 400, -200., 200.);

  m_plots["UnbSTAmuonTime" + sampleTag]        = new TH1F("UnbSTAmuonTime_" + sampleTag ,"UnbSTAmuonTime; time (ns) ; # entries", 400, -200., 200.);
  m_plots["UnbSTAmuonTimeBarrel" + sampleTag]  = new TH1F("UnbSTAmuonTimeBarrel_" + sampleTag ,"UnbSTAmuonTimeBarrel; time (ns) ; # entries", 400, -200., 200.);
  m_plots["UnbSTAmuonTimeEndcap" + sampleTag]  = new TH1F("UnbSTAmuonTimeEndcap_" + sampleTag ,"UnbSTAmuonTimeEndcap; time (ns) ; # entries", 400, -200., 200.);

  // m_plots["GLBmuonTime"]        = new TH1F("GLBmuonTime_" + sampleTag ,"GLBmuonTime; time (ns) ; # entries", 400, -200., 200.);
  // m_plots["GLBmuonTimeBarrel"]  = new TH1F("GLBmuonTimeBarrel_" + sampleTag ,"GLBmuonTimeBarrel; time (ns) ; # entries", 400, -200., 200.);
  // m_plots["GLBmuonTimeEndcap"]  = new TH1F("GLBmuonTimeEndcap_" + sampleTag ,"GLBmuonTimeEndcap; time (ns) ; # entries", 400, -200., 200.);
 
 
  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
    {
         
      TString etaTag = "_fEtaMin" + (*fEtaMinIt) + "_fEtaMax" + (*fEtaMaxIt);

      //Id var
      m_plots["nHitsGLB" + sampleTag + etaTag]        = new TH1F("nHitsGLB_" + sampleTag + etaTag ,"nHitsGLB; # hits ; # entries", 80, 0., 80.);
      m_plots["nHitsTRK" + sampleTag + etaTag]        = new TH1F("nHitsTRK_" + sampleTag + etaTag ,"nHitsTRK; # hits ; # entries", 40, 0., 40.);
      m_plots["nHitsSTA" + sampleTag + etaTag]        = new TH1F("nHitsSTA_" + sampleTag + etaTag ,"nHitsSTA; # hits ; # entries", 60, 0., 60.);

      m_plots["nChi2GLB" + sampleTag + etaTag]        = new TH1F("nChi2GLB_" + sampleTag + etaTag ,"nChi2GLB; chi2/ndof ; # entries", 50, 0., 100.);
      m_plots["nChi2TRK" + sampleTag + etaTag]        = new TH1F("nChi2TRK_" + sampleTag + etaTag ,"nChi2TRK; chi2/ndof ; # entries", 100, 0., 50.);

      m_plots["matchedStation" + sampleTag + etaTag]  = new TH1F("matchedStation_" + sampleTag + etaTag ,"matchedStation; # stations ; # entries", 10, 0., 10.);
      m_plots["muonValidHitsGLB" + sampleTag + etaTag]= new TH1F("muonValidHitsGLB_" + sampleTag + etaTag ,"muonValidHitsGLB; # hits ; # entries", 60, 0., 60.);
      m_plots["PixelHitsTRK" + sampleTag + etaTag]    = new TH1F("PixelHitsTRK_" + sampleTag + etaTag ,"PixelHitsTRK; # hits ; # entries", 10, 0., 10.);
      m_plots["PixelLayersTRK" + sampleTag + etaTag]  = new TH1F("PixelLayersTRK_" + sampleTag + etaTag ,"PixelLayersTRK; # layers ; # entries", 10, 0., 10.);
      m_plots["TrackerLayersTRK" + sampleTag + etaTag]= new TH1F("TrackerLayersTRK_" + sampleTag + etaTag ,"TrackerLayersTRK; # layers ; # entries", 30, 0., 30.);

      m_plots["HitFractionTRK" + sampleTag + etaTag]  = new TH1F("HitFractionTRK_" + sampleTag + etaTag ,"HitFractionTRK; fraction ; # entries", 20, 0., 1.);
      m_plots["TrkStaChi2" + sampleTag + etaTag]      = new TH1F("TrkStaChi2_" + sampleTag + etaTag ,"TrkStaChi2; chi2 ; # entries", 50, 0., 100.);
      m_plots["TrkKink" + sampleTag + etaTag]         = new TH1F("TrkKink_" + sampleTag + etaTag ,"TrkKink", 48, 0., 1200.);
      m_plots["segmentComp" + sampleTag + etaTag]     = new TH1F("segmentComp_" + sampleTag + etaTag ,"segmentComp", 100, 0., 1.);

      m_plots["probePt" + sampleTag + etaTag]  = new TH1F("probePt_" + sampleTag + etaTag," ; muon p_[T] ; # entries", 75,0.,150.);
      m_plots["probeEta" + sampleTag + etaTag] = new TH1F("probeEta_" + sampleTag + etaTag," ; muon #eta ; # entries", 50,-2.5,2.5);
      m_plots["probePhi" + sampleTag + etaTag] = new TH1F("probePhi_" + sampleTag + etaTag," ; muon #phi ; # entries", 50,-TMath::Pi(),TMath::Pi());
      m_plots["probeDxy" + sampleTag + etaTag] = new TH1F("probeDxy_" + sampleTag + etaTag," ; muon dxy ; # entries", 100,-0.5,0.5);
      m_plots["probeDz" + sampleTag + etaTag]  = new TH1F("probeDz_" + sampleTag + etaTag," ;  muon dz ; # entries", 200,-2.5,2.5);
      m_plots["chHadIso" + sampleTag + etaTag]    = new TH1F("chHadIso_" + sampleTag + etaTag," ; muon relative isolation; # entries", 50,0.,5.);
      m_plots["photonIso" + sampleTag + etaTag]   = new TH1F("photonIso_" + sampleTag + etaTag," ; muon relative isolation; # entries", 50,0.,5.);
      m_plots["neutralIso" + sampleTag + etaTag]  = new TH1F("neutralIso_" + sampleTag + etaTag," ; muon relative isolation; # entries", 50,0.,5.);
      m_plots["dBetaRelIso" + sampleTag + etaTag] = new TH1F("dBetaRelIso_" + sampleTag + etaTag," ; muon relative isolation; # entries", 50,0.,2.);

    }

  outFile->mkdir(sampleTag+"/control");
  outFile->cd(sampleTag+"/control");

  m_plots["invMass" + sampleTag] = new TH1F("invMass_" + sampleTag ,"invMass", 100,0.,200.);
  m_plots["dilepPt" + sampleTag] = new TH1F("dilepPt_" + sampleTag ,"dilepPt", 100,0.,200.);
  m_plots["nVertices" + sampleTag] = new TH1F("nVertices_" + sampleTag ,"nVertices", 60,0.,60.);

  m_plots["invMassInRange" + sampleTag] = new TH1F("invMassInRange_" + sampleTag ,"invMass", 100,0.,200.);
  
  m_plots["nProbesVsnTags" + sampleTag] = new TH2F("nProbesVsnTags_" + sampleTag ,"invMass", 10,-0.5,9.,10,-0.5,9.);

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, int nVtx, float weight)
{
  
  TString sampleTag = m_sampleConfig.sampleName;

  /////// for muon timing only
  //std::cout << " Started " << std::endl;
  for (auto & muon : muons)
    {
      if (muon.isStandAlone && ((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) ||
				(fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)))
	m_plots["STAmuonTime" + sampleTag]->Fill(muon.muonTime,weight);
      
      if (muon.isStandAlone && fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20)
	m_plots["STAmuonTimeBarrel" + sampleTag]->Fill(muon.muonTime,weight);
      
      if (muon.isStandAlone && fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11)
	m_plots["STAmuonTimeEndcap" + sampleTag]->Fill(muon.muonTime,weight);
      
      // if (muon.isGlobal && muon.glbNormChi2 > 0 && muon.glbNormChi2 < 10.)
      // 	m_plots["GLBmuonTime"]->Fill(muon.muonTime,weight);

      // if (muon.isGlobal && muon.glbNormChi2 > 0. && muon.glbNormChi2 < 10. && fabs(muon.eta) < 0.9)
      //   m_plots["GLBmuonTimeBarrel"]->Fill(muon.muonTime,weight);

      // if (muon.isGlobal && muon.glbNormChi2 > 0. && muon.glbNormChi2 < 10. && fabs(muon.eta) > 1.2)
      //   m_plots["GLBmuonTimeEndcap"]->Fill(muon.muonTime,weight);

    }
  //  std::cout << "**** 2*********** "<< std::endl; 
	  //std::cout << " Finished " << std::endl;
  ///////////////////////////  
  
  bool pathHasFired = false;

  for (auto path : hlt.triggers)
    {
      if (path.find(m_tnpConfig.hlt_path) != std::string::npos)
	{
	  pathHasFired = true;
	  break;
	}
    }

  if (!pathHasFired) return;

  std::vector<muon_pog::Muon> tagMuons;

  for (auto & muon : muons)
    {
      //std::cout << hasGoodId(muon,"tag") << " " << hasFilterMatch(muon,hlt) << " " << muonTk(muon).Pt() << " " << muon.isoPflow04 << std::endl; 
      if (hasGoodId(muon,"tag") && hasFilterMatch(muon,hlt) &&
	  muonTk(muon).Pt() > m_tnpConfig.tag_minPt   &&
	  muon.isoPflow04 < m_tnpConfig.tag_isoCut)
	tagMuons.push_back(muon);
    }
  
  //if(tagMuons.size()>1) std::cout << " **** # tags: " << tagMuons.size() << " **** " << std::endl; 

  std::vector<muon_pog::Muon> probeMuons;

  for (auto & muon : muons)
    {
      for (auto & tagMuon : tagMuons)
	{

	  if( fabs(tagMuon.eta - muon.eta)>0.001 || 
	      fabs(tagMuon.phi - muon.phi)>0.001 || 
	      fabs(tagMuon.pt  - muon.pt )>0.001    ) // there exists at least another tag in the event  
	    { 

	      if(muon.isStandAlone) {
		if((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) || (fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)) 
		  m_plots["UnbSTAmuonTime" + sampleTag]->Fill(muon.muonTime,weight);
      
		if(fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20) 
		  m_plots["UnbSTAmuonTimeBarrel" + sampleTag]->Fill(muon.muonTime,weight);
      
		if(fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11) 
		  m_plots["UnbSTAmuonTimeEndcap" + sampleTag]->Fill(muon.muonTime,weight);
	      }

	  if ( chargeFromTrk(tagMuon) * chargeFromTrk(muon) == -1 &&    
	       (muon.isGlobal || muon.isTrackerArb)               && 
	       muon.pt > 5                                          ) // CB minimal cuts on potential probe // DT: changed from isTracker to isTrackerArb!  
	    {
	      
	      TLorentzVector tagMuTk = muonTk(tagMuon);
	      TLorentzVector muTk    = muonTk(muon);
	      
	      Float_t mass = (tagMuTk+muTk).M();

	      // CB Fill control plots
	      m_plots["invMass" + sampleTag]->Fill(mass,weight);

	      if ( mass > m_tnpConfig.pair_minInvMass &&
		   mass < m_tnpConfig.pair_maxInvMass )
		{
		  // m_plots["invMassInRange"]->Fill(mass,weight);
	      
		  // Float_t dilepPt = (tagMuTk+muTk).Pt();
		  // m_plots["dilepPt"]->Fill(dilepPt,weight);
		  // //std::cout << " # vertices = " << nVtx << std::endl;
		  // m_plots["nVertices"]->Fill(nVtx,weight);

		  m_plots["nVertices" + sampleTag]->Fill(nVtx, weight);

		  if(hasGoodId(muon,"probe")) { 
		    m_plots["invMassInRange" + sampleTag]->Fill(mass, weight);
		    m_plots["dilepPt" + sampleTag]->Fill((tagMuTk+muTk).Pt(), weight);
		  }

		  std::vector<TString>::const_iterator fEtaMinIt  = m_tnpConfig.probe_fEtaMin.begin();
		  std::vector<TString>::const_iterator fEtaMinEnd = m_tnpConfig.probe_fEtaMin.end();

		  std::vector<TString>::const_iterator fEtaMaxIt  = m_tnpConfig.probe_fEtaMax.begin();
		  std::vector<TString>::const_iterator fEtaMaxEnd = m_tnpConfig.probe_fEtaMax.end();
  
		  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt) { 

		    TString etaTag = "_fEtaMin" + (*fEtaMinIt) + "_fEtaMax" + (*fEtaMaxIt);

		    //id var
		    m_plots["nHitsGLB" + sampleTag + etaTag]->Fill(muon.nHitsGlobal ,weight);
		    m_plots["nHitsTRK" + sampleTag + etaTag]->Fill(muon.nHitsTracker,weight);
		    m_plots["nHitsSTA" + sampleTag + etaTag]->Fill(muon.nHitsStandAlone,weight);
		  
		    m_plots["nChi2GLB" + sampleTag + etaTag]->Fill(muon.glbNormChi2, weight);
		    m_plots["nChi2TRK" + sampleTag + etaTag]->Fill(muon.trkNormChi2, weight);
		  
		    m_plots["matchedStation" + sampleTag + etaTag]->Fill(muon.trkMuonMatchedStations, weight)  ;
		    m_plots["muonValidHitsGLB" + sampleTag + etaTag]->Fill(muon.glbMuonValidHits, weight);
		    m_plots["PixelHitsTRK" + sampleTag + etaTag]->Fill(muon.trkPixelValidHits, weight);
		    m_plots["PixelLayersTRK" + sampleTag + etaTag]->Fill(muon.trkPixelLayersWithMeas, weight); 
		    m_plots["TrackerLayersTRK" + sampleTag + etaTag]->Fill(muon.trkTrackerLayersWithMeas, weight); 
		  
		    m_plots["HitFractionTRK" + sampleTag + etaTag]->Fill(muon.trkValidHitFrac, weight);   
		    m_plots["TrkStaChi2" + sampleTag + etaTag]->Fill(muon.trkStaChi2, weight);       
		    m_plots["TrkKink" + sampleTag + etaTag]->Fill(muon.trkKink, weight);          
		    m_plots["segmentComp" + sampleTag + etaTag]->Fill(muon.muSegmComp, weight);      

		    m_plots["probeDxy" + sampleTag + etaTag]->Fill(muon.dxy,weight);
		    m_plots["probeDz" + sampleTag + etaTag]->Fill(muon.dz,weight);
		  }
	      
		  probeMuons.push_back(muon);
		  //continue; // CB If a muon is already a probe don't look for other tags
		  break; // DT: a continue statement has no effect, it must be a break!  
		}
	    }
	    }
	}
    }
  //std::cout << "**** 4*********** "<< std::endl;

  m_plots["nProbesVsnTags" + sampleTag]->Fill(tagMuons.size(),probeMuons.size());
  
  for (auto & probeMuon : probeMuons)
    {

      if(hasGoodId(probeMuon,"probe")) { 

	std::vector<TString>::const_iterator fEtaMinIt  = m_tnpConfig.probe_fEtaMin.begin();
	std::vector<TString>::const_iterator fEtaMinEnd = m_tnpConfig.probe_fEtaMin.end();

	std::vector<TString>::const_iterator fEtaMaxIt  = m_tnpConfig.probe_fEtaMax.begin();
	std::vector<TString>::const_iterator fEtaMaxEnd = m_tnpConfig.probe_fEtaMax.end();
  
	TLorentzVector probeMuTk = muonTk(probeMuon);
	for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt) { 
	  
	  if (fabs(probeMuTk.Eta()) > fEtaMinIt->Atof() &&
	      fabs(probeMuTk.Eta()) < fEtaMaxIt->Atof() )
	    {
	      
	      TString etaTag = "_fEtaMin" + TString((*fEtaMinIt)) + "_fEtaMax" + TString((*fEtaMaxIt));
	      m_plots["probePt" + sampleTag + etaTag]->Fill(probeMuTk.Pt(),weight);
	      m_plots["probeEta" + sampleTag + etaTag]->Fill(probeMuTk.Eta(),weight);
	      m_plots["probePhi" + sampleTag + etaTag]->Fill(probeMuTk.Phi(),weight);
	      
	      // Fill isolation plots for muons passign a given identification (programmable from cfg)
	      m_plots["photonIso" + sampleTag + etaTag]->Fill(probeMuon.photonIso,weight);
	      m_plots["chHadIso" + sampleTag + etaTag]->Fill(probeMuon.chargedHadronIso,weight);
	      m_plots["neutralIso" + sampleTag + etaTag]->Fill(probeMuon.neutralHadronIso,weight);
	      m_plots["dBetaRelIso" + sampleTag + etaTag]->Fill(probeMuon.isoPflow04,weight);
	    }
	  
	}
      }
    }

  //std::cout << "**** 4B*********** "<< std::endl;  
}

bool muon_pog::Plotter::hasGoodId(const muon_pog::Muon & muon, std::string leg)
{
  std::string & muId = leg == "tag" ? m_tnpConfig.tag_ID : m_tnpConfig.probe_ID ;

  if (muId == "GLOBAL")      return muon.isGlobal == 1 ;
  else if (muId == "TIGHT")  return muon.isTight  == 1;
  else if (muId == "MEDIUM") return muon.isMedium == 1;
  else if (muId == "LOOSE")  return muon.isLoose == 1;
  else if (muId == "HIGHPT") return muon.isHighPt == 1;
  else if (muId == "SOFT")   return muon.isSoft == 1;
  else
    {
      std::cout << "[Plotter::hasGoodId]: Invalid muon id : "
		<< muId << std::endl;
      exit(900);
    }
  
}

bool muon_pog::Plotter::hasFilterMatch(const muon_pog::Muon & muon,
				       const muon_pog::HLT  & hlt )
{
  std::string & filter = m_tnpConfig.tag_hltFilter;
  TLorentzVector muTk = muonTk(muon);

  for (auto object : hlt.objects)
    {
      float Deta = muTk.Eta() - object.eta; 
      float Dphi = muTk.Phi() - object.phi; // Mind phi boundaries! 
      if      (Dphi >   TMath::Pi()) Dphi -= 2*TMath::Pi();
      else if (Dphi <= -TMath::Pi()) Dphi += 2*TMath::Pi();
      if(object.filterTag.find(filter) != std::string::npos && sqrt(Deta*Deta + Dphi*Dphi) < m_tnpConfig.tag_hltDrCut ) 
	return true;
    }

  return false;
  	  
}


Int_t muon_pog::Plotter::chargeFromTrk(const muon_pog::Muon & muon)
{
  std::string & trackType = m_tnpConfig.muon_trackType;

  if (trackType == "PF")         return muon.charge;
  else if (trackType == "TUNEP") return muon.charge_tuneP;
  else if (trackType == "GLB")   return muon.charge_global;
  else if (trackType == "INNER") return muon.charge_tracker;
  else
    {
      std::cout << "[Plotter::chargeFromTrk]: Invalid track type: "
		<< trackType << std::endl;
      exit(900);
    }
  
}

TLorentzVector muon_pog::Plotter::muonTk(const muon_pog::Muon & muon)
{
  std::string & trackType = m_tnpConfig.muon_trackType;

  TLorentzVector result; 
  if (trackType == "PF")
    result.SetPtEtaPhiM(muon.pt,muon.eta,muon.phi,.10565);
  else if (trackType == "TUNEP")
    result.SetPtEtaPhiM(muon.pt_tuneP,muon.eta_tuneP,muon.phi_tuneP,.10565);
  else if (trackType == "GLB")
    result.SetPtEtaPhiM(muon.pt_global,muon.eta_global,muon.phi_global,.10565);
  else if (trackType == "INNER")
    result.SetPtEtaPhiM(muon.pt_tracker,muon.eta_tracker,muon.phi_tracker,.10565);
  else
    {
      std::cout << "[Plotter::muonTk]: Invalid track type: "
		<< trackType << std::endl;
      exit(900);
    }

  return result;
  
}

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

  //std::cout << "**** 4C*********** "<< std::endl;
}

void muon_pog::comparisonPlot(TFile *outFile,TString plotName,
			      std::vector<muon_pog::Plotter> & plotters)
{

  //std::cout << "**** 4D*********** "<< std::endl;

  THStack hMc(plotName,"");
  TH1 * hData = 0;
  //std::cout << "**** 5*********** "<< std::endl; 

  int colorMap[5] { kAzure+7, kGreen+1, kOrange+7, kOrange, kGray+1};

  int iColor = 0;

  float integralData = 1.;
  float integralMC   = 1.;
  //float totalXSec = 0.;
  //float totalEvents = 0.;

  float integralTimeData = 1.;
  float integralTimeMC   = 1.;

  for (auto plotter : plotters)
    {
      TString sampleTag = plotter.m_sampleConfig.sampleName;

      if(std::string(plotter.m_sampleConfig.sampleName.Data()).find("Data") != std::string::npos)
	{
	  integralData = plotter.m_plots["invMassInRange" + sampleTag]->Integral(); // CB Now scales using the inv mass in range integral
	  if(plotName.Contains("muonTime")) { 
	    Int_t firstBin = plotter.m_plots[plotName]->FindBin(-2.5); 
	    Int_t lastBin  = plotter.m_plots[plotName]->FindBin(2.5); 
	    integralTimeData = plotter.m_plots[plotName]->Integral(firstBin, lastBin); 
	  } 
	}
      else
	{
	  integralMC = plotter.m_plots["invMassInRange" + sampleTag]->Integral(); // CB Now scales using the inv mass in range integral
	  //totalXSec += plotter.m_sampleConfig.cSection;
	  //totalEvents += plotter.m_sampleConfig.eventi;
	  if(plotName.Contains("muonTime")) { 
	    Int_t firstBin = plotter.m_plots[plotName]->FindBin(-2.5); 
	    Int_t lastBin  = plotter.m_plots[plotName]->FindBin(2.5); 
	    integralTimeMC = plotter.m_plots[plotName]->Integral(firstBin, lastBin); 
	  } 
	}

    }

  for (auto plotter : plotters)
    {
      if(std::string(plotter.m_sampleConfig.sampleName.Data()).find("Data") != std::string::npos)
	{
	  hData = plotter.m_plots[plotName];
	  hData->Sumw2();
	}
      else
	{
	  plotter.m_plots[plotName]->SetFillColor(colorMap[iColor]);
	  plotter.m_plots[plotName]->SetMarkerStyle(0);
	  TH1* plot = plotter.m_plots[plotName];
	  //float scale = integralData / integralMC; 
	  // float scale = plotter.m_sampleConfig.cSection / totalXSec *
	  //   integralData / integralMC;
	  // float mcfactor = plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.eventi
	  //float mcfactor = totalEvents/plotter.m_sampleConfig.eventi; // a che serve 'sta roba??  
	  //plot->Scale(scale);
	  //plot->Scale(mcfactor);
	  plot->Scale(integralData/integralMC); 
	  if(plotName.Contains("muonTime")) { 
	    plot->Scale((integralMC/integralData)*(integralTimeData/integralTimeMC)); 
	  } 
	  hMc.Add(plot);
	  iColor++;
	}
    }

  TCanvas *canvas = new TCanvas("c"+plotName, "c"+plotName, 500, 500);

  // std::cout << "**** 8*********** "<< std::endl; 

  canvas->cd();
  canvas->SetLogy(1);

  hData->Draw();
  hMc.Draw("samehist");
  hData->Draw("same");

  canvas->Write();
  canvas->SaveAs("c"+plotName+".png");
  //std::cout << "**** 9*********** "<< std::endl; 

  //Canvas with ratio plots
  TCanvas *rcanvas = new TCanvas("rc"+plotName, "rc"+plotName, 500, 700);
  rcanvas->Divide(1,2);
  rcanvas->cd(1);
  TPad *pad1 = (TPad*)rcanvas->GetPad(1);
  pad1->SetPad(0.,0.2,1.,1.);
  pad1->SetLogy(1);
  hData->Draw();
  hMc.Draw("samehist");
  hData->Draw("same");
  
  rcanvas->cd(2);
  TPad *pad2 = (TPad*)rcanvas->GetPad(2);
  pad2->SetPad(0.,0.,1.,0.2);
 
  TH1 *h2 = (TH1*)hData->Clone();
  Int_t xmin  = h2->GetMinimumBin();
  Int_t xmax  = h2->GetMaximumBin();
  Int_t nbins = h2->GetNbinsX();
  TH1 *hs = (TH1*)hMc.GetStack()->Last();  

  //TH1F *h2 = new TH1F("h2","residuals",100,xmin,xmax);
  //h2->GetXaxis()->SetLabelFont(63);
  //h2->GetXaxis()->SetLabelSize(16);
  //h2->GetXaxis()->SetTitle("P_{t} (GeV)");
  //h2->GetYaxis()->SetLabelFont(63);
  //h2->GetYaxis()->SetLabelSize(16);
  h2->GetYaxis()->SetTitle("Data/MC");
  
  //  std::cout << "**** 10*********** "<< std::endl; 
  
  // }
  h2->Divide(hs);
  h2->Draw();
  if(plotName == "nVertices"){
    for (Int_t i=1;i<=nbins;i++) {
      std::cout << " Bin # " << i << " PUweight = " << h2->GetBinContent(i) << std::endl;
    }
  }
    
  rcanvas->Write();
  rcanvas->SaveAs("rc"+plotName+".png");
}
