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
#include "tdrstyle.C"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
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
    Float_t nEvents;
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
    std::vector<std::pair<TString,TString> > probe_fEtaBins;
     
    std::string hlt_path; 
   
    TagAndProbeConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    TagAndProbeConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~TagAndProbeConfig() {};
    
  private:
    std::vector<TString> toArray(const std::string & entries); 
    std::vector<std::pair<TString,TString> > toPairArray(const std::vector<TString> &,
							 const std::vector<TString> &); 
  
  };

  class Observable {

  public :
    
    Observable() {};

    Observable(TString hName, TString sampleTag, TString xTitle, TString yTitle,
	       Int_t nBins, Float_t min, Float_t max, bool kinPlots);
    
    ~Observable() {};//{m_plots.clear(); };

    void fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Float_t mcScale, Int_t PV);
    
    std::vector<TH1 *> & plots() { return m_plots; }
    
  private :
    
    std::vector<TH1 *> m_plots;

  };

  class Plotter {

  public :

    enum HistoType { KIN=0, ID, ISO, TIME, CONT };
      
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, int nVtx, float weight);

    std::map<Plotter::HistoType, std::map<TString, muon_pog::Observable> > m_plots;
    TagAndProbeConfig m_tnpConfig;
    SampleConfig m_sampleConfig;
        
  private :
   
    double deltaR(double eta1, double phi1, double eta2, double phi2);
    
    bool hasGoodId(const muon_pog::Muon & muon,
		   TString leg);
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
	// for (Long64_t iEvent=0; iEvent<1000; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;
	  
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;

	  //weight for Data-Mc comparison config4.ini 2015D vs All MC 
	  float PUweight[60] = {0,6.68472,4.41435,3.49289,3.14212,3.173,2.94319,2.71174,2.44124,2.0717,1.7103,1.39812,1.11666,0.855547,0.618711,0.453196,0.328689,0.228592,0.166069,0.100075,0.070327,0.0567964,0.0405713,0.0308903,0.0186064,0.0208503,0.00264836,0.00779018,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	  //weight for Data-Mc comparison config4.ini 2015D vs All MC  
	  // float PUweight[60] = {0,
	  // 			6.68472,
	  // 			4.41435,
	  // 			3.49289,
	  // 			3.14212,
	  // 			3.173,
	  // 			2.94319,
	  // 			2.71174,
	  // 			2.44124,
	  // 			2.0717,
	  // 			1.7103,
	  // 			1.39812,
	  // 			1.11666,
	  // 			0.855547,
	  // 			0.618711,
	  // 			0.453196,
	  // 			0.328689,
	  // 			0.228592,
	  // 			0.166069,
	  // 			0.100075,
	  // 			0.070327,
	  // 			0.0567964,
	  // 			0.0405713,
	  // 			0.0308903,
	  // 			0.0186064,
	  // 			0.0208503,
	  // 			0.00264836,
	  // 			0.00779018,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0};
	  

	  //weight for Data-Mc comparison config3.ini 2015C_Trocino vs All MC  
	  // float PUweight[60] = {0,
	  // 			2.59946,
	  // 			1.97202,
	  // 			0.879335,
	  // 			0.976222,
	  // 			1.18325,
	  // 			1.18647,
	  // 			1.17652,
	  // 			1.39881,
	  // 			1.40331,
	  // 			1.35794,
	  // 			1.35981,
	  // 			1.24513,
	  // 			1.25133,
	  // 			1.02661,
	  // 			0.990347,
	  // 			0.861541,
	  // 			0.720651,
	  // 			0.661865,
	  // 			0.504051,
	  // 			0.523631,
	  // 			0.37145,
	  // 			0.311188,
	  // 			0.308495,
	  // 			0.273079,
	  // 			0.174498,
	  // 			0.131307,
	  // 			0.12496,
	  // 			0.149624,
	  // 			0,
	  // 			0.194291,
	  // 			0.102388,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0,
	  // 			0};
	 
				 
	  //weight for Data to DY comparison obsoleto
	  //float PUweight[60] = {0, 6.23938, 3.93412, 1.0139, 0.92968, 0.852813, 0.86357, 0.991675, 1.13882, 1.29982, 1.17863, 1.25089, 1.17501, 1.29242, 1.11248, 1.08851, 0.952472, 0.873848, 0.842575, 0.620411, 0.6636, 0.50119, 0.404807, 0.369032, 0.483837, 0.266346, 0.308719, 0.207252, 0.175888, 0.0609816, 0.423387, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	  
	  //weight for Data to Data comparison 2015C 25ns "old" vs 50 ns
	  //float PUweight[60] = {0, 2.59252, 7.0705, 6.74054, 7.77755, 3.26155, 3.42985, 2.89984, 2.28284, 2.13837, 1.70971, 1.50287, 1.20367, 1.17027, 0.956894, 0.881644, 0.7398, 0.635878, 0.595542, 0.494059, 0.499049, 0.408914, 0.344204, 0.286958, 0.388219, 0.244787, 0.265499, 0.216043, 0.172834, 0.0540108, 0.360072, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	 
	  if(plotter.m_sampleConfig.applyReweighting==true)
	    weight *= ev->nVtx < 60 ? PUweight[ev->nVtx] : 0;
	    
	  plotter.fill(ev->muons, ev->hlt, ev->nVtx, weight);
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
      probe_fEtaBins = toPairArray(toArray(vt.second.get<std::string>("probe_fEtaMin")),
				   toArray(vt.second.get<std::string>("probe_fEtaMax")));
  
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
      nEvents = vt.second.get<Float_t>("nEvents"); //CB do we really need this? can't we take nEvents from the file itself?
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

std::vector<std::pair<TString,TString> > muon_pog::TagAndProbeConfig::toPairArray(const std::vector<TString> & fEtaMin,
										  const std::vector<TString> & fEtaMax)
{

  std::vector<std::pair<TString,TString> > result;

  std::vector<TString>::const_iterator fEtaMinIt  = fEtaMin.begin();
  std::vector<TString>::const_iterator fEtaMinEnd = fEtaMin.end();

  std::vector<TString>::const_iterator fEtaMaxIt  = fEtaMax.begin();
  std::vector<TString>::const_iterator fEtaMaxEnd = fEtaMax.end();
  
  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
    {
      result.push_back(std::make_pair((*fEtaMinIt),
				      (*fEtaMaxIt)));
    }

  return result;

}

muon_pog::Observable::Observable(TString hName, TString sampleTag, TString xTitle, TString yTitle,
				 Int_t nBins, Float_t min, Float_t max, bool kinPlots)
{
  m_plots.push_back(new TH1F("h" + hName + "_" + sampleTag, hName + " ;" + xTitle + ";" + yTitle, nBins, min, max));
  if (kinPlots)
    { // CB book here for other plots vs kin variables
      if(sampleTag.Contains("Data")|| sampleTag.Contains("TWantitop")){
	m_plots.push_back(new TProfile("h" + hName + "VsEta_" + sampleTag, hName + " vs #eta;  #eta;"        + xTitle, 24, -2.4, 2.4, min, max));
	m_plots.push_back(new TProfile("h" + hName + "VsPhi_" + sampleTag, hName + " vs #phi;  #phi;"        + xTitle, 25, -TMath::Pi(),TMath::Pi(), min, max));
	m_plots.push_back(new TProfile("h" + hName + "VsPt_"  + sampleTag, hName + " vs p_{T}; p_{T} (GeV);" + xTitle, 50,  0., 150., min, max));
	m_plots.push_back(new TProfile("h" + hName + "VsPV_"  + sampleTag, hName + " vs PV;    # of PV;"     + xTitle, 60,  0., 60., min, max));
      }
           
    }
  
}

void muon_pog::Observable::fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Float_t mcScale, Int_t PV)
{
  
  m_plots.at(0)->Fill(value, weight);
  if (m_plots.size() > 1){
    weight *= mcScale;
    // CB fillhere for other plots vs kin variables
    ((TProfile*)m_plots.at(1))->Fill(muonTk.Eta(), value, weight);
  if (m_plots.size() > 2)
    ((TProfile*)m_plots.at(2))->Fill(muonTk.Phi(), value, weight);
  if (m_plots.size() > 3)
    ((TProfile*)m_plots.at(3))->Fill(muonTk.Pt(),  value, weight);
  if (m_plots.size() > 4)
    ((TProfile*)m_plots.at(4))->Fill(PV,           value, weight);
  }
}

void muon_pog::Plotter::book(TFile *outFile)
{

  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);

  outFile->mkdir(sampleTag+"/timing");
  outFile->mkdir(sampleTag+"/id_variables");
  outFile->mkdir(sampleTag+"/kinematical_variables");
  outFile->mkdir(sampleTag+"/isolation");
  outFile->mkdir(sampleTag+"/control");
      
  outFile->cd(sampleTag+"/timing");
  
  std::cout << sampleTag << "  Begin: " << m_plots.size()<< std::endl;

  m_plots[TIME]["STAmuonTime"]       = muon_pog::Observable("STAmuonTime", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["STAmuonTimeBarrel"] = muon_pog::Observable("STAmuonTimeBarrel", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["STAmuonTimeEndcap"] = muon_pog::Observable("STAmuonTimeEndcap", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);

  m_plots[TIME]["UnbSTAmuonTime"] = muon_pog::Observable("UnbSTAmuonTime", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["UnbSTAmuonTimeBarrel"] = muon_pog::Observable("UnbSTAmuonTimeBarrel", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["UnbSTAmuonTimeEndcap"] = muon_pog::Observable("UnbSTAmuonTimeEndcap", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);

  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
    {
      
      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;

      outFile->cd(sampleTag+"/id_variables");

      //Id var

      bool pippo = false;
      if(sampleTag.Contains("Data")|| sampleTag.Contains("TWantitop"))
	bool pippo = true;

      m_plots[ID]["NHitsGLB" + etaTag] = muon_pog::Observable("NHitsGLB" + etaTag, sampleTag, "# hits", "# entries", 80, 0., 80., pippo);
      m_plots[ID]["NHitsTRK" + etaTag] = muon_pog::Observable("NHitsTRK" + etaTag, sampleTag, "# hits", "# entries", 40, 0., 40., pippo);
      m_plots[ID]["NHitsSTA" + etaTag] = muon_pog::Observable("NHitsSTA" + etaTag, sampleTag, "# hits", "# entries", 60, 0., 60., pippo);
      
      m_plots[ID]["Chi2GLB" + etaTag]  = muon_pog::Observable("Chi2GLB_" + etaTag, sampleTag, "chi2/ndof", "# entries", 50, 0., 100., pippo);
      m_plots[ID]["Chi2TRK" + etaTag]  = muon_pog::Observable("Chi2TRK_" + etaTag, sampleTag, "chi2/ndof", "# entries", 100, 0., 50., pippo);
      
      m_plots[ID]["NMatchedStation" + etaTag]  = muon_pog::Observable("NatchedStation_" + etaTag, sampleTag,"# stations", "# entries", 10, 0., 10., pippo);
      m_plots[ID]["NMuonValidHitsGLB" + etaTag]= muon_pog::Observable("NMuonValidHitsGLB_" + etaTag, sampleTag,"# hits", "# entries", 60, 0., 60., pippo);
      m_plots[ID]["PixelHitsTRK" + etaTag]    = muon_pog::Observable("PixelHitsTRK" + etaTag, sampleTag,"# hits", "# entries", 10, 0., 10., pippo);
      m_plots[ID]["PixelLayersTRK" + etaTag]  = muon_pog::Observable("PixelLayersTRK" + etaTag, sampleTag,"# layers", "# entries", 10, 0., 10., pippo);
      m_plots[ID]["TrackerLayersTRK" + etaTag]= muon_pog::Observable("TrackerLayersTRK" + etaTag, sampleTag,"# layers"," # entries", 30, 0., 30., pippo);
      
      m_plots[ID]["HitFractionTRK" + etaTag]  = muon_pog::Observable("HitFractionTRK" + etaTag, sampleTag, "fraction", " # entries", 20, 0., 1., pippo);
      m_plots[ID]["TrkStaChi2" + etaTag]      = muon_pog::Observable("TrkStaChi2" + etaTag, sampleTag, "chi2", "# entries", 50, 0., 100., pippo);
      m_plots[ID]["TrkKink" + etaTag]         = muon_pog::Observable("TrkKink" + etaTag, sampleTag, "prob.", "# entries", 100, 0., 250., pippo);
      m_plots[ID]["SegmentComp" + etaTag]     = muon_pog::Observable("SegmentComp" + etaTag, sampleTag,"prob.", "segmentComp", 100, 0., 1., pippo);
      
      m_plots[ID]["Dxy" + etaTag]             = muon_pog::Observable("Dxy" + etaTag, sampleTag, "d_{xy} (cm)", "# entries", 100,-0.5,0.5, pippo);
      m_plots[ID]["Dz" + etaTag]              = muon_pog::Observable("Dz" + etaTag, sampleTag, "d_{z} (cm)", "# entries", 200,-2.5,2.5, pippo);    

      outFile->cd(sampleTag+"/kinematical_variables");

      for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  m_plots[KIN]["ProbePt" + etaTag + IDTag]  = muon_pog::Observable("ProbePt" + etaTag + IDTag, sampleTag, "p_{T} (GeV)", "# entries", 75,0.,150., false);
	  m_plots[KIN]["ProbeEta" + etaTag + IDTag] = muon_pog::Observable("hProbeEta_" + etaTag + IDTag, sampleTag, "#eta", "# entries", 50,-2.5,2.5, false);
	  m_plots[KIN]["ProbePhi" + etaTag + IDTag] = muon_pog::Observable("hProbePhi_" + etaTag + IDTag, sampleTag, "#phi", "# entries", 50,-TMath::Pi(),TMath::Pi(), false);      
	}
  
      outFile->cd(sampleTag+"/isolation");

      for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  m_plots[ISO]["ChHadIso" + etaTag + IDTag]    = muon_pog::Observable("ChHadIso_"    + etaTag + IDTag, sampleTag, "charged had. iso 0.4", "# entries", 50,0.,10., pippo);
	  m_plots[ISO]["ChHadIsoPU" + etaTag + IDTag]  = muon_pog::Observable("ChHadIsoPU_"  + etaTag + IDTag, sampleTag, "PU Charged had. iso 0.4", "# entries", 50,0.,10., pippo);
	  m_plots[ISO]["PhotonIso" + etaTag + IDTag]   = muon_pog::Observable("PhotonIso_"   + etaTag + IDTag, sampleTag, "photon iso 0.4", " # entries", 50,0.,10., pippo);
	  m_plots[ISO]["NeutralIso" + etaTag + IDTag]  = muon_pog::Observable("NeutralIso_"  + etaTag + IDTag, sampleTag, "neutral had. iso 0.4", "# entries", 50,0.,10., pippo);
	  m_plots[ISO]["DBetaRelIso" + etaTag + IDTag] = muon_pog::Observable("DBetaRelIso_" + etaTag + IDTag, sampleTag, "PFIso 0.4 (dBeta)", "# entries", 50,0.,2., pippo);
	}
    }

  outFile->cd(sampleTag+"/control");
  
  m_plots[CONT]["01_invMass"] = muon_pog::Observable("invMass", sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
  m_plots[CONT]["02_dilepPt"] = muon_pog::Observable("dilepPt", sampleTag ,"p_{T} (GeV)", "# entries", 100,0.,200., false);
  m_plots[CONT]["03_nVertices"] = muon_pog::Observable("nVertices", sampleTag ,"# vertices", "# entries", 60,0.,60., false);
  
  m_plots[CONT]["99_invMassInRange"] = muon_pog::Observable("invMassInRange", sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);

  
  std::cout << sampleTag << "  End: " << m_plots.size()<< std::endl;
  

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, int nVtx, float weight)
{
  
  TLorentzVector emptyTk;
  float mcScale  = 0;
  if(m_sampleConfig.sampleName.Contains("Data"))
    mcScale = 1.;
  else
    mcScale = m_sampleConfig.cSection/m_sampleConfig.nEvents;
  
  //muon timing only
  for (auto & muon : muons)
    {
      if (muon.isStandAlone && ((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) ||
				(fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)))
	m_plots[TIME]["STAmuonTime"].fill(muon.muonTime,emptyTk,weight, mcScale,nVtx);
      
      if (muon.isStandAlone && fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20)
	m_plots[TIME]["STAmuonTimeBarrel"].fill(muon.muonTime,emptyTk,weight, mcScale,nVtx);
      
      if (muon.isStandAlone && fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11)
	m_plots[TIME]["STAmuonTimeEndcap"].fill(muon.muonTime,emptyTk,weight, mcScale,nVtx);
    }
  
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

  std::vector<const muon_pog::Muon *> tagMuons;
  
  for (auto & muon : muons)
    {
      if (muonTk(muon).Pt() > m_tnpConfig.tag_minPt &&
	  hasFilterMatch(muon,hlt) &&
	  hasGoodId(muon,m_tnpConfig.tag_ID) && 
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
	       chargeFromTrk(tagMuon) * chargeFromTrk(muon) == -1)
	    {
	      if(muon.isStandAlone) {
		if((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) || (fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)) 
		  m_plots[TIME]["UnbSTAmuonTime"].fill(muon.muonTime,emptyTk,weight, mcScale,nVtx);
		
		if(fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20) 
		  m_plots[TIME]["UnbSTAmuonTimeBarrel"].fill(muon.muonTime,emptyTk,weight, mcScale,nVtx);
		
		if(fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11) 
		  m_plots[TIME]["UnbSTAmuonTimeEndcap"].fill(muon.muonTime,emptyTk,weight, mcScale,nVtx);
	      }
	          
	      if((muon.isGlobal || muon.isTracker) &&  // CB minimal cuts on potental probe 
		 muon.pt > m_tnpConfig.probe_minPt)
		{
		  TLorentzVector tagMuTk(muonTk(tagMuon));
		  TLorentzVector muTk(muonTk(muon));
		  
		  Float_t mass = (tagMuTk+muTk).M();
		  
		  // CB Fill control plots
		  m_plots[CONT]["01_invMass"].fill(mass,emptyTk,weight, mcScale,nVtx);
		  if ( mass > m_tnpConfig.pair_minInvMass &&
		       mass < m_tnpConfig.pair_maxInvMass )
		    {
		      m_plots[CONT]["99_invMassInRange"].fill(mass,emptyTk,weight, mcScale,nVtx);
		      
		      Float_t dilepPt = (tagMuTk+muTk).Pt();
		      m_plots[CONT]["02_dilepPt"].fill(dilepPt,emptyTk,weight, mcScale,nVtx);
		      m_plots[CONT]["03_nVertices"].fill(nVtx,emptyTk,weight, mcScale,nVtx);
 
		      probeMuons.push_back(&muon);
		      continue; // CB If a muon is already a probe don't loo on other tags
		    }
		}
	    }
	}
    }
  
  // m_plots[CONT]["04_nProbesVsnTags"]->Fill(tagMuons.size(),probeMuons.size());
 
  
  for (auto probeMuonPointer : probeMuons)
    {
      const muon_pog::Muon & probeMuon = *probeMuonPointer;
      			  
      for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
	{

	  TLorentzVector probeMuTk(muonTk(probeMuon));
	  
	  if (fabs(probeMuTk.Eta()) > fEtaBin.first.Atof() &&
	      fabs(probeMuTk.Eta()) < fEtaBin.second.Atof() )
	    {
	      
	      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;

	      //id var
	      m_plots[ID]["NHitsGLB"  + etaTag].fill(probeMuon.nHitsGlobal, probeMuTk ,weight, mcScale,nVtx);
	      m_plots[ID]["NHitsTRK"  + etaTag].fill(probeMuon.nHitsTracker, probeMuTk, weight, mcScale, nVtx);
	      m_plots[ID]["NHitsSTA"  + etaTag].fill(probeMuon.nHitsStandAlone, probeMuTk, weight, mcScale, nVtx);

	      m_plots[ID]["Chi2GLB"  + etaTag].fill(probeMuon.glbNormChi2, probeMuTk, weight, mcScale, nVtx);
	      m_plots[ID]["Chi2TRK"  + etaTag].fill(probeMuon.trkNormChi2, probeMuTk, weight, mcScale, nVtx);

	      m_plots[ID]["PixelHitsTRK"  + etaTag].fill(probeMuon.trkPixelValidHits, probeMuTk, weight, mcScale, nVtx);
	      m_plots[ID]["PixelLayersTRK"  + etaTag].fill(probeMuon.trkPixelLayersWithMeas, probeMuTk, weight, mcScale, nVtx); 
	      m_plots[ID]["TrackerLayersTRK"  + etaTag].fill(probeMuon.trkTrackerLayersWithMeas, probeMuTk, weight, mcScale, nVtx); 

	      m_plots[ID]["NMatchedStation"  + etaTag].fill(probeMuon.trkMuonMatchedStations, probeMuTk, weight, mcScale, nVtx);
	      m_plots[ID]["NMuonValidHitsGLB"  + etaTag].fill(probeMuon.glbMuonValidHits, probeMuTk, weight, mcScale, nVtx);

	      m_plots[ID]["HitFractionTRK"  + etaTag].fill(probeMuon.trkValidHitFrac, probeMuTk, weight, mcScale, nVtx);   
	      m_plots[ID]["SegmentComp"  + etaTag].fill(probeMuon.muSegmComp, probeMuTk, weight, mcScale, nVtx);      
	      m_plots[ID]["TrkStaChi2"  + etaTag].fill(probeMuon.trkStaChi2, probeMuTk, weight, mcScale, nVtx);       
	      m_plots[ID]["TrkKink"  + etaTag].fill(probeMuon.trkKink, probeMuTk, weight, mcScale, nVtx);          
	      
	      m_plots[ID]["Dxy"  + etaTag].fill(probeMuon.dxy, probeMuTk, weight, mcScale, nVtx);
	      m_plots[ID]["Dz"  + etaTag].fill(probeMuon.dz, probeMuTk, weight, mcScale, nVtx);			

	      for (auto & probe_ID : m_tnpConfig.probe_IDs)
	      	{
	      	  TString IDTag = "_" + probe_ID;
		  
	      	  if(hasGoodId(probeMuon,probe_ID)) 
	      	    {	 
	      	      m_plots[KIN]["ProbePt" + etaTag + IDTag].fill(probeMuTk.Pt(), emptyTk, weight, mcScale, nVtx);
	      	      m_plots[KIN]["ProbeEta" + etaTag + IDTag].fill(probeMuTk.Eta(), emptyTk, weight, mcScale, nVtx);
	      	      m_plots[KIN]["ProbePhi" + etaTag + IDTag].fill(probeMuTk.Phi(), emptyTk, weight, mcScale, nVtx);
		      
	      	      // Fill isolation plots for muons passign a given identification (programmable from cfg)
	      	      m_plots[ISO]["PhotonIso" + etaTag + IDTag].fill(probeMuon.photonIso, probeMuTk, weight, mcScale, nVtx);
	      	      m_plots[ISO]["ChHadIso" + etaTag + IDTag].fill(probeMuon.chargedHadronIso, probeMuTk, weight, mcScale, nVtx);
	      	      m_plots[ISO]["ChHadIsoPU" + etaTag + IDTag].fill(probeMuon.chargedHadronIsoPU, probeMuTk, weight, mcScale, nVtx);
	      	      m_plots[ISO]["NeutralIso" + etaTag + IDTag].fill(probeMuon.neutralHadronIso, probeMuTk, weight, mcScale, nVtx);
	      	      m_plots[ISO]["DBetaRelIso" + etaTag + IDTag].fill(probeMuon.isoPflow04, probeMuTk, weight, mcScale, nVtx);
	      	    }
	      	}
	    }
	}
    }
}


//Functions

double muon_pog::Plotter::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1 - eta2;
  double dphi = phi1 - phi2;
  while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

  return sqrt(deta*deta + dphi*dphi);
}


bool muon_pog::Plotter::hasGoodId(const muon_pog::Muon & muon, TString muId)
{
  //std::string & muId = leg == "tag" ? m_tnpConfig.tag_ID : m_tnpConfig.probe_ID;

  if (muId == "GLOBAL")      return muon.isGlobal == 1;
  else if (muId == "TRACKER")return muon.isTracker  == 1;
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
      if (object.filterTag.find(filter) != std::string::npos &&
	  deltaR(muTk.Eta(), muTk.Phi(), object.eta, object.phi) < m_tnpConfig.tag_hltDrCut)
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


// void muon_pog::setTProfY(TH1 &prof1, TProfile &prof2)
// {

//   Double_t max1 = prof1.GetBinContent(prof1.GetMaximumBin()) - prof1.GetBinError(prof1.GetMaximumBin());
//   Double_t max2 = prof2.GetBinContent(prof2.GetMaximumBin()) - prof2.GetBinError(prof2.GetMaximumBin());
//   Double_t min1 = prof1.GetBinContent(prof1.GetMinimumBin()) - prof1.GetBinError(prof1.GetMinimumBin());
//   Double_t min2 = prof2.GetBinContent(prof2.GetMinimumBin()) - prof2.GetBinError(prof2.GetMinimumBin());
  
//   Double_t max = max1 > max2 ? max1 : max2;
//   Double_t min = min1 < min2 ? min1 : min2;

//   prof1.GetYaxis()->SetRangeUser(min, max);

// }


void muon_pog::addUnderFlow(TH1 &hist)
{
  //Add UF
  hist.SetBinContent(1, hist.GetBinContent(0) + hist.GetBinContent(1));
  hist.SetBinError  (1, sqrt(hist.GetBinError(0)*hist.GetBinError(0) + hist.GetBinError(1)*hist.GetBinError(1)));
  hist.SetBinContent(0, 0); 
  hist.SetBinError  (0, 0);  

}


void muon_pog::addOverFlow(TH1 &hist)
{
  //Add OF		  
  Int_t lastBin = hist.GetNbinsX(); 
  hist.SetBinContent(lastBin, hist.GetBinContent(lastBin) + hist.GetBinContent(lastBin+1));
  hist.SetBinError  (lastBin, sqrt(hist.GetBinError(lastBin)*hist.GetBinError(lastBin) + hist.GetBinError(lastBin+1)*hist.GetBinError(lastBin+1))); 
  hist.SetBinContent(lastBin+1, 0) ; 
  hist.SetBinError  (lastBin+1, 0) ; 

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

  // system("mkdir -p " + outputDir + "/comparison/isolation/profiles");
  // system("mkdir -p " + outputDir + "/comparison/id_variables/profiles");
 
  // system("mkdir -p " + outputDir + "/comparison/isolation/profiles/no_ratio");
  // system("mkdir -p " + outputDir + "/comparison/id_variables/profiles/no_ratio");
 
  system("mkdir -p " + outputDir + "/comparison/timing/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/control/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/isolation/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/id_variables/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/kinematical_variables/no_ratio");

 
  outFile->cd("comparison");

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

      TString outputDirMap[5] {"/comparison/kinematical_variables", "/comparison/id_variables",
	  "/comparison/isolation", "/comparison/timing", "/comparison/control"};
      
      Plotter::HistoType plotType = plotTypeAndName.first;
      TString plotName = plotTypeAndName.second;
      outFile->cd(outputDirMap[plotType]);

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
	  TH1 * hData = 0;	  

	  int colorMap[5] {kGray+1, kRed, kAzure+7, kGreen+1, kOrange};
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
		    integralData = plotter.m_plots[Plotter::CONT]["99_invMassInRange"].plots().at(0)->Integral();
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
		      integralMC += (plotter.m_plots[Plotter::CONT]["99_invMassInRange"].plots().at(0)->Integral() *
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
		  hData->Sumw2();

		  //FP Add OF and UF for iso variables only
		  if (!plot->IsA()->InheritsFrom("TProfile"))
		    {
		      addOverFlow(*hData);
		      if(plotName.Contains("Iso"))
			addUnderFlow(*hData);
		    }
		  
		  Double_t stat = hData->GetEntries();	
		  leg->AddEntry(hData,Form(plotter.m_sampleConfig.sampleName+" [%10.0f  ]",stat),"LP");  
		}
	      else
		{
		  plot->SetFillColor(colorMap[iColor]);
		  plot->SetMarkerStyle(0);
		  plot->Sumw2();
		  
		  float scale = plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents;
		  TString option = "G";

		  if (!plot->IsA()->InheritsFrom("TProfile"))
		    {
		      //FP Add OF and UF for iso variables only
		      addOverFlow(*plot);
		      if(plotName.Contains("Iso"))
			addUnderFlow(*hData);

		      float scaleToData = integralData/integralMC;
		      plot->Scale(scale*scaleToData);
	     
		      hMc.Add(plot);
		      leg->AddEntry(plot, plotter.m_sampleConfig.sampleName, "LF"); 
		      iColor++;
		    }
		  else
		    {
		      // if (!pMc)
		      // 	{
		      pMc = (TProfile*)(plot->Clone("_pClone"));
		      pMc->Clear();
		      pMc->SetFillColor(0);
		      //pMc->SetErrorOption(option);
		      // pMc->Sumw2();
		      //pMc->Add(plot,scale);
		      pMc->Add(plot);
		      leg->AddEntry(pMc, "Weighted sum of MCs", "LP"); 
		    }
		  // else
		  // 	{
		  // 	  // pMc->Add(plot,scale);
		  
		  // 	}
		  //}
		}
	    }
	  
	  TString histoName = hData->GetName();

	  TCanvas *canvas = new TCanvas("c"+histoName, "c"+histoName, 500, 500);
	  canvas->cd();
	  
	  if(pMc){
	    //	    setTProfY(*hData, pMc);
	    
	    Double_t max1 = hData->GetBinContent(hData->GetMaximumBin()) - hData->GetBinError(hData->GetMaximumBin());
	    Double_t max2 = pMc->GetBinContent(pMc->GetMaximumBin()) - pMc->GetBinError(pMc->GetMaximumBin());
	    Double_t min1 = hData->GetBinContent(hData->GetMinimumBin()) - hData->GetBinError(hData->GetMinimumBin());
	    Double_t min2 = pMc->GetBinContent(pMc->GetMinimumBin()) - pMc->GetBinError(pMc->GetMinimumBin());
	    
	    Double_t max = max1 > max2 ? max1 : max2;
	    Double_t min = min1 < min2 ? min1 : min2;
	    
	    hData->GetYaxis()->SetRangeUser(min, max);
	    
	    hData->Draw();
	    pMc->Draw("same");
	  }
	  else{
	    canvas->SetLogy(1);
	    hData->Draw();
	    hMc.Draw("samehist");
	    hData->Draw("same");
	  }  
	 
	  leg->Draw();
	  
	  canvas->Update();
	  canvas->Write();
	  canvas->SaveAs(outputDir + outputDirMap[plotType] + "/no_ratio/c" + histoName + ".png");
	 	  
	  //Canvas with ratio plots
	  TCanvas *ratioCanvas = new TCanvas("rc"+histoName, "rc"+histoName, 500, 700);
	  ratioCanvas->Divide(1,2);
	  ratioCanvas->cd(1);
	  TPad *plotPad = (TPad*)ratioCanvas->GetPad(1);
	  plotPad->SetPad(0.,0.2,1.,1.);
	
	  hData->GetXaxis()->SetLabelSize(0.);
	  hData->GetXaxis()->SetTitleSize(0.);
	 		  
	  if(pMc){
	    //setTProfY(*hData, pMc);

	    Double_t max1 = hData->GetBinContent(hData->GetMaximumBin()) - hData->GetBinError(hData->GetMaximumBin());
	    Double_t max2 = pMc->GetBinContent(pMc->GetMaximumBin()) - pMc->GetBinError(pMc->GetMaximumBin());
	    Double_t min1 = hData->GetBinContent(hData->GetMinimumBin()) - hData->GetBinError(hData->GetMinimumBin());
	    Double_t min2 = pMc->GetBinContent(pMc->GetMinimumBin()) - pMc->GetBinError(pMc->GetMinimumBin());
	    
	    Double_t max = max1 > max2 ? max1 : max2;
	    Double_t min = min1 < min2 ? min1 : min2;
	    
	    hData->GetYaxis()->SetRangeUser(min, max);

	    hData->Draw();
	    pMc->Draw("same");
	  }
	  else{
	    plotPad->SetLogy(1);
	    hData->Draw();
	    hMc.Draw("samehist");
	    hData->Draw("same");
	  }
	  
	  leg->Draw();
	  
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
      
	  if(plotName == "03_nVertices" && !pMc)
	    {
	      Int_t nbins = hRatio->GetNbinsX();
	      
	      for (Int_t i=1;i<=nbins;i++)
		{
		  std::cout << " PUweight: " << hRatio->GetBinContent(i) << ", " << std::endl;
		}
	    }
      
	  ratioCanvas->Write();
	  ratioCanvas->SaveAs(outputDir+ outputDirMap[plotType] + "/rc" + histoName + ".png");
      
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
