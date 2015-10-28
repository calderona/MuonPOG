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
#include <iostream>
#include <fstream> 
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
    
    ~Observable() { m_plots.clear(); };

    // void fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Float_t mcScale, Int_t PV);
    void fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Int_t PV);
    
    std::vector<TH1 *> & plots() { return m_plots; }
    
  private :
    
    std::vector<TH1 *> m_plots;

  };

  class Plotter {

  public :

    enum HistoType { KIN=0, ID, ISO, TIME, CONT, EFF};
      
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
	//for (Long64_t iEvent=0; iEvent<50000; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;
	  
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;

	  //weight for dataD-Startup

	  // float PUweight[60] = {  0,
	  // 			  24.6564,
	  // 			  10.7846,
	  // 			  9.75522,
	  // 			  8.29381,
	  // 			  7.5551,
	  // 			  6.52192,
	  // 			  5.46865,
	  // 			  4.58941,
	  // 			  3.76588,
	  // 			  2.9334,
	  // 			  2.31009,
	  // 			  1.74043,
	  // 			  1.2886,
	  // 			  0.935587,
	  // 			  0.678059,
	  // 			  0.47279,
	  // 			  0.332189,
	  // 			  0.230285,
	  // 			  0.151753,
	  // 			  0.105131,
	  // 			  0.0694238,
	  // 			  0.0463378,
	  // 			  0.0312494,
	  // 			  0.019313,
	  // 			  0.0125001,
	  // 			  0.00766361,
	  // 			  0.00596304,
	  // 			  0.00421311,
	  // 			  0.0012315,
	  // 			  0.00091048,
	  // 			  0.000602633,
	  // 			  0,
	  // 			  0.00110585,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0};



	  //weight for dataD-DYasymptotic comparisons.
	  
	  // float PUweight[60] = {  0,
	  // 			  5.60605,
	  // 			  1.7979,
	  // 			  1.72712,
	  // 			  1.67129,
	  // 			  1.88794,
	  // 			  1.95294,
	  // 			  2.00093,
	  // 			  1.95066,
	  // 			  1.82526,
	  // 			  1.67957,
	  // 			  1.51496,
	  // 			  1.27748,
	  // 			  1.0463,
	  // 			  0.862064,
	  // 			  0.672516,
	  // 			  0.5231,
	  // 			  0.395426,
	  // 			  0.295751,
	  // 			  0.212734,
	  // 			  0.162864,
	  // 			  0.121367,
	  // 			  0.0900407,
	  // 			  0.0711253,
	  // 			  0.0487265,
	  // 			  0.0376877,
	  // 			  0.0286602,
	  // 			  0.0229082,
	  // 			  0.0220411,
	  // 			  0.00732292,
	  // 			  0.00630842,
	  // 			  0.00568248,
	  // 			  0,
	  // 			  0.0192158,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0,
	  // 			  0};			 

	  
	  //weight for dataC-dataD comparisons.
	  float PUweight[60] = {  0,
				  0.818782,
				  0.777041,
				  1.85606,
				  1.59566,
				  1.80562,
				  1.81222,
				  1.79909,
				  1.50801,
				  1.3662,
				  1.29328,
				  1.13194,
				  1.05819,
				  0.821376,
				  0.823951,
				  0.67575,
				  0.569473,
				  0.528223,
				  0.429261,
				  0.386063,
				  0.293675,
				  0.286689,
				  0.245489,
				  0.222443,
				  0.133406,
				  0.194772,
				  0.154812,
				  0.220177,
				  0.106648,
				  0,
				  0.0275221,
				  0.0275221,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0,
				  0};


	  //weight for Data-Mc comparison config4.ini 2015D vs All MC 
	  //float PUweight[60] = {0,6.68472,4.41435,3.49289,3.14212,3.173,2.94319,2.71174,2.44124,2.0717,1.7103,1.39812,1.11666,0.855547,0.618711,0.453196,0.328689,0.228592,0.166069,0.100075,0.070327,0.0567964,0.0405713,0.0308903,0.0186064,0.0208503,0.00264836,0.00779018,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


	  //weight for Data-Mc comparison config5.ini 2015D-v3 FullStat vs All MC 
	  //float PUweight[60] = {0, 2.83379, 1.78476, 1.68828, 1.68078, 1.89932, 1.94938, 1.98453, 1.95908, 1.85889, 1.67011, 1.50816, 1.29412,1.07878, 0.855702, 0.6922, 0.517513, 0.397616, 0.30073, 0.211495, 0.165432, 0.112824, 0.0852806, 0.0690589, 0.0472335, 0.0365517, 0.0262368, 0.0229718, 0.0228724, 0.00656045, 0.00546197,0.00366466, 0, 0.0131781, 0, 0, 0.0596201, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	  //weight for Data-Mc comparison configStartupDY.ini 2015D-v3 FullStat vs DY Startup 
	  // float PUweight[60] = {0, 23.5718, 11.623, 9.75004, 8.38437, 7.49407, 6.5513, 5.50713, 4.56761, 3.75242, 2.94015, 2.3024, 1.73567, 1.28708, 0.940598, 0.677024, 0.473512, 0.331202, 0.23009, 0.152347, 0.104783, 0.0694506, 0.046307, 0.031299, 0.0193667, 0.0124493, 0.00760724, 0.00600702, 0.00425913, 0.00123227, 0.00091714, 0.00060431, 0, 0.00110741, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	  
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
      m_plots.push_back(new TProfile("h" + hName + "VsEta_"      + sampleTag, hName + " vs #eta;    #eta;"        + xTitle, 24, -2.4, 2.4, min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPhi_"      + sampleTag, hName + " vs #phi;    #phi;"        + xTitle, 24, -TMath::Pi(),TMath::Pi(), min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPhiPlus_"  + sampleTag, hName + " vs #phi for #eta +;  #phi;"        + xTitle, 24, -TMath::Pi(),TMath::Pi(), min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPhiMinus_" + sampleTag, hName + " vs #phi for #eta -;  #phi;"        + xTitle, 24, -TMath::Pi(),TMath::Pi(), min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPt_"       + sampleTag, hName + " vs p_{T};   p_{T} (GeV);" + xTitle, 50,  0., 150., min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPV_"       + sampleTag, hName + " vs PV;      # of PV;"     + xTitle, 60,  0., 60., min, max));
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
      if (muonTk.Eta() > 0)
	((TProfile*)m_plots.at(3))->Fill(muonTk.Phi(), value, weight);
      else
	((TProfile*)m_plots.at(4))->Fill(muonTk.Phi(), value, weight);
      ((TProfile*)m_plots.at(5))->Fill(muonTk.Pt(),  value, weight);
      ((TProfile*)m_plots.at(6))->Fill(PV,           value, weight);
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

  m_plots[TIME]["UnbSTAmuonTime"] = muon_pog::Observable("UnbSTAmuonTime", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["UnbSTAmuonTimeBarrel"] = muon_pog::Observable("UnbSTAmuonTimeBarrel", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);
  m_plots[TIME]["UnbSTAmuonTimeEndcap"] = muon_pog::Observable("UnbSTAmuonTimeEndcap", sampleTag, "time (ns)", "# entries", 400, -200., 200., false);

  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
    {
      
      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;

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

      m_plots[ID]["qOverPtTrkSta" + etaTag]      = muon_pog::Observable("qOverPtTrkSta" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.03,0.03, true); 
      m_plots[ID]["qOverPtTrkStaOverPt" + etaTag]      = muon_pog::Observable("qOverPtTrkStaOverPt" + etaTag, sampleTag, "(q/p_{T}^{sta} - q/p_{T}^{trk})/q/p_{T}^{trk}", "# entries", 50,-5.,5., true); 
      m_plots[ID]["qOverPtTrkStaPlus" + etaTag]      = muon_pog::Observable("qOverPtTrkStaPlus" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.03,0.03, true); 
      m_plots[ID]["qOverPtTrkSta200" + etaTag]      = muon_pog::Observable("qOverPtTrkSta200" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 
      m_plots[ID]["qOverPtTrkSta200Plus" + etaTag]      = muon_pog::Observable("qOverPtTrkSta200Plus" + etaTag, sampleTag, "q/p_{T}^{sta} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 

      m_plots[ID]["qOverPtTrkGlb" + etaTag]      = muon_pog::Observable("qOverPtTrkGlb" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.03,0.03, true); 
      m_plots[ID]["qOverPtTrkGlbPlus" + etaTag]      = muon_pog::Observable("qOverPtTrkGlbPlus" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.03,0.03, true); 
      m_plots[ID]["qOverPtTrkGlb200" + etaTag]      = muon_pog::Observable("qOverPtTrkGlb200" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 
      m_plots[ID]["qOverPtTrkGlb200Plus" + etaTag]      = muon_pog::Observable("qOverPtTrkGlb200Plus" + etaTag, sampleTag, "q/p_{T}^{glb} - q/p_{T}^{trk}", "# entries", 50,-0.2,0.2, true); 

      outFile->cd(sampleTag+"/kinematical_variables");

      for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  m_plots[KIN]["ProbePt" + etaTag + IDTag]  = muon_pog::Observable("ProbePt" + etaTag + IDTag, sampleTag, "p_{T} (GeV)", "# entries", 75,0.,150., true);
	  m_plots[KIN]["ProbeEta" + etaTag + IDTag] = muon_pog::Observable("hProbeEta_" + etaTag + IDTag, sampleTag, "#eta", "# entries", 48,-2.4, 2.4, true);
	  m_plots[KIN]["ProbePhi" + etaTag + IDTag] = muon_pog::Observable("hProbePhi_" + etaTag + IDTag, sampleTag, "#phi", "# entries", 48,-TMath::Pi(),TMath::Pi(), true); 
	  m_plots[KIN]["goodMuMass" + etaTag + IDTag]  = muon_pog::Observable("goodMuMass" + etaTag + IDTag, sampleTag, "invariant mass (GeV)", "# entries", 30,85.,115., true);
	  m_plots[KIN]["goodMuMassPlus" + etaTag + IDTag]  = muon_pog::Observable("goodMuMassPlus" + etaTag + IDTag, sampleTag, "invariant mass (GeV)", "# entries", 20 ,86.5,96.5, true);
	  m_plots[KIN]["goodMuMassMinus" + etaTag + IDTag]  = muon_pog::Observable("goodMuMassMinus" + etaTag + IDTag, sampleTag, "invariant mass (GeV)", "# entries", 20,86.5,96.5, true);
   
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
  
  m_plots[CONT]["01_invMass"] = muon_pog::Observable("invMass", sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
  m_plots[CONT]["02_dilepPt"] = muon_pog::Observable("dilepPt", sampleTag ,"p_{T} (GeV)", "# entries", 100,0.,200., false);
  m_plots[CONT]["03_nVertices"] = muon_pog::Observable("nVertices", sampleTag ,"# vertices", "# entries", 60,0.,60., false);
  
  m_plots[CONT]["99_invMassInRange"] = muon_pog::Observable("invMassInRange", sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);

  
  //  std::cout << sampleTag << "  End: " << m_plots.size()<< std::endl;
  

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, int nVtx, float weight)
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
		  m_plots[TIME]["UnbSTAmuonTime"].fill(muon.muonTime,emptyTk,weight,nVtx);
		
		if(fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20) 
		  m_plots[TIME]["UnbSTAmuonTimeBarrel"].fill(muon.muonTime,emptyTk,weight,nVtx);
		
		if(fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11) 
		  m_plots[TIME]["UnbSTAmuonTimeEndcap"].fill(muon.muonTime,emptyTk,weight,nVtx);
	      }
   
	      // General Probe Muons	      
	      if((muon.isGlobal || muon.isTrackerArb) &&  // CB minimal cuts on potental probe 
		 muon.pt > m_tnpConfig.probe_minPt)
		{
		  TLorentzVector tagMuTk(muonTk(tagMuon));
		  TLorentzVector muTk(muonTk(muon));
		  
		  Float_t mass = (tagMuTk+muTk).M();
		  
		  // CB Fill control plots
		  m_plots[CONT]["01_invMass"].fill(mass,emptyTk,weight,nVtx);
		  if ( mass > m_tnpConfig.pair_minInvMass &&
		       mass < m_tnpConfig.pair_maxInvMass )
		    {
		      m_plots[CONT]["99_invMassInRange"].fill(mass,emptyTk,weight,nVtx);
		      
		      Float_t dilepPt = (tagMuTk+muTk).Pt();
		      m_plots[CONT]["02_dilepPt"].fill(dilepPt,emptyTk,weight,nVtx);
		      m_plots[CONT]["03_nVertices"].fill(nVtx,emptyTk,weight,nVtx);
		      
		      for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
			{
			  if (fabs(muTk.Eta()) > fEtaBin.first.Atof() &&
			      fabs(muTk.Eta()) < fEtaBin.second.Atof() )
			    {

			      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;
	
			      for (auto & probe_ID : m_tnpConfig.probe_IDs)
				{
				  TString IDTag = "_" + probe_ID;

				  if(hasGoodId(muon,probe_ID) && muon.isoPflow04 < 0.25 && mass > 60 && mass < 115)
				    {
				    m_plots[KIN]["goodMuMass" + etaTag + IDTag].fill(mass, muTk, weight, nVtx);
				    if (mass > 86.5 && mass < 96.5)
				      {
					if (chargeFromTrk(muon) > 0)
					  m_plots[KIN]["goodMuMassPlus" + etaTag + IDTag].fill(mass, muTk, weight, nVtx);
					else
					  m_plots[KIN]["goodMuMassMinus" + etaTag + IDTag].fill(mass, muTk, weight, nVtx);
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
	      if (probeMuon.isGlobal)
		{
		  m_plots[ID]["NHitsGLB"  + etaTag].fill(probeMuon.nHitsGlobal, probeMuTk ,weight,nVtx);
		  m_plots[ID]["Chi2GLB"  + etaTag].fill(probeMuon.glbNormChi2, probeMuTk, weight, nVtx);
		  m_plots[ID]["NMuonValidHitsGLB"  + etaTag].fill(probeMuon.glbMuonValidHits, probeMuTk, weight, nVtx);
		  m_plots[ID]["TrkStaChi2"  + etaTag].fill(probeMuon.trkStaChi2, probeMuTk, weight, nVtx);       
		  m_plots[ID]["TrkKink"  + etaTag].fill(probeMuon.trkKink, probeMuTk, weight, nVtx);          
		}

	      if (probeMuon.isTracker && probeMuon.isGlobal)    
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

	      if (probeMuon.isTracker)
		{
		  m_plots[ID]["NMatchedStation"  + etaTag].fill(probeMuon.trkMuonMatchedStations, probeMuTk, weight, nVtx);
		}


	      if (probeMuon.isGlobal)
		{
		  Float_t qOverPtTrk = probeMuon.charge_tracker / probeMuon.pt_tracker;
		  Float_t qOverPtSta = probeMuon.charge_standalone / probeMuon.pt_standalone;
		  Float_t qOverPtGlb = probeMuon.charge_global / probeMuon.pt_global;

		  m_plots[ID]["qOverPtTrkSta" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
		  m_plots[ID]["qOverPtTrkStaOverPt" + etaTag].fill((qOverPtSta - qOverPtTrk)/qOverPtTrk, probeMuTk, weight, nVtx);  
		  m_plots[ID]["qOverPtTrkGlb" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  

		  if (probeMuon.charge_tracker > 0 )
		    {
		      m_plots[ID]["qOverPtTrkGlbPlus" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  
		      m_plots[ID]["qOverPtTrkStaPlus" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
		    }
		  if (probeMuon.pt_tracker > 200)
		    {
		      m_plots[ID]["qOverPtTrkSta200" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
		      m_plots[ID]["qOverPtTrkGlb200" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  
		      if (probeMuon.charge_tracker > 0 )
			{
			  m_plots[ID]["qOverPtTrkSta200Plus" + etaTag].fill(qOverPtSta - qOverPtTrk, probeMuTk, weight, nVtx);  
			  m_plots[ID]["qOverPtTrkGlb200Plus" + etaTag].fill(qOverPtGlb - qOverPtTrk, probeMuTk, weight, nVtx);  
			}
		    }
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
  outFile->mkdir("comparison/efficiencies");
  
  system("mkdir -p " + outputDir + "/comparison/control/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/timing/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/isolation/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/id_variables/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/kinematical_variables/no_ratio");
  system("mkdir -p " + outputDir + "/comparison/efficiencies/no_ratio");
  
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

      TString outputDirMap[6] {"/comparison/kinematical_variables", "/comparison/id_variables", "/comparison/isolation", "/comparison/timing", "/comparison/control", "/comparison/efficiencies"};
      
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
	  TH1 *proj = 0;
	  Float_t errors[1000] = {0.};
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
		      if(plotName.Contains("Iso") || plotName.Contains("Dxy") || plotName.Contains("Dz"))
			addUnderFlow(*hData);
		    }
		  
		  Double_t stat = hData->GetEntries();	
		  leg->AddEntry(hData,Form(plotter.m_sampleConfig.sampleName+" [%10.0f  ]",stat),"LP");  
		}
	      else
		{
		  plot->SetFillColor(colorMap[iColor]);
		  plot->SetMarkerStyle(0);
		 
		  float scale = plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents;
		  
		  if (!plot->IsA()->InheritsFrom("TProfile"))
		    {
		      //FP Add OF and UF for iso variables only
		      addOverFlow(*plot);
		      if(plotName.Contains("Iso") || plotName.Contains("Dxy") || plotName.Contains("Dz"))
			addUnderFlow(*plot);

		      float scaleToData = integralData/integralMC;
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
			  
			  pMc->Add(plot,scale);
			  //for(int i = 0; i < pMc->GetNbinsX(); ++i){
			  // errors[i] =  sqrt(errors[i]*errors[i] + pMc->GetBinError(i)*pMc->GetBinError(i)*scale*scale);
			    // errors[i] =  sqrt(errors[i]*errors[i] + pMc->GetBinError(i)*pMc->GetBinError(i));
			  // }
			  leg->AddEntry(pMc, "Weighted sum of MCs", "LP"); 
			}
		      else
		  	{
			  pMc->Add(plot,scale);
			  //  for(int i = 0; i < pMc->GetNbinsX(); ++i){
			  //  errors[i] =  sqrt(errors[i]*errors[i] + plot->GetBinError(i)*plot->GetBinError(i)*scale*scale);
			  //errors[i] =  sqrt(errors[i]*errors[i] + plot->GetBinError(i)*plot->GetBinError(i));
			  //}
			}
		    }
		}
	    }
	  
	  TString histoName = hData->GetName();
	  
	  TCanvas *canvas = new TCanvas("c"+histoName, "c"+histoName, 500, 500);
	  canvas->cd();
	  
	  if(pMc){
	    Double_t max = hData->GetMaximum() > pMc->GetMaximum() ? hData->GetMaximum() : pMc->GetMaximum();
	    Double_t min = hData->GetMinimum() < pMc->GetMinimum() ? hData->GetMinimum() : pMc->GetMinimum();
	    
	    hData->GetYaxis()->SetRangeUser(0.8*min, 1.2*max);
	    hData->Draw();
	    
	    // proj = pMc->ProjectionX();
	    // proj->SetMarkerStyle(26);
	    // for(int i = 0; i < pMc->GetNbinsX(); ++i){
	    //   proj->SetBinError(i, errors[i]);
	    // }
	    // proj->SetOption("E");
	    // proj->Draw("same");
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
	    Double_t max = hData->GetMaximum() > pMc->GetMaximum() ? hData->GetMaximum() : pMc->GetMaximum();
	    Double_t min = hData->GetMinimum() < pMc->GetMinimum() ? hData->GetMinimum() : pMc->GetMinimum();
	    
	    hData->GetYaxis()->SetRangeUser(0.8*min, 1.2*max);
	    hData->Draw();
	      		
	    //proj->Draw("same");
	    pMc->Draw("same");
	  }
	  else{
	    plotPad->SetLogy(1);
	    hData->Draw();
	    hMc.Draw("samehist");
	    hData->Draw("same");
	  }
	  
	  leg->Draw();
	  
	  //Get Integral before/after ID cuts
	  std::vector<TString> fEtaBinMin = {"0.0","2.1"};
	  std::vector<TString> fEtaBinMax = {"2.4","2.4"};
	  
	  for (int i = 0; i < 2; ++i){
	    TString etaTag = "_fEtaMin" + fEtaBinMin[i] + "_fEtaMax" + fEtaBinMax[i];	    

	    if(plotName == "HitFractionTRK" + etaTag && !pMc)
	      {
		//the cut is > 0.8
		double x1 = hData->FindBin(0.8);
		double x2 = hData->GetNbinsX();
		double Data_HitFractionTRK_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
		double Data_HitFractionTRK_b = hData->Integral(1,x1)/hData->Integral(1,x2);
		double MC_HitFractionTRK_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_HitFractionTRK_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
		integrals << "Data: HitFractionTRK_a = " << Data_HitFractionTRK_a << "  HitFractionTRK_b = " << Data_HitFractionTRK_b << "  --> " << Data_HitFractionTRK_a+Data_HitFractionTRK_b << std::endl;
		integrals << "MC  : HitFractionTRK_a = " << MC_HitFractionTRK_a << "  HitFractionTRK_b = " << MC_HitFractionTRK_b << "  --> " << MC_HitFractionTRK_a+MC_HitFractionTRK_b << std::endl;
		integrals << std::endl;
	      }
	    
	    if(plotName == "Chi2GLB" + etaTag && !pMc) // chi2<3
	      {
		double x1 = hData->FindBin(3.);
		double x2 = hData->GetNbinsX();
		double Data_Chi2GLB_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
		double Data_Chi2GLB_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
		double MC_Chi2GLB_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_Chi2GLB_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
		integrals << "Data: Chi2GLB_a = " << Data_Chi2GLB_a << "  Chi2GLB_b = " << Data_Chi2GLB_b << "  --> " << Data_Chi2GLB_a+Data_Chi2GLB_b << std::endl;
		integrals << "MC  : Chi2GLB_a = " << MC_Chi2GLB_a << "  Chi2GLB_b = " << MC_Chi2GLB_b << "  --> " << MC_Chi2GLB_a+MC_Chi2GLB_b << std::endl;
		integrals << std::endl;
	      }
	    
	    
	    if(plotName == "TrkStaChi2" + etaTag && !pMc)
	      {
		double x1 = hData->FindBin(12.);
		double x2 = hData->GetNbinsX();
		double Data_TrkStaChi2_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
		double Data_TrkStaChi2_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
		double MC_TrkStaChi2_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_TrkStaChi2_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		integrals << std::endl;
		integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
		integrals << "Data: TrkStaChi2_a = " << Data_TrkStaChi2_a << "  TrkStaChi2_b = " << Data_TrkStaChi2_b << "  --> " << Data_TrkStaChi2_a+Data_TrkStaChi2_b << std::endl;
		integrals << "MC  : TrkStaChi2_a = " << MC_TrkStaChi2_a << "  TrkStaChi2_b = " << MC_TrkStaChi2_b << "  --> " << MC_TrkStaChi2_a+MC_TrkStaChi2_b << std::endl;
		integrals << std::endl;
		
	      }
	    
	    if(plotName == "TrkKink" + etaTag && !pMc)
	      {
		double x1 = hData->FindBin(20.);
		double x2 = hData->GetNbinsX();
		double Data_TrkKink_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
		double Data_TrkKink_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
		double MC_TrkKink_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_TrkKink_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
		integrals << "Data: TrkKink_a = " << Data_TrkKink_a << "  TrkKink_b = " << Data_TrkKink_b << "  --> " << Data_TrkKink_a+Data_TrkKink_b << std::endl;
		integrals << "MC  : TrkKink_a = " << MC_TrkKink_a << "  TrkKink_b = " << MC_TrkKink_b << "  --> " << MC_TrkKink_a+MC_TrkKink_b << std::endl;
		integrals << std::endl;
		
	      }
	    
	    if(plotName == "SegmentComp" + etaTag && !pMc)
	      {
		double x0  = hData->FindBin(.303);
		double x1  = hData->FindBin(.451);
		double x2  = hData->GetNbinsX();
		
		double Data_SegmentCompLoose_a = hData->Integral(x0+1,x2)/hData->Integral(1,x2);
		double Data_SegmentCompTight_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
		double Data_SegmentCompLoose_b = hData->Integral(1,x0)/hData->Integral(1,x2);
		double Data_SegmentCompTight_b = hData->Integral(1,x1)/hData->Integral(1,x2);
		double MC_SegmentCompLoose_a = ((TH1*)hMc.GetStack()->Last())->Integral(x0+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_SegmentCompTight_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_SegmentCompLoose_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x0)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_SegmentCompTight_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Medium ID " + etaTag + " ********" << std::endl;
		integrals << "Data: SegmentCompLoose_a = " << Data_SegmentCompLoose_a << "  SegmentCompLoose_b = " << Data_SegmentCompLoose_b << "  --> " << Data_SegmentCompLoose_a+Data_SegmentCompLoose_b << std::endl;
		integrals << "MC  : SegmentCompLoose_a = " << MC_SegmentCompLoose_a << "  SegmentCompLoose_b = " << MC_SegmentCompLoose_b << "  --> " << MC_SegmentCompLoose_a+MC_SegmentCompLoose_b << std::endl;
		integrals << "Data: SegmentCompTight_a = " << Data_SegmentCompTight_a << "  SegmentCompTight_b = " << Data_SegmentCompTight_b << "  --> " << Data_SegmentCompLoose_a+Data_SegmentCompLoose_b << std::endl;
		integrals << "MC  : SegmentCompTight_a = " << MC_SegmentCompTight_a << "  SegmentCompTight_b = " << MC_SegmentCompTight_b << "  --> " << MC_SegmentCompLoose_a+MC_SegmentCompLoose_b << std::endl;
		integrals << std::endl;
	      }
	    
	    //Tight ID
	    if(plotName == "Chi2GLB" + etaTag && !pMc)
	      {
		double x1 = hData->FindBin(10.);
		double x2 = hData->GetNbinsX();
		double Data_Chi2GLB_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
		double Data_Chi2GLB_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
		double MC_Chi2GLB_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_Chi2GLB_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
		integrals << "Data: Chi2GLB_a = " << Data_Chi2GLB_a << "  Chi2GLB_b = " << Data_Chi2GLB_b << "  --> " << Data_Chi2GLB_a+Data_Chi2GLB_b << std::endl;
		integrals << "MC  : Chi2GLB_a = " << MC_Chi2GLB_a << "  Chi2GLB_b = " << MC_Chi2GLB_b << "  --> " << MC_Chi2GLB_a+MC_Chi2GLB_b << std::endl;
		integrals << std::endl;
	      }
	    
	    if(plotName == "NMuonValidHitsGLB" + etaTag && !pMc)
	      {
		double x1 = hData->FindBin(0.);
		double x2 = hData->GetNbinsX();
		double Data_NMuonValidHitsGLB_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
		double Data_NMuonValidHitsGLB_b = hData->Integral(1,x1)/hData->Integral(1,x2);
		double MC_NMuonValidHitsGLB_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_NMuonValidHitsGLB_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
		integrals << "Data: NMuonValidHitsGLB_a = " << Data_NMuonValidHitsGLB_a << "  NMuonValidHitsGLB_b = " << Data_NMuonValidHitsGLB_b << "  --> " << Data_NMuonValidHitsGLB_a+Data_NMuonValidHitsGLB_b << std::endl;
		integrals << "MC  : NMuonValidHitsGLB_a = " << MC_NMuonValidHitsGLB_a << "  NMuonValidHitsGLB_b = " << MC_NMuonValidHitsGLB_b << "  --> " << MC_NMuonValidHitsGLB_a+MC_NMuonValidHitsGLB_b << std::endl;
		integrals << std::endl;
	      }
	    
	    if(plotName == "NMatchedStation" + etaTag && !pMc)
	      {
		double x1 = hData->FindBin(1.);
		double x2 = hData->GetNbinsX();
		double Data_NMatchedStation_a = hData->Integral(x1+1,x2)/hData->Integral(1,x2);
		double Data_NMatchedStation_b = hData->Integral(1,x1)/hData->Integral(1,x2);
		double MC_NMatchedStation_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1+1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_NMatchedStation_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
		integrals << "Data: NMatchedStation_a = " << Data_NMatchedStation_a << "  NMatchedStation_b = " << Data_NMatchedStation_b << "  --> " << Data_NMatchedStation_a+Data_NMatchedStation_b << std::endl;
		integrals << "MC  : NMatchedStation_a = " << MC_NMatchedStation_a << "  NMatchedStation_b = " << MC_NMatchedStation_b << "  --> " << MC_NMatchedStation_a+MC_NMatchedStation_b << std::endl;
		integrals << std::endl;
	      }
	    
	    if(plotName == "Dxy" + etaTag && !pMc)
	      {
		double mx1 = hData->FindBin(-0.2);
		double x1  = hData->FindBin(.2);
		double x2  = hData->GetNbinsX();
		double Data_Dxy_a = (hData->Integral(x1,x2) + hData->Integral(1,mx1))/hData->Integral(1,x2);  
		double Data_Dxy_b = hData->Integral(mx1+1,x1-1)/hData->Integral(1,x2);
		double MC_Dxy_a = (((TH1*)hMc.GetStack()->Last())->Integral(x1,x2) + ((TH1*)hMc.GetStack()->Last())->Integral(1,mx1))/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_Dxy_b = ((TH1*)hMc.GetStack()->Last())->Integral(mx1+1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
		integrals << "Data: Dxy_a = " << Data_Dxy_a << "  Dxy_b = " << Data_Dxy_b << "  --> " << Data_Dxy_a+Data_Dxy_b << std::endl;
		integrals << "MC  : Dxy_a = " << MC_Dxy_a << "  Dxy_b = " << MC_Dxy_b << "  --> " << MC_Dxy_a+MC_Dxy_b << std::endl;
		integrals << std::endl;
	      }
	    
	    if(plotName == "Dz" + etaTag && !pMc)
	      {
		double mx1 = hData->FindBin(-0.5);
		double x1 = hData->FindBin(.5);
		double x2 = hData->GetNbinsX();
		double Data_Dz_a = (hData->Integral(x1,x2) + hData->Integral(1,mx1))/hData->Integral(1,x2); //the cut is above th 
		double Data_Dz_b = hData->Integral(mx1+1,x1-1)/hData->Integral(1,x2);
		double MC_Dz_a = (((TH1*)hMc.GetStack()->Last())->Integral(x1,x2) + ((TH1*)hMc.GetStack()->Last())->Integral(1,mx1))/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_Dz_b = ((TH1*)hMc.GetStack()->Last())->Integral(mx1+1,x1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
		integrals << "Data: Dz_a = " << Data_Dz_a << "  Dz_b = " << Data_Dz_b << "  --> " << Data_Dz_a+Data_Dz_b << std::endl;
		integrals << "MC  : Dz_a = " << MC_Dz_a << "  Dz_b = " << MC_Dz_b << "  --> " << MC_Dz_a+MC_Dz_b << std::endl;
		integrals << std::endl;
	      }
	    
	    
	    if(plotName == "PixelHitsTRK" + etaTag && !pMc)
	      {
		double x1 = hData->FindBin(1);
		double x2 = hData->GetNbinsX();
		double Data_PixelHitsTRK_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
		double Data_PixelHitsTRK_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
		double MC_PixelHitsTRK_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_PixelHitsTRK_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
		integrals << "Data: PixelHitsTRK_a = " << Data_PixelHitsTRK_a << "  PixelHitsTRK_b = " << Data_PixelHitsTRK_b << "  --> " << Data_PixelHitsTRK_a+Data_PixelHitsTRK_b << std::endl;
		integrals << "MC  : PixelHitsTRK_a = " << MC_PixelHitsTRK_a << "  PixelHitsTRK_b = " << MC_PixelHitsTRK_b << "  --> " << MC_PixelHitsTRK_a+MC_PixelHitsTRK_b << std::endl;
		integrals << std::endl;
	      }
	    
	    if(plotName == "TrackerLayersTRK" + etaTag && !pMc)
	      {
		double x1 = hData->FindBin(6);
		double x2 = hData->GetNbinsX();
		double Data_TrackerLayersTRK_a = hData->Integral(x1,x2)/hData->Integral(1,x2);
		double Data_TrackerLayersTRK_b = hData->Integral(1,x1-1)/hData->Integral(1,x2);
		double MC_TrackerLayersTRK_a = ((TH1*)hMc.GetStack()->Last())->Integral(x1,x2)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		double MC_TrackerLayersTRK_b = ((TH1*)hMc.GetStack()->Last())->Integral(1,x1-1)/((TH1*)hMc.GetStack()->Last())->Integral(1,x2);
		
		integrals << std::endl;
		integrals << "********* Tight ID " + etaTag + " ********" << std::endl;
		integrals << "Data: TrackerLayersTRK_a = " << Data_TrackerLayersTRK_a << "  TrackerLayersTRK_b = " << Data_TrackerLayersTRK_b << "  --> " << Data_TrackerLayersTRK_a+Data_TrackerLayersTRK_b << std::endl;
		integrals << "MC  : TrackerLayersTRK_a = " << MC_TrackerLayersTRK_a << "  TrackerLayersTRK_b = " << MC_TrackerLayersTRK_b << "  --> " << MC_TrackerLayersTRK_a+MC_TrackerLayersTRK_b << std::endl;
		integrals << std::endl;
	      }
	  }
	  
	  
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
	  hRatio->Sumw2();	  

	  if(pMc){
	    ((TH1*)pMc)->Sumw2();
	    hRatio->Divide((TH1*)pMc);
	  }
	  else{
	    ((TH1*)hMc.GetStack()->Last())->Sumw2();
	    hRatio->Divide((TH1*)hMc.GetStack()->Last());
	  }

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
