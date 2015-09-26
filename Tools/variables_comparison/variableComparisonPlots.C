#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
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


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/include/zip.hpp>

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
  
    std::vector<TString> probe_ID;
    typename boost:result_of::zip< std::vector<TString>,
				   std::vector<TString> >::type probe_fEtaBins;
     
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

    enum HistoType { KIN=0, CONT, ID, ISO, TIME };
      
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, int nVtx, float weight);

    std::map<TString,std::map<int,TH1 *> > m_plots;
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
  
  void comparisonPlot(TFile *outFile, TString plotName,
		      std::vector<Plotter> & plotters,
		      TString &  outputDir);

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
	// for (Long64_t iEvent=0; iEvent<10000; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;
	  
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;
		  
	  //weight for Data to MC(DY_8 + TTbar_1 + tWtop + tWantitop) comparison 
	  float PUweight[60] = {0.,
				6.13641,
				3.96095,
				1.069,
				0.947777,
				0.856153,
				0.880403,
				0.991013,
				1.14579,
				1.29139,
				1.17103,
				1.25331,
				1.16784,
				1.3081,
				1.1052,
				1.09359,
				0.955982,
				0.875554,
				0.833031,
				0.611113,
				0.677309,
				0.496612,
				0.407461,
				0.357879,
				0.485327,
				0.282095,
				0.293168,
				0.203256,
				0.168494,
				0.0594542,
				0.417589,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.};	  
				 
	  //weight for Data to DY comparison 
	  //float PUweight[60] = {0, 6.23938, 3.93412, 1.0139, 0.92968, 0.852813, 0.86357, 0.991675, 1.13882, 1.29982, 1.17863, 1.25089, 1.17501, 1.29242, 1.11248, 1.08851, 0.952472, 0.873848, 0.842575, 0.620411, 0.6636, 0.50119, 0.404807, 0.369032, 0.483837, 0.266346, 0.308719, 0.207252, 0.175888, 0.0609816, 0.423387, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	  
	  //weight for Data to Data comparison
	  //float PUweight[60] = {0, 2.59252, 7.0705, 6.74054, 7.77755, 3.26155, 3.42985, 2.89984, 2.28284, 2.13837, 1.70971, 1.50287, 1.20367, 1.17027, 0.956894, 0.881644, 0.7398, 0.635878, 0.595542, 0.494059, 0.499049, 0.408914, 0.344204, 0.286958, 0.388219, 0.244787, 0.265499, 0.216043, 0.172834, 0.0540108, 0.360072, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	  
	  if(plotter.m_sampleConfig.sampleName.Contains("Data"))
	    weight *= ev->nVtx < 60 ? PUweight[ev->nVtx] : 0;
	  
	  plotter.fill(ev->muons, ev->hlt, ev->nVtx, weight);
	}
      
      delete ev;
      inputFile->Close();
    }
  
  muon_pog::comparisonPlots(plotters,outputFile,dirName);
  
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

      probe_ID = toArray(vt.second.get<std::string>("probe_muonID"));
      probe_fEtaBins = boost::zip(toArray(vt.second.get<std::string>("probe_fEtaMin"))
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

  outFile->mkdir(sampleTag+"timing");
  outFile->cd(sampleTag+"timing");

  m_plots[TIME]["STAmuonTime"]        = new TH1F("hSTAmuonTime_" + sampleTag ,"STAmuonTime; time (ns) ; # entries", 400, -200., 200.);
  m_plots[TIME]["STAmuonTimeBarrel"]  = new TH1F("hSTAmuonTimeBarrel_" + sampleTag ,"STAmuonTimeBarrel; time (ns) ; # entries", 400, -200., 200.);
  m_plots[TIME]["STAmuonTimeEndcap"]  = new TH1F("hSTAmuonTimeEndcap_" + sampleTag ,"STAmuonTimeEndcap; time (ns) ; # entries", 400, -200., 200.);

  m_plots[TIME]["UnbSTAmuonTime"]        = new TH1F("hUnbSTAmuonTime_" + sampleTag ,"UnbSTAmuonTime; time (ns) ; # entries", 400, -200., 200.);
  m_plots[TIME]["UnbSTAmuonTimeBarrel"]  = new TH1F("hUnbSTAmuonTimeBarrel_" + sampleTag ,"UnbSTAmuonTimeBarrel; time (ns) ; # entries", 400, -200., 200.);
  m_plots[TIME]["UnbSTAmuonTimeEndcap"]  = new TH1F("hUnbSTAmuonTimeEndcap_" + sampleTag ,"UnbSTAmuonTimeEndcap; time (ns) ; # entries", 400, -200., 200.);

  outFile->mkdir(sampleTag+"id_variables");
  outFile->cd(sampleTag+"id_variables");
  
  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
    {
      TString etaTag = "_fEtaMin" + fEtaBin.at(0) + "_fEtaMax" + fEtaBin.at(1);
      //Id var
      m_plots[ID]["NHitsGLB" + etaTag]        = new TH1F("hNHitsGLB_" + sampleTag + etaTag,"nHitsGLB; # hits ; # entries", 80, 0., 80.);
      m_plots[ID]["NHitsTRK" + etaTag]        = new TH1F("hNHitsTRK_" + sampleTag + etaTag,"nHitsTRK; # hits ; # entries", 40, 0., 40.);
      m_plots[ID]["NHitsSTA" + etaTag]        = new TH1F("hNHitsSTA_" + sampleTag + etaTag,"nHitsSTA; # hits ; # entries", 60, 0., 60.);
      
      m_plots[ID]["Chi2GLB" + etaTag]        = new TH1F("hChi2GLB_" + sampleTag + etaTag,"Chi2GLB; chi2/ndof ; # entries", 50, 0., 100.);
      m_plots[ID]["Chi2TRK" + etaTag]        = new TH1F("hChi2TRK_" + sampleTag + etaTag,"Chi2TRK; chi2/ndof ; # entries", 100, 0., 50.);
      
      m_plots[ID]["NMatchedStation" + etaTag]  = new TH1F("hNatchedStation_" + sampleTag + etaTag,"matchedStation; # stations ; # entries", 10, 0., 10.);
      m_plots[ID]["NMuonValidHitsGLB" + etaTag]= new TH1F("hNMuonValidHitsGLB_" + sampleTag + etaTag,"muonValidHitsGLB; # hits ; # entries", 60, 0., 60.);
      m_plots[ID]["PixelHitsTRK" + etaTag]    = new TH1F("hPixelHitsTRK_" + sampleTag + etaTag,"PixelHitsTRK; # hits ; # entries", 10, 0., 10.);
      m_plots[ID]["PixelLayersTRK" + etaTag]  = new TH1F("hPixelLayersTRK_" + sampleTag + etaTag,"PixelLayersTRK; # layers ; # entries", 10, 0., 10.);
      m_plots[ID]["TrackerLayersTRK" + etaTag]= new TH1F("hTrackerLayersTRK_" + sampleTag + etaTag,"TrackerLayersTRK; # layers ; # entries", 30, 0., 30.);
      
      m_plots[ID]["HitFractionTRK" + etaTag]  = new TH1F("hHitFractionTRK_" + sampleTag + etaTag,"HitFractionTRK; fraction ; # entries", 20, 0., 1.);
      m_plots[ID]["TrkStaChi2" + etaTag]      = new TH1F("hTrkStaChi2_" + sampleTag + etaTag,"TrkStaChi2; chi2 ; # entries", 50, 0., 100.);
      m_plots[ID]["TrkKink" + etaTag]         = new TH1F("hTrkKink_" + sampleTag + etaTag,"TrkKink; ; # entries", 48, 0., 1200.);
      m_plots[ID]["SegmentComp" + etaTag]     = new TH1F("hSegmentComp_" + sampleTag + etaTag,"segmentComp", 100, 0., 1.);
      
      m_plots[ID]["Dxy" + etaTag]             = new TH1F("Dxy_" + sampleTag + etaTag,"Transverse IP; dxy (cm); # entries", 100,-0.5,0.5);
      m_plots[ID]["Dz" + etaTag]              = new TH1F("Dz_" + sampleTag + etaTag,"Longitudinal IP; dz (cm); # entries", 200,-2.5,2.5);
    

  outFile->mkdir(sampleTag+"kinematical_variables");
  outFile->cd(sampleTag+"kinematical_variables");

  for ( auto & probe_ID : probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	    	  
	  m_plots[KIN]["probePt" + etaTag + IDTag]  = new TH1F("probePt_" + sampleTag + etaTag + IDTag, IDTag+" muons: p_{T} ; p_{T} (GeV) ; # entries", 75,0.,150.);
	  m_plots[KIN]["probeEta" + etaTag + IDTag] = new TH1F("probeEta_" + sampleTag + etaTag + IDTag, IDTag+" muons: #eta ; #eta ; # entries", 50,-2.5,2.5);
	  m_plots[KIN]["probePhi" + etaTag + IDTag] = new TH1F("probePhi_" + sampleTag + etaTag + IDTag, IDTag+" muons: #phi ; #phi ; # entries", 50,-TMath::Pi(),TMath::Pi());      
	}
    }

  outFile->mkdir(sampleTag+"isolation");
  outFile->cd(sampleTag+"isolation");

  for ( auto & probe_ID : probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	    	  
	  m_plots[ISO]["chHadIso" + etaTag + IDTag]    = new TH1F("chHadIso_" + sampleTag + etaTag + IDTag, IDTag+" muons: charged had. iso 0.4; muon relative isolation; # entries", 50,0.,5.);
	  m_plots[ISO]["chHadIsoPU" + etaTag + IDTag]  = new TH1F("chHadIsoPU_" + sampleTag + etaTag + IDTag, IDTag+" muons: PU Charged had. iso 0.4; muon relative isolation;# entries", 50,0.,5.);
	  m_plots[ISO]["photonIso" + etaTag + IDTag]   = new TH1F("photonIso_" + sampleTag + etaTag + IDTag, IDTag+" muons: photon iso 0.4; muon relative isolation; # entries", 50,0.,5.);
	  m_plots[ISO]["neutralIso" + etaTag + IDTag]  = new TH1F("neutralIso_" + sampleTag + etaTag + IDTag, IDTag+" muons: neutral had. iso 0.4; muon relative isolation; # entries", 50,0.,5.);
	  m_plots[ISO]["dBetaRelIso" + etaTag + IDTag] = new TH1F("dBetaRelIso_" + sampleTag + etaTag + IDTag, IDTag+" muons: PFIso 0.4; muon relative isolation; # entries", 50,0.,2.);
	}
    }

  outFile->mkdir(sampleTag+"/control");
  outFile->cd(sampleTag+"/control");
  
  m_plots[CONT]["invMass"] = new TH1F("invMass_" + sampleTag ,"Dimuon mass; mass (GeV); # entries", 100,0.,200.);
  m_plots[CONT]["dilepPt"] = new TH1F("dilepPt_" + sampleTag ,"Dilepton p_{T}; p_{T} (GeV); # entries", 100,0.,200.);
  m_plots[CONT]["nVertices"] = new TH1F("nVertices_" + sampleTag ,"Vertices; # vertices; # entries", 60,0.,60.);
  
  m_plots[CONT]["invMassInRange"] = new TH1F("invMassInRange_" + sampleTag ,"Dimuon mass; mass (GeV); # entries", 100,0.,200.);
  
  m_plots[CONT]["nProbesVsnTags"] = new TH2F("nProbesVsnTags_" + sampleTag ,"Probes vs Tags; # probes; # tags", 10,-0.5,9.,10,-0.5,9.);

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, int nVtx, float weight)
{

  //  TString sampleTag = m_sampleConfig.sampleName;

  //muon timing only
  for (auto & muon : muons)
    {
      if (muon.isStandAlone && ((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) ||
				(fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)))
	m_plots[TIME]["STAmuonTime"]->Fill(muon.muonTime,weight);
      
      if (muon.isStandAlone && fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20)
	m_plots[TIME]["STAmuonTimeBarrel"]->Fill(muon.muonTime,weight);
      
      if (muon.isStandAlone && fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11)
	m_plots[TIME]["STAmuonTimeEndcap"]->Fill(muon.muonTime,weight);
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

  std::vector<muon_pog::Muon *> tagMuons;
  
  for (auto & muon : muons)
    {
      if (muonTk(muon).Pt() > m_tnpConfig.tag_minPt &&
	  hasFilterMatch(muon,hlt) &&
	  hasGoodId(muon,m_tnpConfig.tag_ID) && 
	  muon.isoPflow04 < m_tnpConfig.tag_isoCut)
	tagMuons.push_back(&muon);
    }
  
  std::vector<muon_pog::Muon *> probeMuons;

  for (auto & muon : muons)
    {
      for (auto tagMuonPointer : tagMuons)
	{
	  muon_pog::Muon & tagMuon = *tagMuonPointer;
	  
	  if ( tagMuonPointer != &muon
	       chargeFromTrk(tagMuon) * chargeFromTrk(muon) == -1)
	    {
	      if(muon.isStandAlone) {
		if((fabs(muon.eta) < 1.2  && muon.nHitsStandAlone > 20) || (fabs(muon.eta) >= 1.2 && muon.nHitsStandAlone > 11)) 
		  m_plots[TIME]["UnbSTAmuonTime"]->Fill(muon.muonTime,weight);
		
		if(fabs(muon.eta) < 0.9  && muon.nHitsStandAlone > 20) 
		  m_plots[TIME]["UnbSTAmuonTimeBarrel"]->Fill(muon.muonTime,weight);
		
		if(fabs(muon.eta) > 1.2  && muon.nHitsStandAlone > 11) 
		  m_plots[TIME]["UnbSTAmuonTimeEndcap"]->Fill(muon.muonTime,weight);
	      }
	          
	      if((muon.isGlobal || muon.isTracker) &&  // CB minimal cuts on potental probe 
		 muon.pt > m_tnpConfig.probe_minPt)
		{
		  TLorentzVector tagMuTk(muonTk(tagMuon));
		  TLorentzVector muTk(muonTk(muon));
		  
		  Float_t mass = (tagMuTk+muTk).M();
		  
		  // CB Fill control plots
		  m_plots[CONT]["invMass"]->Fill(mass,weight);
		  if ( mass > m_tnpConfig.pair_minInvMass &&
		       mass < m_tnpConfig.pair_maxInvMass )
		    {
		      m_plots[CONT]["invMassInRange"]->Fill(mass,weight);
		      
		      Float_t dilepPt = (tagMuTk+muTk).Pt();
		      m_plots[CONT]["dilepPt"]->Fill(dilepPt,weight);

		      for (auto & fEtaBin : m_tnpConfig.probe_fEtaBins)
			{
			  TString etaTag = "_fEtaMin" + fEtaBin.at(0) + "_fEtaMax" + fEtaBin.at(1);
			}
		      
		      probeMuons.push_back(&muon);
		      continue; // CB If a muon is already a probe don't loo on other tags
		    }
		}
	    }
	}
    }
  
  m_plots[CONT]["nProbesVsnTags"]->Fill(tagMuons.size(),probeMuons.size());
  m_plots[CONT]["nVertices"]->Fill(nVtx,weight);

  
  for (auto probeMuonPointer : probeMuons)
    {
      muon_pog::Muon & probeMuon = *tagMuonPointer;
      			  
      //id var
      m_plots[ID]["NHitsGLB"  + etaTag]->Fill(probeMuon.nHitsGlobal ,weight);
      m_plots[ID]["NHitsTRK"  + etaTag]->Fill(probeMuon.nHitsTracker,weight);
      m_plots[ID]["NHitsSTA"  + etaTag]->Fill(probeMuon.nHitsStandAlone,weight);
      
      m_plots[ID]["Chi2GLB"  + etaTag]->Fill(probeMuon.glbNormChi2, weight);
      m_plots[ID]["Chi2TRK"  + etaTag]->Fill(probeMuon.trkNormChi2, weight);
      
      m_plots[ID]["MatchedStation"  + etaTag]->Fill(probeMuon.trkMuonMatchedStations, weight)  ;
      m_plots[ID]["MuonValidHitsGLB"  + etaTag]->Fill(probeMuon.glbMuonValidHits, weight);
      m_plots[ID]["PixelHitsTRK"  + etaTag]->Fill(probeMuon.trkPixelValidHits, weight);
      m_plots[ID]["PixelLayersTRK"  + etaTag]->Fill(probeMuon.trkPixelLayersWithMeas, weight); 
      m_plots[ID]["TrackerLayersTRK"  + etaTag]->Fill(probeMuon.trkTrackerLayersWithMeas, weight); 
      
      m_plots[ID]["HitFractionTRK"  + etaTag]->Fill(probeMuon.trkValidHitFrac, weight);   
      m_plots[ID]["TrkStaChi2"  + etaTag]->Fill(probeMuon.trkStaChi2, weight);       
      m_plots[ID]["TrkKink"  + etaTag]->Fill(probeMuon.trkKink, weight);          
      m_plots[ID]["SegmentComp"  + etaTag]->Fill(probeMuon.muSegmComp, weight);      
      
      m_plots[ID]["Dxy"  + etaTag]->Fill(probeMuon.dxy,weight);
      m_plots[ID]["Dz"  + etaTag]->Fill(probeMuon.dz,weight);			

      for (auto & probe_ID : m_tnpConfig.probe_IDs)
	{
	  TString IDTag = "_" + probe_ID;
	  
	  if(hasGoodId(probeMuon,probe_ID)) 
	    {
	      TLorentzVector probeMuTk(muonTk(probeMuon));
	  
	      for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
		{

		  if (fabs(probeMuTk.Eta()) > fEtaBin.at(0).Atof() &&
		      fabs(probeMuTk.Eta()) < fEtaBin.at(1).Atof() )
		    {
		      
		      TString etaTag = "_fEtaMin" + fEtaBin.at(0) + "_fEtaMax" + fEtaBin.at(1);

		      m_plots[KIN]["probePt" + etaTag + IDTag]->Fill(probeMuTk.Pt(),weight);
		      m_plots[KIN]["probeEta" + etaTag + IDTag]->Fill(probeMuTk.Eta(),weight);
		      m_plots[KIN]["probePhi" + etaTag + IDTag]->Fill(probeMuTk.Phi(),weight);
		      
		      // Fill isolation plots for muons passign a given identification (programmable from cfg)
		      m_plots[ISO]["photonIso" + etaTag + IDTag]->Fill(probeMuon.photonIso,weight);
		      m_plots[ISO]["chHadIso" + etaTag + IDTag]->Fill(probeMuon.chargedHadronIso,weight);
		      m_plots[ISO]["chHadIsoPU" + etaTag + IDTag]->Fill(probeMuon.chargedHadronIsoPU,weight);
		      m_plots[ISO]["neutralIso" + etaTag + IDTag]->Fill(probeMuon.neutralHadronIso,weight);
		      m_plots[ISO]["dBetaRelIso" + etaTag + IDTag]->Fill(probeMuon.isoPflow04,weight);
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

  outputFile->cd("/");
  outputFile->mkdir("comparison");
  
  outputFile->mkdir("comparison/control");
  outputFile->mkdir("comparison/isolation");
  outputFile->mkdir("comparison/id_variables");
  outputFile->mkdir("comparison/kinematical_variables");

  outputFile->cd("comparison");

  std::map<int,std::vector<TString> > plotTypesAndNames;
  for (auto & plotPairType : plotter.at(0).m_plots)
    {
      for (auto & plotPairName : plotPairType->second)
	{
	  plotTypesAndNames[plotPairType->first].push_back(plotPairName->first);
	}
 
  for (auto & plotTypeAndName : plotTypesNames)
    {

      TString outputDirMap[4] {"comparison/control", "comparison/isolation",
	  "comparison/id_variables", "comparison/kinematical_variables"};
      
      outputFile->cd(outputDirMap[plotTypeAndName->first]);

      TString plotName = plotTypeAndName->second;
      	   
      THStack hMc(plotName,"");
      TH1 * hData = 0;

      int colorMap[5] {kRed, kAzure+7, kGreen+1, kOrange, kGray+1};
      int iColor = 0;
  
      float integralData = 0;
      float integralMC   = 0;

      for (auto & plotter : plotters)
	{        
	  if(plotter.m_sampleConfig.sampleName.Contains("Data"))
	    {
	      if(plotName.Contains("muonTime"))
		{
		  TH1* plot = plotter.m_plots[plotName];
		  integralData = plot->Integral(plot->FindBin(-2.5),plot->FindBin(2.5)); 
		}
	      else
		integralData = plotter.m_plots["invMassInRange"]->Integral();
	    }
	  else
	    {
	      if(plotName.Contains("muonTime"))
		{
		  TH1* plot = plotter.m_plots[plotName];
		  integralMC += (plot->Integral(plot->FindBin(-2.5),plot->FindBin(2.5)) *
				 plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents);
		}
	      else
		{
		  integralMC += (plotter.m_plots["invMassInRange"]->Integral() *
				 plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents);
		}
	    }
	}
      
      for (auto & plotter : plotters)
	{        
	  if(plotter.m_sampleConfig.sampleName.Contains("Data"))
	    {
	      hData = plotter.m_plots[plotName];
	      hData->Sumw2();
	      leg->AddEntry(hData, plotter.m_sampleConfig.sampleName, "LF"); 
	    }
	  else
	    {
	      plotter.m_plots[plotName]->SetFilegColor(colorMap[iColor]);
	      plotter.m_plots[plotName]->SetMarkerStyle(0);
	      
	      TH1* plot = plotter.m_plots[plotName];
	      
	      float scale = plotter.m_sampleConfig.cSection/plotter.m_sampleConfig.nEvents;
	      plot->Scale(scale);
	      
	      float scaleToData = integralData/integralMC; 
	      plot->Scale(scaleToData);
	      
	      hMc.Add(plot);
	      leg->AddEntry(plot, plotter.m_sampleConfig.sampleName, "LF"); 
	      iColor++;
	    }
	}
      
      TCanvas *canvas = new TCanvas("c"+plotName, "c"+plotName, 500, 500);
      canvas->cd();
      canvas->SetLogy(1);
      
      hData->Draw();
      hMc.Draw("samehist");
      hData->Draw("same");
      
      //leg->AddEntry(hData, "#splitline{SingleMuon}{2015C, 50 ns}", "LF"); 
      //leg->AddEntry(hr, "#splitline{SingleMuon}{2015C, 25 ns}", "PE"); 
      leg->Draw();
      
      canvas->Write();
      canvas->SaveAs(outputDir + "/c" + plotName + ".png");
      
      //Canvas with ratio plots
      TCanvas *ratioCanvas = new TCanvas("rc"+plotName, "rc"+plotName, 500, 700);
      ratioCanvas->Divide(1,2);
      ratioCanvas->cd(1);
      TPad *plotPad = (TPad*)ratioCanvas->GetPad(1);
      plotPad->SetLogy(1);
      plotPad->SetPad(0.,0.2,1.,1.);
      hData->Draw();
      hMc.Draw("samehist");
      hData->Draw("same");
      leg->Draw();
      
      ratioCanvas->cd(2);
      TPad *ratioPad = (TPad*)ratioCanvas->GetPad(2);
      ratioPad->SetPad(0.,0.,1.,0.2);
      
      TH1 *hRatio = (TH1*)hData->Clone();
      
      hRatio->SetTitle(" ");
      hRatio->GetYaxis()->SetTitle("Data/MC");
      hRatio->GetYaxis()->SetRangeUser(-5,5);
      
      hRatio->Divide((TH1*)hMc.GetStack()->Last());
      hRatio->Draw();
      
      if(plotName == "nVertices")
	{
	  Int_t nbins = hRatio->GetNbinsX();
	  
	  for (Int_t i=1;i<=nbins;i++)
	    {
	      std::cout << " PUweight: " << hRatio->GetBinContent(i) << ", " << std::endl;
	    }
	}
      
      ratioCanvas->Write();
      ratioCanvas->SaveAs(outputDir + "/rc" + plotName + ".png");
      
      delete canvas;
      delete ratioCanvas;
    }
      
}
