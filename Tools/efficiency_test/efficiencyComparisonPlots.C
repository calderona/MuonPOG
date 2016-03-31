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
#include "TRegexp.h"
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
    Float_t nEvents;
    Bool_t applyReweighting;
    std::vector<int> runs;
        
    SampleConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    SampleConfig(boost::property_tree::ptree::value_type & vt); 
#endif

    ~SampleConfig() {};

  private:
    std::vector<int> toArray(const std::string & entries); 
    
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

    Float_t     gen_DrCut;
    
    std::string muon_trackType; // applies to both, tag and probe     
  
    Float_t probe_minPt;      
    Float_t probe_isoCut;

    std::vector<TString> probe_IDs;
    std::vector<std::pair<TString,TString> > probe_fEtaBins;
    std::vector<std::pair<TString,TString> > probe_ptBins;
     
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
    void fill(Float_t value, TLorentzVector & muonTk, Float_t weight, const muon_pog::Event & ev);
    
    std::vector<TH1 *> & plots() { return m_plots; }
    
  private :
    
    std::vector<TH1 *> m_plots;

  };


  class EffObservable {

  public :
    
    EffObservable() {};

    EffObservable(TString hName, TString sampleTag);
    
    ~EffObservable() { m_effs.clear(); };

    // void fill(Float_t value, TLorentzVector & muonTk, Float_t weight, Float_t mcScale, Int_t PV);
    void fill(Bool_t pass, TLorentzVector & muonTk, Float_t weight, const muon_pog::Event & ev);
    
    std::vector<TEfficiency *> & effs() { return m_effs; }
    
  private :
    
    std::vector<TEfficiency *> m_effs;

  };

  class Plotter {

  public :

    enum HistoType { VAR=0, CONT, TIGHT };
      
    Plotter(muon_pog::TagAndProbeConfig tnpConfig, muon_pog::SampleConfig & sampleConfig) :
      m_tnpConfig(tnpConfig) , m_sampleConfig(sampleConfig) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight);

    std::map<Plotter::HistoType, std::map<TString, muon_pog::Observable> > m_plots;
    std::map<Plotter::HistoType, std::map<TString, muon_pog::EffObservable> > m_effs;
    TagAndProbeConfig m_tnpConfig;
    SampleConfig m_sampleConfig;
        
  private :
   
    double deltaR(double eta1, double phi1, double eta2, double phi2);
    
    bool hasGoodId(const muon_pog::Muon & muon,
		   TString leg);
    bool hasMother(const muon_pog::GenParticle & gen,
		   Int_t pdgId);
    bool hasFilterMatch(const muon_pog::Muon & muon,
			const muon_pog::HLT  & hlt,
			std::string & filter);
    const muon_pog::GenParticle * hasGenMatch(const muon_pog::Muon & muon,
					const std::vector<muon_pog::GenParticle>  & gens,
					Int_t pdgId = 0);
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
 
  for (auto plotter : plotters)
    {

      TString fileName = plotter.m_sampleConfig.fileName;
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
      int nEntries = plotter.m_sampleConfig.nEvents > 0 ? plotter.m_sampleConfig.nEvents : tree->GetEntriesFast();
      std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

      int nFilteredEvents = 0;

      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;

	  if (iEvent % 25000 == 0 )
	    std::cout << "[" << argv[0] << "] processing event : " << iEvent << "\r" << std::flush;
	    
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;
	 
	  // CB to be fixed when kown what to do for MC !!!

	  //	  if(plotter.m_sampleConfig.applyReweighting==true)
	  //  weight *= ev->nVtx < 60 ? PUweight[ev->nVtx] : 0;
	  
	  plotter.fill(ev->muons, ev->hlt, (*ev), weight);
	}
      
      delete ev;
      inputFile->Close();
      std::cout << std::endl;
	   
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

      gen_DrCut  = vt.second.get<Float_t>("gen_DrCut");
      
      muon_trackType = vt.second.get<std::string>("muon_trackType");

      probe_minPt  = vt.second.get<Float_t>("probe_minPt");
      probe_isoCut = vt.second.get<Float_t>("probe_isoCut"); // CB for now just comb reliso dBeta R04
      probe_IDs    = toArray(vt.second.get<std::string>("probe_muonIDs"));
      probe_fEtaBins = toPairArray(toArray(vt.second.get<std::string>("probe_fEtaMin")),
				   toArray(vt.second.get<std::string>("probe_fEtaMax")));
      probe_ptBins = toPairArray(toArray(vt.second.get<std::string>("probe_ptMin")),
				  toArray(vt.second.get<std::string>("probe_ptMax")));
  
    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
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
      runs = toArray(vt.second.get<std::string>("runs"));
    }
  
  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[TagAndProbeConfig] Can't get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }
  
}

std::vector<int> muon_pog::SampleConfig::toArray(const std::string& entries)
{
  
  std::vector<int> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(atoi(item.c_str()));
  return result;

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
      m_plots.push_back(new TProfile("h" + hName + "VsEta_"      + sampleTag, hName + " vs #eta;    #eta;"          + xTitle, 24, -2.4, 2.4, min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPhi_"      + sampleTag, hName + " vs #phi;    #phi;"          + xTitle, 24, -TMath::Pi(),TMath::Pi(), min, max));     
      m_plots.push_back(new TProfile("h" + hName + "VsPt_"       + sampleTag, hName + " vs p_{T};   p_{T} (GeV);"   + xTitle, 50,  0., 150., min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsPV_"       + sampleTag, hName + " vs PV;      # of PV;"       + xTitle, 60,  0., 60., min, max));
      m_plots.push_back(new TProfile("h" + hName + "VsInstLumi_" + sampleTag, hName + " vs Inst. Lumi. ; Inst. Lumi. [10E30];"       + xTitle, 50,  0., 5000., min, max));
      
    }
}

void muon_pog::Observable::fill(Float_t value, TLorentzVector & muonTk, Float_t weight, const muon_pog::Event & ev)
{
  
  m_plots.at(0)->Fill(value, weight);
  // CB fill here other plots vs kin variables
  if (m_plots.size() > 1)
    {

      Int_t nVtx     = ev.nVtx;
      Int_t instLumi = ev.instLumi;
      
      ((TProfile*)m_plots.at(1))->Fill(muonTk.Eta(), value, weight);
      ((TProfile*)m_plots.at(2))->Fill(muonTk.Phi(), value, weight);
      ((TProfile*)m_plots.at(3))->Fill(muonTk.Pt(),  value, weight);
      ((TProfile*)m_plots.at(4))->Fill(nVtx,         value, weight);
      ((TProfile*)m_plots.at(5))->Fill(instLumi,     value, weight);

    }
}

muon_pog::EffObservable::EffObservable(TString hName, TString sampleTag)
{
  const double etaBins[18] = {-2.4,-2.1,-1.6,-1.2,-1.05,-0.9,-0.6,-0.3,-0.2,
   			      0.2, 0.3, 0.6, 0.9, 1.05, 1.2, 1.6, 2.1, 2.4};

  // const double etaBins[8] = {-2.4, -2.1, -1.2, -0.9, 0.9, 1.2, 2.1, 2.4};

  m_effs.push_back(new TEfficiency("h" + hName + "VsEta_"      + sampleTag, hName + " vs #eta;    #eta;"          , 17, etaBins));
  m_effs.push_back(new TEfficiency("h" + hName + "VsPhi_"      + sampleTag, hName + " vs #phi;    #phi;"          , 24, -TMath::Pi(),TMath::Pi()));     
  m_effs.push_back(new TEfficiency("h" + hName + "VsPt_"       + sampleTag, hName + " vs p_{T};   p_{T} (GeV);"   , 50,  0., 150.));
  m_effs.push_back(new TEfficiency("h" + hName + "VsPV_"       + sampleTag, hName + " vs PV;      # of PV;"       , 60,  0., 60.));
  m_effs.push_back(new TEfficiency("h" + hName + "VsInstLumi_" + sampleTag, hName + " vs Inst. Lumi. ; Inst. Lumi. [10E30];" , 50,  0., 5000.));
  m_effs.push_back(new TEfficiency("h" + hName + "OneBin_"     + sampleTag, hName + " one bin ;;" , 1,  0.5, 1.5));
      
}

void muon_pog::EffObservable::fill(Bool_t pass, TLorentzVector & muonTk, Float_t weight, const muon_pog::Event & ev)
{
  
  Int_t nVtx     = ev.nVtx;
  Int_t instLumi = ev.instLumi;
  
  m_effs.at(0)->FillWeighted(pass, weight, muonTk.Eta());
  m_effs.at(1)->FillWeighted(pass, weight, muonTk.Phi());
  m_effs.at(2)->FillWeighted(pass, weight, muonTk.Pt());
  m_effs.at(3)->FillWeighted(pass, weight, nVtx);
  m_effs.at(4)->FillWeighted(pass, weight, instLumi);
  m_effs.at(5)->FillWeighted(pass, weight, 1.);
  
}


void muon_pog::Plotter::book(TFile *outFile)
{

  TString sampleTag = m_sampleConfig.sampleName;
  
  outFile->cd("/");
  outFile->mkdir(sampleTag);
  outFile->cd(sampleTag);

  outFile->mkdir(sampleTag+"/efficiencies");      
  outFile->mkdir(sampleTag+"/kinematical_variables");
  outFile->mkdir(sampleTag+"/control");
 
  //Tight ID sequence plots

  for (auto ptBin : m_tnpConfig.probe_ptBins)
    {
      for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
	{

	  TString ptTag = "_ptMin" + ptBin.first + "_ptMax" + ptBin.second;
	  TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;
      
	  outFile->cd(sampleTag+"/control");
      
	  m_plots[CONT]["01_invMass" + etaTag + ptTag] = muon_pog::Observable("invMass" + etaTag + ptTag, sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
	  m_plots[CONT]["01_invMassFromTau" + etaTag + ptTag]   = muon_pog::Observable("invMassFromTau" + etaTag + ptTag, sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
	  m_plots[CONT]["01_invMassFromGamma" + etaTag + ptTag] = muon_pog::Observable("invMassFromGamma" + etaTag + ptTag, sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
	  m_plots[CONT]["01_invMassFromKaon" + etaTag + ptTag]  = muon_pog::Observable("invMassFromKaon" + etaTag + ptTag, sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
	  m_plots[CONT]["01_invMassFromMu" + etaTag + ptTag]    = muon_pog::Observable("invMassFromMu" + etaTag + ptTag, sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);
	  m_plots[CONT]["01_invMassFromZ" + etaTag + ptTag]     = muon_pog::Observable("invMassFromZ" + etaTag + ptTag, sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);

	  m_plots[CONT]["10_motherId" + etaTag + ptTag]    = muon_pog::Observable("motherId" + etaTag + ptTag, sampleTag ,"muon mother pdgId", "# entries", 1000,-0.5,999.5, false);
	  m_plots[CONT]["11_motherId50" + etaTag + ptTag]  = muon_pog::Observable("motherId50" + etaTag + ptTag, sampleTag ,"muon mother pdgId", "# entries", 50,-0.5,49.5, false);
	  
	  outFile->cd(sampleTag+"/kinematical_variables");
      
	  for ( auto & probe_ID : m_tnpConfig.probe_IDs)
	    {
	      TString IDTag = "_" + probe_ID;
	  
	      m_plots[VAR]["ProbePt" + etaTag + ptTag + IDTag]  = muon_pog::Observable("ProbePt" + etaTag + ptTag + IDTag, sampleTag, "p_{T} (GeV)", "# entries", 150,0.,150., false);
	      m_plots[VAR]["ProbeEta" + etaTag + ptTag + IDTag] = muon_pog::Observable("hProbeEta_" + etaTag + ptTag + IDTag, sampleTag, "#eta", "# entries", 48,-2.4, 2.4, false);
	      m_plots[VAR]["ProbePhi" + etaTag + ptTag + IDTag] = muon_pog::Observable("hProbePhi_" + etaTag + ptTag + IDTag, sampleTag, "#phi", "# entries", 48,-TMath::Pi(),TMath::Pi(), false); 

	      m_plots[VAR]["ProbeFromTauPt" + etaTag + ptTag + IDTag]  = muon_pog::Observable("ProbeFromTauPt" + etaTag + ptTag + IDTag, sampleTag, "p_{T} (GeV)", "# entries", 150,0.,150., false);
	      m_plots[VAR]["ProbeFromMuPt" + etaTag + ptTag + IDTag]  = muon_pog::Observable("ProbeFromMuPt" + etaTag + ptTag + IDTag, sampleTag, "p_{T} (GeV)", "# entries", 150,0.,150., false);
	      m_plots[VAR]["ProbeFromZPt" + etaTag + ptTag + IDTag]  = muon_pog::Observable("ProbeFromZPt" + etaTag + ptTag + IDTag, sampleTag, "p_{T} (GeV)", "# entries", 150,0.,150., false);

	      m_plots[CONT]["110_iso" + etaTag + ptTag + IDTag] = muon_pog::Observable("iso" + etaTag + ptTag + IDTag, sampleTag ,"PF combined isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["110_isoFromTau" + etaTag + ptTag + IDTag] = muon_pog::Observable("isoFromTau" + etaTag + ptTag + IDTag, sampleTag ,"PF combined isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["110_isoFromMu" + etaTag + ptTag + IDTag]  = muon_pog::Observable("isoFromMu" + etaTag + ptTag + IDTag,  sampleTag ,"PF combined isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["110_isoFromZ" + etaTag + ptTag + IDTag]   = muon_pog::Observable("isoFromZ" + etaTag + ptTag + IDTag,   sampleTag ,"PF combined isolation", "# entries", 100,0.,20., false);

	      m_plots[CONT]["111_chHadIso" + etaTag + ptTag + IDTag] = muon_pog::Observable("chHadIso" + etaTag + ptTag + IDTag, sampleTag ,"Ch. hadron isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["111_chHadIsoFromTau" + etaTag + ptTag + IDTag] = muon_pog::Observable("chHadIsoFromTau" + etaTag + ptTag + IDTag, sampleTag ,"Ch. hadron isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["111_chHadIsoFromMu" + etaTag + ptTag + IDTag]  = muon_pog::Observable("chHadIsoFromMu" + etaTag + ptTag + IDTag,  sampleTag ,"Ch. hadron isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["111_chHadIsoFromZ" + etaTag + ptTag + IDTag]   = muon_pog::Observable("chHadIsoFromZ" + etaTag + ptTag + IDTag,   sampleTag ,"Ch. hadron isolation", "# entries", 100,0.,20., false);

	      m_plots[CONT]["112_neuHadIso" + etaTag + ptTag + IDTag] = muon_pog::Observable("neuHadIso" + etaTag + ptTag + IDTag, sampleTag ,"Neu. hadron isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["112_neuHadIsoFromTau" + etaTag + ptTag + IDTag] = muon_pog::Observable("neuHadIsoFromTau" + etaTag + ptTag + IDTag, sampleTag ,"Neu. hadron isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["112_neuHadIsoFromMu" + etaTag + ptTag + IDTag]  = muon_pog::Observable("neuHadIsoFromMu" + etaTag + ptTag + IDTag,  sampleTag ,"Neu. hadron isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["112_neuHadIsoFromZ" + etaTag + ptTag + IDTag]   = muon_pog::Observable("neuHadIsoFromZ" + etaTag + ptTag + IDTag,   sampleTag ,"Neu. hadron isolation", "# entries", 100,0.,20., false);

	      m_plots[CONT]["113_phIso" + etaTag + ptTag + IDTag] = muon_pog::Observable("phIso" + etaTag + ptTag + IDTag, sampleTag ,"Photon isolation", "# entries", 100,0.,20., false);  
	      m_plots[CONT]["113_phIsoFromTau" + etaTag + ptTag + IDTag] = muon_pog::Observable("phIsoFromTau" + etaTag + ptTag + IDTag, sampleTag ,"Photon isolation", "# entries", 100,0.,20., false); 
	      m_plots[CONT]["113_phIsoFromMu" + etaTag + ptTag + IDTag]  = muon_pog::Observable("phIsoFromMu" + etaTag + ptTag + IDTag,  sampleTag ,"Photon isolation", "# entries", 100,0.,20., false); 
	      m_plots[CONT]["113_phIsoFromZ" + etaTag + ptTag + IDTag]   = muon_pog::Observable("phIsoFromZ" + etaTag + ptTag + IDTag,   sampleTag ,"Photon isolation", "# entries", 100,0.,20., false);

	    }

	  outFile->cd(sampleTag+"/efficiencies");
	  
	  m_effs[TIGHT]["Tight" + etaTag + ptTag] = muon_pog::EffObservable("Tight" + etaTag + ptTag, sampleTag);
	  m_effs[TIGHT]["TightFromTau" + etaTag + ptTag] = muon_pog::EffObservable("TightFromTau" + etaTag + ptTag, sampleTag);
	  m_effs[TIGHT]["TightFromMu" + etaTag + ptTag]  = muon_pog::EffObservable("TightFromMu" + etaTag + ptTag, sampleTag);
	  m_effs[TIGHT]["TightFromZ" + etaTag + ptTag]   = muon_pog::EffObservable("TightFromZ" + etaTag + ptTag, sampleTag);

	}
    }  

  outFile->cd(sampleTag+"/control");
  
  m_plots[CONT]["02_invMassInRange"] = muon_pog::Observable("invMassInRange", sampleTag ,"mass (GeV)", "# entries", 100,0.,200., false);  
  m_plots[CONT]["03_dilepPt"] = muon_pog::Observable("dilepPt", sampleTag ,"p_{T} (GeV)", "# entries", 100,0.,200., false);
  m_plots[CONT]["04_nVertices"] = muon_pog::Observable("nVertices", sampleTag ,"# vertices", "# entries", 60,0.,60., false);

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, const muon_pog::Event & ev, float weight)
{

  bool isGoodRun = false;

  for (auto run :m_sampleConfig.runs)
    {
      if (run == 0 || run == ev.runNumber)
	{
	  isGoodRun = true;
	  break;
	}
    }

  if (!isGoodRun) return;
  
  TLorentzVector emptyTk;
    
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
      if (hasGenMatch(muon,ev.genParticles) &&
	  muonTk(muon).Pt() > m_tnpConfig.tag_minPt &&
	  hasFilterMatch(muon,hlt,m_tnpConfig.tag_hltFilter) &&
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

	      // General Probe Muons	      
	      if((muon.isGlobal || muon.isTrackerArb) &&  // CB : minimal cuts on potental probe 
		 muon.pt > m_tnpConfig.probe_minPt    &&  // FOR EFFICIENCIES INCLUDES ISO!
		 hasGenMatch(muon,ev.genParticles))
		{
		  TLorentzVector tagMuTk(muonTk(tagMuon));
		  TLorentzVector muTk(muonTk(muon));
		  
		  Float_t mass = (tagMuTk+muTk).M();

		  Float_t dilepPt = (tagMuTk+muTk).Pt();
		  m_plots[CONT]["03_dilepPt"].fill(dilepPt,emptyTk,weight,ev);
		  m_plots[CONT]["04_nVertices"].fill(ev.nVtx,emptyTk,weight,ev);

		  const muon_pog::GenParticle * gen = hasGenMatch(muon,ev.genParticles);
		  
		  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
		    {
		      if (fabs(muTk.Eta()) > fEtaBin.first.Atof() &&
			  fabs(muTk.Eta()) < fEtaBin.second.Atof() )
			{

			  for (auto ptBin : m_tnpConfig.probe_ptBins)
			    {
			      if (muTk.Pt() > ptBin.first.Atof() &&
				  muTk.Pt() < ptBin.second.Atof() )
				{
	      
				  TString ptTag = "_ptMin" + ptBin.first + "_ptMax" + ptBin.second;		  	  
				  TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;

				  for (auto motherId : gen->mothers)
				    {
				      m_plots[CONT]["10_motherId" + etaTag + ptTag].fill(fabs(motherId),emptyTk,weight,ev);
				      m_plots[CONT]["11_motherId50" + etaTag + ptTag].fill(fabs(motherId),emptyTk,weight,ev);

				      if ((fabs(motherId)) == 23 && hasGenMatch(muon,ev.genParticles,13))
					std::cout << "CAZZO1" << std::endl;
				      if ((fabs(motherId)) == 13 && hasGenMatch(muon,ev.genParticles,23))
					std::cout << "CAZZO2" << std::endl;
				    }
				  
				  // CB Fill control plots
				  m_plots[CONT]["01_invMass" +  etaTag +  ptTag].fill(mass,emptyTk,weight,ev);
				  if (hasGenMatch(muon,ev.genParticles,15))
				    m_plots[CONT]["01_invMassFromTau" +  etaTag +  ptTag].fill(mass,emptyTk,weight,ev);
				  else if (hasGenMatch(muon,ev.genParticles,22))
				    m_plots[CONT]["01_invMassFromGamma" +  etaTag +  ptTag].fill(mass,emptyTk,weight,ev);
				  else if (hasGenMatch(muon,ev.genParticles,411))
				    m_plots[CONT]["01_invMassFromKaon" +  etaTag +  ptTag].fill(mass,emptyTk,weight,ev);
				  else if (hasGenMatch(muon,ev.genParticles,13))
				    m_plots[CONT]["01_invMassFromMu" +  etaTag +  ptTag].fill(mass,emptyTk,weight,ev);
				  else if (hasGenMatch(muon,ev.genParticles,23))
				    m_plots[CONT]["01_invMassFromZ" +  etaTag +  ptTag].fill(mass,emptyTk,weight,ev);
				  
				}
			    }
			}
		    }

		  if ( mass > m_tnpConfig.pair_minInvMass &&
		       mass < m_tnpConfig.pair_maxInvMass )
		    {
		      m_plots[CONT]["02_invMassInRange"].fill(mass,emptyTk,weight,ev);
		      
		      probeMuons.push_back(&muon);
		      
		      continue; // CB If a muon is already a probe don't loo on other tags
		    }
		}
	    }
	}
    }
  
  for (auto probeMuonPointer : probeMuons)
    {
      const muon_pog::Muon & probeMuon = *probeMuonPointer;
  
      TLorentzVector probeMuTk(muonTk(probeMuon));

      Bool_t tight = probeMuon.isPF && probeMuon.isGlobal && probeMuon.glbMuonValidHits > 0 
	&& probeMuon.trkMuonMatchedStations > 1 && probeMuon.trkTrackerLayersWithMeas > 5 
	&& probeMuon.trkPixelValidHits > 0 && probeMuon.glbNormChi2 < 10. 
	&& fabs(probeMuon.dzBest) < 0.5 && fabs(probeMuon.dxyBest) < 0.2 
	&& probeMuon.isoPflow04 < m_tnpConfig.probe_isoCut;

      for (auto ptBin : m_tnpConfig.probe_ptBins)
	{

	  for (auto fEtaBin : m_tnpConfig.probe_fEtaBins)
	    {
	      
	      if (probeMuTk.Pt() > ptBin.first.Atof() &&
		  probeMuTk.Pt() < ptBin.second.Atof() )
		{
	  
		  if (fabs(probeMuTk.Eta()) > fEtaBin.first.Atof() &&
		      fabs(probeMuTk.Eta()) < fEtaBin.second.Atof() )
		    {
		      
		      TString ptTag = "_ptMin" + ptBin.first + "_ptMax" + ptBin.second;
		      TString etaTag = "_fEtaMin" + fEtaBin.first + "_fEtaMax" + fEtaBin.second;
		      
		      m_effs[TIGHT]["Tight" + etaTag + ptTag].fill(tight,probeMuTk,weight,ev);
		      if (hasGenMatch(probeMuon,ev.genParticles,15))
			{
			  m_effs[TIGHT]["TightFromTau" + etaTag + ptTag].fill(tight,probeMuTk,weight,ev);
			}
		      else if (hasGenMatch(probeMuon,ev.genParticles,13))
			{
			  m_effs[TIGHT]["TightFromMu" + etaTag + ptTag].fill(tight,probeMuTk,weight,ev);
			}
		      else if (hasGenMatch(probeMuon,ev.genParticles,23))
			{
			  m_effs[TIGHT]["TightFromZ" + etaTag + ptTag].fill(tight,probeMuTk,weight,ev);
			}

		      for (auto & probe_ID : m_tnpConfig.probe_IDs)
			{
			  TString IDTag = "_" + probe_ID;
			  
			  Float_t pfIso = probeMuon.photonIso + probeMuon.chargedHadronIso + probeMuon.neutralHadronIso;
			  
			  if(hasGoodId(probeMuon,probe_ID)) 
			    {	 
			      m_plots[VAR]["ProbePt" + etaTag + ptTag + IDTag].fill(probeMuTk.Pt(), probeMuTk, weight, ev);
			      m_plots[CONT]["110_iso" + etaTag + ptTag + IDTag].fill(pfIso,emptyTk,weight, ev);
			      m_plots[CONT]["111_chHadIso" + etaTag + ptTag + IDTag].fill(probeMuon.chargedHadronIso,emptyTk,weight, ev);
			      m_plots[CONT]["112_neuHadIso" + etaTag + ptTag + IDTag].fill(probeMuon.neutralHadronIso,emptyTk,weight,ev);
			      m_plots[CONT]["113_phIso" + etaTag + ptTag + IDTag].fill(probeMuon.photonIso,emptyTk,weight,ev);

			      if (hasGenMatch(probeMuon,ev.genParticles,15))
				{
				  m_plots[VAR]["ProbeFromTauPt" + etaTag + ptTag + IDTag].fill(probeMuTk.Pt(), probeMuTk, weight, ev);
				  m_plots[CONT]["110_isoFromTau" + etaTag + ptTag + IDTag].fill(pfIso,emptyTk,weight,ev);
				  m_plots[CONT]["111_chHadIsoFromTau" + etaTag + ptTag + IDTag].fill(probeMuon.chargedHadronIso,emptyTk,weight, ev);
				  m_plots[CONT]["112_neuHadIsoFromTau" + etaTag + ptTag + IDTag].fill(probeMuon.neutralHadronIso,emptyTk,weight,ev);
				  m_plots[CONT]["113_phIsoFromTau" + etaTag + ptTag + IDTag].fill(probeMuon.photonIso,emptyTk,weight,ev);
				}
			      else if (hasGenMatch(probeMuon,ev.genParticles,13))
				{
				  m_plots[VAR]["ProbeFromMuPt" + etaTag + ptTag + IDTag].fill(probeMuTk.Pt(), probeMuTk, weight, ev);
				  m_plots[CONT]["110_isoFromMu" + etaTag + ptTag + IDTag].fill(pfIso,emptyTk,weight,ev);
				  m_plots[CONT]["111_chHadIsoFromMu" + etaTag + ptTag + IDTag].fill(probeMuon.chargedHadronIso,emptyTk,weight, ev);
				  m_plots[CONT]["112_neuHadIsoFromMu" + etaTag + ptTag + IDTag].fill(probeMuon.neutralHadronIso,emptyTk,weight,ev);
				  m_plots[CONT]["113_phIsoFromMu" + etaTag + ptTag + IDTag].fill(probeMuon.photonIso,emptyTk,weight,ev);
				}
			      else if (hasGenMatch(probeMuon,ev.genParticles,23))
				{
				  m_plots[VAR]["ProbeFromZPt" + etaTag + ptTag + IDTag].fill(probeMuTk.Pt(), probeMuTk, weight, ev);
				  m_plots[CONT]["110_isoFromZ" + etaTag + ptTag + IDTag].fill(pfIso,emptyTk,weight,ev);
				  m_plots[CONT]["111_chHadIsoFromZ" + etaTag + ptTag + IDTag].fill(probeMuon.chargedHadronIso,emptyTk,weight, ev);
				  m_plots[CONT]["112_neuHadIsoFromZ" + etaTag + ptTag + IDTag].fill(probeMuon.neutralHadronIso,emptyTk,weight,ev);
				  m_plots[CONT]["113_phIsoFromZ" + etaTag + ptTag + IDTag].fill(probeMuon.photonIso,emptyTk,weight,ev);

				}
			      m_plots[VAR]["ProbeEta" + etaTag + ptTag + IDTag].fill(probeMuTk.Eta(), probeMuTk, weight, ev);
			      m_plots[VAR]["ProbePhi" + etaTag + ptTag + IDTag].fill(probeMuTk.Phi(), probeMuTk, weight, ev);
			    }
			}
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

bool muon_pog::Plotter::hasMother(const muon_pog::GenParticle & gen, Int_t pdgId)
{

  for (auto motherId : gen.mothers)
    {
      if (abs(motherId) == pdgId)
	{
	  return true;
	}
    }

  return false;
 
}


bool muon_pog::Plotter::hasFilterMatch(const muon_pog::Muon & muon,
				       const muon_pog::HLT  & hlt,
				       std::string & filter)
{
  TLorentzVector muTk = muonTk(muon);

  for (auto object : hlt.objects)
    {
      if (object.filterTag.find(filter) != std::string::npos &&
	  deltaR(muTk.Eta(), muTk.Phi(), object.eta, object.phi) < m_tnpConfig.tag_hltDrCut)
	return true;
    }

  return false;
}

const muon_pog::GenParticle * muon_pog::Plotter::hasGenMatch(const muon_pog::Muon & muon,
						       const std::vector<muon_pog::GenParticle>  & gens,
						       Int_t pdgId)
{
  TLorentzVector muTk = muonTk(muon);


  const muon_pog::GenParticle * bestGen = 0;
  Float_t bestDr = 999.;

  for (auto & gen : gens)
    {
      if (fabs(gen.pdgId) == 13)
	{
	  Float_t dr = deltaR(muTk.Eta(), muTk.Phi(), gen.eta, gen.phi);
	  if (dr < m_tnpConfig.gen_DrCut && dr< bestDr)
	    {
	      bestGen = &gen;
	      bestDr = dr;
	    }
	}
    }

  return ((pdgId == 0) || hasMother(*bestGen,pdgId)) ? bestGen : 0;
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
  outFile->mkdir("comparison/variables");
  outFile->mkdir("comparison/efficiencies_tight");
  
  outFile->cd("comparison");

  system("mkdir -p " + outputDir + "/comparison/control/");
  system("mkdir -p " + outputDir + "/comparison/variables/");
  system("mkdir -p " + outputDir + "/comparison/efficiencies_tight/");
  

  
  std::vector<std::pair<Plotter::HistoType,TString> > plotTypesAndNames;
  for (auto & plotPairType : plotters.at(0).m_plots)
    {
      for (auto & plotPairName : plotPairType.second)
	{
	  if ( !plotPairName.first.Contains("From")) 
	    plotTypesAndNames.push_back(std::make_pair(plotPairType.first, plotPairName.first));
	}
    }

  for (auto & plotPairType : plotters.at(0).m_effs)
    {
      for (auto & plotPairName : plotPairType.second)
	{
	  if ( !plotPairName.first.Contains("From")) 
	    plotTypesAndNames.push_back(std::make_pair(plotPairType.first, plotPairName.first));
	}
    }

  for (auto & plotTypeAndName : plotTypesAndNames)
    {

      TString outputDirMap[3] {"/comparison/variables/", "/comparison/control/",
	  "/comparison/efficiencies_tight/"};
      
      Plotter::HistoType plotType = plotTypeAndName.first;
      TString observableName = plotTypeAndName.second;
      std::string base(observableName.Data());
      
      std::regex binRegEx("fEtaMin[0-9\\.]+\\_fEtaMax[0-9\\.]+\\_ptMin[0-9\\.]+\\_ptMax[0-9\\.]+");
      std::smatch binMatch;

      std::regex_search(base ,binMatch, binRegEx);
      std::string binTag = binMatch.str();

      outFile->cd(outputDirMap[plotType]);

      Bool_t isEff = (plotType != Plotter::HistoType::VAR &&
		      plotType != Plotter::HistoType::CONT);
      
      std::vector<TH1 *>::size_type nPlots =  isEff ?
	plotters.at(0).m_effs[plotType][observableName].effs().size() :
	plotters.at(0).m_plots[plotType][observableName].plots().size();;
      
      for (std::vector<TH1 *>::size_type iPlot=0; iPlot < nPlots; ++iPlot)
	{
	  TString plotName = isEff ? 
	    plotters.at(0).m_effs[plotType][observableName].effs().at(iPlot)->GetName() :
	    plotters.at(0).m_plots[plotType][observableName].plots().at(iPlot)->GetName();

	  TLegend *leg = new TLegend(0.45,0.7,0.95,0.95);
	  leg->SetBorderSize(0);
	  leg->SetLineWidth(0);
	  leg->SetFillColor(0);
	  leg->SetFillStyle(0);	  

	  int colorMap[5]  {kGray+2, kGreen+2, kOrange +7, kAzure+2, kRed};
	  int markerMap[5] {20, 21, 22, 23, 24};

	  TCanvas *canvas = new TCanvas("c"+plotName, "c"+plotName, 500, 500);

	  int iColor = 0;

	  auto plotter = plotters.begin();
	  
	  for (;plotter != plotters.end(); 
	       ++plotter, ++iColor )
	    {
	      if (isEff)
		{

		  TEfficiency * plot = plotter->m_effs[plotType][observableName].effs().at(iPlot);
	      
		  leg->AddEntry(plot,Form(plotter->m_sampleConfig.sampleName+" [%10.0f  ]",plot->GetTotalHistogram()->GetEntries()),"LP");  

		  plot->SetLineColor(colorMap[iColor]);
		  plot->SetFillColor(colorMap[iColor]);
		  plot->SetMarkerColor(colorMap[iColor]);
		  plot->SetMarkerStyle(markerMap[iColor]);
		  canvas->cd();

		  plot->Draw(iColor == 0 ? "" : "same");

		  // TEfficiency * effFrom = plotter->m_effs[plotType][observableName.ReplaceAll("Barrel","BarrelFromTau")].effs().at(iPlot);
		  // effFrom->SetLineColor(colorMap[iColor+1]);
		  // effFrom->SetFillColor(colorMap[iColor+1]);
		  // effFrom->SetMarkerColor(colorMap[iColor+1]);
		  // effFrom->SetMarkerStyle(markerMap[iColor+1]);
		  // effFrom->Draw("same");

		  TEfficiency * effFrom = plotter->m_effs[plotType][observableName.ReplaceAll("Tight","TightFromMu")].effs().at(iPlot);
		  effFrom->SetLineColor(colorMap[iColor+2]);
		  effFrom->SetFillColor(colorMap[iColor+2]);
		  effFrom->SetMarkerColor(colorMap[iColor+2]);
		  effFrom->SetMarkerStyle(markerMap[iColor+2]);
		  effFrom->Draw("same");
		  
		  effFrom = plotter->m_effs[plotType][observableName.ReplaceAll("TightFromMu","TightFromZ")].effs().at(iPlot);
		  effFrom->SetLineColor(colorMap[iColor+3]);
		  effFrom->SetFillColor(colorMap[iColor+3]);
		  effFrom->SetMarkerColor(colorMap[iColor+3]);
		  effFrom->SetMarkerStyle(markerMap[iColor+3]);
		  effFrom->Draw("same");
		  
		  observableName.ReplaceAll("TightFromZ","Tight");
		      
		  canvas->Update();
		  if (plot->GetPaintedGraph())
		    plot->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.5,1.05);

		}
	      else
		{
		  TH1* plot = plotter->m_plots[plotType][observableName].plots().at(iPlot);
	      
		  leg->AddEntry(plot,Form(plotter->m_sampleConfig.sampleName+" [%10.0f  ]",plot->GetEntries()),"LP");  

		  if(!plot->IsA()->InheritsFrom("TProfile"))
		    {
		      plot->Sumw2();		  
		      if (observableName.Contains("iso") ||
			  observableName.Contains("Iso")  )
			{ 
			  plot->Scale(1./plot->Integral());
			  canvas->SetLogy();
			  addOverFlow(*plot);		     
			}
		      //addUnderFlow(*plot);

		    }
		  
		  plot->SetLineColor(colorMap[iColor]);
		  plot->SetFillColor(colorMap[iColor]);
		  plot->SetMarkerColor(colorMap[iColor]);
		  plot->SetMarkerStyle(markerMap[iColor]);
		  
		  canvas->cd();
		  plot->Draw(iColor == 0 ? "" : "same");
		  
		  if (observableName.Contains("invMass_"))
		    {
		      
		      THStack *hs = new THStack("hs" + observableName,"Stacked invMass");
		      TH1* histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("invMass","invMassFromTau")].plots().at(iPlot);
		      histoFrom->SetFillColor(colorMap[iColor+1]);
		      hs->Add(histoFrom);
		      
		      histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("invMassFromTau","invMassFromMu")].plots().at(iPlot);
		      
		      histoFrom->SetFillColor(colorMap[iColor+2]);
		      hs->Add(histoFrom);
		      
		      histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("invMassFromMu","invMassFromZ")].plots().at(iPlot);
		      
		      histoFrom->SetFillColor(colorMap[iColor+3]);
		      hs->Add(histoFrom);
		      
		      hs->Draw("same");
		      
		    }
		  if (observableName.Contains("iso") ||
		      observableName.Contains("Iso") )
		    {
		      
		      TH1 * histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("so","soFromTau")].plots().at(iPlot);
		      histoFrom->Scale(1./ histoFrom->Integral());
		      histoFrom->SetLineColor(colorMap[iColor+1]);
		      histoFrom->SetMarkerColor(colorMap[iColor+1]);
		      histoFrom->SetMarkerStyle(markerMap[iColor+1]);
		      histoFrom->Draw("samep");

		      histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("soFromTau","soFromMu")].plots().at(iPlot);
		      histoFrom->Scale(1./ histoFrom->Integral());
		      histoFrom->SetLineColor(colorMap[iColor+2]);
		      histoFrom->SetMarkerColor(colorMap[iColor+2]);
		      histoFrom->SetMarkerStyle(markerMap[iColor+2]);
		      histoFrom->Draw("samep");

		      histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("soFromMu","soFromZ")].plots().at(iPlot);
		      histoFrom->Scale(1./ histoFrom->Integral());
		      histoFrom->SetLineColor(colorMap[iColor+3]);
		      histoFrom->SetMarkerColor(colorMap[iColor+3]);
		      histoFrom->SetMarkerStyle(markerMap[iColor+3]);
		      histoFrom->Draw("samep");
		      
		    }

		  if (observableName.Contains("ProbePt"))
		    {
		     
		      THStack *hs = new THStack("hs" + observableName,"Stacked ProbePt");
		      TH1* histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("ProbePt","ProbeFromTauPt")].plots().at(iPlot);
		      histoFrom->SetFillColor(colorMap[iColor+1]);
		      hs->Add(histoFrom);
		      
		      histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("ProbeFromTauPt","ProbeFromMuPt")].plots().at(iPlot);
		      
		      histoFrom->SetFillColor(colorMap[iColor+2]);
		      hs->Add(histoFrom);
		      
		      histoFrom = plotter->m_plots[plotType][observableName.ReplaceAll("ProbeFromMuPt","ProbeFromZPt")].plots().at(iPlot);
		      
		      histoFrom->SetFillColor(colorMap[iColor+3]);
		      hs->Add(histoFrom);
		      
		      hs->Draw("same");
		      
		    }
		  
		}
	    }
      
	  canvas->cd();	  
	  leg->Draw();
	  
	  canvas->Update();
	  canvas->Write();
	  system("mkdir -p " + outputDir + outputDirMap[plotType] + binTag);
	  canvas->SaveAs(outputDir + outputDirMap[plotType] + binTag + "/c" + plotName + ".png");

	  delete canvas;

	  // if(plotName == "05_nVertices" && !pMc)
	  //   {
	  //     Int_t nbins = hRatio->GetNbinsX();
	      
	  //     for (Int_t i=1;i<=nbins;i++)
	  // 	{
	  // 	  std::cout << " PUweight: " << hRatio->GetBinContent(i) << ", " << std::endl;
	  // 	}
	  //   }
	        
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
