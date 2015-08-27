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
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt, float weight);

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
	{
	  if (tree->LoadTree(iEvent)<0) break;
	  
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;

	  plotter.fill(ev->muons, ev->hlt, weight);
	  
	}
      
      delete ev;
      delete evBranch;
      
      inputFile->Close();
      
    }

  outputFile->cd("/");
  outputFile->mkdir("comparison");
  outputFile->cd("comparison");

  muon_pog::comparisonPlot(outputFile,"invMass",plotters);
  muon_pog::comparisonPlot(outputFile,"dilepPt",plotters);

  std::vector<TString>::const_iterator fEtaMinIt  = tnpConfig.probe_fEtaMin.begin();
  std::vector<TString>::const_iterator fEtaMinEnd = tnpConfig.probe_fEtaMin.end();

  std::vector<TString>::const_iterator fEtaMaxIt  = tnpConfig.probe_fEtaMax.begin();
  std::vector<TString>::const_iterator fEtaMaxEnd = tnpConfig.probe_fEtaMax.end();
  
  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
    {
      TString etaTag = "_fEtaMin" + (*fEtaMinIt) + "_fEtaMax" + (*fEtaMaxIt);
      muon_pog::comparisonPlot(outputFile,"probePt" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"probeEta" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"probePhi" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"probeDxy" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"probeDz" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"chHadIso" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"photonIso" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"neutralIso" + etaTag,plotters);
      muon_pog::comparisonPlot(outputFile,"dBetaRelIso" + etaTag,plotters);
    }
  
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
  
  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
    {
         
      TString etaTag = "_fEtaMin" + (*fEtaMinIt) + "_fEtaMax" + (*fEtaMaxIt);
      m_plots["probePt" + etaTag]  = new TH1F("probePt_" + sampleTag + etaTag," ; # entries; muon p_[T] ", 75,0.,150.);
      m_plots["probeEta" + etaTag] = new TH1F("probeEta_" + sampleTag + etaTag," ; # entries; muon #eta ", 50,-2.5,2.5);
      m_plots["probePhi" + etaTag] = new TH1F("probePhi_" + sampleTag + etaTag," ; # entries; muon #phi ", 50,-TMath::Pi(),TMath::Pi());
      m_plots["probeDxy" + etaTag] = new TH1F("probeDxy_" + sampleTag + etaTag," ; # entries; muon #phi ", 100,-0.5,0.5);
      m_plots["probeDz" + etaTag]  = new TH1F("probeDz_" + sampleTag + etaTag," ; # entries; muon #phi ", 200,-2.5,2.5);
      m_plots["chHadIso" + etaTag]    = new TH1F("chHadIso_" + sampleTag + etaTag," ; # entries; muon relative isolation", 50,0.,5.);
      m_plots["photonIso" + etaTag]   = new TH1F("photonIso_" + sampleTag + etaTag," ; # entries; muon relative isolation", 50,0.,5.);
      m_plots["neutralIso" + etaTag]  = new TH1F("neutralIso_" + sampleTag + etaTag," ; # entries; muon relative isolation", 50,0.,5.);
      m_plots["dBetaRelIso" + etaTag] = new TH1F("dBetaRelIso_" + sampleTag + etaTag," ; # entries; muon relative isolation", 50,0.,2.);

    }

  outFile->mkdir(sampleTag+"/control");
  outFile->cd(sampleTag+"/control");

  m_plots["invMass"] = new TH1F("invMass_" + sampleTag ,"invMass", 100,0.,200.);
  m_plots["dilepPt"] = new TH1F("dilepPt_" + sampleTag ,"dilepPt", 100,0.,200.);

  m_plots["nProbesVsnTags"] = new TH2F("nProbesVsnTags_" + sampleTag ,"invMass", 10,-0.5,9.,10,-0.5,9.);

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			     const muon_pog::HLT & hlt, float weight)
{

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
      if (hasGoodId(muon,"tag") && hasFilterMatch(muon,hlt) &&
	  muonTk(muon).Pt() > m_tnpConfig.tag_minPt   &&
	  muon.isoPflow04 < m_tnpConfig.tag_isoCut)
	tagMuons.push_back(muon);
    }
  
  std::vector<muon_pog::Muon> probeMuons;

  for (auto & muon : muons)
    {
      for (auto & tagMuon : tagMuons)
	{
	  if ( tagMuon.eta != muon.eta &&
	       tagMuon.phi != muon.phi &&
	       tagMuon.pt  != muon.pt  &&
	       chargeFromTrk(tagMuon) * chargeFromTrk(muon) == -1 &&    
	       (muon.isGlobal || muon.isTracker) ) // CB minimal cuts on potental probe 
	    {
	      
	      TLorentzVector tagMuTk = muonTk(tagMuon);
	      TLorentzVector muTk    = muonTk(muon);
	      
	      Float_t mass = (tagMuTk+muTk).M();

	      // CB Fill control plots
	      m_plots["invMass"]->Fill(mass,weight);
	      if ( mass > m_tnpConfig.pair_minInvMass &&
		   mass < m_tnpConfig.pair_maxInvMass )
		{
		  Float_t dilepPt = (tagMuTk+muTk).Pt();
		  m_plots["dilepPt"]->Fill(dilepPt,weight);
	      	  probeMuons.push_back(muon);
		  continue; // CB If a muon is already a probe don't loo on other tags
		}
	    }
	}
    }

  m_plots["nProbesVsnTags"]->Fill(tagMuons.size(),probeMuons.size());
  
  for (auto & probeMuon : probeMuons)
    {

      TLorentzVector probeMuTk = muonTk(probeMuon);

      std::vector<TString>::const_iterator fEtaMinIt  = m_tnpConfig.probe_fEtaMin.begin();
      std::vector<TString>::const_iterator fEtaMinEnd = m_tnpConfig.probe_fEtaMin.end();

      std::vector<TString>::const_iterator fEtaMaxIt  = m_tnpConfig.probe_fEtaMax.begin();
      std::vector<TString>::const_iterator fEtaMaxEnd = m_tnpConfig.probe_fEtaMax.end();
  
      for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
	{
	  
	  if (fabs(probeMuTk.Eta()) > fEtaMinIt->Atof() &&
	      fabs(probeMuTk.Eta()) < fEtaMaxIt->Atof() )
	    {
	      
	      TString etaTag = "_fEtaMin" + TString((*fEtaMinIt)) + "_fEtaMax" + TString((*fEtaMaxIt));
	      m_plots["probePt" + etaTag]->Fill(probeMuTk.Pt(),weight);
	      m_plots["probeEta" + etaTag]->Fill(probeMuTk.Eta(),weight);
	      m_plots["probePhi" + etaTag]->Fill(probeMuTk.Phi(),weight);
	      
	      m_plots["probeDxy" + etaTag]->Fill(probeMuon.dxy,weight);
	      m_plots["probeDz" + etaTag]->Fill(probeMuon.dz,weight);
	      
	      if(hasGoodId(probeMuon,"probe")) 
		{
		  // Fill isolation plots for muons passign a given identification (programmable from cfg)
		  m_plots["chHadIso" + etaTag]->Fill(probeMuon.isoPflow04,weight);
		  m_plots["photonIso" + etaTag]->Fill(probeMuon.isoPflow04,weight);
		  m_plots["neutralIso" + etaTag]->Fill(probeMuon.isoPflow04,weight);
		  m_plots["dBetaRelIso" + etaTag]->Fill(probeMuon.isoPflow04,weight);
		}
								    
	    }
	}
    }

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
      if (object.filterTag.find(filter) != std::string::npos &&
	  sqrt((muTk.Eta() - object.eta) * (muTk.Eta() - object.eta) +
	       (muTk.Phi() - object.phi) * (muTk.Phi() - object.phi))
	  < m_tnpConfig.tag_hltDrCut )
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

void muon_pog::comparisonPlot(TFile *outFile,TString plotName,
			      std::vector<muon_pog::Plotter> & plotters)
{

  THStack hMc(plotName,"");
  TH1 * hData = 0;

  int colorMap[5] { kOrange+7, kAzure+7, kGreen+1, kOrange, kGray+1};

  int iColor = 0;

  float integral  = 0;
  float totalXSec = 0;
  
  for (auto plotter : plotters)
    {
      
      if(std::string(plotter.m_sampleConfig.sampleName.Data()).find("Data") != std::string::npos)
	{
	  integral = plotter.m_plots[plotName]->Integral();
	}
      else
	{
	  totalXSec += plotter.m_sampleConfig.cSection;
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
	  float scale = plotter.m_sampleConfig.cSection / totalXSec *
	    integral / plot->Integral();
	  plot->Scale(scale);
	  hMc.Add(plot);
	  iColor++;
	}

    }

  TCanvas *canvas = new TCanvas("c"+plotName, "c"+plotName, 500, 500);

  canvas->cd();

  hData->Draw();
  hMc.Draw("samehist");
  hData->Draw("same");

  canvas->Write();

}
