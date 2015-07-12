#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
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

namespace muon_pog {

  class PlotterConfig {

  public :

    // config parameters (public for direct access)

    TString general_title;

    Float_t plot_minInvMass;
    Float_t plot_maxInvMass;
    
    std::vector<TString> plot_fRapidityMin;
    std::vector<TString> plot_fRapidityMax;

    TString plot_signFitFunc;
    TString plot_bkgFitFunc;
    
    Float_t muon_minPt;   

    std::string muon_ID;  
    std::string muon_trackType;  

    Float_t muon_isoCut;
  
    std::vector<TString> muon_fEtaMin;
    std::vector<TString> muon_fEtaMax;
    
    std::string hlt_path; 
    
   
    PlotterConfig() {};
    
#ifndef __MAKECINT__ // CB CINT doesn't like boost :'-(    
    PlotterConfig(std::string & configFile); 
#endif

    ~PlotterConfig() {};

  private:
    std::vector<TString> toArray(const std::string& entries); 
  
  };

  class Plotter {

  public :
    
    Plotter(std::string config) : m_config(config) {};
    ~Plotter() {};
    
    void book(TFile *outFile);
    void fill(const std::vector<muon_pog::Muon> & muons, const muon_pog::HLT & hlt);
    void fit() {}; //CB empty before roofit 	  
    
  private :

    bool hasGoodId(const muon_pog::Muon & muon);
    Int_t chargeFromTrk(const muon_pog::Muon & muon);
    TLorentzVector muonTk(const muon_pog::Muon & muon);
    
    PlotterConfig m_config;
    std::map<TString,TH1 *> m_plots;
    
  };

}



int main(int argc, char* argv[]){
  using namespace muon_pog;


  if (argc < 3) 
    {
      std::cout << "Usage : "
		<< argv[0] << " PATH_TO_INPUT_FILE PAT_TO_CONFIG_FILE(s)\n";
      exit(100);
    }

  // Input root file
  TString fileName = argv[1];

  std::cout << "[" << argv[0] << "] Processing file " << fileName.Data() << std::endl;
  
  std::vector<Plotter> plotters;
  for (int iConfig = 2; iConfig < argc; ++iConfig)
    {
        std::cout << "[" << argv[0] << "] Using config file " << argv[iConfig] << std::endl;
	plotters.push_back(std::string(argv[iConfig]));
    }
  
  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  //setTDRStyle(); what to do here?
   
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

  system("mkdir -p results");
  
  TFile* outputFile = TFile::Open("results/results.root","RECREATE"); // CB find a better name for output file  

  for (auto & plotter : plotters)
    plotter.book(outputFile);
      
  // Watch number of entries
  int nEntries = tree->GetEntriesFast();
  std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

  int nFilteredEvents = 0;
  
  for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
    {
      if (tree->LoadTree(iEvent)<0) break;

      evBranch->GetEntry(iEvent);

      for (auto & plotter : plotters)
	plotter.fill(ev->muons, ev->hlt);

    }

  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;
}

// CB Helpers: the configuration class!
muon_pog::PlotterConfig::PlotterConfig (std::string & configFile)
{

  boost::property_tree::ptree pt;
  
  
  try
    {
      boost::property_tree::ini_parser::read_ini(configFile, pt);
    }
  catch (boost::property_tree::ini_parser::ini_parser_error iniParseErr)
    {
      std::cout << "[PlotterConfig] Can't open : " << iniParseErr.filename()
		<< "\n\tin line : " << iniParseErr.line()
		<< "\n\thas error :" << iniParseErr.message()
		<< std::endl;
      throw std::runtime_error("Bad INI parsing");
    }

  try
    {

      general_title = TString(pt.get<std::string>("general.title"));
 
      plot_minInvMass = pt.get<Float_t>("plot.minInvMass");
      plot_maxInvMass = pt.get<Float_t>("plot.maxInvMass");

      plot_fRapidityMin = toArray(pt.get<std::string>("plot.fRapidityMin"));
      plot_fRapidityMax = toArray(pt.get<std::string>("plot.fRapidityMax"));

      plot_signFitFunc = TString(pt.get<std::string>("plot.signFitFunc").c_str()); // CB is the casting needed?
      plot_bkgFitFunc  = TString(pt.get<std::string>("plot.bkgFitFunc").c_str());
      
      
      muon_minPt     = pt.get<Float_t>("muon.minPt");
      muon_ID        = pt.get<std::string>("muon.muonID");
      muon_trackType = pt.get<std::string>("muon.trackType");
      muon_isoCut    = pt.get<Float_t>("muon.isoCut"); // CB for now just dBeta R04

      muon_fEtaMin = toArray(pt.get<std::string>("muon.fEtaMin"));
      muon_fEtaMax = toArray(pt.get<std::string>("muon.fEtaMax"));
      
     hlt_path = pt.get<std::string>("hlt.path");

    }

  catch (boost::property_tree::ptree_bad_data bd)
    {
      std::cout << "[PlotterConfig] Cant'g get data : has error : "
		<< bd.what() << std::endl;
      throw std::runtime_error("Bad INI variables");
    }

}

std::vector<TString> muon_pog::PlotterConfig::toArray(const std::string& entries)
{
  std::vector<TString> result;
  std::stringstream sentries(entries);
  std::string item;
  while(std::getline(sentries, item, ','))
    result.push_back(TString(item));
  return result;
}

std::ostream&
operator<<(std::ostream& out, const muon_pog::PlotterConfig &config)
{

  out << "PlotterConfig Configuration :" // CB to be filled 
      << std::endl;
   
  return out;

}

void muon_pog::Plotter::book(TFile *outFile)
{

  TString titleTag = m_config.general_title;
  
  outFile->cd("/");
  outFile->mkdir(titleTag);
  outFile->cd(titleTag);

  std::vector<TString>::const_iterator rMinIt  = m_config.plot_fRapidityMin.begin();
  std::vector<TString>::const_iterator rMinEnd = m_config.plot_fRapidityMin.end();

  std::vector<TString>::const_iterator rMaxIt  = m_config.plot_fRapidityMax.begin();
  std::vector<TString>::const_iterator rMaxEnd = m_config.plot_fRapidityMax.end();
  
  for (; rMinIt != rMinEnd || rMaxIt != rMaxEnd; ++rMinIt, ++rMaxIt)
    {
         
      TString hName = TString("hInvMass") + "_rMin" + TString((*rMinIt))
	              + "_rMax" + TString((*rMaxIt));  
      m_plots[hName] = new TH1F(hName + "_" + titleTag,hName,100, m_config.plot_minInvMass,
				m_config.plot_maxInvMass);
    }

  std::vector<TString>::const_iterator fEtaMinIt  = m_config.muon_fEtaMin.begin();
  std::vector<TString>::const_iterator fEtaMinEnd = m_config.muon_fEtaMin.end();

  std::vector<TString>::const_iterator fEtaMaxIt  = m_config.muon_fEtaMax.begin();
  std::vector<TString>::const_iterator fEtaMaxEnd = m_config.muon_fEtaMax.end();
  
  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
    {
         
      TString hName = TString("hInvMass") + "_fEtaMin" + (*fEtaMinIt)
	              + "_fEtaMax" + (*fEtaMaxIt);  
      m_plots[hName] = new TH1F(hName + "_" + titleTag,hName,100, m_config.plot_minInvMass,
				m_config.plot_maxInvMass);
    }

  m_plots["mu1Pt"] = new TH1F("hMu1Pt_" + titleTag,"mu1Pt",200,0.,200.);
  m_plots["mu2Pt"] = new TH1F("hMu2Pt_" + titleTag,"mu2Pt",200,0.,200.);

  m_plots["mu1EtaPhi"] = new TH2F("hMu1EtaPhi_" + titleTag,"mu1EtaPhi",100,-2.5,2.5,100,-TMath::Pi(),TMath::Pi());
  m_plots["mu2EtaPhi"] = new TH2F("hMu2EtaPhi_" + titleTag,"mu2EtaPhi",100,-2.5,2.5,100,-TMath::Pi(),TMath::Pi());

}

void muon_pog::Plotter::fill(const std::vector<muon_pog::Muon> & muons,
			    const muon_pog::HLT & hlt)
{

  bool pathHasFired = false;

  for (auto path : hlt.triggers)
    {
      if (path.find(m_config.hlt_path) != std::string::npos)
	{
	  pathHasFired = true;
	  break;
	}
    }

  if (!pathHasFired) return;

  std::vector<muon_pog::Muon> goodMuons;

  for (auto & muon : muons)
    {
      if (hasGoodId(muon) &&
	  muonTk(muon).Pt() > m_config.muon_minPt &&
	  muon.isoPflow04 < m_config.muon_isoCut)
	goodMuons.push_back(muon);
    }

  std::vector<muon_pog::Muon>::const_iterator goodMu1It  = goodMuons.begin();
  std::vector<muon_pog::Muon>::const_iterator goodMuEnd  = goodMuons.end();

  for (; goodMu1It != goodMuEnd; ++goodMu1It)
    {

      std::vector<muon_pog::Muon>::const_iterator goodMu2It  = goodMu1It;
      
      for (goodMu2It++; goodMu2It != goodMuEnd; ++goodMu2It)
	{
	  
	  if (chargeFromTrk((*goodMu1It)) * chargeFromTrk((*goodMu2It)) != -1)
	    continue;

	  TLorentzVector mu1Tk = muonTk((*goodMu1It));
	  TLorentzVector mu2Tk = muonTk((*goodMu2It));

	  Float_t mass = (mu1Tk+mu2Tk).M();

	  m_plots["mu1Pt"]->Fill(mu1Tk.Pt());
	  m_plots["mu2Pt"]->Fill(mu2Tk.Pt());
	  
	  m_plots["mu1EtaPhi"]->Fill(mu1Tk.Eta(),mu1Tk.Phi());
	  m_plots["mu2EtaPhi"]->Fill(mu2Tk.Eta(),mu2Tk.Phi());
	  
	  std::vector<TString>::const_iterator rMinIt  = m_config.plot_fRapidityMin.begin();
	  std::vector<TString>::const_iterator rMinEnd = m_config.plot_fRapidityMin.end();

	  std::vector<TString>::const_iterator rMaxIt  = m_config.plot_fRapidityMax.begin();
	  std::vector<TString>::const_iterator rMaxEnd = m_config.plot_fRapidityMax.end();
  
	  for (; rMinIt != rMinEnd || rMaxIt != rMaxEnd; ++rMinIt, ++rMaxIt)
	    {

	      Float_t rapidity = (mu1Tk+mu2Tk).Rapidity();

	      if (fabs(rapidity) > rMinIt->Atof() &&
		  fabs(rapidity) < rMaxIt->Atof() )
		{
		  
		  TString hName = TString("hInvMass")
		    + "_rMin" + TString((*rMinIt))
		    + "_rMax" + TString((*rMaxIt));  
		  m_plots[hName]->Fill(mass);

		}
	      
	    }

	  std::vector<TString>::const_iterator fEtaMinIt  = m_config.muon_fEtaMin.begin();
	  std::vector<TString>::const_iterator fEtaMinEnd = m_config.muon_fEtaMin.end();

	  std::vector<TString>::const_iterator fEtaMaxIt  = m_config.muon_fEtaMax.begin();
	  std::vector<TString>::const_iterator fEtaMaxEnd = m_config.muon_fEtaMax.end();
  
	  for (; fEtaMinIt != fEtaMinEnd || fEtaMaxIt != fEtaMaxEnd; ++fEtaMinIt, ++fEtaMaxIt)
	    {
         
	      if (fabs(mu1Tk.Eta()) > fEtaMinIt->Atof() &&
		  fabs(mu1Tk.Eta()) < fEtaMaxIt->Atof() &&
		  fabs(mu2Tk.Eta()) > fEtaMinIt->Atof() &&
		  fabs(mu2Tk.Eta()) < fEtaMaxIt->Atof())
		{
		  
		  TString hName = TString("hInvMass")
		    + "_fEtaMin" + TString((*fEtaMinIt))
		    + "_fEtaMax" + TString((*fEtaMaxIt));  
		  m_plots[hName]->Fill(mass);

		}
	    }
	}
    }
      
}

bool muon_pog::Plotter::hasGoodId(const muon_pog::Muon & muon)
{
  std::string & muId = m_config.muon_ID;

  if (muId == "GLOBAL")      return muon.isGlobal == 1 ;
  else if (muId == "TIGHT")  return muon.isTight == 1;
  else if (muId == "MEDIUM") return muon.isMedium == 1;
  else if (muId == "LOOSE")  return muon.isLoose == 1;
  else if (muId == "HIGHPT") return muon.isHighPt == 1;

  else
    {
      std::cout << "[Plotter::hasGoodId]: Invalid muon id : "
		<< muId << std::endl;
      exit(900);
    }
  
}

Int_t muon_pog::Plotter::chargeFromTrk(const muon_pog::Muon & muon)
{
  std::string & trackType = m_config.muon_trackType;

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
  std::string & trackType = m_config.muon_trackType;

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


  
  
