#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"

#include "../src/MuonPogTree.h"
#include "tdrstyle.C"

#include "rochcor2015.h"

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
    
    std::string muon_trackType; // applies to both, tag and probe     
  
    Float_t probe_minPt;      
    Float_t probe_isoCut;

    TString probe_ID;
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

}

// Helper classes defintion *****
// 1. parseConfig : parse the full cfg file
// ******************************

namespace muon_pog {
  void parseConfig(const std::string configFile, TagAndProbeConfig & tpConfig,
		   std::vector<SampleConfig> & sampleConfigs);

  double deltaR(double eta1, double phi1, double eta2, double phi2);
  bool hasGoodId(const muon_pog::Muon & muon, TString muId);
  bool hasFilterMatch(const TLorentzVector & muonTk, 
		      const muon_pog::HLT  & hlt,
		      std::string & filter, 
		      double dRCut);

  Int_t chargeFromTrk(const muon_pog::Muon & muon, std::string & trackType);
  TLorentzVector muonTk(const muon_pog::Muon & muon, std::string & trackType);

}



// The main program******** *****
// 1. Get configuration file and produces configurations
// 2. Create Plotters and loop on the event to fill them
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

  //const double etaBins[18] = {-2.4,-2.1,-1.6,-1.2,-1.05,-0.9,-0.6,-0.3,-0.2,
  // 			      0.2, 0.3, 0.6, 0.9, 1.05, 1.2, 1.6, 2.1, 2.4};

  const double etaBins[7] = {-2.4, -1.2, -0.9, 0., 0.9, 1.2, 2.4};

  // const double etaBins[8] = {-2.4, -2.1, -1.2, -0.9, 0.9, 1.2, 2.1, 2.4};

  double binW = TMath::Pi()/6.;

  const double phiBins[13] = {-6*binW, -5*binW, -4*binW, -3*binW, -2*binW, -binW, 0., binW, 2*binW, 3*binW, 4*binW, 5*binW, 6*binW};

  TProfile2D * pScaleMC   = new TProfile2D("pScaleMC",  "Average scale correction in montecarlo", 6, etaBins, 12, phiBins);
  TProfile2D * pScaleDATA = new TProfile2D("pScaleDATA","Average scale correction in real data",  6, etaBins, 12, phiBins);
  TH2F * hResolMC = new TH2F("hResolMC",  "Average resol correction in montecarlo", 6, etaBins, 500,-.1, .1);


  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();

  rochcor2015 *rochCor = new rochcor2015();
 
  TagAndProbeConfig tnpConfig;
  std::vector<SampleConfig> sampleConfigs;
  
  parseConfig(configFile,tnpConfig,sampleConfigs);

  for (auto sampleConfig : sampleConfigs)
    {

      bool runOnMC = sampleConfig.sampleName.Contains("MC");

      TString fileName = sampleConfig.fileName;
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
      int nEntries = 100000.; //tree->GetEntriesFast();
      std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

      int nFilteredEvents = 0;

      for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
	{
	  if (tree->LoadTree(iEvent)<0) break;

	  if (iEvent % 25000 == 0)
	    std::cout << "[" << argv[0] << "] processing event : " << iEvent << "\r" << std::flush;
	    
	  evBranch->GetEntry(iEvent);
	  float weight = ev->genInfos.size() > 0 ?
	    ev->genInfos[0].genWeight/fabs(ev->genInfos[0].genWeight) : 1.;	 
	  

	  bool pathHasFired = false;
	  
	  for (auto path : ev->hlt.triggers)
	    {
	      if (path.find(tnpConfig.hlt_path) != std::string::npos)
		{
		  pathHasFired = true;
		  break;
		}
	    }
  
	  if (!pathHasFired) continue;

	  std::vector<const muon_pog::Muon *> tagMuons;

	  for (auto & muon : ev->muons)
	    {
	      if (muon_pog::muonTk(muon,tnpConfig.muon_trackType).Pt() > tnpConfig.tag_minPt &&
		  muon_pog::hasFilterMatch(muon_pog::muonTk(muon,tnpConfig.muon_trackType),ev->hlt,tnpConfig.tag_hltFilter,tnpConfig.tag_hltDrCut) &&
		  muon_pog::hasGoodId(muon,tnpConfig.tag_ID) && 
		  muon.isoPflow04 < tnpConfig.tag_isoCut)
		tagMuons.push_back(&muon);
	    }
  
	  std::vector<const muon_pog::Muon *> probeMuons;
	  
	  for (auto & muon : ev->muons)
	    {
	      for (auto tagMuonPointer : tagMuons)
		{
		  const muon_pog::Muon & tagMuon = *tagMuonPointer;
		  
		  if ( tagMuonPointer != &muon && 
		       muon_pog::chargeFromTrk(tagMuon,tnpConfig.muon_trackType) * 
		       muon_pog::chargeFromTrk(muon,tnpConfig.muon_trackType) == -1)
		    {

		      // General Probe Muons	      
		      if(muon_pog::hasGoodId(muon,tnpConfig.probe_ID)   && 
			 muon.pt > tnpConfig.probe_minPt    &&
			 muon.isoPflow04 < tnpConfig.probe_isoCut)
			{
			  TLorentzVector tagMuTk(muon_pog::muonTk(tagMuon,tnpConfig.muon_trackType));
			  TLorentzVector muTk(muon_pog::muonTk(muon,tnpConfig.muon_trackType));
		  
			  //Float_t mass = (tagMuTk+muTk).M();
			  //Float_t dilepPt = (tagMuTk+muTk).Pt();
			  
			  Float_t origPtTag   = tagMuTk.Pt();
			  Float_t origPtProbe = muTk.Pt();
			      
			  Int_t chTag   = muon_pog::chargeFromTrk(tagMuon,tnpConfig.muon_trackType);
			  Int_t chProbe = muon_pog::chargeFromTrk(muon,tnpConfig.muon_trackType);

			  Float_t q = 1.;

			  if (runOnMC)
			    {

			      TLorentzVector tagMuTkResol = tagMuTk;
			      TLorentzVector muTkResol = muTk;
 
			      rochCor->momcor_mc(tagMuTk,chTag,0,q,false);
			      rochCor->momcor_mc(muTk,chProbe,0,q,false);

			      rochCor->momcor_mc(tagMuTkResol,chTag,0,q,true);
			      rochCor->momcor_mc(muTkResol,chProbe,0,q,true);

			      Float_t corrPtTag   = tagMuTk.Pt();
			      Float_t corrPtProbe = muTk.Pt();

			      Float_t resolPtTag   = tagMuTkResol.Pt();
			      Float_t resolPtProbe = muTkResol.Pt();

			      pScaleMC->Fill(tagMuTk.Eta(),tagMuTk.Phi(),chTag*(corrPtTag-origPtTag)/corrPtTag);
			      pScaleMC->Fill(muTk.Eta(),muTk.Phi(),chProbe*(corrPtProbe-origPtProbe)/corrPtProbe);
			      
			      // why need to apply cut on resol ?????
			      if(fabs(resolPtTag-corrPtTag)>0.0001)
				hResolMC->Fill(tagMuTk.Eta(),(resolPtTag-corrPtTag)/corrPtTag);
			      if(fabs(resolPtProbe-corrPtProbe)>0.0001)
				hResolMC->Fill(muTk.Eta(),(resolPtProbe-corrPtProbe)/corrPtProbe);
			    }
			  else
			    {
			      
			      rochCor->momcor_data(tagMuTk,chTag,0,q);
			      rochCor->momcor_data(muTk,chProbe,0,q);

			      Float_t corrPtTag   = tagMuTk.Pt();
			      Float_t corrPtProbe = muTk.Pt();

			      pScaleDATA->Fill(tagMuTk.Eta(),tagMuTk.Phi(),chTag*(corrPtTag-origPtTag)/corrPtTag);
			      pScaleDATA->Fill(muTk.Eta(),muTk.Phi(),chProbe*(corrPtProbe-origPtProbe)/corrPtProbe);
			      
			    }
			}
		    }
		}
	    }
	  
	}
      
      delete ev;

      outputFile->cd();

      std::cout << std::endl;

      if (runOnMC)
	{
	  for (Int_t iBin = 0; iBin < hResolMC->GetNbinsX(); ++iBin)
	    {
	      TH1D *proj = hResolMC->ProjectionY(TString("resolProjectionBin_")+TString(iBin),iBin,iBin+1);
	      proj->Fit("gaus","Q","",-3*proj->GetRMS(),3*proj->GetRMS());
	      std::cout << "Resolution result for [" << etaBins[iBin] << "," << etaBins[iBin+1] << "] = rms "
			<< proj->GetRMS() << "\t fit "<< proj->GetFunction("gaus")->GetParameter(2) << std::endl;
	      proj->Write();
	    }
	}

      inputFile->Close();
      std::cout << std::endl;
	   
    }

  TProfile2D pScaleDiff(*pScaleDATA);
  
  pScaleDiff.Add(pScaleMC,-1);

  for (Int_t iBinX = 1; iBinX <= hResolMC->GetNbinsX(); ++iBinX)
    {

      Float_t max = 0.;

      for (Int_t iBinY = 1; iBinY <= hResolMC->GetNbinsY(); ++iBinY)
	{
	  Float_t fDiff = fabs(pScaleDiff.GetBinContent(iBinX,iBinY));
	    if (fDiff>max)
	      max = fDiff;	  
	}
      
      std::cout << "Scale result for [" << etaBins[iBinX-1] << "," << etaBins[iBinX] << 
	"] = max scale deviation in phi " << " " << max << std::endl;
      
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

      probe_minPt  = vt.second.get<Float_t>("probe_minPt");
      probe_isoCut = vt.second.get<Float_t>("probe_isoCut"); // CB for now just comb reliso dBeta R04
      probe_ID     = vt.second.get<std::string>("probe_muonID");
      probe_fEtaBins = toPairArray(toArray(vt.second.get<std::string>("probe_fEtaMin")),
				   toArray(vt.second.get<std::string>("probe_fEtaMax")));
  
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

//Functions

double muon_pog::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1 - eta2;
  double dphi = phi1 - phi2;
  while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

  return sqrt(deta*deta + dphi*dphi);
}


bool muon_pog::hasGoodId(const muon_pog::Muon & muon, TString muId)
{

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


bool muon_pog::hasFilterMatch(const TLorentzVector & muTk,
			      const muon_pog::HLT  & hlt,
			      std::string & filter, 
			      double dRCut)
{

  for (auto object : hlt.objects)
    {
      if (object.filterTag.find(filter) != std::string::npos &&
	  muon_pog::deltaR(muTk.Eta(), muTk.Phi(), object.eta, object.phi) < dRCut)
	return true;
    }

  return false;
}


Int_t muon_pog::chargeFromTrk(const muon_pog::Muon & muon, std::string & trackType)
{

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

TLorentzVector muon_pog::muonTk(const muon_pog::Muon & muon, std::string & trackType)
{

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

