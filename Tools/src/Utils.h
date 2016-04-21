#ifndef Utils_h__
#define Utils_h__

#include "MuonPogTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <string>

namespace muon_pog
{

  // The names says it all
  double deltaR(double eta1, double phi1, double eta2, double phi2)
  {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

    return sqrt(deta*deta + dphi*dphi);
  }

    
  // Check if a muon_pog::Muon passes a given ID (embedded using official selectors)
  // Valid IDs are GLOBAL, TRACKER, SOFT, LOOSE, MEDIUM, TIGHT, HIGHPT
  bool hasGoodId(const muon_pog::Muon & muon, TString muId)
  {
    
    if (muId == "GLOBAL")      return muon.isGlobal  == 1;
    else if (muId == "TRACKER")return muon.isTracker == 1;
    else if (muId == "TIGHT")  return muon.isTight   == 1;
    else if (muId == "MEDIUM") return muon.isMedium  == 1;
    else if (muId == "LOOSE")  return muon.isLoose   == 1;
    else if (muId == "HIGHPT") return muon.isHighPt  == 1;
    else if (muId == "SOFT")   return muon.isSoft == 1;
    else
      {
	std::cout << "[Plotter::hasGoodId]: Invalid muon id : "
		  << muId << std::endl;
	exit(900);
      }
    
    return 0;
 
  };


  // Returns the charge muon_pog::Muon for a given fit 
  // Valid track fits are: PF, TUNEP, GLB, INNER 
  Int_t chargeFromTrk(const muon_pog::Muon & muon, 
		      const std::string & trackType)
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

    return 999;
    
  }


  // Return a TLorentz vector out of a given fit from muon_pog::Muon 
  // Valid track fits are: PF, TUNEP, GLB, INNER 
  TLorentzVector muonTk(const muon_pog::Muon & muon, 
			const std::string & trackType)
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


  // Checks if a trigger path fired using muon_pog::HLT
  // if the path name is "none" returns always true
  bool pathHasFired(const muon_pog::HLT  & hlt, std::string pathName)
  {
    
    if (pathName == "none")
      return true;
    
    for (auto path : hlt.triggers)
      {
	if (path.find(pathName) != std::string::npos)
	  {
	    return true;
	  }
      }
  }
  

  // Checks if the iner track of muon_pog::Muon matches geometrically 
  // in dR with a given muon_pog::HLT object filter
  // if the filter name is "none" returns always true
  bool hasFilterMatch(const muon_pog::Muon & muon,
		      const muon_pog::HLT  & hlt,
		      std::string & filter, Float_t dR)
  {

    if (filter == "none")
      return true;
    
    TLorentzVector muTk = muonTk(muon,std::string("INNER"));
    
    for (auto object : hlt.objects)
      {
	if (object.filterTag.find(filter) != std::string::npos &&
	    deltaR(muTk.Eta(), muTk.Phi(), object.eta, object.phi) < dR)
	  return true;
      }
    
    return false;
    
  }

  bool hasMother(const muon_pog::GenParticle & gen, Int_t pdgId)
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


  const muon_pog::GenParticle * hasGenMatch(const muon_pog::Muon & muon,
					    const std::vector<muon_pog::GenParticle>  & gens,
					    Float_t dRCut, Int_t motherPdgId = 0, Int_t vetoPdgId = 0)
  {
   
    TLorentzVector muTk = muonTk(muon,std::string("INNER"));
    
    const muon_pog::GenParticle * bestGen = 0;
    Float_t bestDr = 999.;
    
    for (auto & gen : gens)
      {
	if (fabs(gen.pdgId) == 13)
	  {
	    Float_t dr = deltaR(muTk.Eta(), muTk.Phi(), gen.eta, gen.phi);
	    if (dr < dRCut && dr < bestDr && 
		((vetoPdgId == 0) || !hasMother(gen,vetoPdgId)))
	      {
		bestGen = &gen;
		bestDr = dr;
	      }
	  }
      }

    return (bestGen && ((motherPdgId == 0) || hasMother(*bestGen,motherPdgId))) ? bestGen : 0;
  }


  //From Fede, the function name says it all
  void addUnderFlow(TH1 &hist)
  {
    hist.SetBinContent(1, hist.GetBinContent(0) + 
		          hist.GetBinContent(1));
    hist.SetBinError  (1, sqrt(hist.GetBinError(0)*hist.GetBinError(0) + 
			       hist.GetBinError(1)*hist.GetBinError(1)));
    hist.SetBinContent(0, 0); 
    hist.SetBinError  (0, 0);  
  }



  //From Fede, the function name says it all
  void addOverFlow(TH1 &hist)
  {
    Int_t lastBin = hist.GetNbinsX(); 
    hist.SetBinContent(lastBin, hist.GetBinContent(lastBin) + 
		                hist.GetBinContent(lastBin+1));
    hist.SetBinError  (lastBin, sqrt(hist.GetBinError(lastBin)*hist.GetBinError(lastBin) + 
				     hist.GetBinError(lastBin+1)*hist.GetBinError(lastBin+1))); 
    hist.SetBinContent(lastBin+1, 0) ; 
    hist.SetBinError  (lastBin+1, 0) ; 
    
  }
  
}

#endif


