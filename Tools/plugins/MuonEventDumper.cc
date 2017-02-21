///////////////////////////////////////
// Dump muon related information for a specific event
//////////////////////////////////////

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include <DataFormats/PatCandidates/interface/MET.h>


#include <TH1.h>
#include <TFile.h>

#include <algorithm>
#include <iostream>

template<typename T> edm::Handle<T> conditionalGet(const edm::Event & ev,
						   const edm::EDGetTokenT<T> & token,
						   const std::string & collectionName)
{

  edm::Handle<T> collection ;

  if (!token.isUninitialized()) 
    {
      if (!ev.getByToken(token, collection)) 
	edm::LogError("") << "[MuonEventDumper]::conditionalGet: "
			  << collectionName << " collection does not exist !!!";
    }

  return collection;

}
  
class MuonEventDumper : public edm::EDAnalyzer 
{
public:

  MuonEventDumper(const edm::ParameterSet &);
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
  
  static int counterAll; 
  static  int counterDuplicated;
  static int counterFilter;
  static int counter1; 
  static int counter2; 

private:
  
  void printGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &) const ;

  void printGenParticles(const edm::Handle<reco::GenParticleCollection> &) const ;

  void printHlt(const edm::Handle<edm::TriggerResults> &, 
		const edm::Handle<trigger::TriggerEvent> &,
		const edm::TriggerNames &) const ;
  
  void printPV(const edm::Handle<std::vector<reco::Vertex> > &) const ;
  
  void printMuons(const edm::Handle<edm::View<reco::Muon> > &,
		  const edm::Handle<std::vector<reco::Vertex> > &,
		  const edm::Handle<reco::BeamSpot> &) const ;

  void printMediumIDCuts(const reco::Muon &) const ;

  void printTrack(const reco::Track * track,
		  const std::string & trackType) const;

		  // CB copy of PF muon selections as they are in 740
  bool isIsolatedMuon(const reco::Muon& muon) const ;
  bool isGlobalTightMuon( const reco::Muon& muon ) const ; 
  bool isTrackerTightMuon( const reco::Muon& muon ) const ;

  void printIsolation( const reco::MuonIsolation & iso, const std::string & cone) const;
  void printPFIsolation( const reco::MuonPFIsolation & iso, const std::string & cone) const;
  
  void printMET(const edm::Handle<std::vector<pat::MET> > & pfmet_) const;

  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigSummaryToken_;

  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > primaryVertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileUpInfoToken_;

  edm::EDGetTokenT<edm::View<reco::Track> > adHocTrackToken_;
  std::string adHocTrackLabel_;  

  edm::EDGetTokenT<std::vector<pat::MET> > metToken_;


  bool isMINIAOD_;


  TFile * file;

  TH1F * h_met;
  TH1F * h_bad_met;
  TH1F * h_good_met;

  TH1F * h_bad_segmentComp;  
  TH1F * h_bad_glbChi2;  
  TH1F * h_bad_staChi2;  
  TH1F * h_bad_TkStaChi2;
  TH1F * h_bad_type;
  TH1F * h_bad_nMatches; 
  TH1F * h_bad_MuIdAlgo;

  TH1F * h_good_segmentComp;  
  TH1F * h_good_glbChi2; 
  TH1F * h_good_staChi2;
  TH1F * h_good_TkStaChi2;
  TH1F * h_good_type;
  TH1F * h_good_nMatches;
  TH1F * h_good_MuIdAlgo;
 
};


int MuonEventDumper::counterAll = 0; 
int MuonEventDumper::counterDuplicated = 0;
int MuonEventDumper::counter1 = 0;
int MuonEventDumper::counter2 = 0;
int MuonEventDumper::counterFilter = 0;

MuonEventDumper::MuonEventDumper( const edm::ParameterSet & cfg )
{

  // Input collections
  edm::InputTag tag = cfg.getUntrackedParameter<edm::InputTag>("TrigResultsTag", edm::InputTag("TriggerResults::HLT"));
  if (tag.label() != "none") trigResultsToken_ = consumes<edm::TriggerResults>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("TrigSummaryTag", edm::InputTag("hltTriggerSummaryAOD::HLT")); 
  if (tag.label() != "none") trigSummaryToken_ =consumes<trigger::TriggerEvent>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"));
  if (tag.label() != "none") muonToken_ = consumes<edm::View<reco::Muon> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag", edm::InputTag("offlinePrimaryVertices"));
  if (tag.label() != "none") primaryVertexToken_ = consumes<std::vector<reco::Vertex> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag", edm::InputTag("offlineBeamSpot"));
  if (tag.label() != "none") beamSpotToken_ = consumes<reco::BeamSpot>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("GenTag", edm::InputTag("prunedGenParticles"));
  if (tag.label() != "none") genToken_ = consumes<reco::GenParticleCollection>(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("PileUpInfoTag", edm::InputTag("pileupInfo"));
  if (tag.label() != "none") pileUpInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("AdHocTrackTag", edm::InputTag("AdHocTracks"));
  if (tag.label() != "none") adHocTrackToken_ = consumes<edm::View<reco::Track> >(tag);

  tag = cfg.getUntrackedParameter<edm::InputTag>("METTag", edm::InputTag("mets"));
  if (tag.label() != "none") metToken_ = consumes<std::vector<pat::MET> >(tag);

  isMINIAOD_= cfg.getUntrackedParameter<bool>("isMINIAOD");
	     //if (tag.label() != "none") isMINIAOD_ = tag;

}


void MuonEventDumper::beginJob() 
{  

  file = new TFile("outfile.root","recreate");
  //histo = new TH1F("pt","pt",1000,0,100);  

  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  h_met = new TH1F("met", "met", 100, 0, 400); 
  h_good_met = new TH1F("good_met", "met", 100, 0, 400); 
  h_bad_met = new TH1F("bad_met", "met", 100, 0, 400); 

  h_bad_segmentComp = new TH1F("bad_segComp", "segComp", 50, 0, 1);  
  h_bad_glbChi2 = new TH1F("bad_glbChi2", "glbChi2", 50, 0, 20);  
  h_bad_staChi2 = new TH1F("bad_staChi2", "staChi2", 50, 0, 20); 
  h_bad_TkStaChi2 = new TH1F("bad_TkStaChi2", "TkStaChi2", 50, 0, 20); 
  h_bad_type = new TH1F("bad_type", "type", 4, 1, 5);
  h_bad_nMatches = new TH1F("bad_nMatches", "nMatches", 5, 0, 5);
  h_bad_MuIdAlgo = new TH1F("bad_MuIdAlgo", "MuIdAlgo", 14, 0, 14);
  
  h_good_segmentComp = new TH1F("good_segComp", "segComp", 50, 0, 1);  
  h_good_glbChi2 = new TH1F("good_glbChi2", "glbChi2", 50, 0, 20);
  h_good_staChi2 = new TH1F("good_staChi2", "staChi2", 50, 0, 20);   
  h_good_TkStaChi2 = new TH1F("good_TkStaChi2", "TkStaChi2", 50, 0, 20); 
  h_good_type = new TH1F("good_type", "type", 4, 1, 5);
  h_good_nMatches = new TH1F("good_nMatches", "nMatches", 5, 0, 5);
  h_good_MuIdAlgo = new TH1F("good_MuIdAlgo", "MuIdAlgo", 14, 0, 14);


  TH1::AddDirectory(oldAddDir);
}


void MuonEventDumper::beginRun(const edm::Run & run, const edm::EventSetup & config )
{
  
} 


void MuonEventDumper::endJob() 
{

  std::cout << "PRINT COUNTER: " << counterAll << std::endl;
  std::cout << "PRINT COUNTER DUPLICATED: " << counterDuplicated << std::endl;
  std::cout << "PRINT STA&GLB&TK&PF: " << counter1 << std::endl;
  std::cout << "PRINT TK&PF: " << counter2 << std::endl;
  std::cout << "PRINT COUNTER FILTER: " << counterFilter << std::endl;


  file->Write();
  file->Close();

}

void MuonEventDumper::analyze (const edm::Event & ev, const edm::EventSetup &)
{

  counterAll++;

  std::cout << "[MuonEventDumper::analyze]: "
	    << "\nProcessing Run: " << ev.id().run()
	    << "\tLumiBlock: " << ev.id().luminosityBlock()
	    << "\tevent: " << ev.id().event() << std::endl;


  //  if ( ev.id().run() == 276831  &&  ev.id().luminosityBlock() == 141 &&  ev.id().event() == 161758936                                                                              ) {
  // Print GEN pile up information
  if (!ev.isRealData()) 
    {
      auto puInfo = conditionalGet<std::vector<PileupSummaryInfo> >(ev,pileUpInfoToken_,"Pile-Up Info");
      //if (puInfo.isValid()) printGenInfo(puInfo);
    }
  

  // Print GEN particles information
  if (!ev.isRealData()) 
    {
      auto genParticles= conditionalGet<reco::GenParticleCollection>(ev,genToken_,"GenParticle Collection");
      //if (genParticles.isValid()) printGenParticles(genParticles);
    } 
  
  // Print trigger information
  auto triggerResults = conditionalGet<edm::TriggerResults>(ev,trigResultsToken_,"TriggerResults");
  auto triggerEvent   = conditionalGet<trigger::TriggerEvent>(ev,trigSummaryToken_,"TriggerEvent");
      
  //if (triggerResults.isValid() && triggerEvent.isValid()) 
    //printHlt(triggerResults, triggerEvent,ev.triggerNames(*triggerResults));
  
  // Print vertex information
  auto vertexes = conditionalGet<std::vector<reco::Vertex> >(ev,primaryVertexToken_,"Vertex");
  //if (vertexes.isValid()) printPV(vertexes);
  
  // Get beam spot for muons
  auto beamSpot = conditionalGet<reco::BeamSpot>(ev,beamSpotToken_,"BeamSpot");

  auto mets = conditionalGet<std::vector<pat::MET> >(ev,metToken_,"MET");
  // if (mets.isValid()) printMET(mets);

  // Get muons  
  auto muons = conditionalGet<edm::View<reco::Muon> >(ev,muonToken_,"MuonCollection");

  // Print muon information
  //  if (muons.isValid() && vertexes.isValid() && beamSpot.isValid()) 
    //printMuons(muons,vertexes,beamSpot);

  if ( !isMINIAOD_) {

    auto tracks = conditionalGet<edm::View<reco::Track> >(ev,adHocTrackToken_, "AdHocTracks");

    if (tracks.isValid() )
      {
      std::cout << "Ad-hoc trackcollection is " + adHocTrackLabel_ + " and has size: " << tracks->size() << std::endl;

      edm::View<reco::Track> ::const_iterator trackIt  = tracks->begin();
      edm::View<reco::Track> ::const_iterator trackEnd = tracks->end();

      for (; trackIt != trackEnd; ++trackIt) 
	{
	  printTrack(&(*trackIt),adHocTrackLabel_);
	}      
      }
  }
  //  


  //// DO MUON ANALISIS

  
  edm::View<reco::Muon> ::const_iterator muonIt1  = muons->begin();
  edm::View<reco::Muon> ::const_iterator muonIt2  = muonIt1+1;

  edm::View<reco::Muon> ::const_iterator muonEnd = muons->end();

  float minDR = 999;
   edm::View<reco::Muon> ::const_iterator muA = muonIt1;
   edm::View<reco::Muon> ::const_iterator muB = muonIt1;

  if (muons->size() > 1) {

    for (; muonIt1 != muonEnd; ++muonIt1) {
      for (; muonIt2 != muonEnd; ++muonIt2) {
 
	if (muonIt1 == muonIt2) continue;

	const reco::Muon& mu1 = (*(muonIt1));
	const reco::Muon& mu2 = (*(muonIt2));

	float deta = mu1.eta() - mu2.eta();
	float dphi = mu1.phi() - mu2.phi();
	
	float dR = sqrt(deta*deta + dphi*dphi);
      
	if (dR < minDR) {
	  minDR = dR; 
	  muA = muonIt1;
	  muB = muonIt2;
	} 
      
      }
    }  
  }


 
 bool goodEvent = false;

  if ( true) { //minDR < 0.2) {

  

  edm::View<reco::Muon> ::const_iterator muonIt  = muons->begin();

  



  // check overlap between muons   
 
  // shared segmets 
  int nSharedSeg = muon::sharedSegments(*muA, *muB);

  for (; muonIt != muonEnd; ++muonIt) 
    {

    
      //      if ( muonIt != muA && muonIt != muB ) continue;

      const reco::Muon& mu = (*muonIt);
      const reco::Vertex & vertex = vertexes->at(0); // CB for now vertex is always valid, but add a protection	    
      
      //histo->Fill(mu.innerTrack()->pt()); 

      //if ( !mu.isPFMuon() ) continue; 
      //if ( mu.muonBestTrack()->pt() < 100 ) continue;
      
      //if (muon::segmentCompatibility(mu) < 0.3) continue;
      //if ( mu.numberOfMatches() == 0 ) continue;
      //if (!muon::isGoodMuon(mu, muon::TM2DCompatibilityLoose)) continue;

      if ( !muon::isLooseMuon(mu)) continue;
      

      goodEvent = true;

      if ( mu.muonBestTrack()->pt() >= 95 ){

	
	if (mu.isStandAloneMuon() &&  mu.isGlobalMuon() && mu.isTrackerMuon()) 
	  h_bad_type->Fill(1.0);
	if (!mu.isGlobalMuon() && mu.isTrackerMuon()) 
	  h_bad_type->Fill(2.0);
	if (!mu.isStandAloneMuon() && mu.isGlobalMuon() && mu.isTrackerMuon()) 
	  h_bad_type->Fill(3.0);
	//h_bad_type->Fill(4.0);
	

	h_bad_segmentComp->Fill(muon::segmentCompatibility(mu));  
	if( mu.isGlobalMuon()) {
	  h_bad_glbChi2->Fill( mu.globalTrack()->normalizedChi2());
	  h_bad_TkStaChi2->Fill(mu.combinedQuality().chi2LocalPosition);
	}
	if(mu.isStandAloneMuon())h_bad_staChi2->Fill( mu.outerTrack()->normalizedChi2());

	h_bad_nMatches->Fill(mu.numberOfMatches());

	// arbitration
	//muon::isGoodMuon(recoMu, TMOneStationTight) 
	//if(muon::isGoodMuon(*muon, muon::TrackerMuonArbitrated)) // resolve ambiguity of sharing segments

	// Muon ID selection. As described in AN-2008/098
	if(muon::isGoodMuon(mu, muon::All))                                // dummy options - always true
      h_bad_MuIdAlgo->Fill(0);
    if(muon::isGoodMuon(mu, muon::AllStandAloneMuons))                 // checks isStandAloneMuon flag
      h_bad_MuIdAlgo->Fill(1);
    if(muon::isGoodMuon(mu, muon::AllTrackerMuons))                    // checks isTrackerMuon flag
      h_bad_MuIdAlgo->Fill(2);
    if(muon::isGoodMuon(mu, muon::TrackerMuonArbitrated))              // resolve ambiguity of sharing segments
      h_bad_MuIdAlgo->Fill(3);
    if(muon::isGoodMuon(mu, muon::AllArbitrated))                      // all muons with the tracker muon arbitrated
      h_bad_MuIdAlgo->Fill(4);
    if(muon::isGoodMuon(mu, muon::GlobalMuonPromptTight))              // global muons with tighter fit requirements
      h_bad_MuIdAlgo->Fill(5);
    if(muon::isGoodMuon(mu, muon::TMLastStationLoose))                 // penetration depth loose selector
      h_bad_MuIdAlgo->Fill(6);
    if(muon::isGoodMuon(mu, muon::TMLastStationTight))                 // penetration depth tight selector
      h_bad_MuIdAlgo->Fill(7);
    if(muon::isGoodMuon(mu, muon::TM2DCompatibilityLoose))             // likelihood based loose selector
      h_bad_MuIdAlgo->Fill(8);
    if(muon::isGoodMuon(mu, muon::TM2DCompatibilityTight))             // likelihood based tight selector
      h_bad_MuIdAlgo->Fill(9);
    if(muon::isGoodMuon(mu, muon::TMOneStationLoose))                  // require one well matched segment
      h_bad_MuIdAlgo->Fill(10);
    if(muon::isGoodMuon(mu, muon::TMOneStationTight))                  // require one well matched segment
      h_bad_MuIdAlgo->Fill(11);
    if(muon::isGoodMuon(mu, muon::TMLastStationOptimizedLowPtLoose))   // combination of TMLastStation and TMOneStation
      h_bad_MuIdAlgo->Fill(12);
    if(muon::isGoodMuon(mu, muon::TMLastStationOptimizedLowPtTight))   // combination of TMLastStation and TMOneStation
      h_bad_MuIdAlgo->Fill(13);
	//calo compatibility 
	// hMuCaloCompatibility->Fill(muon->caloCompatibility());

      
      } else{

	
	if (mu.isStandAloneMuon() &&  mu.isGlobalMuon() && mu.isTrackerMuon()) 
	  h_good_type->Fill(1.0);
	if (!mu.isGlobalMuon() && mu.isTrackerMuon()) 
	  h_good_type->Fill(2.0);
	if (!mu.isStandAloneMuon() && mu.isGlobalMuon() && mu.isTrackerMuon()) 
	  h_good_type->Fill(3.0);
	//h_bad_type->Fill(4.0);
	

	h_good_segmentComp->Fill(muon::segmentCompatibility(mu));  
	if(mu.isGlobalMuon()) { 
	  h_good_glbChi2->Fill( mu.globalTrack()->normalizedChi2());
	  h_good_TkStaChi2->Fill(mu.combinedQuality().chi2LocalPosition);
	}
	if(mu.isStandAloneMuon())h_good_staChi2->Fill( mu.outerTrack()->normalizedChi2());
	h_good_nMatches->Fill(mu.numberOfMatches());
      
      	// Muon ID selection. As described in AN-2008/098
	if(muon::isGoodMuon(mu, muon::All))                                // dummy options - always true
      h_good_MuIdAlgo->Fill(0);
    if(muon::isGoodMuon(mu, muon::AllStandAloneMuons))                 // checks isStandAloneMuon flag
      h_good_MuIdAlgo->Fill(1);
    if(muon::isGoodMuon(mu, muon::AllTrackerMuons))                    // checks isTrackerMuon flag
      h_good_MuIdAlgo->Fill(2);
    if(muon::isGoodMuon(mu, muon::TrackerMuonArbitrated))              // resolve ambiguity of sharing segments
      h_good_MuIdAlgo->Fill(3);
    if(muon::isGoodMuon(mu, muon::AllArbitrated))                      // all muons with the tracker muon arbitrated
      h_good_MuIdAlgo->Fill(4);
    if(muon::isGoodMuon(mu, muon::GlobalMuonPromptTight))              // global muons with tighter fit requirements
      h_good_MuIdAlgo->Fill(5);
    if(muon::isGoodMuon(mu, muon::TMLastStationLoose))                 // penetration depth loose selector
      h_good_MuIdAlgo->Fill(6);
    if(muon::isGoodMuon(mu, muon::TMLastStationTight))                 // penetration depth tight selector
      h_good_MuIdAlgo->Fill(7);
    if(muon::isGoodMuon(mu, muon::TM2DCompatibilityLoose))             // likelihood based loose selector
      h_good_MuIdAlgo->Fill(8);
    if(muon::isGoodMuon(mu, muon::TM2DCompatibilityTight))             // likelihood based tight selector
      h_good_MuIdAlgo->Fill(9);
    if(muon::isGoodMuon(mu, muon::TMOneStationLoose))                  // require one well matched segment
      h_good_MuIdAlgo->Fill(10);
    if(muon::isGoodMuon(mu, muon::TMOneStationTight))                  // require one well matched segment
      h_good_MuIdAlgo->Fill(11);
    if(muon::isGoodMuon(mu, muon::TMLastStationOptimizedLowPtLoose))   // combination of TMLastStation and TMOneStation
      h_good_MuIdAlgo->Fill(12);
    if(muon::isGoodMuon(mu, muon::TMLastStationOptimizedLowPtTight))   // combination of TMLastStation and TMOneStation
      h_good_MuIdAlgo->Fill(13);
      }

      //if( muon::segmentCompatibility(mu) < 0.3 ) continue;
     //     if( mu.globalTrack()->normalizedChi2()>10) continue;
    
    }
                                       
  }


 const pat::MET & pfmet = mets->front();   

 h_met->Fill(pfmet.pt());

 if (goodEvent == true) h_good_met->Fill(pfmet.pt());
 if (goodEvent == false) h_bad_met->Fill(pfmet.pt());





}


void MuonEventDumper::printGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo) const
{

  // muon_hlt::GenInfo genInfo;
  
  // genInfo.trueNumberOfInteractions   = -1.;
  // genInfo.actualNumberOfInteractions = -1 ;    
  
  // std::vector<PileupSummaryInfo>::const_iterator puInfoIt  = puInfo->begin();
  // std::vector<PileupSummaryInfo>::const_iterator puInfoEnd = puInfo->end();

  // for(; puInfoIt != puInfoEnd; ++puInfoIt) 
  //   {
    
  //     int bx = puInfoIt->getBunchCrossing();
	  
  //     if(bx == 0) 
  // 	{ 
  // 	  genInfo.trueNumberOfInteractions   = puInfoIt->getTrueNumInteractions();
  // 	  genInfo.actualNumberOfInteractions = puInfoIt->getPU_NumInteractions();
  // 	  continue;
  // 	}
  //   }
  
  // event_.genInfos.push_back(genInfo);

  return; // to be implementded
  
}


void MuonEventDumper::printGenParticles(const edm::Handle<reco::GenParticleCollection> & genParticles) const
{
  
  // unsigned int gensize = genParticles->size();
  
  // // Do not record the initial protons
  // for (unsigned int i=0; i<gensize; ++i) 
  //   {

  //     const reco::GenParticle& part = genParticles->at(i);
    
  //     muon_hlt::GenParticle gensel;
  //     gensel.pdgId = part.pdgId();
  //     gensel.status = part.status();
  //     gensel.energy = part.energy();
  //     gensel.pt = part.pt();
  //     gensel.eta = part.eta();
  //     gensel.phi = part.phi();
  //     gensel.vx = part.vx();
  //     gensel.vy = part.vy();
  //     gensel.vz = part.vz();

  //     gensel.mothers.clear();
  //     unsigned int nMothers = part.numberOfMothers();

  //     for (unsigned int iMother=0; iMother<nMothers; ++iMother) 
  // 	{
  // 	  gensel.mothers.push_back(part.motherRef(iMother)->pdgId());
  // 	}

  //     // Protect agains bug in genParticles (missing mother => first proton)
  //     if (i>=2 && nMothers==0) gensel.mothers.push_back(0);
      
  //     event_.genParticles.push_back(gensel);
  //   }

  return; // to be implementded
}


void MuonEventDumper::printHlt(const edm::Handle<edm::TriggerResults> & triggerResults, 
				  const edm::Handle<trigger::TriggerEvent> & triggerEvent,
			       const edm::TriggerNames & triggerNames) const
{    

  std::cout << "[MuonEventDumper::printHlt]: "
	    << "\nPatsh from trigger results : " << std::endl;

    for (unsigned int iTrig=0; iTrig<triggerNames.size(); ++iTrig) 
    {
      
      if (triggerResults->accept(iTrig)) 
	{
	  std::string pathName = triggerNames.triggerName(iTrig);
 	  std::cout << "[PATH]: " << pathName << std::endl; 
	}
    }
      
  std::cout << "[MuonEventDumper::printHlt]: "
	    << "\nTrigger Objects : " << std::endl;

    const trigger::size_type nFilters(triggerEvent->sizeFilters());

  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
    {
	
      std::string filterTag = triggerEvent->filterTag(iFilter).encode();

      trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
      const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
	
      for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
	{  
	  trigger::size_type objKey = objectKeys.at(iKey);
	  const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
	  
	  float trigObjPt = triggerObj.pt();
	  float trigObjEta = triggerObj.eta();
	  float trigObjPhi = triggerObj.phi();	  
	  
 	  std::cout << "[OBJECT]: (pt, eta, phi) ("
		    << trigObjPt << ", "
		    << trigObjEta << ", "
		    << trigObjPhi << ") "
		    <<filterTag << std::endl; 

	}       
    }

}


void MuonEventDumper::printPV(const edm::Handle<std::vector<reco::Vertex> > & vertexes) const
{
      
  std::cout << "[MuonEventDumper::printPV]: "
	    << "\nFirst Vertex : " << std::endl;

    int nVtx = 0;

  std::vector<reco::Vertex>::const_iterator vertexIt  = vertexes->begin();
  std::vector<reco::Vertex>::const_iterator vertexEnd = vertexes->end();

  for (; vertexIt != vertexEnd; ++vertexIt) 
    {

      const reco::Vertex& vertex = *vertexIt;

      if (!vertex.isValid()) continue;
      ++nVtx;

      if (vertexIt == vertexes->begin()) 
	{
	  std::cout << "Coordinates (x, y, z) : ("
		    << vertex.x() << ", "
		    << vertex.y() << ", "
		    << vertex.z() << " )"
		    <<  std::endl;
	}
    }
  
  std::cout << "Number of vertices : " << nVtx << std::endl;
  
}

void MuonEventDumper::printTrack(const reco::Track * track, const std::string & trackType) const
{
  
  std::cout << "[From " << trackType << "] : pT :" << track->pt()
	    << " , sigma(pT)/pT :" << track->ptError()/track->pt()
	    << " , eta :" << track->eta()
	    << " , phi :" << track->phi()
	    << " , charge :" << track->charge()
	    << " , n.d.o.f. :" << track->ndof()
	    << " , chi2/n.d.o.f. :" << track->normalizedChi2()
	    << " , # pixel hits :" << track->hitPattern().numberOfValidPixelHits()
	    << " , # trk layers :" << track->hitPattern().trackerLayersWithMeasurement()
	    << " , # muon hits :" << track->hitPattern().numberOfMuonHits()
	    << " , # muon valid hits :" << track->hitPattern().numberOfValidMuonHits()
	    << " , # muon bad hits :" << track->hitPattern().numberOfBadMuonHits()
	    << " , # DT valid hits :" << track->hitPattern().numberOfValidMuonDTHits()
	    << " , # CSC valid hits :" << track->hitPattern().numberOfValidMuonCSCHits()
	    << " , # RPC valid hits :" << track->hitPattern().numberOfValidMuonRPCHits()
	    << " , # DT bad hits :" << track->hitPattern().numberOfBadMuonDTHits()
	    << " , # CSC bad hits :" << track->hitPattern().numberOfBadMuonCSCHits()
	    << " , # RPC bad hits :" << track->hitPattern().numberOfBadMuonRPCHits()
	    << " , # DT lost hits :" << track->hitPattern().numberOfLostMuonDTHits()
	    << " , # CSC lost hits :" << track->hitPattern().numberOfLostMuonCSCHits()
	    << " , # RPC lost hits :" << track->hitPattern().numberOfLostMuonRPCHits()
	    << " , # muon stations with valid hits :" << track->hitPattern().muonStationsWithValidHits()
	    << " , # DT stations with valid hits :" << track->hitPattern().dtStationsWithValidHits()
	    << " , # CSC stations with valid hits :" << track->hitPattern().cscStationsWithValidHits()
	    << " , # RPC stations with valid hits :" << track->hitPattern().rpcStationsWithValidHits()
	    << " , # muon stations with valid hits :" << track->hitPattern().muonStationsWithValidHits()
	    << " , # DT stations with any hits :" << track->hitPattern().dtStationsWithAnyHits()
	    << " , # CSC stations with any hits :" << track->hitPattern().cscStationsWithAnyHits()
	    << " , # RPC stations with any hits :" << track->hitPattern().rpcStationsWithAnyHits()
	    << std::endl;  

  // CB compte print of the hit pattern
  //track->hitPattern().print(reco::HitPattern::HitCategory::TRACK_HITS);
  
}

void MuonEventDumper::printMediumIDCuts(const reco::Muon & muon) const
{

  std::cout << "[MUON MEDIUM ID CUTS]: " 
	    << "\n\tIs Loose: " << muon::isLooseMuon(muon);

  if(muon::isLooseMuon(muon))
    {
      std::cout << "\n\tValid Trk Hits Fraction: " << muon.innerTrack()->validFraction()
		<< "\n\tIs Global: " << muon.isGlobalMuon();

      if (muon.isGlobalMuon())
	std::cout << "\n\tNorm. Glb. Chi2: " << muon.globalTrack()->normalizedChi2()
		  << "\n\tTrk-Sta position match: " << muon.combinedQuality().chi2LocalPosition
		  << "\n\tKink: " << muon.combinedQuality().trkKink;
      
      std::cout << "\n\tSegm Comp: " << muon::segmentCompatibility(muon);
    }
  
  std::cout << std::endl;

}


void MuonEventDumper::printMET(const edm::Handle<std::vector<pat::MET> > & mets) const
{

 const pat::MET & pfmet = mets->front();

 std::cout << "MET:: " << pfmet.pt() << std::endl;

}


void MuonEventDumper::printMuons(const edm::Handle<edm::View<reco::Muon> > & muons,
				    const edm::Handle<std::vector<reco::Vertex> > & vertexes,
				    const edm::Handle<reco::BeamSpot> & beamSpot) const
{

  std::cout << "[MuonEventDumper::printMuons]: " << std::endl;
  std::cout << "[MUON COLLECTION SIZE]: " << muons->size() << std::endl;

  edm::View<reco::Muon> ::const_iterator muonIt1  = muons->begin();
  edm::View<reco::Muon> ::const_iterator muonIt2  = muonIt1+1;

  edm::View<reco::Muon> ::const_iterator muonEnd = muons->end();

  float minDR = 999;
   edm::View<reco::Muon> ::const_iterator muA = muonIt1;
   edm::View<reco::Muon> ::const_iterator muB = muonIt1;

  if (muons->size() > 1) {

    for (; muonIt1 != muonEnd; ++muonIt1) {
      for (; muonIt2 != muonEnd; ++muonIt2) {
 
	if (muonIt1 == muonIt2) continue;

	const reco::Muon& mu1 = (*(muonIt1));
	const reco::Muon& mu2 = (*(muonIt2));

	float deta = mu1.eta() - mu2.eta();
	float dphi = mu1.phi() - mu2.phi();
	
	float dR = sqrt(deta*deta + dphi*dphi);
      
	if (dR < minDR) {
	  minDR = dR; 
	  muA = muonIt1;
	  muB = muonIt2;
	} 
      
      }
    }  
  }



  std::cout << " Min DR:: " << minDR << std::endl;
 

  if ( minDR < 0.2) {

    counterDuplicated++;

  edm::View<reco::Muon> ::const_iterator muonIt  = muons->begin();

  int idx = 0; 



  // check overlap between muons 

 
  // shared segmets 
  int nSharedSeg = muon::sharedSegments(*muA, *muB);

  std::cout << "Number of segments shared: " << nSharedSeg << std::endl;


  for (; muonIt != muonEnd; ++muonIt) 
    {


      if ( muonIt != muA && muonIt != muB ) continue;

      const reco::Muon& mu = (*muonIt);
      const reco::Vertex & vertex = vertexes->at(0); // CB for now vertex is always valid, but add a protection	    
      

     std::cout << "[MUON DETAILS]: " << idx << std::endl;

     int nSegments = mu.numberOfMatches(reco::Muon::SegmentArbitration); 

     std::cout << "How many segments?: " << nSegments << std::endl;

      //if ( !mu.isPFMuon() ) continue; 
      //if ( mu.muonBestTrack()->pt() < 100 ) continue;

     //     if( muon::segmentCompatibility(mu) < 0.3 ) continue;
     //if(nSegments==0) continue;
 
     if ( mu.innerTrack()->pt() >= 95 ) counter1++; 
     if ( mu.innerTrack()->pt() < 95 ) counter2++; 

      
      std::cout << "It is a : "
		<< (mu.isStandAloneMuon() ? "STANDALONE " : "")
		<< (mu.isGlobalMuon() ? "GLOBAL " : "")
	<< (mu.isTrackerMuon() ? "TRACKER " : "")
		<< (mu.isPFMuon() ? "PF " : "")
		<< "muon." << std::endl;
      

      std::cout << "With eta : " << mu.eta()
		<< " , phi : " << mu.phi()
		<< " , charge : " << mu.charge()
		<< std::endl ;

      std::cout << "With  " << mu.numberOfMatchedStations()
		<< " matched stations" << std::endl ; 

      if (mu.isGlobalMuon())
	std::cout << "With (dYX, dZ) for global track : ("
		  << mu.globalTrack()->dxy() << ", "
		  << mu.globalTrack()->dz() << ")"
		  << std::endl;
      
      if (!mu.innerTrack().isNull())
	std::cout << "With (dYX, dZ) for inner track : ("
		  << mu.innerTrack()->dxy() << ", "
		  << mu.innerTrack()->dz() << ")"
		  << std::endl;

      if ( mu.innerTrack()->pt() >= 100 ) { 
	if (!mu.outerTrack().isNull()) {
	 if(mu.outerTrack()->normalizedChi2()>10) counterFilter++;
	}}

      std::cout << "PF id variables: "
		<< "GLOBAL_TIGHT : "  << (isGlobalTightMuon(mu)  ? "TRUE  " : "FALSE  ")
		<< "TRACKER_TIGHT : " << (isTrackerTightMuon(mu) ? "TRUE  " : "FALSE  ")
		<< "ISOLATED : "      << (isIsolatedMuon(mu)     ? "TRUE  " : "FALSE  ")
		<< std::endl;


      std::cout << "With track parameters : \n";

      if (mu.isGlobalMuon())
	printTrack(mu.tunePMuonBestTrack().get(),"TUNEP");
      
      if (mu.isPFMuon())
	printTrack(mu.muonBestTrack().get(),"PF");

      if (!mu.innerTrack().isNull())
	printTrack(mu.innerTrack().get(),"INNER");

      if (!mu.globalTrack().isNull())
	printTrack(mu.globalTrack().get(),"GLB");

      if (!mu.pickyTrack().isNull() && mu.pickyTrack().isAvailable())
       	printTrack(mu.pickyTrack().get(),"PICKY");

      if (!mu.tpfmsTrack().isNull() && mu.tpfmsTrack().isAvailable())
	printTrack(mu.tpfmsTrack().get(),"TPFMS");

      if (!mu.dytTrack().isNull() && mu.dytTrack().isAvailable())
       	printTrack(mu.dytTrack().get(),"DYT");
      
      if (!mu.outerTrack().isNull())
	printTrack(mu.outerTrack().get(),"STANDALONE");
     
      if (mu.isTimeValid())
	{
	  const reco::MuonTime time = mu.time();
	  const reco::MuonTime rpcTime = mu.rpcTime();
                                                                                                               
          // Implement algorithm from Piotr on OOT muon rejection                                              
          bool veto = 0;                                                                                       
          bool cmbok =(time.nDof > 7);                                                                         
          bool rpcok =(rpcTime.nDof > 1 && rpcTime.timeAtIpInOutErr==0);                                       
          if (rpcok) {                                                                                         
            if ((fabs(rpcTime.timeAtIpInOut) > 10) && !(cmbok && fabs(time.timeAtIpInOut) < 10))              
              veto=1;                                                                                          
            else                                                                                               
              if (cmbok && (time.timeAtIpInOut > 20 || time.timeAtIpInOut < -45))                              
                veto=1;                                                                                        
          }        


	  std::cout << "[Muon TIME] : n.d.o.f. :" << time.nDof                                          
                      << " , time at IP IN->OUT :" << time.timeAtIpInOut                                       
                      << " , error on time at IP IN->OUT :" << time.timeAtIpInOutErr                           
                      << " , RPC time at IP IN->OUT:" << rpcTime.timeAtIpInOut                                 
                      << " , OOT veto :" << veto                                                               
                      << std::endl;

	  /*	    std::cout << "[Muon TIME] : n.d.o.f. :" << time.nDof
		      << " , time at IP IN->OUT :" << time.timeAtIpInOut
		      << " , error on time at IP IN->OUT :" << time.timeAtIpInOutErr
		      << std::endl;*/
	}      

      std::cout << "It passes the : "
		<< (muon::isSoftMuon(mu,vertex)     ? "SOFT " : "")
		<< (muon::isLooseMuon(mu)           ? "LOOSE " : "")
		<< (muon::isMediumMuon(mu)         ? "MEDIUM " : "")
		<< (muon::isTightMuon(mu,vertex)    ? "TIGHT " : "")
		<< (muon::isHighPtMuon(mu,vertex) ? "HIGH_PT " : "")
		<< " IDs." << std::endl;      

      printMediumIDCuts(mu);

      reco::MuonIsolation iso03 = mu.isolationR03();
      printIsolation(iso03,"0.3");
      
      reco::MuonPFIsolation isoPF03 = mu.pfIsolationR04();
      reco::MuonPFIsolation isoPF04 = mu.pfIsolationR04();

      printPFIsolation(isoPF03,"0.3");
      printPFIsolation(isoPF04,"0.4");

     idx++;

    }

  }
}

void MuonEventDumper::printIsolation( const reco::MuonIsolation & iso, const std::string & cone) const
{

  double sumPt   = iso.sumPt;
  double emEt    = iso.emEt;
  double hadEt   = iso.hadEt;
  int nTracks = iso.nTracks;
  int nJets   = iso.nJets;

  std::cout << "Detector Based Isolation variables are (R=" << cone << ") : "
	    << "\nsumPt: " << sumPt 
	    << "\nemEt : " << emEt
	    << "\nhadEt : " << hadEt
	    << "\nnTracks : " << nTracks
	    << "\nnJets: " << nJets
	    << std::endl;

}

void MuonEventDumper::printPFIsolation( const reco::MuonPFIsolation & iso, const std::string & cone) const
{

  double chargedHadronIso = iso.sumChargedHadronPt;
  double neutralHadronIso = iso.sumNeutralHadronEt;
  double photonIso        = iso.sumPhotonEt;
  double pu               = iso.sumPUPt;
      
  std::cout << "PF Isolation variables are (R=" << cone << ") : "
	    << "\nCh. Had (from vtx) : " << chargedHadronIso
	    << "\nNeu. Had : " << neutralHadronIso
	    << "\nPhot. : " << photonIso
	    << "\nPU. : " << pu
	    << std::endl;

}

bool MuonEventDumper::isGlobalTightMuon( const reco::Muon& muon ) const
{

 if ( !muon.isGlobalMuon() ) return false;
 if ( !muon.isStandAloneMuon() ) return false;
 
 
 if ( muon.isTrackerMuon() ) { 
   
   bool result = muon::isGoodMuon(muon,muon::GlobalMuonPromptTight);
 
   std::cout << "TEST: IS GLB prompt tight? " << result << std::endl;
  
   bool isTM2DCompatibilityTight =  muon::isGoodMuon(muon,muon::TM2DCompatibilityTight);   

   std::cout << "TEST: IS TM2DCompatibilityTight? " << isTM2DCompatibilityTight << std::endl;

   int nMatches = muon.numberOfMatches();

   bool quality = nMatches > 2 || isTM2DCompatibilityTight;

   std::cout << "TEST: nMatches" << nMatches  << std::endl;

   return result && quality;
   
 } else {
 
   reco::TrackRef standAloneMu = muon.standAloneMuon();
   
    // No tracker muon -> Request a perfect stand-alone muon, or an even better global muon
    bool result = false;
      
    // Check the quality of the stand-alone muon : 
    // good chi**2 and large number of hits and good pt error
    if ( ( standAloneMu->hitPattern().numberOfValidMuonDTHits() < 22 &&
	   standAloneMu->hitPattern().numberOfValidMuonCSCHits() < 15 ) ||
	 standAloneMu->normalizedChi2() > 10. || 
	 standAloneMu->ptError()/standAloneMu->pt() > 0.20 ) {
      result = false;

      std::cout << "TEST: STA quality?" << result  << std::endl;
   
    } else { 
      
      reco::TrackRef combinedMu = muon.combinedMuon();
      reco::TrackRef trackerMu = muon.track();
            
      // If the stand-alone muon is good, check the global muon
      if ( combinedMu->normalizedChi2() > standAloneMu->normalizedChi2() ) {
	// If the combined muon is worse than the stand-alone, it 
	// means that either the corresponding tracker track was not 
	// reconstructed, or that the sta muon comes from a late 
	// pion decay (hence with a momentum smaller than the track)
	// Take the stand-alone muon only if its momentum is larger
	// than that of the track
	result = standAloneMu->pt() > trackerMu->pt() ;
	std::cout << "TEST: check pt?" << result  << std::endl;
     } else { 
	// If the combined muon is better (and good enough), take the 
	// global muon
	result = 
	  combinedMu->ptError()/combinedMu->pt() < 
	  std::min(0.20,standAloneMu->ptError()/standAloneMu->pt());

	std::cout << "TEST: check pt error?" << result  << std::endl;  
    }
    }      

    return result;    
  }

  return false;

}

bool MuonEventDumper::isTrackerTightMuon( const reco::Muon& muon ) const
{

  if(!muon.isTrackerMuon()) return false;
  
  reco::TrackRef trackerMu = muon.track();
  const reco::Track& track = *trackerMu;
  
  unsigned nTrackerHits =  track.hitPattern().numberOfValidTrackerHits();
  
  if(nTrackerHits<=12) return false;
  
  bool isAllArbitrated = muon::isGoodMuon(muon,muon::AllArbitrated);
  
  bool isTM2DCompatibilityTight = muon::isGoodMuon(muon,muon::TM2DCompatibilityTight);
  
  if(!isAllArbitrated || !isTM2DCompatibilityTight)  return false;

  if((trackerMu->ptError()/trackerMu->pt() > 0.10)){
    //std::cout<<" PT ERROR > 10 % "<< trackerMu->pt() <<std::endl;
    return false;
  }
  return true;
  
}

bool MuonEventDumper::isIsolatedMuon( const reco::Muon& muon ) const
{

  if ( !muon.isIsolationValid() ) return false;
  
  // Isolated Muons which are missed by standard cuts are nearly always global+tracker
  if ( !muon.isGlobalMuon() ) return false;

  // If it's not a tracker muon, only take it if there are valid muon hits

  reco::TrackRef standAloneMu = muon.standAloneMuon();

  if ( !muon.isTrackerMuon() ){
    if(standAloneMu->hitPattern().numberOfValidMuonDTHits() == 0 &&
       standAloneMu->hitPattern().numberOfValidMuonCSCHits() ==0) return false;
  }
  
  // for isolation, take the smallest pt available to reject fakes

  reco::TrackRef combinedMu = muon.combinedMuon();
  double smallestMuPt = combinedMu->pt();
  
  if(standAloneMu->pt()<smallestMuPt) smallestMuPt = standAloneMu->pt();
  
  if(muon.isTrackerMuon())
    {
      reco::TrackRef trackerMu = muon.track();
      if(trackerMu->pt() < smallestMuPt) smallestMuPt= trackerMu->pt();
    }
     
  double sumPtR03 = muon.isolationR03().sumPt;
  
  double relIso = sumPtR03/smallestMuPt;

  if(relIso<0.1) return true;
  else return false;

}



#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuonEventDumper);
