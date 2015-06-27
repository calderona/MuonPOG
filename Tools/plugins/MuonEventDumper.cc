//////////////////////////////////////
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

#include <algorithm>
#include <iostream>

template<typename T> edm::Handle<T> conditionalGet(const edm::Event & ev,
						   const edm::InputTag & inputTag,
						   const std::string & collectionName)
{

  edm::Handle<T> collection ;

  if (inputTag.label() != "none") 
    {
      if (!ev.getByLabel(inputTag, collection)) 
	edm::LogError("") << "[MuonEventDumper]::conditionalGet: "
			  << collectionName << " " << inputTag.label() << " collection does not exist !!!";
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
  
private:
  
  void printGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &) const ;

  void printGenParticles(const edm::Handle<reco::GenParticleCollection> &) const ;

  void printHlt(const edm::Handle<edm::TriggerResults> &, 
		const edm::Handle<trigger::TriggerEvent> &,
		const edm::TriggerNames &) const ;
  
  void printPV(const edm::Handle<std::vector<reco::Vertex> > &) const ;
  
  void printMuons(const edm::Handle<reco::MuonCollection> &,
		  const edm::Handle<std::vector<reco::Vertex> > &,
		  const edm::Handle<reco::BeamSpot> &) const ;

  void printTrack(const reco::Track * track,
		  const std::string & trackType) const;

		  // CB copy of PF muon selections as they are in 740
  bool isIsolatedMuon(const reco::Muon& muon) const ;
  bool isGlobalTightMuon( const reco::Muon& muon ) const ; 
  bool isTrackerTightMuon( const reco::Muon& muon ) const ;
  
  edm::InputTag trigResultsTag_;
  edm::InputTag trigSummaryTag_;

  edm::InputTag muonTag_;
  edm::InputTag primaryVertexTag_;
  edm::InputTag beamSpotTag_;

  edm::InputTag genTag_;
  edm::InputTag pileUpInfoTag_;

};


MuonEventDumper::MuonEventDumper( const edm::ParameterSet & cfg ) :
  // Input collections
  trigResultsTag_(cfg.getUntrackedParameter<edm::InputTag>("TrigResultsTag", edm::InputTag("TriggerResults::HLT"))),
  trigSummaryTag_(cfg.getUntrackedParameter<edm::InputTag>("TrigSummaryTag", edm::InputTag("hltTriggerSummaryAOD::HLT"))),

  muonTag_(cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"))),
  primaryVertexTag_(cfg.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag", edm::InputTag("offlinePrimaryVertices"))),
  beamSpotTag_(cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag", edm::InputTag("offlineBeamSpot"))),

  genTag_(cfg.getUntrackedParameter<edm::InputTag>("GenTag", edm::InputTag("prunedGenParticles"))),
  pileUpInfoTag_(cfg.getUntrackedParameter<edm::InputTag>("PileUpInfoTag", edm::InputTag("pileupInfo")))
{

}


void MuonEventDumper::beginJob() 
{  

}


void MuonEventDumper::beginRun(const edm::Run & run, const edm::EventSetup & config )
{
  
}


void MuonEventDumper::endJob() 
{

}

void MuonEventDumper::analyze (const edm::Event & ev, const edm::EventSetup &)
{

  std::cout << "[MuonEventDumper::analyze]: "
	    << "\nProcessing Run: " << ev.id().run()
	    << "\tLumiBlock: " << ev.id().luminosityBlock()
	    << "\tevent: " << ev.id().event() << std::endl;

  // Print GEN pile up information
  if (!ev.isRealData()) 
    {
      auto puInfo = conditionalGet<std::vector<PileupSummaryInfo> >(ev,pileUpInfoTag_,"Pile-Up Info");
      if (puInfo.isValid()) printGenInfo(puInfo);
    }
  

  // Print GEN particles information
  if (!ev.isRealData()) 
    {
      auto genParticles= conditionalGet<reco::GenParticleCollection>(ev,genTag_,"GenParticle Collection");
      if (genParticles.isValid()) printGenParticles(genParticles);
    } 
  
  // Print trigger information
  auto triggerResults = conditionalGet<edm::TriggerResults>(ev,trigResultsTag_,"TriggerResults");
  auto triggerEvent   = conditionalGet<trigger::TriggerEvent>(ev,trigSummaryTag_,"TriggerEvent");
      
  if (triggerResults.isValid() && triggerEvent.isValid()) 
    printHlt(triggerResults, triggerEvent,ev.triggerNames(*triggerResults));
  
  // Print vertex information
  auto vertexes = conditionalGet<std::vector<reco::Vertex> >(ev,primaryVertexTag_,"Vertex");
  if (vertexes.isValid()) printPV(vertexes);
  
  // Get beam spot for muons
  auto beamSpot = conditionalGet<reco::BeamSpot>(ev,beamSpotTag_,"BeamSpot");

  // Get muons  
  auto muons = conditionalGet<reco::MuonCollection>(ev,muonTag_,"MuonCollection");

  // Print muon information
  if (muons.isValid() && vertexes.isValid() && beamSpot.isValid()) 
    printMuons(muons,vertexes,beamSpot);
  
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
		    << vertex.x() << " )"
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
	    << " , # muon hits :" << track->hitPattern().numberOfValidMuonHits()
	    << " , # muon valid hits :" << track->hitPattern().numberOfMuonHits()
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
  
}

void MuonEventDumper::printMuons(const edm::Handle<reco::MuonCollection> & muons,
				    const edm::Handle<std::vector<reco::Vertex> > & vertexes,
				    const edm::Handle<reco::BeamSpot> & beamSpot) const
{

  std::cout << "[MuonEventDumper::printMuons]: " << std::endl;


  reco::MuonCollection::const_iterator muonIt  = muons->begin();
  reco::MuonCollection::const_iterator muonEnd = muons->end();

  for (; muonIt != muonEnd; ++muonIt) 
    {

      const reco::Muon& mu = (*muonIt);
      const reco::Vertex & vertex = vertexes->at(0); // CB for now vertex is always valid, but add a protection	    

      if (!(mu.isStandAloneMuon() && 
	    mu.outerTrack().get()->pt()>2000. && 
	    (mu.isTrackerMuon() || mu.isGlobalMuon() )
	    )
	  )
	continue;

      std::cout << "[MUON DETAILS]: " << std::endl;
      
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

      if (!mu.pickyTrack().isNull())
	printTrack(mu.pickyTrack().get(),"PICKY");

      if (!mu.tpfmsTrack().isNull())
	printTrack(mu.tpfmsTrack().get(),"TPFMS");

      if (!mu.dytTrack().isNull())
	printTrack(mu.dytTrack().get(),"DYT");
      
      if (!mu.outerTrack().isNull())
	printTrack(mu.outerTrack().get(),"STANDALONE");
     
      if (mu.isTimeValid())
	{
	  const reco::MuonTime time = mu.time();
	    std::cout << "[Muon TIME] : n.d.o.f. :" << time.nDof
		      << " , time at IP IN->OUT :" << time.timeAtIpInOut
		      << " , error on time at IP IN->OUT :" << time.timeAtIpInOutErr
		      << std::endl;
	}
      
      std::cout << "It passes the : "
		<< (muon::isSoftMuon(mu,vertex)     ? "SOFT " : "")
		<< (muon::isLooseMuon(mu)           ? "LOOSE " : "")
		<< (muon::isTightMuon(mu,vertex)    ? "TIGHT " : "")
		<< (muon::isHighPtMuon(mu,vertex) ? "HIGH_PT " : "")
		<< " IDs." << std::endl;      

      reco::MuonIsolation iso03 = mu.isolationR03();
      double sumPt   = iso03.sumPt;
      double emEt    = iso03.emEt;
      double hadEt   = iso03.hadEt;
      int nTracks = iso03.nTracks;
      int nJets   = iso03.nJets;

      std::cout << "Detector Based Isolation variables are (R=0.3) : "
		<< "\nsumPt: " << sumPt 
		<< "\nemEt : " << emEt
		<< "\nhadEt : " << hadEt
		<< "\nnTracks : " << nTracks
		<< "\nnJets: " << nJets
		<< std::endl;
      
      reco::MuonPFIsolation iso04 = mu.pfIsolationR04();
      double chargedHadronIso = iso04.sumChargedHadronPt;
      double neutralHadronIso = iso04.sumNeutralHadronEt;
      double photonIso        = iso04.sumPhotonEt;
      double pu               = iso04.sumPUPt;
      
      double relIsodBeta = (iso04.sumChargedHadronPt+ std::max(0.,iso04.sumPhotonEt+iso04.sumNeutralHadronEt - 0.5*iso04.sumPUPt)) / mu.pt();

      std::cout << "PF Isolation variables are (R=0.4) : "
		<< "\nCh. Had (from vtx) : " << chargedHadronIso
		<< "\nNeu. Had : " << neutralHadronIso
		<< "\nPhot. : " << photonIso
		<< "\nPU. : " << pu
		<< "\nRelative isolation (combined, dBeta) : " << relIsodBeta
		<< std::endl;

    }

}

bool MuonEventDumper::isGlobalTightMuon( const reco::Muon& muon ) const
{

 if ( !muon.isGlobalMuon() ) return false;
 if ( !muon.isStandAloneMuon() ) return false;
 
 
 if ( muon.isTrackerMuon() ) { 
   
   bool result = muon::isGoodMuon(muon,muon::GlobalMuonPromptTight);
   
   bool isTM2DCompatibilityTight =  muon::isGoodMuon(muon,muon::TM2DCompatibilityTight);   
   int nMatches = muon.numberOfMatches();
   bool quality = nMatches > 2 || isTM2DCompatibilityTight;
   
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
     } else { 
	// If the combined muon is better (and good enough), take the 
	// global muon
	result = 
	  combinedMu->ptError()/combinedMu->pt() < 
	  std::min(0.20,standAloneMu->ptError()/standAloneMu->pt());
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
