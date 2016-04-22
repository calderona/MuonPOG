#ifndef MuonPOG_Tools_MuonPogTree_H
#define MuonPOG_Tools_MuonPogTree_H

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


namespace muon_pog {


  class GenInfo {
  public:
    Float_t trueNumberOfInteractions;   // Number of simultaneous interactions generated (before poissonian ev by ev smearing)
    Int_t   actualNumberOfInteractions; // Number of simultaneous interactions generated (after poissonian ev by ev smearing)
    Float_t genWeight;
    GenInfo(){};
    virtual ~GenInfo(){};
    
    ClassDef(GenInfo,1)
  };
  
  class GenParticle {
  public:
    
    GenParticle(){};
    virtual ~GenParticle(){};
    
    Int_t pdgId;  // PDG identifier
    Int_t status; // MC status
    Float_t energy; // energy [GeV]
    Float_t pt; // pt [GeV]
    Float_t eta; // eta
    Float_t phi; // phi
    Float_t vx; // x coordinate of production vertex [cm]
    Float_t vy; // y coordinate of production vertex [cm]
    Float_t vz; // z coordinate of production vertex [cm]
    std::vector<Int_t> mothers; // vector of indices of mothers
    std::vector<bool>  flags;   // vector of flags, in the same order of
                                //  of "DataFormats/HepMCCandidate/interface/GenStatusFlag.h"

    void setFlags(const reco::GenStatusFlags & statusFlags)
    {
      flags.clear();
      for (unsigned int flag = 0; flag < 15; flag++) 
	flags.push_back(statusFlags.flags_[flag]);      
    }
    
    
  private:
    
    ClassDef(GenParticle,1)
  };

  class METs {
  public:
    Float_t pfMet;   // raw PF MET [GeV]
    Float_t pfChMet; // raw PF charged MET [GeV]
    Float_t caloMet; // raw Calo MET [GeV]

    METs(){};
    virtual ~METs(){};
    
    ClassDef(METs,1)
  };

  enum MuonDetType { DT=0, CSC, RPC };

  class ChambMatch {
  public:
    Int_t r;   // station/disk
    Int_t phi; // sector
    Int_t eta; // ring/wheel
    
    MuonDetType type;
    
    Float_t dx;  // 999999 if not matched with a segment (I think) 
    Float_t dy;  // 999999 if not matched with a segment (I think)
    
    Float_t errxTk; 
    Float_t erryTk; 
    
    Float_t errxSeg;  // 999999 if not matched with a segment (I think)
    Float_t errySeg;  // 999999 if not matched with a segment (I think) 
    
    ChambMatch(){};
    virtual ~ChambMatch(){};
    
    ClassDef(ChambMatch,1)
  };

  class HitInfo {
  public:
    Int_t r; // station/disk
    Int_t phi; // sector
    Int_t eta;   // ring/wheel
    
    MuonDetType type;

    Int_t nHits; 
    Int_t nHitsPhi; 
    Int_t nHitsTheta; 

    HitInfo(){};
    virtual ~HitInfo(){};
    
    ClassDef(HitInfo,1)
  };

  class Muon {
  public:

    Float_t pt;  // pt [GeV]   
    Float_t eta; // eta
    Float_t phi; // phi

    Int_t   charge;    // charge

    Float_t pt_tuneP;  // pt [GeV]
    Float_t eta_tuneP; // eta
    Float_t phi_tuneP; // phi

    Int_t   charge_tuneP;    // charge

    Float_t pt_global;  // pt [GeV]
    Float_t eta_global; // eta
    Float_t phi_global; // phi

    Int_t   charge_global;    // charge

    Float_t pt_tracker;  // pt [GeV]
    Float_t eta_tracker; // eta
    Float_t phi_tracker; // phi

    Int_t   charge_tracker;    // charge

    Float_t pt_standalone;  // pt [GeV]
    Float_t eta_standalone; // eta
    Float_t phi_standalone; // phi

    Int_t   charge_standalone;    // charge

    Int_t   isGlobal;
    Int_t   isTracker;
    Int_t   isTrackerArb;
    Int_t   isRPC;
    Int_t   isStandAlone;
    Int_t   isPF;

    Int_t   isSoft;
    Int_t   isLoose;
    Int_t   isTight;
    Int_t   isMedium;
    Int_t   isHighPt;
    
    //Detector Based Isolation
    Float_t trackerIso;
    Float_t EMCalIso;
    Float_t HCalIso;

    // PF Isolation
    Float_t chargedHadronIso;
    Float_t chargedHadronIsoPU;
    Float_t photonIso;
    Float_t neutralHadronIso;


    Float_t isoPflow04; // PF isolation in dR<0.4 cone dBeta
    Float_t isoPflow03; // PF isolation in dR<0.3 cone dBeta

    Float_t dxy;       // signed transverse distance to primary vertex [cm]
    Float_t dz;        // signed longitudinal distance to primary vertex at min. transv. distance [cm]
    Float_t edxy;      // uncertainty on dxy [cm]
    Float_t edz;       // uncertainty on dz [cm]
    Float_t dxybs;     // signed transverse distance to beamspot [cm]
    Float_t dzbs;      // signed longitudinal distance to beamspot [cm]

    Int_t   nHitsGlobal;
    Int_t   nHitsTracker;
    Int_t   nHitsStandAlone; 

    // Variables for ID 
    //  - General (Tight, HighPt, Soft) 
    Float_t glbNormChi2; 
    Float_t trkNormChi2; 
    Int_t   trkMuonMatchedStations; 
    Int_t   glbMuonValidHits; 
    Int_t   trkPixelValidHits; 
    Int_t   trkPixelLayersWithMeas; 
    Int_t   trkTrackerLayersWithMeas; 

    //  - HighPt 
    Float_t bestMuPtErr; 

    //  - Medium 
    Float_t trkValidHitFrac; 
    Float_t trkStaChi2; 
    Float_t trkKink; 
    Float_t muSegmComp; 

    //  - Soft 
    Int_t   isTrkMuOST; 
    Int_t   isTrkHP; 

    Float_t dxyBest; 
    Float_t dzBest; 
    Float_t dxyInner; 
    Float_t dzInner; 

    // Muon time 
    Float_t muonTimeDof; 
    Float_t muonTime; 
    Float_t muonTimeErr;

    std::vector<HitInfo> hits;
    std::vector<ChambMatch> matches;

    Muon(){};
    virtual ~Muon(){};

    ClassDef(Muon,2)
  };

  class HLTObject {
  public:

    std::string filterTag; // name of filter passed by the object
    Float_t pt;            // pt of the object passing the filter [GeV]
    Float_t eta;           // eta of the object passing the filter
    Float_t phi;           // phi of the object passing the filter
    
    HLTObject(){};
    virtual ~HLTObject(){};

    ClassDef(HLTObject,1)

  };

  class HLT {
  public:
    std::vector<std::string> triggers; // vector of strings with HLT paths
    std::vector<muon_pog::HLTObject>   objects;  // vector of hlt objects assing filters

    HLT(){};
    virtual ~HLT(){};
    bool match( const std::string & path ) {
      if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )
	return true;
      
      return false;
    }

    bool find( const std::string & path ) {
      for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
	if ( it->find ( path ) != std::string::npos ) return true;
      }
      return false;
    }

    ClassDef(HLT,1)

  };

  class EventId {
  public:

    Int_t runNumber;             // run number
    Int_t luminosityBlockNumber; // luminosity block number
    Int_t eventNumber;           // event number

    EventId(){};
    virtual ~EventId(){};

    ClassDef(EventId,1)
  };

  class Event {
  public:

    Int_t runNumber;             // run number
    Int_t luminosityBlockNumber; // luminosity block number
    Int_t eventNumber;           // event number

    Int_t bxId;                  // bunch crossing number
    unsigned long long orbit;    // orbit number
    Float_t instLumi;            // inst lumi from scalers [10E30]

    Int_t nVtx;                      // number of valid reconstructed primary vertices 
    Float_t primaryVertex[3];        // 3d coordinates of PV [cm]
    Float_t cov_primaryVertex[3][3]; // 3x3 covariance matrix of PV estimation [cm*cm]

    std::vector <muon_pog::GenInfo> genInfos;        // venctor of genInfos; size=0 in data
    std::vector<muon_pog::GenParticle> genParticles; // venctor of genParticles size=0 in data
    std::vector<muon_pog::Muon> muons; // vector of muons
    muon_pog::METs mets;  // vector of different MET definitions 
    muon_pog::HLT hlt;                 // HLT objects

    Event(){};
    virtual ~Event(){};

    ClassDef(Event,3)
  };

}
#endif
