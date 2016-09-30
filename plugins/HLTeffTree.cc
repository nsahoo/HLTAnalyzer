// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TNtuple.h>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtTotal.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416
using namespace std;

class HLTTree : public edm::EDAnalyzer {
public:
  explicit HLTTree(const edm::ParameterSet&);
  ~HLTTree();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
    
    TTree* HLTeffTree;
   
    //Tree info
    int Ntrkoffline;
    int NtrkFull;
    int NtrkPixel;
    int NvtxTrkOffline;
    int NvtxTrkFull;
    double ETT;
    double OfflineVtxX;
    double OfflineVtxY;
    double OfflineVtxZ;
    double HLTVtxX;
    double HLTVtxY;
    double HLTVtxZ;
    int NvtxOffline;
    int NvtxHLT;
    double OfflineLeadingPt;
    double HLTLeadingPt;
    
    //cuts for HighMultiplicity
    double min_Pt_; // min pt cut for number of tracks at HLT
    double max_Pt_; // max pt cut for number of tracks at HLT
    double max_Eta_; // max eta cut for number of tracks at HLT
    double max_Vz_; // max vz cut for number of tracks at HLT
    double min_sep_; // minimum separation in z of track-vertex in phi-eta for pixel tracks
    double min_sep_full_; // minimum separation in z of track-vertex in phi-eta for full tracks
    
    //cuts for HighPt track
    double pTerr_offline_;
    double pTerr_HLT_;
    double DCA_offline_;
    double DCA_HLT_;
    int nHit_offline_;
    int nHit_HLT_;
    double eta_offline_;
    double eta_HLT_;
    
    edm::InputTag pixelVertices_;
    edm::InputTag pixelTracks_;
    edm::InputTag fullVertices_;
    edm::InputTag fullTracks_;
    edm::InputTag HLTTrack_;
    
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    
    edm::EDGetTokenT<reco::VertexCollection> tok_onlinePixelVtx_;
    edm::EDGetTokenT<reco::RecoChargedCandidateCollection> tok_onlinePixelTrk_;
    
    edm::EDGetTokenT<reco::VertexCollection> tok_onlineFullVtx_;
    edm::EDGetTokenT<reco::RecoChargedCandidateCollection> tok_onlineFullTrk_;
    edm::EDGetTokenT<reco::TrackCollection> tok_HLTTrk_;
    
    edm::EDGetTokenT<l1t::EtSumBxCollection> tok_EtSum_Stage2_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HLTTree::HLTTree(const edm::ParameterSet& iConfig)
{

  //now do what ever initialization is needed
    min_Pt_ = iConfig.getUntrackedParameter<double>("min_Pt");
    max_Pt_ = iConfig.getUntrackedParameter<double>("max_Pt");
    max_Eta_ = iConfig.getUntrackedParameter<double>("max_Eta");
    max_Vz_ = iConfig.getUntrackedParameter<double>("max_Vz");
    min_sep_ = iConfig.getUntrackedParameter<double>("min_sep");
    min_sep_full_ = iConfig.getUntrackedParameter<double>("min_sep_full");
    
    pTerr_offline_ = iConfig.getUntrackedParameter<double>("pTerr_offline");
    pTerr_HLT_ = iConfig.getUntrackedParameter<double>("pTerr_HLT");
    eta_offline_ = iConfig.getUntrackedParameter<double>("eta_offline");
    eta_HLT_ = iConfig.getUntrackedParameter<double>("eta_HLT");
    DCA_offline_ = iConfig.getUntrackedParameter<double>("DCA_offline");
    DCA_HLT_ = iConfig.getUntrackedParameter<double>("DCA_HLT");
    nHit_offline_ = iConfig.getUntrackedParameter<int>("nHit_offline");
    nHit_HLT_ = iConfig.getUntrackedParameter<int>("nHit_HLT");

    
    double pTerr_offline_;
    double pTerr_HLT_;
    double DCA_offline_;
    double DCA_HLT_;
    int nHit_offline_;
    int nHit_HLT_;
    double eta_offline_;
    double eta_online_;

    
    pixelVertices_ = iConfig.getParameter<edm::InputTag>("PixelVertices");
    pixelTracks_ = iConfig.getParameter<edm::InputTag>("PixelTracks");
    fullVertices_ = iConfig.getParameter<edm::InputTag>("FullVertices");
    fullTracks_ = iConfig.getParameter<edm::InputTag>("FullTracks");
    HLTTrack_ = iConfig.getParameter<edm::InputTag>("HLTTrack");

    tok_offlinePV_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
    
    tok_onlinePixelVtx_ = consumes<reco::VertexCollection>(pixelVertices_);
    tok_onlinePixelTrk_ = consumes<reco::RecoChargedCandidateCollection>(pixelTracks_);
    
    tok_onlineFullVtx_ = consumes<reco::VertexCollection>(fullVertices_);
    tok_onlineFullTrk_ = consumes<reco::RecoChargedCandidateCollection>(fullTracks_);
    
    tok_HLTTrk_ = consumes<reco::TrackCollection>(HLTTrack_);
    
    tok_EtSum_Stage2_ = consumes<l1t::EtSumBxCollection>(edm::InputTag("hltCaloStage2Digis","EtSum"));
}


HLTTree::~HLTTree()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HLTTree::analyze(const edm::Event& iEvent, const edm::EventSetup& 
iSetup)
{
    using std::vector;
    using namespace edm;
    
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);

    edm::Handle<reco::TrackCollection> recotracks;
    iEvent.getByToken(tok_generalTrk_, recotracks);
    
    edm::Handle<reco::VertexCollection> vertexCollection;
    iEvent.getByToken( tok_onlinePixelVtx_, vertexCollection );
    
    edm::Handle<reco::RecoChargedCandidateCollection> trackCollection;
    iEvent.getByToken(tok_onlinePixelTrk_,trackCollection);
    
    edm::Handle<reco::VertexCollection> vertexCollectionFull;
    iEvent.getByToken( tok_onlineFullVtx_, vertexCollectionFull );
    
    edm::Handle<reco::RecoChargedCandidateCollection> trackCollectionFull;
    iEvent.getByToken(tok_onlineFullTrk_,trackCollectionFull);

    edm::Handle<reco::TrackCollection> tracksHLT;
    iEvent.getByToken(tok_HLTTrk_, tracksHLT);

    edm::Handle<l1t::EtSumBxCollection> L1EtSum;
    iEvent.getByToken(tok_EtSum_Stage2_, L1EtSum);
    
    //reco vtx
    double bestvz=0, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    NvtxOffline = -1;
    NvtxTrkOffline = -1;

    if(vertices.isValid())
    {
        const reco::Vertex & vtx = (*vertices)[0];
        bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
        bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
        NvtxOffline = vertices->size();
        NvtxTrkOffline = vtx.nTracks();
    }
    
    OfflineVtxX = bestvx;
    OfflineVtxY = bestvy;
    OfflineVtxZ = bestvz;
    
    //if(bestvz < -15.0 || bestvz>15.0) return;
    
    //ntrk offline
    Ntrkoffline = -1;
    OfflineLeadingPt = -1;
    
    if(recotracks.isValid())
    {
        for(unsigned it=0; it<recotracks->size(); ++it){
            
            const reco::Track & trk = (*recotracks)[it];
            
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzvtx = trk.dz(bestvtx);
            double dxyvtx = trk.dxy(bestvtx);
            double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
            double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
            
            if(!trk.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk.ptError())/trk.pt()>0.10) continue;
            if(fabs(dzvtx/dzerror) > 3) continue;
            if(fabs(dxyvtx/dxyerror) > 3) continue;
            
            double eta = trk.eta();
            double pt  = trk.pt();
            
            if(fabs(eta)>2.4) continue;
            if(pt<=0.4) continue;
            Ntrkoffline++;
        }
        
        for(unsigned i=0; i<recotracks->size(); ++i){
            const reco::Track & trk = (*recotracks)[i];
            
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzvtx = trk.dz(bestvtx);
            double dxyvtx = trk.dxy(bestvtx);
            double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
            double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
            
            if(!trk.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk.ptError())/trk.pt()>pTerr_offline_) continue;
            if(fabs(dzvtx/dzerror) > DCA_offline_) continue;
            if(fabs(dxyvtx/dxyerror) > DCA_offline_) continue;
            
            if(fabs(trk.eta())>eta_offline_) continue;
            
            if(trk.numberOfValidHits()<nHit_offline_) continue;
            
            if(!(trk.algo()==4 || trk.algo()==5 || trk.algo()==6 || trk.algo()==7 || trk.algo()==8 || trk.algo()==9 || trk.algo()==10 || trk.algo()==11)) continue;
            
            if(trk.pt()>OfflineLeadingPt) OfflineLeadingPt = trk.pt();
        }
    }

    //ntrk online pixel
    NtrkPixel = -1;
    double vzmax = -999;
    int nmax = 0;
    
    if(vertexCollection.isValid())
    {
        const reco::VertexCollection * vertices = vertexCollection.product();
        int npixelvertices = vertices->size();
        if (npixelvertices!=0)
        {
            reco::VertexCollection::const_iterator verticesItr;
            for (verticesItr=vertices->begin(); verticesItr!=vertices->end(); ++verticesItr)
            {
                int ntracksize = verticesItr->tracksSize();
                double vz = verticesItr->z();
                if(fabs(vz) > max_Vz_) continue;
                if( ntracksize > nmax)
                {
                    vzmax = vz;
                    nmax = ntracksize;
                }
            }

            if(trackCollection.isValid())
            {
                const reco::RecoChargedCandidateCollection * tracks = trackCollection.product();
                reco::RecoChargedCandidateCollection::const_iterator tracksItr;
                for (tracksItr=tracks->begin(); tracksItr!=tracks->end(); ++tracksItr)
                {
                    double eta = tracksItr->eta();
                    if(fabs(eta) > max_Eta_) continue;
                    double pt  = tracksItr->pt();
                    if(pt < min_Pt_ || pt > max_Pt_) continue;
                    double vz = tracksItr->vz();
                    if(fabs(vz-vzmax) > min_sep_) continue;
                    
                    NtrkPixel++;
                }
                
            }
        }
    }
    
    //ntrk online full
    NtrkFull = -1;
    NvtxTrkFull = -1;
    HLTLeadingPt = -1;
    HLTVtxX = -999;
    HLTVtxY = -999;
    HLTVtxZ = -999;
    NvtxHLT = -1;
    int   nmax_full = 0;
    double HLTVtxXerr = -999;
    double HLTVtxYerr = -999;
    double HLTVtxZerr = -999;
    double vxmax_full = -999;
    double vymax_full = -999;
    double vzmax_full = -999;
    
    if(vertexCollectionFull.isValid())
    {
        const reco::VertexCollection * vertices = vertexCollectionFull.product();
        int nvertices_full = vertices->size();
        NvtxHLT = nvertices_full;
        if (nvertices_full!=0)
        {
            reco::VertexCollection::const_iterator verticesItr;
            for (verticesItr=vertices->begin(); verticesItr!=vertices->end(); ++verticesItr)
            {
                int ntracksize = verticesItr->tracksSize();
                double vz = verticesItr->z();
                if(fabs(vz) > max_Vz_) continue;
                if( ntracksize > nmax_full)
                {
                    vzmax_full = vz;
                    nmax_full = ntracksize;
                    vxmax_full = verticesItr->x();
                    vymax_full = verticesItr->y();
                    
                    HLTVtxXerr = verticesItr->xError();
                    HLTVtxYerr = verticesItr->yError();
                    HLTVtxZerr = verticesItr->zError();
                }
            }
            
            NvtxTrkFull = nmax_full;
            HLTVtxX = vxmax_full;
            HLTVtxY = vymax_full;
            HLTVtxZ = vzmax_full;
            
            if(trackCollectionFull.isValid())
            {
                const reco::RecoChargedCandidateCollection * tracks = trackCollectionFull.product();
                reco::RecoChargedCandidateCollection::const_iterator tracksItr;
                for (tracksItr=tracks->begin(); tracksItr!=tracks->end(); ++tracksItr)
                {
                    double eta = tracksItr->eta();
                    if(fabs(eta) > max_Eta_) continue;
                    double pt  = tracksItr->pt();
                    if(pt < min_Pt_ || pt > max_Pt_) continue;
                    double vz = tracksItr->vz();
                    if(fabs(vz-vzmax) > min_sep_full_) continue;
                    
                    NtrkFull++;
                }
                
            }
            
            if(tracksHLT.isValid())
            {
                for(unsigned i=0; i<tracksHLT->size(); ++i){
                    const reco::Track & trk = (*tracksHLT)[i];
                    
                    math::XYZPoint bestvtx(HLTVtxX,HLTVtxY,HLTVtxZ);
                    
                    double dzvtx = trk.dz(bestvtx);
                    double dxyvtx = trk.dxy(bestvtx);
                    double dzerror = sqrt(trk.dzError()*trk.dzError()+HLTVtxZerr*HLTVtxZerr);
                    double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+HLTVtxXerr*HLTVtxYerr);
                    
                    if(!trk.quality(reco::TrackBase::highPurity)) continue;
                    if(fabs(trk.ptError())/trk.pt()>pTerr_HLT_) continue;
                    if(fabs(dzvtx/dzerror) > DCA_HLT_) continue;
                    if(fabs(dxyvtx/dxyerror) > DCA_HLT_) continue;
                    
                    if(fabs(trk.eta())>eta_HLT_) continue;
                    
                    if(trk.numberOfValidHits()<nHit_HLT_) continue;
                    
                    if(trk.pt()>HLTLeadingPt) HLTLeadingPt = trk.pt();
                }
            }
        }
    }
    
    //ETT
    ETT = -1;
    
    if(L1EtSum.isValid())
    {
        ETT = L1EtSum->begin()->pt();
    }

    HLTeffTree->Fill();
}


// ------------ method called once each job just before starting event
//loop  ------------
void 
HLTTree::beginJob()
{
    edm::Service<TFileService> fs;
        
    TH1D::SetDefaultSumw2();

    HLTeffTree = fs->make< TTree>("HLTeff","HLTeff");
    
    HLTeffTree->Branch("Ntrkoffline",&Ntrkoffline);
    HLTeffTree->Branch("NtrkFull",&NtrkFull);
    HLTeffTree->Branch("NtrkPixel",&NtrkPixel);
    HLTeffTree->Branch("NvtxTrkOffline",&NvtxTrkOffline);
    HLTeffTree->Branch("NvtxTrkFull",&NvtxTrkFull);
    HLTeffTree->Branch("TowerCount",&ETT);
    HLTeffTree->Branch("OfflineVtxX",&OfflineVtxX);
    HLTeffTree->Branch("OfflineVtxY",&OfflineVtxY);
    HLTeffTree->Branch("OfflineVtxZ",&OfflineVtxZ);
    HLTeffTree->Branch("HLTVtxX",&HLTVtxX);
    HLTeffTree->Branch("HLTVtxY",&HLTVtxY);
    HLTeffTree->Branch("HLTVtxZ",&HLTVtxZ);
    HLTeffTree->Branch("NvtxOffline",&NvtxOffline);
    HLTeffTree->Branch("NvtxHLT",&NvtxHLT);
    HLTeffTree->Branch("OfflineLeadingPt",&OfflineLeadingPt);
    HLTeffTree->Branch("HLTLeadingPt",&HLTLeadingPt);

}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
HLTTree::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTTree);



