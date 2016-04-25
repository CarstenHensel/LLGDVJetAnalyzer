// -*- C++ -*-
//
// Package:    LLGDVJetAnalyzer/LLGDVJetAnalyzer
// Class:      LLGDVJetAnalyzer
// 
/**\class LLGDVJetAnalyzer LLGDVJetAnalyzer.cc LLGDVJetAnalyzer/LLGDVJetAnalyzer/plugins/LLGDVJetAnalyzer.cc

 Description: 
 Simple analysis class to dump output in a plain rootfile

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Hamer
//         Created:  Tue, 13 Jan 2015 17:58:12 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TFile.h"
#include "TTree.h"

std::vector<double> CalculateVertex( std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> weight, std::vector<int> charge, std::vector<double> distance, int &nConsidered, double &weightednConsidered, std::vector<double> &error, double &maxScore );
//
//
// class declaration
//

class LLGDVJetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit LLGDVJetAnalyzer(const edm::ParameterSet&);
      ~LLGDVJetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::JetCollection> jetTokennoCHS_;
      edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> secVtxToken_;
      
      // the output file and tree
      TFile *fOutput = new TFile("RecoOutput.root", "RECREATE");
      TTree *tOutput = new TTree("RecoData", "RecoData");


      // missing transverse energy
      
      // the CHS jet variables
      std::vector<double> *jet_eta = new std::vector<double>;
      std::vector<double> *jet_phi = new std::vector<double>;
      std::vector<double> *jet_pt = new std::vector<double>;

      // the CHS jet constituents 
      std::vector<std::vector<double> >* jet_constVertex_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertex_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertex_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertexRef_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertexRef_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertexRef_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pt      = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_eta     = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_phi     = new std::vector<std::vector<double> >;
      std::vector<std::vector<int> >*    jet_const_charge  = new std::vector<std::vector<int> >;
      std::vector<std::vector<double> >* jet_const_pca0_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pca0_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pca0_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_dxy = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_dz = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_d = new std::vector<std::vector<double> >;
      
      // the nonCHS jet variables
      std::vector<double> *jetnoCHS_eta = new std::vector<double>;
      std::vector<double> *jetnoCHS_phi = new std::vector<double>;
      std::vector<double> *jetnoCHS_pt = new std::vector<double>;

      // the noCHS jet constituents 
      std::vector<std::vector<double> >* jetnoCHS_constVertex_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_constVertex_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_constVertex_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_constVertexRef_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_constVertexRef_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_constVertexRef_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_pt      = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_eta     = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_phi     = new std::vector<std::vector<double> >;
      std::vector<std::vector<int> >*    jetnoCHS_const_charge  = new std::vector<std::vector<int> >;
      std::vector<std::vector<double> >* jetnoCHS_const_pca0_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_pca0_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_pca0_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_closestVertex_dxy = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_closestVertex_dz = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jetnoCHS_const_closestVertex_d = new std::vector<std::vector<double> >;

      // the primary vertex information
      std::vector<double> *vertex_x = new std::vector<double>;
      std::vector<double> *vertex_y = new std::vector<double>;
      std::vector<double> *vertex_z = new std::vector<double>;
      std::vector<double> *vertex_nTracks = new std::vector<double>;
      std::vector<double> *vertex_pt = new std::vector<double>;
      std::vector<double> *vertex_ndof = new std::vector<double>;
      std::vector<double> *vertex_d0 = new std::vector<double>;
      std::vector<double> *vertex_dx = new std::vector<double>;
      std::vector<double> *vertex_dy = new std::vector<double>;
      std::vector<double> *vertex_dz = new std::vector<double>;

      // the secondary vertices
      std::vector<double> *secVertex_x = new std::vector<double>;
      std::vector<double> *secVertex_y = new std::vector<double>;
      std::vector<double> *secVertex_z = new std::vector<double>;
      std::vector<double> *secVertex_pt = new std::vector<double>;
      std::vector<double> *secVertex_ndof = new std::vector<double>;
      std::vector<double> *secVertex_chi2 = new std::vector<double>;
      std::vector<double> *secVertex_dx = new std::vector<double>;
      std::vector<double> *secVertex_dy = new std::vector<double>;
      std::vector<double> *secVertex_dz = new std::vector<double>;


      // event metadata
      unsigned int RunNumber = 0;
      unsigned long long EventNumber = 0;
      unsigned int LuminosityBlock = 0;
      
      double generatorWeight = 0.;

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
LLGDVJetAnalyzer::LLGDVJetAnalyzer(const edm::ParameterSet& iConfig):
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetTokennoCHS_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsnochs"))),
  genEvtInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GenEventInfo") )),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  secVtxToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices")))
  {
   //now do what ever initialization is needed
   
   // define the triggers we (might) want to use
   // these seem to be interesting for us
   /* triggerNames->push_back( "HLT_PFJet260_v1" );
    * triggerNames->push_back( "HLT_JetE30_NoBPTX_v1" );
    * triggerNames->push_back( "HLT_JetE30_NoBPTX3BX_NoHalo_v1" );
    * triggerNames->push_back( "HLT_JetE50_NoBPTX3BX_NoHalo_v1" );
    * triggerNames->push_back( "HLT_JetE70_NoBPTX3BX_NoHalo_v1" );
   */
   
   // set the output branches for the tree
   tOutput -> Branch("RecoCHSJet_eta", &jet_eta );
   tOutput -> Branch("RecoCHSJet_phi", &jet_phi );
   tOutput -> Branch("RecoCHSJet_pt", &jet_pt );
   
   tOutput -> Branch("RecoCHSJet_constVertex_x", &jet_constVertex_x );
   tOutput -> Branch("RecoCHSJet_constVertex_y", &jet_constVertex_y );
   tOutput -> Branch("RecoCHSJet_constVertex_z", &jet_constVertex_z );
   tOutput -> Branch("RecoCHSJet_const_pt", &jet_const_pt );
   tOutput -> Branch("RecoCHSJet_const_charge", &jet_const_charge );
   tOutput -> Branch("RecoCHSJet_const_pca0_x", &jet_const_pca0_x );
   tOutput -> Branch("RecoCHSJet_const_pca0_y", &jet_const_pca0_y );
   tOutput -> Branch("RecoCHSJet_const_pca0_z", &jet_const_pca0_z );
   tOutput -> Branch("RecoCHSJet_const_closestVertex_dxy", &jet_const_closestVertex_dxy );
   tOutput -> Branch("RecoCHSJet_const_closestVertex_dz", &jet_const_closestVertex_dz );
   tOutput -> Branch("RecoCHSJet_const_closestVertex_d", &jet_const_closestVertex_d );
   tOutput -> Branch("RecoCHSJet_const_eta", &jet_const_eta );
   tOutput -> Branch("RecoCHSJet_const_phi", &jet_const_phi );
   
   tOutput -> Branch("RecoNoCHSJet_eta", &jetnoCHS_eta );
   tOutput -> Branch("RecoNoCHSJet_phi", &jetnoCHS_phi );
   tOutput -> Branch("RecoNoCHSJet_pt", &jetnoCHS_pt );
   
   tOutput -> Branch("RecoNoCHSJet_constVertex_x", &jetnoCHS_constVertex_x );
   tOutput -> Branch("RecoNoCHSJet_constVertex_y", &jetnoCHS_constVertex_y );
   tOutput -> Branch("RecoNoCHSJet_constVertex_z", &jetnoCHS_constVertex_z );
   tOutput -> Branch("RecoNoCHSJet_const_pt", &jetnoCHS_const_pt );
   tOutput -> Branch("RecoNoCHSJet_const_charge", &jetnoCHS_const_charge );
   tOutput -> Branch("RecoNoCHSJet_const_pca0_x", &jetnoCHS_const_pca0_x );
   tOutput -> Branch("RecoNoCHSJet_const_pca0_y", &jetnoCHS_const_pca0_y );
   tOutput -> Branch("RecoNoCHSJet_const_pca0_z", &jetnoCHS_const_pca0_z );
   tOutput -> Branch("RecoNoCHSJet_const_closestVertex_dxy", &jetnoCHS_const_closestVertex_dxy );
   tOutput -> Branch("RecoNoCHSJet_const_closestVertex_dz", &jetnoCHS_const_closestVertex_dz );
   tOutput -> Branch("RecoNoCHSJet_const_closestVertex_d", &jetnoCHS_const_closestVertex_d );
   tOutput -> Branch("RecoNoCHSJet_const_eta", &jetnoCHS_const_eta );
   tOutput -> Branch("RecoNoCHSJet_const_phi", &jetnoCHS_const_phi );
    
    

   tOutput -> Branch("RecoVertex_x", &vertex_x );
   tOutput -> Branch("RecoVertex_y", &vertex_y );
   tOutput -> Branch("RecoVertex_z", &vertex_z );
   tOutput -> Branch("RecoVertex_ndof", &vertex_ndof ); 
   tOutput -> Branch("RecoVertex_d0", &vertex_d0 );
   tOutput -> Branch("RecoVertex_xError", &vertex_dx );
   tOutput -> Branch("RecoVertex_yError", &vertex_dy );
   tOutput -> Branch("RecoVertex_zError", &vertex_dz );
   tOutput -> Branch("RecoVertex_nTracks", &vertex_nTracks );
   tOutput -> Branch("RecoVertex_pt", &vertex_pt );
   
   tOutput -> Branch("RecoSecVertex_x", &secVertex_x );
   tOutput -> Branch("RecoSecVertex_y", &secVertex_y );
   tOutput -> Branch("RecoSecVertex_z", &secVertex_z );
   tOutput -> Branch("RecoSecVertex_ndof", &secVertex_ndof );
   tOutput -> Branch("RecoSecVertex_chi2", &secVertex_chi2 );
   tOutput -> Branch("RecoSecVertex_pt", &secVertex_pt );
   tOutput -> Branch("RecoSecVertex_xError", &secVertex_dx );
   tOutput -> Branch("RecoSecVertex_yError", &secVertex_dy );
   tOutput -> Branch("RecoSecVertex_zError", &secVertex_dz );
 
   tOutput -> Branch("EventNumber", &EventNumber );
   tOutput -> Branch("RunNumber", &RunNumber );
   tOutput -> Branch("LuminosityBlock", &LuminosityBlock );
   tOutput -> Branch("GeneratorWeight", &generatorWeight );

}



LLGDVJetAnalyzer::~LLGDVJetAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
LLGDVJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  
   // clear all variables
   
   jet_eta->clear();
   jet_phi->clear();
   jet_pt->clear();
   jet_constVertex_x->clear();
   jet_constVertex_y->clear();
   jet_constVertex_z->clear();
   jet_const_pt->clear();
   jet_const_eta->clear();
   jet_const_phi->clear();
   jet_const_charge->clear();
   jet_const_pca0_x->clear();
   jet_const_pca0_y->clear();
   jet_const_pca0_z->clear();
   jet_const_closestVertex_dxy->clear();
   jet_const_closestVertex_dz->clear();
   jet_const_closestVertex_d->clear();
   
   jetnoCHS_eta->clear();
   jetnoCHS_phi->clear();
   jetnoCHS_pt->clear();
   jetnoCHS_constVertex_x->clear();
   jetnoCHS_constVertex_y->clear();
   jetnoCHS_constVertex_z->clear();
   jetnoCHS_const_pt->clear();
   jetnoCHS_const_eta->clear();
   jetnoCHS_const_phi->clear();
   jetnoCHS_const_charge->clear();
   jetnoCHS_const_pca0_x->clear();
   jetnoCHS_const_pca0_y->clear();
   jetnoCHS_const_pca0_z->clear();
   jetnoCHS_const_closestVertex_dxy->clear();
   jetnoCHS_const_closestVertex_dz->clear();
   jetnoCHS_const_closestVertex_d->clear();
   
   vertex_x -> clear();
   vertex_y -> clear();
   vertex_z -> clear();
   vertex_d0 -> clear();
   vertex_dx -> clear();
   vertex_dy -> clear();
   vertex_dz -> clear();
   vertex_nTracks -> clear();
   vertex_pt -> clear();
   vertex_ndof -> clear();
   secVertex_x -> clear();
   secVertex_y -> clear();
   secVertex_z -> clear();
   secVertex_ndof -> clear();
   secVertex_chi2 -> clear();
   secVertex_pt -> clear();
   secVertex_dx -> clear();
   secVertex_dy -> clear();
   secVertex_dz -> clear();
   
   generatorWeight = 0.;


   //edm::Handle<LHEEventProduct> lheEventProduct;
   //iEvent.getByToken( lheEventToken_, lheEventProduct );
   //std::cout << "compare " << generatorWeight << " with " << lheEventProduct->hepeup().XWGTUP << std::endl;
     
   edm::Handle<GenEventInfoProduct> genEvtInfo; 
   iEvent.getByToken( genEvtInfoToken_, genEvtInfo );

   generatorWeight = genEvtInfo->weight();


   // get all the physics objects from the miniAOD
   Handle<pat::JetCollection> jets;
   iEvent.getByToken( jetToken_, jets );
   
   Handle<pat::JetCollection> jetsnoCHS;
   iEvent.getByToken( jetTokennoCHS_, jetsnoCHS );
   
   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken( vtxToken_, vertices );
   
   Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
   iEvent.getByToken( secVtxToken_, secVertices );

   edm::EventAuxiliary aux = iEvent.eventAuxiliary();
   edm::EventID id = aux.id();
  
   EventNumber = id.event();
   RunNumber = id.run();
   LuminosityBlock = id.luminosityBlock();
   
   for( const reco::Vertex &v : *vertices ) {
      bool isFake = v.isFake();
      if( !isFake && v.ndof() >= 4. && v.position().Rho() <= 2.0 && fabs(v.position().z()) <= 24. ) {
        vertex_x -> push_back( v.x() );
        vertex_y -> push_back( v.y() );
        vertex_z -> push_back( v.z() );
        vertex_dx -> push_back( v.xError() );
        vertex_dy -> push_back( v.yError() );
        vertex_dz -> push_back( v.zError() );
        vertex_ndof -> push_back( v.ndof() );
        vertex_d0 -> push_back( v.position().rho() );
        vertex_nTracks -> push_back( v.nTracks() );
        vertex_pt -> push_back( v.p4().pt() );
      }
   }
   
   // jets
   // using the tight selection from:
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
   int ctrJet = -1;
   for( const pat::Jet &j : *jets ) {
     
     if( j.neutralHadronEnergyFraction() >= 0.90 ) continue;
     if( j.neutralEmEnergyFraction() >= 0.90 ) continue;
     if( j.numberOfDaughters() <= 1 ) continue;
     if( j.muonEnergyFraction() >= 0.8 ) continue;
     
     if( fabs(j.eta()) < 2.4 ) {
        if( j.chargedEmEnergyFraction() >= 0.9 ) continue;
        if( j.chargedHadronEnergyFraction() <= 0. ) continue;
        if( j.chargedMultiplicity() <= 0. ) continue;
     }
     if( fabs(j.eta()) > 3.0 ) {
        if( j.neutralMultiplicity() <= 10 ) continue;
     }
     if( j.pt() < 10. ) continue;
     
     ctrJet += 1;
     
     std::vector<double> constVert_x;
     std::vector<double> constVert_y;
     std::vector<double> constVert_z;
     std::vector<double> const_pt;
     std::vector<double> const_eta;
     std::vector<double> const_phi;
     std::vector<int> const_charge;
     std::vector<double> const_pca0_x;
     std::vector<double> const_pca0_y;
     std::vector<double> const_pca0_z;
     std::vector<double> constVert_closestVertex_dxy;
     std::vector<double> constVert_closestVertex_dz;
     std::vector<double> constVert_closestVertex_d;
     std::vector<double> average_distance( vertices->size(), 0. );
     
     // loop over the jet constituents to find the vertex closest to each:
     for( unsigned int iD = 0; iD < j.numberOfDaughters(); ++iD ) {
        
        const pat::PackedCandidate &dau1 = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(iD));
        pat::PackedCandidate dau(dau1);
      
        const_pca0_x.push_back( dau.vertex().x() );
        const_pca0_y.push_back( dau.vertex().y() );
        const_pca0_z.push_back( dau.vertex().z() );

        // minimum distances (total, xy and z) for the jet constituent to any vertex
        double dMin = 100000.;
        double dxyMin = 10000.;
        double dzMin = 10000.;
        int ctr = -1;
       
        // get the unity vector pointing in the direction of the momentum and a reference point to build a self-made pseudo track
        // the 'track' is then (x(t), y(t), z(t) ) = (x,y,z) + t*(px,py,pz);
        // tmin, the time parameter for the point of closest approach is determined by minimising d = sqrt( (vx - x(t))^2 + (vy - y(t))^2 + (vz - z(t))^2);

        double jetVertex_x = -10000.;
        double jetVertex_y = -10000.;
        double jetVertex_z = -10000.;
        // first loop over the primary vertices
        for( const reco::Vertex &v : *vertices ) {
            ctr += 1;
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.position().x();
            double pv_y = v.position().y();
            double pv_z = v.position().z();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);

            // if the vertex is closer than the current reference vertex, set dMin, dxyMin, dzMin, and also change the vertex of reference
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              jetVertex_x = v.position().x();
              jetVertex_y = v.position().y();
              jetVertex_z = v.position().z();
             
            }
        }
  
        // now do the same for the secondary vertices
        // however, for some reason, I have to use additional variable here, jet constitutent won't accept a VertexCompositePtrCandidate as a new reference
        for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.vx();
            double pv_y = v.vy();
            double pv_z = v.vz();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);
            
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              ctr = -1;  
              jetVertex_x = v.vx();
              jetVertex_y = v.vy();
              jetVertex_z = v.vz();
          }
        }
        
        // now fill all the variables for the jet constituent
        constVert_closestVertex_dxy.push_back( dxyMin );
        constVert_closestVertex_dz.push_back( dzMin );
        constVert_closestVertex_d.push_back( dMin );
        constVert_x.push_back( jetVertex_x );
        constVert_y.push_back( jetVertex_y );
        constVert_z.push_back( jetVertex_z );
        const_pt.push_back( dau.pt() );
        const_eta.push_back( dau.eta() );
        const_phi.push_back( dau.phi() );
        const_charge.push_back( dau.charge() );
     }

     jet_pt->push_back( j.pt() );
     jet_eta->push_back( j.eta() );
     jet_phi->push_back( j.phi() );
     
     jet_constVertex_x->push_back( constVert_x ); 
     jet_constVertex_y->push_back( constVert_y ); 
     jet_constVertex_z->push_back( constVert_z ); 
     jet_const_closestVertex_dxy->push_back(constVert_closestVertex_dxy);
     jet_const_closestVertex_dz->push_back(constVert_closestVertex_dz);
     jet_const_closestVertex_d->push_back(constVert_closestVertex_d);
     jet_const_pt->push_back( const_pt );
     jet_const_eta->push_back( const_eta );
     jet_const_phi->push_back( const_phi );
     jet_const_charge->push_back( const_charge );
     jet_const_pca0_x->push_back( const_pca0_x );
     jet_const_pca0_y->push_back( const_pca0_y );
     jet_const_pca0_z->push_back( const_pca0_z );
      
   }
   for( const pat::Jet &j : *jetsnoCHS ) {
     
     if( j.neutralHadronEnergyFraction() >= 0.90 ) continue;
     if( j.neutralEmEnergyFraction() >= 0.90 ) continue;
     if( j.numberOfDaughters() <= 1 ) continue;
     if( j.muonEnergyFraction() >= 0.8 ) continue;
     
     if( fabs(j.eta()) < 2.4 ) {
        if( j.chargedEmEnergyFraction() >= 0.9 ) continue;
        if( j.chargedHadronEnergyFraction() <= 0. ) continue;
        if( j.chargedMultiplicity() <= 0. ) continue;
     }
     if( fabs(j.eta()) > 3.0 ) {
        if( j.neutralMultiplicity() <= 10 ) continue;
     }
     if( j.pt() < 10. ) continue;
     
     ctrJet += 1;
     
     std::vector<double> constVert_x;
     std::vector<double> constVert_y;
     std::vector<double> constVert_z;
     std::vector<double> const_pt;
     std::vector<double> const_eta;
     std::vector<double> const_phi;
     std::vector<int> const_charge;
     std::vector<double> const_pca0_x;
     std::vector<double> const_pca0_y;
     std::vector<double> const_pca0_z;
     std::vector<double> constVert_closestVertex_dxy;
     std::vector<double> constVert_closestVertex_dz;
     std::vector<double> constVert_closestVertex_d;
     std::vector<double> average_distance( vertices->size(), 0. );
     
     // loop over the jet constituents to find the vertex closest to each:
     for( unsigned int iD = 0; iD < j.numberOfDaughters(); ++iD ) {
        
        const pat::PackedCandidate &dau1 = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(iD));
        pat::PackedCandidate dau(dau1);
      
        const_pca0_x.push_back( dau.vertex().x() );
        const_pca0_y.push_back( dau.vertex().y() );
        const_pca0_z.push_back( dau.vertex().z() );

        // minimum distances (total, xy and z) for the jet constituent to any vertex
        double dMin = 100000.;
        double dxyMin = 10000.;
        double dzMin = 10000.;
        int ctr = -1;
       
        // get the unity vector pointing in the direction of the momentum and a reference point to build a self-made pseudo track
        // the 'track' is then (x(t), y(t), z(t) ) = (x,y,z) + t*(px,py,pz);
        // tmin, the time parameter for the point of closest approach is determined by minimising d = sqrt( (vx - x(t))^2 + (vy - y(t))^2 + (vz - z(t))^2);

        double jetVertex_x = -10000.;
        double jetVertex_y = -10000.;
        double jetVertex_z = -10000.;
        // first loop over the primary vertices
        for( const reco::Vertex &v : *vertices ) {
            ctr += 1;
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.position().x();
            double pv_y = v.position().y();
            double pv_z = v.position().z();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);

            // if the vertex is closer than the current reference vertex, set dMin, dxyMin, dzMin, and also change the vertex of reference
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              jetVertex_x = v.position().x();
              jetVertex_y = v.position().y();
              jetVertex_z = v.position().z();
             
            }
        }
  
        // now do the same for the secondary vertices
        // however, for some reason, I have to use additional variable here, jet constitutent won't accept a VertexCompositePtrCandidate as a new reference
        for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.vx();
            double pv_y = v.vy();
            double pv_z = v.vz();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);
            
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              ctr = -1;  
              jetVertex_x = v.vx();
              jetVertex_y = v.vy();
              jetVertex_z = v.vz();
          }
        }
        
        // now fill all the variables for the jet constituent
        constVert_closestVertex_dxy.push_back( dxyMin );
        constVert_closestVertex_dz.push_back( dzMin );
        constVert_closestVertex_d.push_back( dMin );
        constVert_x.push_back( jetVertex_x );
        constVert_y.push_back( jetVertex_y );
        constVert_z.push_back( jetVertex_z );
        const_pt.push_back( dau.pt() );
        const_eta.push_back( dau.eta() );
        const_phi.push_back( dau.phi() );
        const_charge.push_back( dau.charge() );
     }

     jetnoCHS_pt->push_back( j.pt() );
     jetnoCHS_eta->push_back( j.eta() );
     jetnoCHS_phi->push_back( j.phi() );
     
     jetnoCHS_constVertex_x->push_back( constVert_x ); 
     jetnoCHS_constVertex_y->push_back( constVert_y ); 
     jetnoCHS_constVertex_z->push_back( constVert_z ); 
     jetnoCHS_const_closestVertex_dxy->push_back(constVert_closestVertex_dxy);
     jetnoCHS_const_closestVertex_dz->push_back(constVert_closestVertex_dz);
     jetnoCHS_const_closestVertex_d->push_back(constVert_closestVertex_d);
     jetnoCHS_const_pt->push_back( const_pt );
     jetnoCHS_const_eta->push_back( const_eta );
     jetnoCHS_const_phi->push_back( const_phi );
     jetnoCHS_const_charge->push_back( const_charge );
     jetnoCHS_const_pca0_x->push_back( const_pca0_x );
     jetnoCHS_const_pca0_y->push_back( const_pca0_y );
     jetnoCHS_const_pca0_z->push_back( const_pca0_z );
      
  }
 

  
   // now fill the secondary vertex information
   for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
      secVertex_x -> push_back( v.vx() );
      secVertex_y -> push_back( v.vy() );
      secVertex_z -> push_back( v.vz() );
      secVertex_ndof -> push_back( v.vertexNdof() );
      secVertex_chi2 -> push_back( v.vertexChi2() );
      secVertex_pt -> push_back( v.pt() );
      secVertex_dx -> push_back( v.vertexCovariance(0,0) );
      secVertex_dy -> push_back( v.vertexCovariance(1,1) );
      secVertex_dz -> push_back( v.vertexCovariance(2,2) );
   }

   // finally, write it to the tree
   tOutput->Fill(); 

}

std::vector<double> CalculateVertex( std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> weight, std::vector<int> charge, std::vector<double> distance, int &nConsidered, double &weightednConsidered, std::vector<double> &error, double &maxScore ) {

   nConsidered = 0;
   std::vector<double> diff_x;
   std::vector<double> diff_y;
   std::vector<double> diff_z;
   std::vector<double> score;


   for( unsigned int i = 0; i < x.size(); ++i ) {
      if( charge.at(i) == 0 ) continue;
      nConsidered += 1;
      bool knownPoint = false;
      int iKnown = -1;
      for( unsigned int i2 = 0; i2 < diff_x.size(); ++i2 ) {
        if( fabs( diff_x.at(i2) - x.at(i) ) < 1.e-10 && fabs( diff_y.at(i2) - y.at(i) ) < 1.e-10 && fabs( diff_z.at(i2) - z.at(i) ) < 1.e-10 ) {
            knownPoint = true;
            iKnown = i2;
        }
      }

      if( knownPoint ) {
        if( distance.at(i) == 0. ) score.at(iKnown) += 1.e12;
        else                       score.at(iKnown) += weight.at(i)/distance.at(i);
      }
      else {
        diff_x.push_back( x.at(i) );
        diff_y.push_back( y.at(i) );
        diff_z.push_back( z.at(i) );
        if( distance.at(i) == 0. ) score.push_back( 1.e12 );
        else                       score.push_back( weight.at(i)/distance.at(i) );
      }
   }

   double scoreMax = 0.;
   std::vector<double> mean(3, -10000.);
   for( unsigned int i = 0; i < diff_x.size(); ++i ) {
        if ( score.at(i) > scoreMax ) {
            scoreMax = score.at(i);
            maxScore = scoreMax;
            mean.at(0) = diff_x.at(i);
            mean.at(1) = diff_y.at(i);
            mean.at(2) = diff_z.at(i);
        }
    }
    return mean;
}


// ------------ method called once each job just before starting event loop  ------------
void 
LLGDVJetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LLGDVJetAnalyzer::endJob() 
{

  // save the output tree and file
  gDirectory = fOutput;
  tOutput->Write();
  fOutput->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
LLGDVJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
LLGDVJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
LLGDVJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
LLGDVJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LLGDVJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LLGDVJetAnalyzer);
