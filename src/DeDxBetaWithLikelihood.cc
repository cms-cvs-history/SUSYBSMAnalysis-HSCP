// -*- C++ -*-
//
// Package:    DeDxBetaWithLikelihood
// Class:      DeDxBetaWithLikelihood
// 
/**\class DeDxBetaWithLikelihood DeDxBetaWithLikelihood.cc SUSYBSMAnalysis/DeDxBetaWithLikelihood/src/DeDxBetaWithLikelihood.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rizzi Andrea
//         Created:  Wed Oct 10 12:01:28 CEST 2007
// $Id: DeDxBetaWithLikelihood.cc,v 1.3 2007/10/15 09:32:52 arizzi Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "RecoTracker/DeDx/interface/DeDxEstimatorProducer.h"
#include "DataFormats/TrackReco/interface/TrackDeDxEstimate.h"
#include "DataFormats/TrackReco/interface/TrackDeDxHits.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "DataFormats/TrackReco/interface/Track.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


#include <vector>
#include <TNtuple.h>
#include <TF1.h>
#include <iostream>
//
// class decleration
//

class DeDxBetaWithLikelihood : public edm::EDProducer {
   public:
      explicit DeDxBetaWithLikelihood(const edm::ParameterSet&);
      ~DeDxBetaWithLikelihood();

   private:
      float fit(const reco::TrackDeDxHits &);
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      edm::InputTag m_trackDeDxHitsTag;
      TNtuple * tmpNt;
      TF1 * f1;
      // ----------member data ---------------------------
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
DeDxBetaWithLikelihood::DeDxBetaWithLikelihood(const edm::ParameterSet& iConfig)
{
  edm::Service<TFileService> fs;
  TFileDirectory subDir = fs->mkdir( "Temp" );
  tmpNt =  subDir.make<TNtuple>( "dedx","dedx","dedx");
  f1 =  subDir.make<TF1>( "f1", "landaun" );

   m_trackDeDxHitsTag = iConfig.getParameter<edm::InputTag>("trackDeDxHits");
   produces<std::vector<float> >();

}


DeDxBetaWithLikelihood::~DeDxBetaWithLikelihood()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DeDxBetaWithLikelihood::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;

   edm::Handle<reco::TrackDeDxHitsCollection> trackDeDxHitsCollectionHandle;
   iEvent.getByLabel(m_trackDeDxHitsTag,trackDeDxHitsCollectionHandle);
   const reco::TrackDeDxHitsCollection & hits = *trackDeDxHitsCollectionHandle.product();
  vector<float> * outputCollection = new vector<float>;

   reco::TrackDeDxHitsCollection::const_iterator it= hits.begin();
   for(int j=0;it!=hits.end();++it,j++)
   {
      float val=fit(*it);
      outputCollection->push_back(val);
   }

    std::auto_ptr<vector<float> > estimator(outputCollection);
    iEvent.put(estimator);

}

float DeDxBetaWithLikelihood::fit(const reco::TrackDeDxHits & dedxVec)
{
 using namespace std;
 
 double mpv=-5.;
 double chi=-5.;
 int entries = 0; 
 if(dedxVec.second.size() < 1) return 0;
 tmpNt->Reset();
 // copy data into a tree:
 for (unsigned int i=0; i<dedxVec.second.size();i++) {
// cdout << dedxVec.second[i].charge() << endl;
   if(dedxVec.second[i].subDet() != 1 && dedxVec.second[i].subDet() != 2 )
   {
     tmpNt->Fill(dedxVec.second[i].charge());
     entries++; 
   }
 }
 if( entries < 1) return 0;

 // fit:
 f1->SetParameters(1, 3.0 , 0.3);
 f1->SetParLimits(0, 1, 1); // fix the normalization parameter to 1
 int status = tmpNt->UnbinnedFit("f1","dedx","","Q");
 mpv = f1->GetParameter(1);
 if (status<=0) {
   cout << "(AnalyzeTracks::LandauFit) no convergence!   status = " << status << endl;
 tmpNt->Scan();
    return 0;
 }
 return mpv;
// return 0.123123; 
}

// ------------ method called once each job just before starting event loop  ------------
void 
DeDxBetaWithLikelihood::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DeDxBetaWithLikelihood::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeDxBetaWithLikelihood);
