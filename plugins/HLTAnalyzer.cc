// -*- C++ -*-
//
// Package:    HLTtest/HLTAnalyzer
// Class:      HLTAnalyzer
// 
/**\class HLTAnalyzer HLTAnalyzer.cc HLTtest/HLTAnalyzer/plugins/HLTAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Mon, 14 Apr 2014 09:14:11 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//
// class declaration
//


class HLTAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HLTAnalyzer(const edm::ParameterSet&);
      ~HLTAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
HLTAnalyzer::HLTAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


HLTAnalyzer::~HLTAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HLTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

  edm::Handle<edm::TriggerResults> triggerResultsH_;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjectStandAlonesH_;
  edm::TriggerNames triggerNames_;
  // Get TriggerResults and triggerNames
  iEvent.getByLabel(edm::InputTag("TriggerResults",       "", "HLT"), triggerResultsH_);
  if (! triggerResultsH_.isValid())
    throw cms::Exception("HLTTest::produce: error getting TriggerResults product from Event!");
  triggerNames_ = iEvent.triggerNames(*triggerResultsH_);

  // Get TriggerObjectStandAlones
  iEvent.getByLabel("selectedPatTrigger", triggerObjectStandAlonesH_);
  if (! triggerObjectStandAlonesH_.isValid())
    throw cms::Exception("HLTTest::produce: error getting TriggerObjectsStandAlone product from Event!");

  // Get PackedTriggerPrescales
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesH_; 
  iEvent.getByLabel( "patTrigger", triggerPrescalesH_);
  if (! triggerPrescalesH_.isValid())
    throw cms::Exception("HLTTest::produce: error getting PackedTriggerPrescales product from Event!");
 /////////////////////////////////////////////////////
  // Print all TriggerObjectStandAlones from miniAOD //
  // including the paths and filters they passed.    //
  //
  if ( triggerObjectStandAlonesH_.isValid()) 
    cout<<"Got triggerObjectStandAlonesHandle with size "<<triggerObjectStandAlonesH_->size()<<endl;
  else cout<<"Couldn't find triggerObjectStandAlonesH"<<endl;
  for ( uint i = 0; i < triggerObjectStandAlonesH_->size(); i++ ) {
    pat::TriggerObjectStandAlone TO = triggerObjectStandAlonesH_->at(i);
    TO.unpackPathNames( triggerNames_ );
    cout<<"Trigger Object "<<i<<" has pt, eta, phi = "<< TO.pt() <<" "<<TO.eta()<<" "<<TO.phi()<<endl;
  
    std::vector< std::string > path_namesPASS = TO.pathNames(true); 
    cout<<"Trigger Object "<<i<<" passed "<< path_namesPASS.size() <<" paths: ";
    for (uint j  = 0; j < path_namesPASS.size(); j++) cout<<path_namesPASS[j]<<" ";
    cout<<endl;
  
   std::vector< std::string > filter_labels = TO.filterLabels();
   cout<<"Trigger Object "<<i<<" passed "<< filter_labels.size() <<" filters: ";
   for (uint j  = 0; j < filter_labels.size(); j++) cout<<filter_labels[j]<<" ";
   cout<<endl;
  }
  //
  // END printing miniAOD content                    //
  /////////////////////////////////////////////////////

  // Loop over triggers
  unsigned int nTriggers = triggerResultsH_->size();
  for(unsigned int i = 0; i < nTriggers; ++i)
    {

      // Name
      const string& name = triggerNames_.triggerName(i);
	
      // Prescale
      int prescale = triggerPrescalesH_.isValid() ? triggerPrescalesH_->getPrescaleForIndex(i) : -1;
	
      // Passed...
      if (triggerResultsH_->accept(i)) {
	
	cout<<"Passed trigger path "<<name<<" with prescale "<< prescale <<endl;
	cout<<"Associated trigger objects are: "<<endl;
	for ( uint i = 0; i < triggerObjectStandAlonesH_->size(); i++ ) {
	  pat::TriggerObjectStandAlone TO = triggerObjectStandAlonesH_->at(i);
	  TO.unpackPathNames( triggerNames_ ); 
	  // TO.hasPathName(name, true): true if TO passed the LAST filter in the HLT path.
	  // TO.hasPathName(name, false): true if TO passed ANY filter in the path.
	  if ( TO.hasPathName(name, true) ) { 
	    cout<<"pt, eta, phi: "<<TO.pt()<<", "<<TO.eta()<<", "<<TO.phi()<<endl;
	  // To check whether a specific HLT filter is applied, use: 
	  // TO.hasFilterLabel( const std::string & filterLabel )
	  }
	} 
      }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void 
HLTAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HLTAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HLTAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HLTAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HLTAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HLTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTAnalyzer);
