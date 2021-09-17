#ifndef Validation_MuonCSCDigis_CSCStubResolutionValidation_H
#define Validation_MuonCSCDigis_CSCStubResolutionValidation_H

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "Validation/MuonCSCDigis/interface/CSCBaseValidation.h"
#include "Validation/MuonCSCDigis/interface/CSCStubMatcher.h"


#include "Validation/MuonCSCDigis/interface/CSCBaseValidation.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCLUTReader.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/ComparatorCodeLUT.h"

#include <map>
#include <string>
#include <tuple>

class CSCStubResolutionValidation : public CSCBaseValidation {
public:
  CSCStubResolutionValidation(const edm::ParameterSet &pset, edm::ConsumesCollector &&iC);
  ~CSCStubResolutionValidation() override;

  void bookHistograms(DQMStore::IBooker &);
  void analyze(const edm::Event &, const edm::EventSetup &) override;

  // access to the matcher
  std::shared_ptr<CSCStubMatcher> cscStubMatcher() { return cscStubMatcher_; }
  void setCSCStubMatcher(std::shared_ptr<CSCStubMatcher> s) { cscStubMatcher_ = s; }

  
private:
  bool isSimTrackGood(const SimTrack &t);

  edm::InputTag inputTag_;

  std::shared_ptr<CSCStubMatcher> cscStubMatcher_;
  

  // resolution for each CSC TP; 10 CSC stations;
  MonitorElement *posresCLCT_hs[10];
  MonitorElement *posresCLCT_qs[10];
  MonitorElement *posresCLCT_es[10];

  MonitorElement *bendresCLCT[10];

  edm::EDGetTokenT<edm::SimVertexContainer> simVertexInput_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackInput_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clcts_Token_;

  double simTrackMinPt_;
  double simTrackMinEta_;
  double simTrackMaxEta_;

};

#endif
