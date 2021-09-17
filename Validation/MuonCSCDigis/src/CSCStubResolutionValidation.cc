#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Validation/MuonCSCDigis/interface/CSCStubResolutionValidation.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"

#include "Validation/MuonHits/interface/CSCSimHitMatcher.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"


CSCStubResolutionValidation::CSCStubResolutionValidation(const edm::ParameterSet& pset, edm::ConsumesCollector&& iC)
    : CSCBaseValidation(pset) {
  const auto& simVertex = pset.getParameter<edm::ParameterSet>("simVertex");
  simVertexInput_ = iC.consumes<edm::SimVertexContainer>(simVertex.getParameter<edm::InputTag>("inputTag"));
  const auto& simTrack = pset.getParameter<edm::ParameterSet>("simTrack");
  simTrackInput_ = iC.consumes<edm::SimTrackContainer>(simTrack.getParameter<edm::InputTag>("inputTag"));
  simTrackMinPt_ = simTrack.getParameter<double>("minPt");
  simTrackMinEta_ = simTrack.getParameter<double>("minEta");
  simTrackMaxEta_ = simTrack.getParameter<double>("maxEta");

  // all CSC TPs have the same label
  const auto& stubConfig = pset.getParameterSet("cscALCT");
  inputTag_ = stubConfig.getParameter<edm::InputTag>("inputTag");
  clcts_Token_ = iC.consumes<CSCCLCTDigiCollection>(inputTag_);
  
  // Initialize stub matcher
  cscStubMatcher_.reset(new CSCStubMatcher(pset, std::move(iC)));

}

CSCStubResolutionValidation::~CSCStubResolutionValidation() {}

//create folder for resolution histograms and book them
void CSCStubResolutionValidation::bookHistograms(DQMStore::IBooker& iBooker) {
  iBooker.setCurrentFolder("MuonCSCDigisV/CSCDigiTask/Stub/Resolution/");
  
  for (int i = 1; i <= 10; ++i) {
    int j = i - 1;
    const std::string cn(CSCDetId::chamberName(i));

    //Position resolution; CLCT
    std::string t1 = "CLCTPosRes_hs_" + cn;
    std::string t2 = "CLCTPosRes_qs_" + cn;
    std::string t3 = "CLCTPosRes_es_" + cn;

    posresCLCT_hs[j] = iBooker.book1D(t1, cn + " CLCT Position Resolution (1/2-strip prec.); Strip_{L1T} - Strip_{SIM}; Entries", 50, -1, 1);
    posresCLCT_qs[j] = iBooker.book1D(t2, cn + " CLCT Position Resolution (1/4-strip prec.); Strip_{L1T} - Strip_{SIM}; Entries", 50, -1, 1);
    posresCLCT_es[j] = iBooker.book1D(t3, cn + " CLCT Position Resolution (1/8-strip prec.); Strip_{L1T} - Strip_{SIM}; Entries", 50, -1, 1);

    //Slope resolution; CLCT
    std::string t4 = "CLCTBendRes_" + cn;
   
    bendresCLCT[j] = iBooker.book1D(t4, cn + " CLCT Bend Resolution; Slope_{L1T} - Slope_{SIM}; Entries", 50, -0.5, 0.5);  
  }
}


void CSCStubResolutionValidation::analyze(const edm::Event& e, const edm::EventSetup& eventSetup) {
  // Define handles
  edm::Handle<edm::SimTrackContainer> sim_tracks;
  edm::Handle<edm::SimVertexContainer> sim_vertices;
  edm::Handle<CSCCLCTDigiCollection> clcts;
  
  // Use token to retreive event information
  e.getByToken(simTrackInput_, sim_tracks);
  e.getByToken(simVertexInput_, sim_vertices);
  e.getByToken(clcts_Token_, clcts);
  
  // Initialize StubMatcher
  cscStubMatcher_->init(e, eventSetup);

  
  const edm::SimTrackContainer& sim_track = *sim_tracks.product();
  const edm::SimVertexContainer& sim_vert = *sim_vertices.product();
  
  if (!clcts.isValid()) {
    edm::LogError("CSCStubResolutionValidation") << "Cannot get CLCTs by label " << inputTag_.encode();
  }
  
  // select simtracks for true muons
  edm::SimTrackContainer sim_track_selected;
  for (const auto& t : sim_track) {
    if (!isSimTrackGood(t))
      continue;
    sim_track_selected.push_back(t);
  }

  // Skip events with no selected simtracks
  if (sim_track_selected.empty())
    return;

  // Loop through good tracks, use corresponding vertex to match stubs, then fill hists of chambers where the stub appears.
  for (const auto& t : sim_track_selected) {
    std::vector<bool> hitCLCT(10);
    
    std::vector<float> delta_fhs_clct(10);
    std::vector<float> delta_fqs_clct(10);
    std::vector<float> delta_fes_clct(10);

    std::vector<float> dslope_clct(10);

    // Match track to stubs with appropriate vertex
    cscStubMatcher_->match(t, sim_vert[t.vertIndex()]);

    // Store matched stubs.
    // Key: ChamberID, Value : CSCStubDigiContainer
    const auto& clcts = cscStubMatcher_->clcts();
    

    //CLCTs
    for (auto& [id, container] : clcts) {
      const CSCDetId cscId(id);
      const unsigned chamberType(cscId.iChamberType());

      //get the best clct in chamber
      const auto& clct = cscStubMatcher_->bestClctInChamber(id);
      if (!clct.isValid()) continue;
      
      hitCLCT[chamberType - 1] = true;

      //calculate deltastrip
      int deltaStrip = 0;
      if (cscId.station() == 1 and cscId.ring() == 4 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B)
      	deltaStrip = CSCConstants::NUM_HALF_STRIPS_ME1B;

      //get the matched stub's keystrip 
      // const int hs_clct = clct.getKeyStrip(2);
      // const int qs_clct = clct.getKeyStrip(4);
      // const int es_clct = clct.getKeyStrip(8);

      // fractional strip
       const float fhs_clct = clct.getFractionalStrip(2);
       const float fqs_clct = clct.getFractionalStrip(4);
       const float fes_clct = clct.getFractionalStrip(8);

      // in half-strips per layer
      const float slopeHalfStrip(clct.getFractionalSlope());
      const float slopeStrip(slopeHalfStrip / 2.);
	
      //get the fit hits in chamber for true value
      float stripIntercept, stripSlope;
      cscStubMatcher_->cscDigiMatcher()->muonSimHitMatcher()->fitHitsInChamber(id, stripIntercept, stripSlope);

      //add offset of +0.25 strips for non-ME1/1 chambers
      const bool isME11(cscId.station()==1 and (cscId.ring()==4 or cscId.ring()==1));
      if (!isME11){
      	stripIntercept -= 0.25;
      }

      const float strip_csc_sh = stripIntercept;
      const float bend_csc_sh = stripSlope;
      
      delta_fhs_clct[chamberType - 1] = fhs_clct - deltaStrip - strip_csc_sh;
      delta_fqs_clct[chamberType - 1] = fqs_clct - deltaStrip - strip_csc_sh;
      delta_fes_clct[chamberType - 1] = fes_clct - deltaStrip - strip_csc_sh;

      dslope_clct[chamberType - 1] = slopeStrip - bend_csc_sh;

      //print statements
      if (chamberType -1 ==0){
	std::cout << "clct" << clct << std::endl;
	std::cout << "fhs_clct " << fhs_clct << std::endl;
	std::cout << "deltaStrip " << deltaStrip << std::endl;
	std::cout << "strip_csc_sh " << strip_csc_sh << std::endl;
      }

    }

    for (int i = 0; i < 10; i++) {
      if (hitCLCT[i]) {

	//fill histograms
	posresCLCT_hs[i]->Fill(delta_fhs_clct[i]);
	posresCLCT_qs[i]->Fill(delta_fqs_clct[i]);
	posresCLCT_es[i]->Fill(delta_fes_clct[i]);
	
	bendresCLCT[i]->Fill(dslope_clct[i]);
      }
    }
  }
}

bool CSCStubResolutionValidation::isSimTrackGood(const SimTrack& t) {
  // SimTrack selection
  if (t.noVertex())
    return false;
  if (t.noGenpart())
    return false;
  // only muons
  if (std::abs(t.type()) != 13)
    return false;
  // pt selection
  if (t.momentum().pt() < simTrackMinPt_)
    return false;
  // eta selection
  const float eta(std::abs(t.momentum().eta()));
  if (eta > simTrackMaxEta_ || eta < simTrackMinEta_)
    return false;
  return true;
}


