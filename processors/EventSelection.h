#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

// Includes from ROOT:
#include "TString.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TParameter.h"

// Includes for Criteria Selection (Taken from local classes directory)
#include "../classes/mTowerHit.h"
#include "../classes/mTowerClusterRobbie.h"
#include "../classes/mTowerEvent.h"
#include "../classes/mTowerChipRobbie.h"

// Includes for kT-Algorithm Selection
#include "fastjet/ClusterSequence.hh"                                                   // Requires installation of Fastjet.
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/util.hh"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/mTowerUtil.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_hot_pixels.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_module_thickness.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_inclination_parameters.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_alignment_parameters.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Nov2020/11_id_pileup/cellinfo.h"

bool IsLeftChip(int lane);

bool EventSelectionA(bool CheckRejects, bool CheckAccepts, bool CheckThirdLayer, bool C2, bool C4, bool C6, int nPixelRadiusC2, int nPixelRadiusC4, int nPixelBorderC6, const int laneNumber[], const int laneOffset, const int columnsPerChip, const int rowsPerChip, double nPixelsGap, mTowerChipRobbie hitsInChip[], int eventID, int nHits, int eventIndex);

bool EventSelectionB(mTowerLayer* clayers[], mTowerChip* cchips[], Double_t nHitsTotNew[], Double_t nClustersTotNew[], Double_t dlim, Double_t dcenter, Int_t minNlayerBulk, Double_t W0Bulk, Double_t thDistBulk, Double_t wBinBulk, TH2D* htmpBulk, TH2D* htmpLim, std::map<Int_t, cellinfo> mcellsBulk, std::map<Int_t, cellinfo> mcellsLim, Int_t nla, Int_t nLayerBulk, Int_t max_layer[]);

#endif
