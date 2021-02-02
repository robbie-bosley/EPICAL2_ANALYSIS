#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

bool IsLeftChip(int lane);

bool EventSelectionA(bool DB, bool CheckRejects, bool CheckAccepts, bool CheckThirdLayer, bool C2, bool C4, bool C6, int nPixelRadiusC2, int nPixelRadiusC4, int nPixelBorderC6, const int laneNumber[], const int laneOffset, const int columnsPerChip, const int rowsPerChip, double nPixelsGap, mTowerChipRobbie* hitsInChip[], int eventID);

bool EventSelectionB(mTowerLayer* clayers[], mTowerChip* cchips[], Double_t nHitsTotNew[], Double_t nClustersTotNew[], Double_t dlim, Double_t dcenter, Int_t minNlayerBulk, Double_t W0Bulk, Double_t thDistBulk, Double_t wBinBulk, TH2D* htmpBulk, TH2D* htmpLim, std::map<Int_t, cellinfo> mcellsBulk, std::map<Int_t, cellinfo> mcellsLim);

#endif
