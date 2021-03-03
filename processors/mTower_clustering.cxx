#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <sys/stat.h>
#include <stdio.h>

#include "../classes/mTowerHit.h"
#include "../classes/mTowerClusterRobbie.h"
#include "../classes/mTowerEvent.h"
#include "../classes/mTowerChipRobbie.h"

using namespace std;

const std::map< Int_t, Int_t > lane2chipid_qasim_lut = { 
{ 32,0 },{ 59,2 },{ 34,4 },{ 57,6 }, { 36,8 },{ 61,10 },{ 38, 12 },{ 79, 14 },{ 54, 16 },{ 77, 18 },{ 52, 20 },{ 75,22 },
{ 35,1 }, { 56,3 },{ 33,5 },{ 58,7 },{ 37,9 },{ 60,11 },{ 55, 13 },{ 62, 15 },{ 53, 17 },{ 78, 19 },{ 51, 21 },{ 76,23 },

{ 50,24 },{ 73,26 },{ 48,28 },{ 71,30 }, { 46,32 },{ 69,34 },{ 44, 36 },{ 67, 38 },{ 42, 40 },{ 65, 42 },{ 40, 44 },{ 63, 46 },
{ 49,25 },{ 74,27 },{ 47,29 },{ 72,31 },{ 45,33 },{ 70,35 },{ 43, 37 },{ 68, 39 },{ 41, 41 },{ 66, 43 },{ 39, 45 },{ 64, 47 } };

long find_neigh(int chip, int n_entries, const std::vector< std::vector<pair<int,int> > > &address)
{					
  for(long i=(address[chip].size()-1); i>((address[chip].size()-1) - n_entries); i--) 
    {		
      if(((abs((address[chip][address[chip].size()-1].first)-(address[chip][i-1].first)))<clust_neigh_pix_lim) && ((abs((address[chip][address[chip].size()-1].second)-(address[chip][i-1].second)))<clust_neigh_pix_lim))
	{
	  return i-1;
	}
    }			
  return -1;
}


bool lane2isleft( Int_t lane) {
  Int_t chip = lane2chipid_qasim_lut.at(lane);
  if (chip % 2 == 0) {
    return false;
  } else {
    return true;
  }
}

double row2x( Int_t row, Int_t lane) {
  if (lane2isleft(lane)) {
    double x = -0.02688*row;
  }
  else {
    double x = 0.02688*row;
  }
  return x;
}

double column2y( Int_t col, Int_t lane) {
  double offset = 512*0.02924;
 if (lane2isleft(lane)) {
    double y = 0.02924*col - offset;
  }
  else {
    double y = offset - 0.02924*col;
  }
 return y;
} 


void QasimClustering(TObjArray* hitList, Int_t cur_evt_ch_[], int event, vector<Int_t>& clusterLane, vector<Double_t>& clusterX, vector<Double_t>& clusterY, vector<Int_t>& clusterSize, vector<vector<Int_t> >& clusterHitIndex) {
  

  for (i = 0; i < hitList->GetEntries(); i++) {
    mTowerHit* currentHit = (mTowerHit*)hitList->At(hit);
    Int_t lane = currentHit -> getLane();
    Int_t row = currentHit -> getRow();
    Int_t col = currentHit -> getColumn();
    Int_t chip = lane2chipid_qasim_lut.at( lane );

    

    if( chip < 0 ) continue;
    int n_entries=-1;
    int clust_index=-1;
    if(cur_evt_ch_[chip]!=event) // First instance of hit in this chip, in this event.
      {
	chip_w_hits.push_back(chip);
	ind_nHits[chip].push_back(1);
	clust_size[chip].push_back(1);
	part_per_evt[chip].push_back(1);
	address[chip].push_back( make_pair(col,row) );
	clust_no[chip].push_back(clust_size[chip].size()-1);
	chip_lane[chip] = lane;
	cur_evt_ch_[chip]=event;
      }
    else if (cur_evt_ch_[chip]==event)
      {
	address[chip].push_back( make_pair(col,row));
	n_entries=ind_nHits[chip][ind_nHits[chip].size()-1];
	ind_nHits[chip][ind_nHits[chip].size()-1]+=1;
	clust_index=find_neigh(chip,n_entries,address);

	if(clust_index!=-1) // Some already existing cluster was found.
	  {
	    clust_size[chip][clust_no[chip][clust_index]]+=1;
	    clust_no[chip].push_back(clust_no[chip][clust_index]);
	    clusterHitIndex[chip][clust_no[chip][clust_index]].push_back(i);
	  }
	else // No cluster was found, create a new cluster of size 1.
	  {
	    clust_size[chip].push_back(1);
	    clust_no[chip].push_back(clust_size[chip].size()-1);
	    clusterHitIndex[chip][clust_size[chip].size()-1].push_back(i);
	    chip_lane[chip] = lane;
	    part_per_evt[chip][part_per_evt[chip].size()-1]+=1;
	  }
      }
  }//iHits for loop ends here

  

  for (int chip=0; chip < 48; chip++) {
    for (int clust=0; clust < clust_size[chip].size(); clust++) {
      clusterLane.push_back(chip_lane[chip]);
      clusterSize.push_back(clust_size[chip][clust]);
      vector<Int_t> clushitindex;
      vector<double> clushitx;
      vector<double> clushity;
      if (clust_no[chip] == clust) {
	clushitx.push_back(row2x(address[chip][cluster].first));
	clushity.push_back(column2y(address[chip][cluster].second));
      }
      
      int sumx, sumy;
      for (int i = 0; i < clushitx.size(); i++) {
	sumx += clushitx[i];
	sumy += clushity[i];
      }
      
      double meanx = sumx/clushitx.size();
      double meany = sumy/clushitx.size();
      clusterX.push_back(meanx);
      clusterY.push_back(meany);

    }
  }
  

  chip_w_hits.clear();


}
