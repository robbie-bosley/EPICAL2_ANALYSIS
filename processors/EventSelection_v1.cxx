/////////////////////////////////////////////////////////////////////////////////////
// 1)  INCLUDES

#include "processors/EventSelection_v1.h"

// Including some standard libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <tuple>
#include <fstream>
#include <string>

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
#include "classes/mTowerHit.h"
#include "classes/mTowerClusterRobbie.h"
#include "classes/mTowerEvent.h"
#include "classes/mTowerChipRobbie.h"



// Includes for kT-Algorithm Selection
#include "fastjet/ClusterSequence.hh"                                                   // Requires installation of Fastjet.
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/util.hh"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/mTowerUtil.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_hot_pixels.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_module_thickness.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_inclination_parameters.h"
#include "/eos/project/m/mtower/public/hiroki/event_selection_Feb2021/02_convert_v1/00_util/load_alignment_parameters.h"
//#include "/eos/project/m/mtower/public/hiroki/11_id_pileup/cellinfo.h"

// Adding in namespaces
using namespace std;
using namespace fastjet;
using namespace ROOT::Math;

////////////////////////////////////////////////////////////////////////////
// 2) LOOKUP TABLES

const std::map< Int_t, Int_t > lane2chipidrobbie_lut_ES = {
  {40, 0},{39, 1},{42, 2},{41, 3},{44, 4},{43, 5},{46, 6},{45, 7},
  {48, 8},{47, 9},{50,10},{49,11},{52,12},{51,13},{54,14},{53,15},
  {38,16},{55,17},{36,18},{37,19},{32,20},{35,21},{34,22},{33,23},
  {64,24},{63,25},{66,26},{65,27},{68,28},{67,29},{70,30},{69,31},
  {72,32},{71,33},{74,34},{73,35},{76,36},{75,37},{78,38},{77,39},
  {62,40},{79,41},{60,42},{61,43},{56,44},{59,45},{58,46},{57,47}
};

const std::map< Int_t, Int_t > lane2layerrobbie_lut_ES = {
  {40,22},{39,22},{42,20},{41,20},{44,18},{43,18},{46,16},{45,16},
  {48,14},{47,14},{50,12},{49,12},{52,10},{51,10},{54, 8},{53, 8},
  {38, 6},{55, 6},{36, 4},{37, 4},{32, 0},{35, 0},{34, 2},{33, 2},
  {64,23},{63,23},{66,21},{65,21},{68,19},{67,19},{70,17},{69,17},
  {72,15},{71,15},{74,13},{73,13},{76,11},{75,11},{78, 9},{77, 9},
  {62, 7},{79, 7},{60, 5},{61, 5},{56, 1},{59, 1},{58, 3},{57, 3}
};

const std::map<int,bool> layer2isInvrobbie_lut_ES = {
  { 0, kFALSE}, { 1, kTRUE}, { 2, kFALSE}, { 3, kTRUE}, 
  { 4, kFALSE}, { 5, kTRUE}, { 6, kFALSE}, { 7, kTRUE}, 
  { 8, kFALSE}, { 9, kTRUE}, {10, kFALSE}, {11, kTRUE}, 
  {12, kFALSE}, {13, kTRUE}, {14, kFALSE}, {15, kTRUE}, 
  {16, kFALSE}, {17, kTRUE}, {18, kFALSE}, {19, kTRUE}, 
  {20, kFALSE}, {21, kTRUE}, {22, kFALSE}, {23, kTRUE} 
};

///////////////////////////////////////////////////////////
// 3) USEFUL FUNCTIONS

// Chip-based Functions
bool IsLeftChip(int lane){
  int layerNr = lane2layerrobbie_lut_ES.at(lane);
  bool isInv  = layer2isInvrobbie_lut_ES.at(layerNr);
  int chipid = lane2chipidrobbie_lut_ES.at(lane); 
  bool isOdd;
  if (chipid%2 == 1){ isOdd = kTRUE;}
  else {isOdd = kFALSE;}
  bool isLeft = (bool)(isOdd != isInv);
  return isLeft;
}

///////////////////////////////////////////////////////////
// 4) CRITERIA-BASED EVENT SELECTION
bool EventSelectionA_v1(bool CheckRejects, bool CheckAccepts, bool CheckThirdLayer, bool C2, bool C4, bool C6, int nPixelRadiusC2, int nPixelRadiusC4, int nPixelBorderC6, const int laneNumber[], const int laneOffset, const int columnsPerChip, const int rowsPerChip, double nPixelsGap, mTowerChipRobbie* hitsInChip[], int eventID, int nHits, int eventIndex) {

  std::cout << "SELECTIONA" << std::endl;
  
  /////////////////////////////////////////////////////////
  // 4A) CLUSTERING FIRST THREE LAYERS
  vector<vector<int>> vClusters_layer0(4, vector<int> {}); //Vector with for every cluster: lanecode, meanX, meanY, clustersize
  vector<vector<int>> vClusters_layer1(4, vector<int> {}); //Vector with for every cluster: lanecode, meanX, meanY, clustersize
  vector<vector<int>> vClusters_layer2(4, vector<int> {}); //Vector with for every cluster: lanecode, meanX, meanY, clustersize
  int nClusters = 0;
  int nClusters_layer1 = 0;
  int nClusters_layer2 = 0;
  double minsep = 0;
    
  for (int lanecode = 0; lanecode < 2; lanecode++) //loop over chips in first layer to find all clusters
    {
      int l = laneNumber[lanecode];
      if (hitsInChip[l]->getNHits()>0)
	{
	  //get the array of clusters
	  hitsInChip[l]->Clusterize(); // This is the main clustering workhorse. You can find it in mTowerChipRobbie.cxx in the ./classes/ folder
	  TObjArray* clusterlist = hitsInChip[l]->getClusters();
	  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
	    {
	      mTowerClusterRobbie* cluster = (mTowerClusterRobbie*) clusterlist->At(c);
	      if (cluster)
		{
		  nClusters++;
		  
		  //get properties of cluster
		  vClusters_layer0[0].push_back(lanecode); //lanecode
		  
		  //get mean coordinate of cluster.
		  int meanRow = round(cluster -> getMeanRow());
		  int meanColumn = round(cluster -> getMeanColumn());
		  
		  //change chip dependent coordinates to absolute coordinates   
		  int meanX = meanRow;
		  int meanY = meanColumn;
		  if (IsLeftChip(l+laneOffset)) //cluster is in the left chip
		    {
		      meanX = -meanX - 1 - nPixelsGap;
		      meanY = columnsPerChip - 1 - meanY;
		    }
		  meanX += rowsPerChip + nPixelsGap;
		  
		  vClusters_layer0[1].push_back(meanX); //meanX
		  vClusters_layer0[2].push_back(meanY); //meanY
		  
		  vClusters_layer0[3].push_back(cluster->getNHits()); //clustersize
		}
	    }
	}
    }

  //Clustering in second layer
  for (int lanecode = 2; lanecode < 4; lanecode++) //loop over chips in second layer to find all clusters
    {
      int l = laneNumber[lanecode];
      if (hitsInChip[l]->getNHits()>0)
	{
	  //get the array of clusters
	  hitsInChip[l]->Clusterize();
	  TObjArray* clusterlist = hitsInChip[l]->getClusters();
	  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
	    {
	      mTowerClusterRobbie* cluster = (mTowerClusterRobbie*) clusterlist->At(c);
	      if (cluster)
		{
		  nClusters_layer1++;
		  
		  //get properties of cluster
		  vClusters_layer1[0].push_back(lanecode); //lanecode
		  
		  //get mean coordinate of cluster.
		  int meanRow = round(cluster -> getMeanRow());
		  int meanColumn = round(cluster -> getMeanColumn());
		  
		  //change chip dependent coordinates to absolute coordinates   
		  int meanX = meanRow;
		  int meanY = meanColumn;
		  if (IsLeftChip(l+laneOffset)) //cluster is in the left chip
		    {
		      meanX = -meanX - 1 - nPixelsGap;
		      meanY = columnsPerChip - 1 - meanY;
		    }
		  meanX += rowsPerChip + nPixelsGap;
		  
		  vClusters_layer1[1].push_back(meanX); //meanX
		  vClusters_layer1[2].push_back(meanY); //meanY
		  vClusters_layer1[3].push_back(cluster->getNHits()); //clustersize
		}
	    }
	}
    }
  
  for (int lanecode = 4; lanecode < 6; lanecode++) //loop over chips in second layer to find all clusters
    {
      int l = laneNumber[lanecode];
      if (hitsInChip[l]->getNHits()>0)
	{
	  //get the array of clusters
	  hitsInChip[l]->Clusterize();
	  TObjArray* clusterlist = hitsInChip[l]->getClusters();
	  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
	    {
	      mTowerClusterRobbie* cluster = (mTowerClusterRobbie*) clusterlist->At(c);
	      if (cluster)
		{
		  nClusters_layer2++;
		  
		  //get properties of cluster
		  vClusters_layer2[0].push_back(lanecode); //lanecode
		  
		  //get mean coordinate of cluster.
		  int meanRow = round(cluster -> getMeanRow());
		  int meanColumn = round(cluster -> getMeanColumn());
		  
		  //change chip dependent coordinates to absolute coordinates   
		  int meanX = meanRow;
		  int meanY = meanColumn;
		  if (IsLeftChip(l+laneOffset)) //cluster is in the left chip
		    {
		      meanX = -meanX - 1 - nPixelsGap;
		      meanY = columnsPerChip - 1 - meanY;
		    }
		  meanX += rowsPerChip + nPixelsGap;
		  
		  vClusters_layer2[1].push_back(meanX); //meanX
		  vClusters_layer2[2].push_back(meanY); //meanY
		  vClusters_layer2[3].push_back(cluster->getNHits()); //clustersize
		}
	    }
	}
    }

  /////////////////////////////////////////////////////////////////
  // 4B) EVENT SELECTION
  
  bool eventRejected = false; //variable for event criteria
  int acceptedCluster = -1; //the index of the accepted cluster - for now no cluster accepted
  int meanXAC; //the mean x of the accepted cluster
  int meanYAC; //the mean y of the accepted cluster
  int meanXASC = 0; // the mean x of the layer1 cluster used to verify the particle in layer0.
  int meanYASC = 0; // the mean y of the layer1 cluster used to verify the particle in layer0.
  int meanXATC = 0; // the mean x of the layer2 cluster used to verify the particle in layer0.
  int meanYATC = 0; // the mean y of the layer2 cluster used to verify the particle in layer0.
  
  ////////////////////////////////////////////////////////////////
  // 4Bi) CRITERIA C2: IDENTIFY PARTICLES IN LAYER0 USING TRANSVERSE POSITION OF CLUSTERS IN LAYERS 1 AND 2
  
  if (C2 && !(eventRejected))
    {
      int nparticles = 0;
      for (int c = 0; c < nClusters; c++) //loop over clusters for criteria C2
	{
	  if (eventRejected) break;
	  if (vClusters_layer0[0][c] < 2) //if cluster is in first layer
	    {
	      int meanXC = vClusters_layer0[1][c]; //mean x of this cluster
	      int meanYC = vClusters_layer0[2][c]; //mean y of this cluster
	      bool particle = false; // Flag for if this event contains a particle or not.
	      bool doubled = false; // Flag for if this event contains two 'particle' clusters that were identified by the same layer1 cluster.
	      for (int c2 = 0; c2 < nClusters_layer1; c2++)
		{
		  int meanXC2 = vClusters_layer1[1][c2];
		  int meanYC2 = vClusters_layer1[2][c2];
		  if (vClusters_layer1[3][c2] > 1) { //Require that the layer 1 cluster is more than 1 pixel large
		    if (pow(pow(meanXC-meanXC2,2)+pow(meanYC-meanYC2,2),0.5)<nPixelRadiusC2) //Is mean of the layer 1 cluster within shadow of layer 0 cluster?
		      {
			if ((meanXC2 == meanXASC) && (meanYC2 == meanYASC)) { // If true, then we assume this cluster is part of the same particle track as the already accepted particle, so we skip right over it as we already know about this particle!
			  doubled = true;
			  continue;
			}
			particle = true;
			if (CheckRejects || CheckAccepts) {
			  std::cout << "EVENT " << eventID << "  Found a particle in layer0 at (" << meanXC << "," << meanYC << "), from a layer1 cluster at (" << meanXC2 << "," << meanYC2 << ")" << endl;
			}
			if (acceptedCluster != -1) //If there is already an accepted cluster
			  {
			    eventRejected = true;
			    //std::cout << "check1 eventRejected = " << eventRejected << std::endl;
			    if (eventRejected && CheckRejects) {
			      std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
			      std::cout << "EVENT NUMBER " << eventID << " WAS REJECTED BY C2, as there was more than one particle. There were " << nHits << " hits in the event." << endl;
			    }
			  }
			acceptedCluster = c;
			nparticles++;
			meanXAC = meanXC;
			meanYAC = meanYC;
			meanXASC = meanXC2;
			meanYASC = meanYC2;
			break; //stop loop over clusters
		      }
		  }
		}
	      if ((!particle) && (!eventRejected) && (!doubled)) { // Loop over clusters in third layer, do we find any evidence here of a throughgoing particle?
		for (int c3 = 0; c3 < nClusters_layer2; c3++)
		  {
		    int meanXC3 = vClusters_layer2[1][c3];
		    int meanYC3 = vClusters_layer2[2][c3];
		    if ((meanXC3 != meanXATC) && (meanYC3 != meanXATC)) { // This statement avoids doubling up particles being identified by the same layer2 cluster.
		      if (vClusters_layer2[3][c3] > 1) { //Require that the layer 2 cluster is more than 1 pixel large
			if (pow(pow(meanXC-meanXC3,2)+pow(meanYC-meanYC3,2),0.5)<(nPixelRadiusC2)) //Is mean of the layer 2 cluster within shadow of layer 0 cluster?
			  {
			    if (CheckRejects || CheckThirdLayer || CheckAccepts) {
			      std::cout << "EVENT " << eventID << "  Found a particle in layer0 at (" << meanXC << "," << meanYC << "), from a layer2 cluster at (" << meanXC3 << "," << meanYC3 << ")" << endl;
			    }
			    if (acceptedCluster != -1) //If there is already an accepted cluster
			      {
				eventRejected = true;
				//std::cout << "check3 eventRejected = " << eventRejected << std::endl;
				if (eventRejected && CheckRejects) {
				  std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
				  std::cout << "EVENT NUMBER " << eventID << " WAS REJECTED BY C2, as there was more than one particle. There were " << nHits << " hits in the event." << endl;
				}
			      }
			    acceptedCluster = c;
			    nparticles++;
			    meanXAC = meanXC;
			    meanYAC = meanYC;
			    meanXATC = meanXC3;
			    meanYATC = meanYC3;
			    break; //stop loop over clusters
			  }
		      }
		    }
		  }
	      } else if (!doubled) { // This loop will just find a corresponding cluster in layer2, it won't accept or reject as we already have a particle we are happy with!
		for (int c3 = 0; c3 < nClusters_layer2; c3++)
		  {
		    int meanXC3 = vClusters_layer2[1][c3];
		    int meanYC3 = vClusters_layer2[2][c3];
		    if ((meanXC3 != meanXATC) && (meanYC3 != meanXATC)) {
		      if (vClusters_layer2[3][c3] > 1) { //Require that the layer 2 cluster is more than 1 pixel large
			if (pow(pow(meanXC-meanXC3,2)+pow(meanYC-meanYC3,2),0.5)<(nPixelRadiusC2)) //Is mean of the layer 2 cluster within shadow of layer 0 cluster?
			  {
			    if (CheckRejects || CheckThirdLayer) {
			      std::cout << "EVENT " << eventID << "  Particle in layer0 at (" << meanXC << "," << meanYC << "), verified by a layer2 cluster at (" << meanXC3 << "," << meanYC3 << ")" << endl;
			    }
			    meanXATC = meanXC3;
			    meanYATC = meanYC3;
			    break; //stop loop over clusters
			  }
		      }
		    }
		  }
	      } 
	    }//if cluster in first layer
	} //loop over layer0 clusters for criterion C2

      meanXATC = 0;
      meanYATC = 0;
      
      bool caughtaparticle = true;
      if (acceptedCluster == -1) {
	caughtaparticle = false;
      }

      for (int c = 0; c < nClusters_layer1; c++) //loop over clusters for criteria C2
	{
	  if (eventRejected) break;
	  if ((vClusters_layer1[0][c] == 2) || (vClusters_layer1[0][c] == 3)) { //if cluster is in second layer
	    int meanXC2 = vClusters_layer1[1][c];
	    int meanYC2 = vClusters_layer1[2][c];
	    if (vClusters_layer1[3][c] > 1) { //Require that the layer 1 cluster is more than 1 pixel large
	      //something.
	      
	      for (int c3 = 0; c3 < nClusters_layer2; c3++)
		{
		  int meanXC3 = vClusters_layer2[1][c3];
		  int meanYC3 = vClusters_layer2[2][c3];
		  if ((meanXC3 != meanXATC) && (meanYC3 != meanXATC)) { // This statement avoids doubling up particles being identified by the same layer2 cluster.
		    if (vClusters_layer2[3][c3] > 1) { //Require that the layer 2 cluster is more than 1 pixel large
		      if (pow(pow(meanXC2-meanXC3,2)+pow(meanYC2-meanYC3,2),0.5)<(nPixelRadiusC2)) //Is mean of the layer 2 cluster within shadow of layer 1 cluster?
			{
			  if (CheckRejects || CheckThirdLayer || CheckAccepts) {
			    std::cout << "EVENT " << eventID << "  Found a particle in layer1 at (" << meanXC2 << "," << meanYC2 << "), from a layer2 cluster at (" << meanXC3 << "," << meanYC3 << ")" << endl;
			  }
			  if (acceptedCluster != -1) //If there is already an accepted cluster
			    {
			      if (pow(pow(meanXC2-meanXAC,2)+pow(meanYC2-meanYAC,2),0.5)<(2*nPixelRadiusC2)) // Is mean of layer 1 cluster within shadow of accepted cluster?
				{ // If yes, then we can ignore this cluster. It's pretty similar to either an already-found layer 0 cluster, or an already-found layer 1 cluster.
				  if (CheckRejects || CheckThirdLayer || CheckAccepts) {
				    std::cout << "EVENT " << eventID << "This layer 1 particle candidate at (" << meanXC2 << "," << meanYC2 << ") is close enough to already found particle at (" << meanXAC << "," << meanYAC << ") that we can treat them as part of the same shower." << std::endl;
				  }
				} else { // If not, then this is a whole separate track we didn't find before. We now have more than one particle candidate, so the event is rejected.
				eventRejected = true;
				if (eventRejected && CheckRejects) {
				  std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
				  std::cout << "EVENT NUMBER " << eventID << " WAS REJECTED BY C2, as there was more than one particle. There were " << nHits << " hits in the event." << endl;
				}
			      }
			    } else {
			    acceptedCluster = c;
			    nparticles++;
			    meanXAC = meanXC2;
			    meanYAC = meanYC2;
			    meanXATC = meanXC3;
			    meanYATC = meanYC3;
			    break; //stop loop over clusters
			  }
			}
		    }
		  }
		}	      
	    }
	  }
	}

      if ((acceptedCluster != -1) && (!caughtaparticle)) {
	std::cout << "FLAGFLAGFLAG we found a particle from the layer 2 information that we otherwise wouldn't have done!" << std::endl;
      }
      
      if (acceptedCluster == -1) //If there is no cluster accepted
	{
	  eventRejected = true;
	  if (eventRejected && CheckRejects) {
	    std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
	    std::cout << "EVENT NUMBER " << eventID << " WAS REJECTED BY C2, as there were no particles. There were " << nHits << " hits in the event." << endl;
	  }
	}
      if (nparticles == 0) {
	eventRejected = true;
	if (eventRejected && CheckRejects) {
	  std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
	  std::cout << "EVENT NUMBER " << eventID << " WAS REJECTED BY C2, as there were no particles. There were " << nHits << " hits in the event." << endl;
	}
      }
    } //if C2

  ////////////////////////////////////////////////////////////////
  // 4Bii) CRITERIA C6: FIDUCIAL CUT, REMOVING EVENTS WHERE THE 'PARTICLE' CLUSTER IS WITHIN RANGE (USUALLY 50 PIXELS) OF THE EDGE OF THE DETECTOR IN THE TRANSVERSE DIRECTION

  if ((C6) && !(eventRejected)) 
    {
      if ((meanXAC < nPixelBorderC6) || (meanXAC > (columnsPerChip - nPixelBorderC6)) || (meanYAC < nPixelBorderC6) || (meanYAC > (columnsPerChip - nPixelBorderC6)))
	{
	  eventRejected = true;
	}
      if (eventRejected && CheckRejects) {
	std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
	std::cout << "EVENT NUMBER " << eventID << " WAS REJECTED BY C6, AS IT FELL TOO CLOSE TO THE BORDER OF THE CHIP. There were " << nHits << " hits in the event." << endl;
      }
    }

  ////////////////////////////////////////////////////////////////
  // 4Biii) CRITERIA C4: REMOVE EVENTS WHERE A LARGE (NPIXELS > 1) CLUSTER IN LAYER1 IS WITHIN 120 PIXELS IN THE TRANSVERSE DIRECTION OF THE 'PARTICLE' CLUSTER IN LAYER0.
  
  if ((C4) && !(eventRejected))
    {
      for (int c2 = 0; c2 < nClusters_layer1; c2++)
	{
	  int meanXC2 = vClusters_layer1[1][c2]; //mean x of this cluster
	  int meanYC2 = vClusters_layer1[2][c2]; //mean y of this cluster
	  if ((pow(pow(meanXAC-meanXC2,2)+pow(meanYAC-meanYC2,2),0.5) > nPixelRadiusC4) && (vClusters_layer1[3][c2] > 1)) { //Is mean of the layer 1 cluster within shadow of layer 0 cluster? Is the cluster size greater than 1?
	    eventRejected = true;
	    if (eventRejected && CheckRejects && (nHits > 2000)) {
	      std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
	      std::cout << "EVENT NUMBER " << eventIndex << " WAS REJECTED BY C4. There were " << nHits << " hits in the event." << endl;
	    }
	  }
	} //Loop over layer 1 clusters
    } //C4 loop

  /////////////////////////////////////////////////////////////////
  // 4C) SELECTION COMPLETE. RETURN TRUE IS EVENT IS NOT REJECTED, ELSE RETURN FALSE.
  if (eventRejected) {
    return false;
  } else {
    return true;
  }
}

////////////////////////////////////////////////////////////////////
// 5) KT-ALGORITHM BASED SELECTION

bool EventSelectionB_v1(mTowerLayer* clayers[], mTowerChip* cchips[], Double_t nHitsTotNew[], Double_t nClustersTotNew[], Double_t dlim, Double_t dcenter, Int_t minNlayerBulk, Double_t W0Bulk, Double_t thDistBulk, Double_t wBinBulk, TH2D* htmpBulk, TH2D* htmpLim, std::map<Int_t, cellinfo> mcellsBulk, std::map<Int_t, cellinfo> mcellsLim, Int_t nla, Int_t nLayerBulk, Int_t max_layer[]) {

  ///////////////////////////////////////////////////////////////////
  // 5A) FILL CHIPS INTO LAYER
  for(Int_t ich=0; ich<nchip ; ich++){
    auto      ilayer = get<0>( chipid2layer_lut.at(ich) );
    clayers[ilayer]->add_chip( cchips[ich] );
  }

  ////////////////////////////////////////////////////////////////
  // 5B) TOTAL NUMBER OF HITS/CLUSTERS
  for(Int_t i=0; i<nla; i++){
    nClustersTotNew[i] = 0;
    nHitsTotNew    [i] = 0;
  }
  for(Int_t ich=0; ich<nchip ; ich++){
    auto lane  = cchips[ich]->get_nrlane();
    auto nclus = cchips[ich]->get_nClusters();
    auto nhits = cchips[ich]->get_nPixels();
    auto scale = cchips[ich]->get_scale();
    for(Int_t i=0; i<nla; i++)
      if( get<0>( lane2layer_lut.at(lane)) < max_layer[i] ){
	nClustersTotNew[i] += (Double_t)nclus         ;
	nHitsTotNew    [i] += (Double_t)nhits * scale ;
      }
  }

  /////////////////////////////////////////////////////////////
  // 5C) FILL HIT POSITIONS
  for(Int_t ila=0; ila<nlayer; ila++){
    for(UInt_t ic=0; ic<clayers[ila]->get_nChips(); ic++){
      auto chip   = clayers[ila]->get_chip(ic);
      
      for (UInt_t iclus=0; iclus<chip->get_nClusters(); iclus++){
	auto xy   = chip->get_cluster(iclus)->get_pos_glob_in_mm();
	Double_t x = get<0>(xy) ;
	Double_t y = get<1>(xy) ;
	
	//Bulk                                                                                                                                                                                                                        
	if( ila < nLayerBulk ){
	  if( TMath::Abs(x) < dlim && TMath::Abs(y) < dlim ){
	    Int_t    icell = htmpBulk->FindBin(x,y);
	    Double_t xcell = htmpBulk->GetXaxis()->GetBinCenter( htmpBulk->GetXaxis()->FindBin(x) );
	    Double_t ycell = htmpBulk->GetYaxis()->GetBinCenter( htmpBulk->GetYaxis()->FindBin(y) );
	    if( mcellsBulk.find(icell) == mcellsBulk.end() ){
	      mcellsBulk.emplace( icell, cellinfo(xcell, ycell) );
	      mcellsBulk.at(icell).Add(x, y, ila);
	    }else{
	      mcellsBulk.at(icell).Add(x, y, ila);
	    }
	  }
	  else {
	    Int_t    icell = htmpLim->FindBin(x,y);
	    Double_t xcell = htmpLim->GetXaxis()->GetBinCenter( htmpLim->GetXaxis()->FindBin(x) );
	    Double_t ycell = htmpLim->GetYaxis()->GetBinCenter( htmpLim->GetYaxis()->FindBin(y) );
	    if( mcellsLim.find(icell) == mcellsLim.end() ){
	      mcellsLim.emplace( icell, cellinfo(xcell, ycell) );
	      mcellsLim.at(icell).Add(x, y, ila);
	    }else{
	      mcellsLim.at(icell).Add(x, y, ila);
	    }
	  }
	}	
      }//cluster
    }//chip
  }//layer          

  ///////////////////////////////////////////////////////////////////
  // 5D) CELL REMOVAL
  
  //Bulk
  vector<Int_t> vCellRemoveBulk;
  vCellRemoveBulk.clear();
  for ( auto it = mcellsBulk.begin(); it != mcellsBulk.end(); it++ ) {
    auto icell = it->first  ;
    auto info  = it->second ;
    if( info.GetNlayer() < minNlayerBulk || info.GetN() <= W0Bulk )
      vCellRemoveBulk.push_back( icell );
  }
  for( UInt_t i=0; i<vCellRemoveBulk.size(); i++ )
    mcellsBulk.erase(vCellRemoveBulk.at(i));
  
  //Lim
  vector<Int_t> vCellRemoveLim;
  vCellRemoveLim.clear();
  for ( auto it = mcellsLim.begin(); it != mcellsLim.end(); it++ ) {
    auto icell = it->first  ;
    auto info  = it->second ;
    if( info.GetNlayer() < minNlayerBulk || info.GetN() <= W0Bulk )
      vCellRemoveLim.push_back( icell );
  }
  for( UInt_t i=0; i<vCellRemoveLim.size(); i++ )
    mcellsLim.erase(vCellRemoveLim.at(i));
  
  ///////////////////////////////////////////////////////////////
  // 5Ei) PILEUP SUSPICIOUS EVENT (STAGE 1: BULK)

  vector<tuple<Int_t, Double_t, Double_t, Int_t>> vShowerBulk ;
  vShowerBulk.clear();
  {//FJ start
    vector<PseudoJet> particles;
    particles.clear();
    //Bulk
    for ( auto it1 = mcellsBulk.begin(); it1 != mcellsBulk.end(); it1++ ) {
      //cell center
      TVector3 vtmp;
      vtmp.SetPtEtaPhi( (it1->second).GetN(), (it1->second).GetX()/10., (it1->second).GetY()/10. );
      particles.push_back( PseudoJet( vtmp.Px(), vtmp.Py(), vtmp.Pz(), vtmp.Mag()) );
    }
    //Lim
    for ( auto it1 = mcellsLim.begin(); it1 != mcellsLim.end(); it1++ ) {
      //cell center
      TVector3 vtmp;
      vtmp.SetPtEtaPhi( (it1->second).GetN(), (it1->second).GetX()/10., (it1->second).GetY()/10. );
      particles.push_back( PseudoJet( vtmp.Px(), vtmp.Py(), vtmp.Pz(), vtmp.Mag()) );
    }
    
    Double_t R = thDistBulk / 10.;
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    
    // print the jets
    for (unsigned i = 0; i < jets.size(); i++) {
      vector<PseudoJet> constituents = jets[i].constituents();
    }
    for (unsigned i = 0; i < jets.size(); i++) {
      if( jets[i].constituents().size() < 2 ) return false;
      
      vShowerBulk.push_back(make_tuple(
				       (Int_t)(jets[i].pt()),
				       10.*jets[i].rap(),
				       10.*jets[i].phi_std(),
				       (Int_t)(jets[i].constituents().size())
				       ));
    }
  }//FJ end
  
    ///////////////////////////////////////////////
  // 5Eii) EXTRACT STAGE 1 STATUS

  Bool_t status0 = kFALSE;
  Bool_t status1 = kFALSE;
  {
    Int_t nCenter = 0;
    Int_t nLim    = 0;
    Int_t nAll    = 0;
    for( UInt_t i=0; i<vShowerBulk.size(); i++){
      Int_t    a = get<0>(vShowerBulk.at(i)) ;
      Int_t    n = get<3>(vShowerBulk.at(i)) ;
      Double_t x = get<1>(vShowerBulk.at(i)) ;
      Double_t y = get<2>(vShowerBulk.at(i)) ;
      
      nAll++;
      if( TMath::Abs(x) < dcenter && TMath::Abs(y) < dcenter ) nCenter++;
      if( TMath::Abs(x) < dlim    && TMath::Abs(y) < dlim    ) nLim++   ;
    }
    if( nLim   ==nAll && nLim    == 1 ) status0 = kTRUE;
    if( nCenter==nAll && nCenter == 1 ) status1 = kTRUE;
  }

  if( !status1 ){ return false; }

  ////////////////////////////////////////////////////
  // 5Eiii) PILEUP SUSPICIOUS EVENT (STAGE 2: CORE)
  
  Bool_t isSelect = (vShowerBulk.size() < 1 )? kFALSE : kTRUE ;
  
  vector<pair<Double_t, Double_t>> vclus; vclus.clear();
  
  for( UInt_t i=0; i<vShowerBulk.size(); i++){
    Double_t xs = get<1>(vShowerBulk.at(i)) ;
    Double_t ys = get<2>(vShowerBulk.at(i)) ;
    
    Int_t nclus = 0;
    for(Int_t ila=0; ila<2; ila++){
      for(UInt_t ic=0; ic<clayers[ila]->get_nChips(); ic++){
	auto chip   = clayers[ila]->get_chip(ic);
	for (UInt_t iclus=0; iclus<chip->get_nClusters(); iclus++){
	  auto xy   = chip->get_cluster(iclus)->get_pos_glob_in_mm();
	  Double_t xc = get<0>(xy) ;
	  Double_t yc = get<1>(xy) ;
	  Double_t rc = 0.;
	  rc += TMath::Power(xc-xs,2);
	  rc += TMath::Power(yc-ys,2);
	  rc =  TMath::Sqrt(rc);
	  if( rc > thDistBulk ) continue;
	  nclus++;
	  
	  isSelect &= ( rc < 2.*wBinBulk );
	  isSelect &= (( ila==0 && nclus>1 )? kFALSE : kTRUE );
	  vclus.push_back( xy );
	}
      }
    }
    isSelect &= ( nclus >= 2 );
    if(!isSelect) {
      return false;
    }
    
    for (UInt_t iclus=0    ; iclus<vclus.size(); iclus++){
      for (UInt_t jclus=iclus; jclus<vclus.size(); jclus++){
	Double_t x0 = get<0>(vclus.at(iclus)) ;
	Double_t y0 = get<1>(vclus.at(iclus)) ;
	Double_t x1 = get<0>(vclus.at(jclus)) ;
	Double_t y1 = get<1>(vclus.at(jclus)) ;
	Double_t r = 0.;
	r += TMath::Power(x1-x0,2);
	r += TMath::Power(y1-y0,2);
	r =  TMath::Sqrt(r);
	isSelect &= ( r < 0.5 );
	if(!isSelect) return false;
      }
      if(!isSelect) return false;
    }
  }

  /////////////////////////////////////////////////////////////////
  // 5F) IF EVENT HAS NOT ALREADY BEEN REJECTED, CHECK TO SEE IF WE SELECT EVENT, AND RETURN RELEVANT BOOL.
  if( isSelect ){
    return true;
  } else {
    return false;
  }
}
