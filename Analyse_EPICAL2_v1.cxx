///////////////////////////////////////////////////////////////////////////////////
// 1)  INCLUDES

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

// Includes for analysis processors (Taken from local processors directory)
#include "processors/EventSelection_v1.h"

// Includes for kT-Algorithm Selection
/*#include "fastjet/ClusterSequence.hh"                                                   // Requires installation of Fastjet.
#include "/eos/project/m/mtower/public/hiroki/00_util/util.hh"
#include "/eos/project/m/mtower/public/hiroki/00_util/mTowerUtil.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_hot_pixels.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_module_thickness.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_inclination_parameters.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_alignment_parameters.h"
#include "/eos/project/m/mtower/public/hiroki/11_id_pileup/cellinfo.h"
*/
// Adding in namespaces
using namespace std;
using namespace fastjet;
using namespace ROOT::Math;

// Define variables for use with loading/saving event selection data
#define SQS  // true if intending to save results of event selection to a local file
//#define LQS  // true if intending to use local file to perform quick event selection


// Including Processors:
// None yet.

////////////////////////////////////////////////////////////////////////////
// 2) LOOKUP TABLES

// Lookup tables for criteria selection
const std::map< Int_t, Int_t > chipid2lanerobbie_lut = {
  { 0,40},{ 1,39},{ 2,42},{ 3,41},{ 4,44},{ 5,43},{ 6,46},{ 7,45},
  { 8,48},{ 9,47},{10,50},{11,49},{12,52},{13,51},{14,54},{15,53},
  {16,38},{17,55},{18,36},{19,37},{20,32},{21,35},{22,34},{23,33},
  {24,64},{25,63},{26,66},{27,65},{28,68},{29,67},{30,70},{31,69},
  {32,72},{33,71},{34,74},{35,73},{36,76},{37,75},{38,78},{39,77},
  {40,62},{41,79},{42,60},{43,61},{44,56},{45,59},{46,58},{47,57}
};

const std::map< Int_t, Int_t > lane2chipidrobbie_lut = {
  {40, 0},{39, 1},{42, 2},{41, 3},{44, 4},{43, 5},{46, 6},{45, 7},
  {48, 8},{47, 9},{50,10},{49,11},{52,12},{51,13},{54,14},{53,15},
  {38,16},{55,17},{36,18},{37,19},{32,20},{35,21},{34,22},{33,23},
  {64,24},{63,25},{66,26},{65,27},{68,28},{67,29},{70,30},{69,31},
  {72,32},{71,33},{74,34},{73,35},{76,36},{75,37},{78,38},{77,39},
  {62,40},{79,41},{60,42},{61,43},{56,44},{59,45},{58,46},{57,47}
};

const std::map< Int_t, Int_t > chipid2layerrobbie_lut = {
  { 0,22},{ 1,22},{ 2,20},{ 3,20},{ 4,18},{ 5,18},{ 6,16},{ 7,16},
  { 8,14},{ 9,14},{10,12},{11,12},{12,10},{13,10},{14, 8},{15, 8},
  {16, 6},{17, 6},{18, 4},{19, 4},{20, 0},{21, 0},{22, 2},{23, 2},
  {24,23},{25,23},{26,21},{27,21},{28,19},{29,19},{30,17},{31,17},
  {32,15},{33,15},{34,13},{35,13},{36,11},{37,11},{38, 9},{39, 9},
  {40, 7},{41, 7},{42, 5},{43, 5},{44, 1},{45, 1},{46, 3},{47, 3}
};

const std::map< Int_t, Int_t > lane2layerrobbie_lut = {
  {40,22},{39,22},{42,20},{41,20},{44,18},{43,18},{46,16},{45,16},
  {48,14},{47,14},{50,12},{49,12},{52,10},{51,10},{54, 8},{53, 8},
  {38, 6},{55, 6},{36, 4},{37, 4},{32, 0},{35, 0},{34, 2},{33, 2},
  {64,23},{63,23},{66,21},{65,21},{68,19},{67,19},{70,17},{69,17},
  {72,15},{71,15},{74,13},{73,13},{76,11},{75,11},{78, 9},{77, 9},
  {62, 7},{79, 7},{60, 5},{61, 5},{56, 1},{59, 1},{58, 3},{57, 3}
};

const std::map<int,bool> layer2isInvrobbie_lut = {
  { 0, kFALSE}, { 1, kTRUE}, { 2, kFALSE}, { 3, kTRUE}, 
  { 4, kFALSE}, { 5, kTRUE}, { 6, kFALSE}, { 7, kTRUE}, 
  { 8, kFALSE}, { 9, kTRUE}, {10, kFALSE}, {11, kTRUE}, 
  {12, kFALSE}, {13, kTRUE}, {14, kFALSE}, {15, kTRUE}, 
  {16, kFALSE}, {17, kTRUE}, {18, kFALSE}, {19, kTRUE}, 
  {20, kFALSE}, {21, kTRUE}, {22, kFALSE}, {23, kTRUE} 
};

///////////////////////////////////////////////////////////
// 3) USEFUL FUNCTIONS

// Chip-based Functions included from EventSelection.h


// Functions for Hit Maps
int lane2padhitmap(int lane){
  int layerNr = lane2layerrobbie_lut.at(lane);
  bool isLeft = IsLeftChip(lane);
  
  int padid;
  if (layerNr < 6) {padid = layerNr+1;}
  else if (layerNr < 12) {padid = layerNr+7;}
  else if (layerNr < 18) {padid = layerNr+13;}
  else {padid = layerNr+19;}
  if (isLeft) {padid += 6;}
 
  return padid; 
}

int lane2padoccupancy(int lane){
  int layerNr = lane2layerrobbie_lut.at(lane);
  int padid = layerNr+1;
  return padid; 
}


////////////////////////////////////////////////////////////
// 4) Run Analysis
void Analyse_mTower(int run, Double_t energy) // This is the main workhorse. Any added processors will find their place somewhere within here. 'energy_in' should be given in GeV.
{
  //////////////////////////////////////////////////////////
  // 4A) Set hard variables
  // Hard Variables for Using Quick Selection From Local File
  TString fileLocationQuickSelection = ""; // File being used to load quick selection information

  // Hard Variable for choosing selection(s) to use
  Int_t selection = 1;
  /*
    If selection = 0, then no selection is performed. All events are treated as acceptable.
    If selection = 1, then only events which pass the criteria-based cuts are accepted.
    If selection = 2, then only events which pass the kT-algorithm-based cuts are accepted.
    If selection = 3, then only events which pass the criteria-based cuts AND the kT-algorithm-based cuts are accepted.
   */
  
  // Hard Variables for kT-Algorithm Selection
  Double_t  wBinBulk      = 0.50; //[mm]
  Int_t     nLayerBulk    = 24;
  Int_t     minNlayerBulk = 3;
  Double_t  W0Bulk        = 2.;
  Double_t  thDistBulk    = 5.;  
  Int_t     ith           = 0;
  Int_t     ndiv          = 1;
  Double_t dcenter =  5.;
  Double_t dlim    = 12.;
  Double_t wBinLim = 2.*wBinBulk;
  TString s_thickness 	= "/eos/project/m/mtower/public/hiroki/event_selection_Nov2020/00_util/setting/module_thickness_Feb_2020.txt"     ;
  TString s_alignment 	= "/eos/project/m/mtower/public/hiroki/event_selection_Nov2020/00_util/setting/chip_correction_Oct_2020.txt"      ;
  TString s_inclination	= "/eos/project/m/mtower/public/hiroki/event_selection_Nov2020/00_util/setting/detector_inclination_Oct_2020.txt" ;

  // Hard Variables for Criteria Analysis
  TString fileLocationOutputFile = "./";
  TString fileLocation = "/eos/project/m/mtower/Data/Data_TB_February_2020/mTower_Data_DESY_Feb_2020_raw1/"; //The location of the selected data
  TString maskingFileLocation = "/eos/project/m/mtower/public/analysis_fw_sample_HY/data/hotpixel_TB/"; //The location of the masking .txt files
  bool testing = false; //Testing outputfile, run over small amount of events.
  bool DB = false; //debug print statements
  bool HP = false; //find hot pixels
  bool CT = true; //Create TTree
  bool C2 = true; //Criterion 2. Every cluster in the first layer with hits behind it in layer 2 is accepted. Only events with 1 accepted cluster are accepted for analysis.
  bool C4 = false; //Criterion 4. Only events where there are no hits outside of a certain area in the second layer centrated around the accepted cluster (C2) in the first layer are accepted for analysis. NEEDS C2
  bool C6 = true; // Criterion 6. Only events with clusters which are not within a certain number of pixels of the layer border are accepted for analysis. NEEDS C2
  // Activate some debug tools?
  bool CheckRejects = false;
  bool CheckThirdLayer = false;
  bool CheckAccepts = false;
  int nPixelRadiusC2 = 10; //The search area behind a cluster is a circle centered around the average coordinate of the cluster. nPixelRadiusC2 is the radius of that circle.
  int nPixelBorderC6 = 50; //The width of the border in which to reject clusters. If a 'particle' cluster is found in layer0 within this many pixels of the layer's border, that event will be rejected. 
  int nPixelRadiusC4 = 120; //Radius of circle for criterion C4
  int nParts = 1; //In how many parts is this run done?
  int part = 1; //Which part is this?
  const int maxNChips = 48; //48 chips is the maximum size of mTower
  const int laneOffset = 32; 
  std::cout << "*** lane offset used is "<<laneOffset<<endl;
  //for the mTower cosmics test and 2019 and 2020 test beam the offset is 32, for the testbeam data of 2018 there is no offset
  const int rowsPerChip = 512; //512 for ALPIDE chip
  const int columnsPerChip = 1024; //1024 for ALPIDE chip
  double nPixelsGap = 5; //Width of gap between chips in units of pixels
  const int laneNumber[6] = {0,3,27,24,2,1}; //corresponds 'lanecode' with laneNumber-laneOffset. Element i is a chip in the (i/2)th layer. Lane 32, 35 in layer 0, 59, 56 in l1, 34, 33 in l2

  /////////////////////////////////////////////////////////////
  // 4B) Calculate Other Variables

  // Some Standard Variables
  std::ios oldState(nullptr);
  oldState.copyfmt( std::cout );
  TH1::StatOverflows(kTRUE);
  gStyle->SetPalette(kRainBow);
  gStyle->SetLineScalePS(1.0);
  
  // Setting Up Variables For Criteria Selection
  if (!(C2)) {C4 = false; C6 = false;} //If C2 is not applied, all of these criteria cannot be applied
  TString beamEnergy; //for masking inputfile
  if (run == 1309 || run == 1310 || run == 1346 || run == 1375 || run == 1376) {
    beamEnergy = "5R8";
    energy = 5.8;
  }
  else if (run == 1245 || run == 1250 || run == 1261 || run == 1308 || run == 1333 || run == 1339 || run == 1413) {
    beamEnergy = "5R0";
    energy = 5.0;
  }
  else if (run == 1257 || run == 1272 || run == 1274 || run == 1275 || run == 1338 || run == 1345) {
    beamEnergy = "4R0";
    energy = 4.0;
  }
  else if (run == 1335 || run == 1341 || run == 1262) {
    beamEnergy = "3R0";
    energy = 3.0;
  }
  else if (run == 1276 || run == 1337 || run == 1344) {
    beamEnergy = "2R0";
    energy = 2.0;
  }
  else if (run == 1336 || run == 1343 || run == 1263) {
    beamEnergy = "1R0";
    energy = 1.0;
  }
  else {
    beamEnergy = "0";
    energy = 0.0;
  }
  int maxX = 2*rowsPerChip + nPixelsGap - 1;
  int maxY = columnsPerChip - 1;
  //set plain style for histograms
  gROOT->SetStyle("Plain");
  //extra masking
  TH3F* hMaskPtn = new TH3F("hMaskPtn","Hit mask",maxNChips, 0, maxNChips, columnsPerChip, 0, columnsPerChip, rowsPerChip, 0, rowsPerChip); //Histogram with pixels to mask: lane, column, row
  if (beamEnergy != "0")
    {
      ifstream in;
      TString maskingFileName = "outMaskPtn_TB_";
      maskingFileName += beamEnergy;
      maskingFileName += "_GeV.txt";
      in.open(Form(maskingFileName,maskingFileLocation.Data()));
      Float_t chip_id, nr_lane, hot_pixel_column, hot_pixel_row, pixel_entry;
      Double_t difference, average, std_dev;
      Int_t nlines = 0;
      while (1)
	{
	  in >> chip_id >> nr_lane >> hot_pixel_column >> hot_pixel_row >> pixel_entry;
	  if (!in.good()) break;
	  hMaskPtn->Fill(nr_lane-laneOffset,hot_pixel_column,hot_pixel_row);
	  nlines++;
	}
      in.close();
    }
  
  // Setting Up Variables For kT-Algorithm Selection
  TString  suffix0 = Form("%02d_GeV", (Int_t)( 10.*energy ));
  TString  suffix  = Form("%02d_GeV_b%03d_n%02d_m%02d_w%03d_d%03d",
      (Int_t)( 10.*energy       ),
      (Int_t)(100.*wBinBulk     ),
      (Int_t)(     nLayerBulk   ),
      (Int_t)(     minNlayerBulk),
      (Int_t)( 10.*W0Bulk       ),
      (Int_t)( 10.*thDistBulk   )
      );
  auto aposz    = mTowerThickness::get_sensor_position( s_thickness.Data() );
  auto paralign = mTowerAlignment::get_sensor_position( s_alignment.Data() );
  auto parinc   = mTowerInclination::get_inclination_function( s_inclination.Data(), energy );
  const Int_t nla = 4;
  Int_t max_layer[nla] = {6,12,18,24};
  
  ////////////////////////////////////////////////////////////
  // 4C) Set Up Histograms, Input and Output Files

  // Create Output File For TTrees
  TString baseName = "Run_";
  TFile* outputFile;
  fileLocationOutputFile += baseName;
  fileLocationOutputFile += run;
  fileLocationOutputFile += "_PostAnalysis_v1.root";
  if (!(CT)) fileLocationOutputFile = "tobedeleted"; //Otherwise a previously made file might be deleted later
  outputFile = new TFile(fileLocationOutputFile,"recreate");

  // Find and Open the Input Root File (from eos)
  fileLocation += baseName;
  fileLocation += run;
  fileLocation += "/rootout_raw/conv_Run_";
  fileLocation += run;
  fileLocation +=".root";
  std::cout<<endl<<"*** Reading file: "<<fileLocation<<endl<<endl<<"Ignore the following warning if there is one."<<endl;
  TFile* inputFile = TFile::Open(fileLocation);
  std::cout<<endl;

  // Get pixels to be masked
  TString s_file, s_ext, s_path ; s_file = GetFileName (fileLocation, s_path, s_ext );
  TString stmp      = Form("outMaskPtn_TB_%dR%d_GeV.txt", (Int_t)energy, (Int_t)(energy*10)%10 );
  TString s_maskptn = Form("%s/hotpixel_TB/%s",s_path.Data(),stmp.Data());
  vector<pair<Int_t,Int_t>> vmask[nchip];
  cout<< "Load hot pixels from "<< s_maskptn.Data() <<"." <<endl;
  for(Int_t ich=0; ich<nchip; ich++) vmask[ich] = mTowerHotPixels::get_hot_pixels( s_maskptn.Data(), chipid2lane_lut.at( ich ) );

  // Read the Data from Input Root File
  inputFile->cd(); //set to current directory
  TTree* frames = (TTree*)inputFile->Get("Frames");
  int runNumber, fileNumber, eventIndex, nHits, dataSize, eventNumberOriginal;
  UInt_t eventID;
  vector<Int_t>* vlane = new vector<Int_t>();
  vector<Int_t>* vcolumn = new vector<Int_t>();
  vector<Int_t>* vrow = new vector<Int_t>();
  Long64_t packetState  ;
  vector<Int_t>*vst_lane   = new vector<Int_t>();
  vector<Int_t>*vst_error  = new vector<Int_t>();
  vector<Int_t>*vst_roflag = new vector<Int_t>();
  frames->SetBranchAddress("runNumber",&runNumber);
  frames->SetBranchAddress("fileNumber",&fileNumber);
  frames->SetBranchAddress("eventNumber",&eventIndex); //Note that this does not run from 0 to nEvents but resets to 0 sometimes (usually after 63000 events)
  frames->SetBranchAddress("eventID",&eventID); //Runs from 0 to nEvents  
  frames->SetBranchAddress("nHits",&nHits);
  //frames->SetBranchAddress("dataSize",&dataSize);
  frames->SetBranchAddress("lane",&vlane);
  frames->SetBranchAddress("column",&vcolumn);
  frames->SetBranchAddress("row",&vrow);
  frames->SetBranchAddress("packetState" , &packetState );
  frames->SetBranchAddress("st_lane"     , &vst_lane    );
  frames->SetBranchAddress("st_error"    , &vst_error   );
  frames->SetBranchAddress("st_roflag"   , &vst_roflag  );
  int nEvents = frames->GetEntries();
  int nEventsOriginal = nEvents;
  
  // If using preprepared selection data from a local file, load in this data.
#ifdef LQS 
  TFile* inputSelectionFile = TFile::Open(fileLocationQuickSelection);
  inputSelectionFile->cd();
  TTree* inselectiontree = (TTree*)inputSelectionFile->Get("prepreparedselections");
  Int_t inselectionEvent, inselectionRunNumber, inselectionFileNumber;
  Bool_t inselectionHirokiStatus, inselectionRobbieStatus;
  inselectiontree->SetBranchAddress("eventID",&inselectionEvent);
  inselectiontree->SetBranchAddress("runNumber",&inselectionRunNumber);
  inselectiontree->SetBranchAddress("fileNumber",&inselectionFileNumber);
  inselectiontree->SetBranchAddress("status_selection_B",&inselectionHirokiStatus);
  inselectiontree->SetBranchAddress("status_selection_A",&inselectionRobbieStatus);
#endif

  // Create Output Tree(s) and Add Branches
  outputFile->cd();
  TTree* outtree = new TTree("tree","Output Tree");  
  outtree->Branch("runNumber",&runNumber,"r/I");
  outtree->Branch("fileNumber",&fileNumber,"f/I");
  outtree->Branch("eventNumber",&eventIndex,"e/I");
  outtree->Branch("eventID",&eventID,"eid/I"); //This number does not reset to 0 sometimes, it goes from 0 to nEvents
  outtree->Branch("nHits",&nHits,"n/I");
  //outtree->Branch("dataSize",&dataSize,"d/I");
  outtree->Branch("lane",&vlane);
  outtree->Branch("column",&vcolumn);
  outtree->Branch("row",&vrow);
  outtree->Branch("packetState",&packetState);
  outtree->Branch("st_lane"     , &vst_lane    );
  outtree->Branch("st_error"    , &vst_error   );
  outtree->Branch("st_roflag"   , &vst_roflag  );

  // Create Tree for Storing Quick Pre-Prepared Selections, if required.
#ifdef SQS
  TTree* prepreparedtree = new TTree("prepreparedselections","Pre-prepared Event Selections");
  Int_t  prepreparedEvent;
  Int_t  prepreparedRunNumber;
  Int_t  prepreparedFileNumber;
  Bool_t  prepreparedHiroki;
  Bool_t  prepreparedRobbie; 
  prepreparedtree->Branch("eventID", &prepreparedEvent, "eventID/I");
  prepreparedtree->Branch("status_selection_B", &prepreparedHiroki, "status_selection_B/O");
  prepreparedtree->Branch("status_selection_A", &prepreparedRobbie, "status_selection_A/O");
  prepreparedtree->Branch("runNumber", &prepreparedRunNumber, "runNumber/I");
  prepreparedtree->Branch("fileNumber", &prepreparedFileNumber, "fileNumber/I");
#endif

  ////////////////////////////////////////////////////////////
  // 4D) Pre-prepare Some Objects For Loop Over Events

  // Initialise layer/chip objects, temporary-use variables and temporary-use histograms for use in kT-algorithm-based selection
  mTowerLayer* clayers [nlayer];
  mTowerChip*  cchips  [nchip ] ;
  for(Int_t ila=0; ila<nlayer; ila++) clayers[ila] = new mTowerLayer( ila );
  for(Int_t ich=0; ich<nchip ; ich++) cchips [ich] = new mTowerChip ( ich );
  for(Int_t ich=0; ich<nchip ; ich++) cchips [ich]->set_method( mTowerChip::k8N ) ;
  for(Int_t ila=0; ila<nlayer; ila++) clayers[ila]->set_zpos( aposz[ila] );
  Int_t nbin = 1;
  nbin = (Int_t)(16./wBinBulk);
  TH2D* htmpBulk = new TH2D("htmpBulk", "htmpBulk", 2*nbin, -1.*wBinBulk*nbin, wBinBulk*nbin, 2*nbin, -1.*wBinBulk*nbin, wBinBulk*nbin);
  nbin = (Int_t)(16./wBinLim );
  TH2D* htmpLim  = new TH2D("htmpLim" , "htmpLim ", 2*nbin, -1.*wBinLim *nbin, wBinLim *nbin, 2*nbin, -1.*wBinLim *nbin, wBinLim *nbin);
  std::map<Int_t, cellinfo> mcellsBulk{};
  std::map<Int_t, cellinfo> mcellsLim {};  
  Int_t eff_den1 = 0; Int_t eff_neu1 = 0;
  
  // Prepare Useful General Objects For Loop Over Events
  std::cout<<"*** Loop over events in input file"<<endl;
  int nReadEvents = 0; //count how many events the loop has read
  int nEmptyEvents = 0; //count the number of "events" that have less than 1 hits    
  int minEvent = nEvents/nParts*(part-1);
  int maxEvent = nEvents/nParts*part;  
  if (testing)
    {
      minEvent = 0;
      maxEvent = 1000; // If testing, use a small range of events
    }
  Bool_t isFlushTree = kFALSE ;

  ////////////////////////////////////////////////////////////
  // 4E) Loop Over All Events
  for (int event = minEvent; event < maxEvent; event++)
    {
      frames->GetEntry(event);
#ifdef LQS
      inselectiontree->GetEntry(event);
#endif
      prepreparedEvent = event;
      
      if(( event %  10000 == 0) && DB) cout<<"event = "<< event <<endl;
      
      ////////////////////////////////////////////////////////////
      // 4Ei) Quick Initialisation of a few Necessary steps for kT-algorithm Selection
      
      // A few extra parameters, and a reset for kT-algorithm selection temporary histograms
      Double_t  nHitsTotNew    [nla] = {0};
      Double_t  nClustersTotNew[nla] = {0};
      mcellsBulk.clear(); htmpBulk->Reset("ICESM");
      mcellsLim.clear() ; htmpLim ->Reset("ICESM");

      // Check Inclination Parameter Exists
      UInt_t runPeriod = ( runNumber < 1413 )? 0 : 1 ;
      if(runPeriod >= parinc.size() ){ 
	cerr<<" === wrong run number ===  or  === failure on loading inclination parameter ===" <<endl; 
	return; 
      }
      
      // Check Event Status
      Bool_t   isgood = check_event_quality( vlane, vst_lane, vst_error, vst_roflag );
      if( !isgood ) continue;

      // Set Transform Parameters
      Double_t ptmp[4] = {0};
      for(Int_t ich=0; ich<nchip ; ich++){
	auto  nrlayer = chipid2layer_lut.at( ich );
	Int_t layer   = get<0>(nrlayer);
	ptmp[0] = get<0>(paralign.at(ich)) - get<0>(parinc[runPeriod])->Eval( aposz[layer] );
	ptmp[1] = get<1>(paralign.at(ich)) - get<1>(parinc[runPeriod])->Eval( aposz[layer] );
	ptmp[2] = get<2>(paralign.at(ich)) ;
	ptmp[3] = get<3>(paralign.at(ich)) ;
	cchips [ich] ->set_trans( ptmp[0], ptmp[1], ptmp[2], ptmp[3] );
      }

      // Reset chip/layer Objects
      for(Int_t ich=0; ich<nchip ; ich++) cchips [ich]->clear_chip() ;
      for(Int_t ila=0; ila<nlayer; ila++) clayers[ila]->clear_layer();

      ////////////////////////////////////////////////////////////
      // 4Eii) Create mTowerEvent and mTowerHit objects, and add Pixel Mask
      mTowerEvent* currentEvent = new mTowerEvent(runNumber,eventIndex);
      currentEvent->setNHits(nHits);
      currentEvent->setNChips(maxNChips);
      TObjArray* hitList = currentEvent->getHits();
      if (DB)
	{
	  std::cout<<endl<<"(DB) Run: "<<runNumber<<", event: "<<event<<"/"<<nEvents-1<<", number of hits: "<<nHits;
	  if (DB) std::cout<<", original event number: "<<eventNumberOriginal;
	  std::cout<<endl<<"(DB) Loop over hits in event "<<event<<" to add hits to hitlist and apply extra mask"<<endl;
	}
      int nHitsMasked = 0; //number of extra hits masked
      for (int hit=0; hit<nHits; hit++)
	{
	  mTowerHit* currentHit = new mTowerHit();
	  Int_t lane = vlane->at(hit);
	  Int_t col = vcolumn->at(hit);
	  Int_t row = vrow->at(hit);
	  Int_t chipid  = lane2chipidrobbie_lut.at( lane );
	  if (hMaskPtn->GetBinContent(lane-laneOffset+1,col+1,row+1) == 0) //extra masking
	    {
 	      currentHit->setCoordinates(lane, col, row);
	      hitList->Add(currentHit);
	      cchips[chipid]->add_hit( col, row );
	    }
	  else nHitsMasked++;
	}
      if (DB) std::cout<<"(DB) "<<nHitsMasked<<" extra hits masked in this event"<<endl;

      
      ////////////////////////////////////////////////////////////
      // 4Eiii) Read the Event, & Begin Analysis

      nReadEvents++;
      int entries = hitList->GetEntries();
      if (entries < 1) nEmptyEvents++;
      else
      	{
	  //Initialising variables
	  double meanCol[maxNChips] = {0.0};
	  double mean2Col[maxNChips] = {0.0};
	  double meanRow[maxNChips] = {0.0};
	  double mean2Row[maxNChips] = {0.0};
	  int nHitsPerLane[maxNChips]={0};
	  mTowerChipRobbie* hitsInChip[maxNChips];
	  vector<vector<int>> hitsInLayer(maxNChips/2, vector<int>{});
	  for (int l = 0; l<maxNChips ; l++)
	    {
	      hitsInChip[l] = new mTowerChipRobbie(l);
	      hitsInChip[l]->setLane(l+laneOffset);
	    }

	  if (DB) std::cout<<"(DB) Loop over hits in event "<<event<<" to get properties"<<endl;
	  for (int hit = 0; hit < entries; hit++)
	    {
	      mTowerHit* currentHit = (mTowerHit*)hitList->At(hit);
	      int lane = currentHit->getLane();
	      int row = currentHit->getRow();
	      int column = currentHit->getColumn();
	      
	      if ((lane-laneOffset) <maxNChips &&(lane-laneOffset)>-1 )
		{	  
		  // Fill any relevant histograms, if required, and add hits to chip object.
		  if (DB) std::cout << "Early Histograms filled." << std::endl;
		  hitsInLayer[lane2layerrobbie_lut.at(lane)].push_back(hit);
		  hitsInChip[lane-laneOffset]->AddHit(currentHit);
		}
	      else
		{
		  std::cout<<"lane number of hit "<<hit<<" out of range: "<<lane<<endl;
		}
	      
	    } //loop over entries

	  ////////////////////////////////////////////////////////////
	  // 4Eiv) Perform Event Selection
	  bool IsGood = false;
	  Int_t selection_choice;
	  Int_t selection_status = 0;
#ifdef LQS
	  if ((runNumber != inselectionRunNumber) || (fileNumber != inselectionFileNumber) || (eventID != inselectionEvent)) {
	    std::cerr << "Problem! Selection event/run/file number doesn't match up with data. Aborting!" << std::endl;
	    break;
	  }
	  else {
	    switch (selection) {
	    case 0:  // No Selection Performed.
	      IsGood = true;
	      break;
	    case 1:
	      IsGood = inselectionRobbieStatus;
	      break;
	    case 2:
	      IsGood = inselectionHirokiStatus;
	      break;
	    case 3:
	      if (inselectionRobbieStatus && inselectionHirokiStatus) {
		IsGood = true;
	      }
	      break;
	    default:
	      std::cerr << "Problem with 'selection' variable." << std::endl;
	      break;
	    }
	  }
	  
#else 
#ifdef SQS
	  selection_choice = 4;
#else 
	  selection_choice = selection;
#endif
	  std::cout << "selection_choice = " << selection_choice << std::endl;
	  switch (selection_choice) {
	  case 0:  // No Selection Performed.
	    IsGood = true;
	    break;
	  case 1:
	    IsGood = EventSelectionA_v1(CheckRejects, CheckAccepts, CheckThirdLayer, C2, C4, C6, nPixelRadiusC2, nPixelRadiusC4, nPixelBorderC6, laneNumber, laneOffset, columnsPerChip, rowsPerChip, nPixelsGap, hitsInChip, eventID, nHits, eventIndex);
	    break;
	  case 2:
	    IsGood = EventSelectionB_v1(clayers, cchips, nHitsTotNew, nClustersTotNew, dlim, dcenter, minNlayerBulk, W0Bulk, thDistBulk, wBinBulk, htmpBulk, htmpLim, mcellsBulk, mcellsLim, nla, nLayerBulk, max_layer);
	    break;
	  case 3:
	    if (EventSelectionA_v1(CheckRejects, CheckAccepts, CheckThirdLayer, C2, C4, C6, nPixelRadiusC2, nPixelRadiusC4, nPixelBorderC6, laneNumber, laneOffset, columnsPerChip, rowsPerChip, nPixelsGap, hitsInChip, eventID, nHits, eventIndex) && EventSelectionB_v1(clayers, cchips, nHitsTotNew, nClustersTotNew, dlim, dcenter, minNlayerBulk, W0Bulk, thDistBulk, wBinBulk, htmpBulk, htmpLim, mcellsBulk, mcellsLim, nla, nLayerBulk, max_layer)) {
	      IsGood = true;
	    }
	    break;
	  case 4:
	    prepreparedRobbie = EventSelectionA_v1(CheckRejects, CheckAccepts, CheckThirdLayer, C2, C4, C6, nPixelRadiusC2, nPixelRadiusC4, nPixelBorderC6, laneNumber, laneOffset, columnsPerChip, rowsPerChip, nPixelsGap, hitsInChip, eventID, nHits, eventIndex);
	    prepreparedHiroki = EventSelectionB_v1(clayers, cchips, nHitsTotNew, nClustersTotNew, dlim, dcenter, minNlayerBulk, W0Bulk, thDistBulk, wBinBulk, htmpBulk, htmpLim, mcellsBulk, mcellsLim, nla, nLayerBulk, max_layer);
	    std::cout << "selection = " << selection << std::endl;
	    switch (selection) {
	    case 0:
	      IsGood = true;
	      break;
	    case 1:
	      IsGood = prepreparedRobbie;
	      break;
	    case 2:
	      IsGood = prepreparedHiroki;
	      break;
	    case 3:
	      if (prepreparedRobbie && prepreparedHiroki) {
		IsGood = true;
	      }
	      break;
	    default:
	      std::cerr << "Problem with 'selection' variable." << std::endl;
	      break;
	    }
	    break;
	  default:
	    std::cerr << "Problem with selection_choice." << std::endl;
	    break;
	  }
#endif
	  // Check if event is accepted by chosen selection. If not, save any information required, and move on to the next event.
	  if( !IsGood ){
#ifdef SQS
	    prepreparedtree->Fill();
#endif
	    if (DB) std::cout << "Event number " << event << "was rejected." << std::endl;
	    continue;
	  }
	  
	  ////////////////////////////////////////////////////////////
	  // 4Ev) Perform Clustering
	  // To be added by Hiroki/Qasim at a later date.
	  
	  ////////////////////////////////////////////////////////////
	  // 4Evi) Fill Trees
	  
	  outtree->Fill();
#ifdef SQS 
	  prepreparedtree->Fill();
#endif
	  if( isFlushTree ) {
	    outtree->AutoSave("FlushBaskets");
#ifdef SQS
	    prepreparedtree->AutoSave("FlushBaskets");
#endif
	    isFlushTree = kFALSE;
	  }

	  ////////////////////////////////////////////////////////////
	  // 4Evii) Delete any necessary variables, and end loop over events
	  delete currentEvent;
	}
    }

  //////////////////////////////////////////////////////////
  // 5) Write the tree(s) to the disk
  
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  

    
}
