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

#define nchip 48
#define nrepeat_chk 30
#define clust_neigh_pix_lim 2

TString file_name;
TString file_out_path;

const std::map< Int_t, Int_t > lane2chipid_lut = { 
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

void mTower_clustering( const char* filename = "conv_Run_0058.root" ) 
	{
	  if (!(gInterpreter->IsLoaded("vector")))
    	gInterpreter->ProcessLine("#include <vector>");
  	gSystem->Exec("rm -f AutoDict*vector*vector*int*");
	gInterpreter->GenerateDictionary("vector<vector<int> >","vector");
	bool check_repeatetion= false;
	bool chip_0_hit=false;
	bool chip_1_hit=false;
	bool chip_0_cond=false;
	bool chip_1_cond=false;
	TFile*          filein      ;
	TTree*          tree        ;
	Int_t           nlayer      ;
	Int_t           runNumber   ;
	Int_t           fileNumber  ;
	Int_t           eventNumber ;
	Int_t           nHits       ;
	Int_t           dataSize    ;
	vector<Int_t>*  vlane       ;
	vector<Int_t>*  vcolumn     ;
	vector<Int_t>*  vrow        ;
	TH1F *clust_hist[nchip];
	TH1F *nHits_hist[nchip];
	TH1F *ppe_hist[nchip];
		
	TH1F *clust_hist_filtered[nchip];
	TH1F *nHits_hist_filtered[nchip];
	TH1F *ppe_hist_filtered[nchip];

	file_out_path = filename;
	int namelen = strlen(filename);
	file_out_path.Resize(namelen - 5);

	mkdir(file_out_path,S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH);
	std::cout << "current directory " << file_out_path << std::endl;
	TFile f_out_raw (Form("%s/half_mTower_test_filter.root", file_out_path.Data()), "RECREATE");

	TTree *hit_tree=new TTree("hit_tree", "hit tree");

	std::vector< std::vector<pair<int,int> > > address;
	std::vector<std::vector<int> > clust_size;
	std::vector<std::vector<int> > clust_no;
	std::vector<std::vector<int> > part_per_evt;
	std::vector<std::vector<int> > ind_nHits;
	std::vector<int> chip_w_hits;
	for (int i=0; i<nchip; i++)
		{
			address.push_back(vector<pair<int,int> >());
			clust_size.push_back(std::vector<int>());
			clust_no.push_back(std::vector<int>());
			part_per_evt.push_back(std::vector<int>());
			ind_nHits.push_back(std::vector<int>());
			clust_hist[i] = new TH1F(Form("clust_Chip_%d",i), Form("clust_Chip_%d",i), 100, 0, 100);
			nHits_hist[i] = new TH1F(Form("nHits_Chip_%d",i), Form("nHits_Chip_%d",i), 500, 0, 500);
			ppe_hist[i] = new TH1F(Form("ppe_Chip_%d",i), Form("ppe_Chip_%d",i), 100, 0, 100);

			clust_hist_filtered[i] = new TH1F(Form("clust_Chip_filtered_%d",i), Form("clust_Chip_filtered_%d",i), 100, 0, 100);
			nHits_hist_filtered[i] = new TH1F(Form("nHits_Chip_filtered_%d",i), Form("nHits_Chip_filtered_%d",i), 500, 0, 500);
			ppe_hist_filtered[i] = new TH1F(Form("ppe_Chip_filtered_%d",i), Form("ppe_Chip_filtered_%d",i), 100, 0, 100);
	
		}
	hit_tree->Branch("ind_nHits", &ind_nHits);
	hit_tree->Branch("clust_size", &clust_size);
	hit_tree->Branch("part_per_evt", &part_per_evt);

	Int_t cur_evt_ch_[nchip];
	for (int k=0; k <nchip; k++)
		{cur_evt_ch_[k]=-1;} 

	nlayer  = nchip/2;
	filein  = new TFile( filename, "read" );
	tree    = (TTree*)filein->Get("Frames");
	vlane   = new vector<Int_t>();
	vcolumn = new vector<Int_t>();
	vrow    = new vector<Int_t>();
	tree->SetBranchAddress("runNumber"   , &runNumber   );
	tree->SetBranchAddress("fileNumber"  , &fileNumber  );
	tree->SetBranchAddress("eventNumber" , &eventNumber );
	tree->SetBranchAddress("nHits"       , &nHits       );
	tree->SetBranchAddress("dataSize"    , &dataSize    );
	tree->SetBranchAddress("lane"        , &vlane       );
	tree->SetBranchAddress("column"      , &vcolumn     );
	tree->SetBranchAddress("row"         , &vrow        );
	cout<<"Processing starts "  <<endl;

  	Int_t ev_start, ev_end;
    	ev_start = 0 ;
    	ev_end   = tree->GetEntries() ;
	//ev_end   = 1000;

  	for(Int_t iev = ev_start; iev < ev_end; iev++ ){ 
    		tree->GetEntry( iev ); 
		cout<<"Processing Event no = " <<iev <<endl;
    		for (Int_t ihit = 0; ihit < nHits; ihit++){
	      		Int_t lane = vlane  ->at(ihit) ;
			Int_t col  = vcolumn->at(ihit) ;
			Int_t row  = vrow   ->at(ihit) ;
			Int_t chip = lane2chipid_lut.at( lane ); 
			if( chip < 0 ) continue;
		//cout<<"Event = " << iev <<" nHits = " << nHits <<" Chip Id = " << chip << " Col = " << col << " row = " << row  <<endl;
			int n_entries=-1;
			int clust_index=-1;
			 if(cur_evt_ch_[chip]!=iev)
				{
				if(chip==0)	{chip_0_hit=true;}
				if(chip==1)	{chip_1_hit=true;}
				chip_w_hits.push_back(chip);
				ind_nHits[chip].push_back(1);
				clust_size[chip].push_back(1);
				part_per_evt[chip].push_back(1);
				address[chip].push_back( make_pair(col,row) );
				clust_no[chip].push_back(clust_size[chip].size()-1);
				cur_evt_ch_[chip]=iev;
				}
			else if (cur_evt_ch_[chip]==iev)
				{
				address[chip].push_back( make_pair(col,row));
				n_entries=ind_nHits[chip][ind_nHits[chip].size()-1];
				ind_nHits[chip][ind_nHits[chip].size()-1]+=1;
				clust_index=find_neigh(chip,n_entries,address);
				//clust_index=find_neigh(chip,n_entries,address);
				if(clust_index!=-1)
					{
					clust_size[chip][clust_no[chip][clust_index]]+=1;
					clust_no[chip].push_back(clust_no[chip][clust_index]);
					}
				else
					{clust_size[chip].push_back(1);
					clust_no[chip].push_back(clust_size[chip].size()-1);
					part_per_evt[chip][part_per_evt[chip].size()-1]+=1;
					}
				}
		}//iHits for loop ends here
	if ((chip_0_hit)^(chip_1_hit)){
	//if ((chip_0_hit)||(chip_1_hit)){
		if(part_per_evt[0].size()>0)
	 		{chip_0_cond=((part_per_evt[0][part_per_evt[0].size()-1]==1)&&(chip_0_hit));}
		if(part_per_evt[1].size()>0)
			{chip_1_cond=((part_per_evt[1][part_per_evt[1].size()-1]==1)&&(chip_1_hit));}

		if(chip_0_cond||chip_1_cond)

			{//cout<<"first layer chip 1 detects a particle"<<endl;
			for (int p=0;p<chip_w_hits.size();p++){
			nHits_hist_filtered[chip_w_hits[p]]->Fill(ind_nHits[chip_w_hits[p]][ind_nHits[chip_w_hits[p]].size()-1]);
			ppe_hist_filtered[chip_w_hits[p]]->Fill(part_per_evt[chip_w_hits[p]][part_per_evt[chip_w_hits[p]].size()-1]);

			int ppe_last_evt= part_per_evt[chip_w_hits[p]][part_per_evt[chip_w_hits[p]].size()-1];
			int last_clust_address= clust_size[chip_w_hits[p]].size()-1; 
			for(int z=(last_clust_address-ppe_last_evt);z<last_clust_address;z++)
				{
				//cout<<"cluster size is "<<clust_size[chip_w_hits[p]][z]<<"in chip "<<chip_w_hits[p]<<endl;
				clust_hist_filtered[chip_w_hits[p]]->Fill(clust_size[chip_w_hits[p]][z]);}			
				}
			}
		}
	chip_0_hit=false;
	chip_1_hit=false;
	chip_w_hits.clear();
	}//event loop ends here	

	for (int k=0; k<nchip;k++)
		{ cout <<"number of nHits in chip "<<k<<" = "<< address[k].size()<<endl;
		for (int l=0; l<clust_size[k].size();l++)
			{
			clust_hist[k]->Fill(clust_size[k][l]);	
			}
		for (int l=0; l<ind_nHits[k].size();l++)
			{
			nHits_hist[k]->Fill(ind_nHits[k][l]);			
			}
		for (int l=0; l<part_per_evt[k].size();l++)
			{
			ppe_hist[k]->Fill(part_per_evt[k][l]);			
			}
		}
	hit_tree->Fill();
	f_out_raw.Write();

	ind_nHits.clear();
	clust_size.clear();
	clust_no.clear();
	address.clear();
	part_per_evt.clear();
}

