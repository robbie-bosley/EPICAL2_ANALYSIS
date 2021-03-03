#ifndef MTOWERCLUSTERING_H
#define MTOWERCLUSTERING_H

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

long find_neigh(int chip, int n_entries, const std::vector< std::vector<pair<int,int> > > &address);

bool lane2isleft( Int_t lane);

double row2x( Int_t row, Int_t lane);

double column2y( Int_t col, Int_t lane);
