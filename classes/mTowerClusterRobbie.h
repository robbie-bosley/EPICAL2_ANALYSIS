#ifndef MTOWERCLUSTERROBBIE_H
#define MTOWERCLUSTERROBBIE_H

#include "TObjArray.h"
#include "mTowerHit.h"

//Class cluster for the mTower

class mTowerClusterRobbie : public TObject
{

 protected:
  int id; //cluster id number
  int lane; //lane of the chip
  TObjArray* hits; //of type mTowerHit
  

 public:
  mTowerClusterRobbie();
  mTowerClusterRobbie(int i);
  virtual ~mTowerClusterRobbie();

  int getNHits() {return hits->GetEntries();}
  TObjArray* getHits() {return hits;}
  void AddHit(mTowerHit* hit);
  void AddCluster(mTowerClusterRobbie* cluster);
  void setId(int i) {id = i;}
  int getId() {return id;}
  void setLane(int l) {lane = l;}
  int getLane() {return lane;}
  double getMeanRow();
  double getMeanColumn();

  ClassDef(mTowerClusterRobbie,1) //mTower cluster
    
};

#endif

