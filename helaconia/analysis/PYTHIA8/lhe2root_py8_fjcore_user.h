#ifndef HEPMC2ROOT_H
#define HEPMC2ROOT_H

#include <TTree.h>
#include <errno.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <limits>
#include "fjcore.hh"

using namespace std;

vector<int>     PID_v;
vector<double>  P_X_v;
vector<double>  P_Y_v;
vector<double>  P_Z_v;
vector<double>  P_T_v;
vector<double>  E_v;
vector<double>  M_v;
vector<double>  Eta_v;
vector<double>  Phi_v;

vector<double>  Jet_pt_v;
vector<double>  Jet_E_v;
vector<double>  Jet_eta_v;
vector<double>  Jet_phi_v;

fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, 0.4 );
vector<fjcore::PseudoJet> jetParticles;
vector<fjcore::PseudoJet> jets;
fjcore::PseudoJet MET_p;

double weight;
double cross_section; // in unit of nb
double evt_weight; // weight for each event in unit of nb  
int Jet_n = 0;
int n_particles = 0;
double MET_et = 0;
double MET_phi = 0;

double JetPtCut=20;
double ParticlePtCut=0;

std::vector<std::string> files;

void ResetTreeVars();

void SetupTree( TTree* tree ){
  tree->Branch("n_particles", &(n_particles));
  tree->Branch("PID", &(PID_v) );
  tree->Branch("P_X", &(P_X_v) );
  tree->Branch("P_Y", &(P_Y_v) );
  tree->Branch("P_Z", &(P_Z_v) );
  tree->Branch("P_T", &(P_T_v) );
  tree->Branch("E", &(E_v) );
  tree->Branch("M", &(M_v) );
  tree->Branch("Eta", &(Eta_v) );
  tree->Branch("Phi", &(Phi_v) );
  tree->Branch("Jet_n", &(Jet_n));
  tree->Branch("Jet_E", &(Jet_E_v) );
  tree->Branch("Jet_pt", &(Jet_pt_v) );
  tree->Branch("Jet_eta", &(Jet_eta_v) );
  tree->Branch("Jet_phi", &(Jet_phi_v) );
  tree->Branch("MET_et", &(MET_et));
  tree->Branch("MET_phi", &(MET_phi));
  tree->Branch("weight", &(weight));
  tree->Branch("cross_section",&(cross_section));
  tree->Branch("evt_weight",&(evt_weight));
}


int getdir (std::string dir, std::vector<std::string> &files)
{
  DIR *dp;
  struct dirent *dirp;
  if((dp = opendir(dir.c_str())) == NULL) {
    cout << "Error(" << errno << ") opening " << dir << endl;
    return errno;
  }

  while ((dirp = readdir(dp)) != NULL) {
    files.push_back(string(dirp->d_name));
  }
  closedir(dp);
  sort( files.begin(), files.end() );
  return 0;
}



#endif
