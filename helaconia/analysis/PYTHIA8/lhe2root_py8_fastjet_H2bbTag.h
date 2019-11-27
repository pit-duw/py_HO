#ifndef HEPMC2ROOT_H
#define HEPMC2ROOT_H

#include <cmath>
#include <math.h>
#include <TTree.h>
#include <errno.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <limits>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
#include "HEPTopTagger.hh"

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

fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 1.2 );
vector<fastjet::PseudoJet> jetParticles;
vector<fastjet::PseudoJet> jets;
fastjet::PseudoJet MET_p;
fastjet::MassDropTagger htagger(0.67, 0.15 );
double hmasslow=0.;
double hmassup=1e30;

double weight;
double cross_section; // in unit of nb
double evt_weight; // weight for each event in unit of nb  
int Jet_n = 0;
int n_particles = 0;
double MET_et = 0;
double MET_phi = 0;
// using massdrop for higgs tag (without b tagging)
double htag_mu = 0;
double htag_y = 0;
double htag_mass = 0;
double htag_mass_afterfilter = 0;
// using massdrop for higgs tag (with two b tagging)
double hbbtag_mu = 0;
double hbbtag_y = 0;
double hbbtag_mass = 0;
double hbbtag_mass_afterfilter = 0;

double JetPtCut=200;
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
  tree->Branch("htag_mu", &(htag_mu));
  tree->Branch("htag_y",&(htag_y));
  tree->Branch("htag_mass",&(htag_mass));
  tree->Branch("htag_mass_afterfilter",&(htag_mass_afterfilter));
  tree->Branch("hbbtag_mu", &(hbbtag_mu));
  tree->Branch("hbbtag_y",&(hbbtag_y));
  tree->Branch("hbbtag_mass",&(hbbtag_mass));
  tree->Branch("hbbtag_mass_afterfilter",&(hbbtag_mass_afterfilter));
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

vector<fastjet::PseudoJet> gran_jets ( vector<fastjet::PseudoJet> & ori_hadrons,const double & eta_cell, const double & phi_cell, const double & pt_cutoff)
{
  double pi = 3.142592654;
  vector<fastjet::PseudoJet> granulated_jets;
  granulated_jets.clear();

  ori_hadrons = sorted_by_pt(ori_hadrons);
  granulated_jets.push_back(ori_hadrons[0]);
  for (unsigned i= 1; i < ori_hadrons.size(); i++)
    {
      int new_jet= 0;
      for (unsigned j=0; j < granulated_jets.size(); j++)
	{
	  double eta_cell_diff = abs(ori_hadrons[i].pseudorapidity() -
				     granulated_jets[j].pseudorapidity())/eta_cell;
	  double phi_cell_diff = abs(ori_hadrons[i].phi() -
				     granulated_jets[j].phi());
	  if(phi_cell_diff > pi)  phi_cell_diff = 2*pi - phi_cell_diff;
	  phi_cell_diff = phi_cell_diff/phi_cell;

	  if( eta_cell_diff < 1 && phi_cell_diff < 1)
	    {
	      new_jet = 1;

	      double total_energy  = ori_hadrons[i].e() + granulated_jets[j].e();
	      double rescale_factor = sqrt( pow(ori_hadrons[i].px()+granulated_jets[j].px(),2) +
					    pow(ori_hadrons[i].py()+granulated_jets[j].py(),2) +
					    pow(ori_hadrons[i].pz()+granulated_jets[j].pz(),2) );
	      double rescaled_px = total_energy*(ori_hadrons[i].px()+granulated_jets[j].px()) /rescale_factor ;
	      double rescaled_py = total_energy*(ori_hadrons[i].py()+granulated_jets[j].py()) /rescale_factor ;
	      double rescaled_pz = total_energy*(ori_hadrons[i].pz()+granulated_jets[j].pz()) /rescale_factor ;

	      fastjet::PseudoJet comb_jet( rescaled_px, rescaled_py,rescaled_pz, total_energy);
	      comb_jet.set_user_index(ori_hadrons[i].user_index()+granulated_jets[j].user_index());
	      granulated_jets.erase(granulated_jets.begin()+j);
	      granulated_jets.push_back(comb_jet);
	      break;
	    }
	}
      if(new_jet != 1)
	{
	  granulated_jets.push_back(ori_hadrons[i]);
	  granulated_jets = sorted_by_pt(granulated_jets);
	}
    }

  for(unsigned ii=0; ii <granulated_jets.size();ii++)
    {

      if((granulated_jets[ii].perp() < pt_cutoff))
	{
	  granulated_jets.erase(granulated_jets.begin()+ii);
	  ii--;
	}
    }

  return(granulated_jets);

}

#endif
