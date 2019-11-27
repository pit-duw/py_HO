#include "lhe2root_py8_fjcore_user.h"
#include <iostream>
#include "Pythia8/Pythia.h"
#include "fstream"
#include <TFile.h>
#include <math.h>
#include "sstream"

#define PI 3.14159265358979323846264338328

using namespace Pythia8;
using namespace std;

double GetPseudoRapidity(double Px, double Py, double Pz);
double GetPt(double Px, double Py);
double GetPhi(double Px, double Py, double Pz);
double GetMomentum(double Px, double Py, double Pz);
double GetInvMass(double E1, double Px1, double Py1, double Pz1,
		  double E2, double Px2, double Py2, double Pz2);

double muonPionMisidentification(float mom);
double muonKaonMisidentification(float mom);

int main() {
  TFile outFile( "helaconia_hepmc.root", "RECREATE" );
  TTree *output_tree = new TTree("helaconia", "helaconia");
  SetupTree( output_tree );
  
  
  Pythia pythia;

  Event& event = pythia.event;

  string inputname="Pythia8_lhe.cmnd";

  pythia.readFile(inputname.c_str());

  pythia.init();
  string filename = pythia.word("Beams:LHEF");
  // read the total cross section information
  ifstream lhefile(filename.c_str());
  string str;
  cross_section=1.;
  while (getline(lhefile,str)){
    if (str.find("<init>") == std::string::npos ) continue;
    // find the init block
    double temp;
    for ( int i = 0; i < 10; ++i ){
      lhefile>>temp;
    };
    lhefile>>cross_section;
    break;
  };
  lhefile.close();
  //from unit of pb to nb
  cross_section = cross_section*1e-3;

  int nAbort=10;
  int nPrintLHA=1;
  int iAbort=0;
  int iPrintLHA=0;
  int iEventtot=pythia.mode("Main:numberOfEvents");
  int iEventshower=pythia.mode("Main:spareMode1");
  

  double nSelected;
  double norm=1.;

  // Cross section
  double sigmaTotal  = 0.;

  for (int iEvent = 0; ; ++iEvent) {
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      if (++iAbort < nAbort) continue;
      break;
    };
    // the number of events read by Pythia so far
    nSelected=double(pythia.info.nSelected());

    if (nSelected >= iEventshower) break;
    if (pythia.info.isLHA() && iPrintLHA < nPrintLHA) {
      pythia.LHAeventList();
      pythia.info.list();
      pythia.process.list();
      pythia.event.list();
      ++iPrintLHA;
    };

    // Reset output variables 
    ResetTreeVars();

    double evtweight = pythia.info.weight();
    double normhepmc;
    // Add the weight of the current event to the cross section.
    normhepmc = 1. / double(iEventshower);
    sigmaTotal += evtweight*normhepmc;
    weight = evtweight;
    // assume there is no negative weight
    evt_weight = normhepmc*cross_section;

    // Loop over all particles in the event
    for (int i=1; i < event.size(); ++i){
      int pdgId = event[i].id();
      if ( fabs(pdgId) == 12 || fabs(pdgId) == 14 || fabs(pdgId) == 16 ){
	// Add ALL neutrinos to missing eT
	MET_p += fjcore::PseudoJet( event[i].px(), event[i].py(), event[i].pz(), event[i].e() );
      } else {
	// Fill all OTHER final state particles into a vector of (pseudo)jets for fastjet to cluster later
	jetParticles.push_back( fjcore::PseudoJet( event[i].px(), event[i].py(), event[i].pz(), event[i].e() ) );
      }
      double perp=GetPt(event[i].px(),event[i].py());
      // Only store particles with pT above ParticlePtCut 
      // Change ParticlePtCut in the head file
      if ( perp > ParticlePtCut ){
	double eta=GetPseudoRapidity(event[i].px(),event[i].py(),event[i].pz());
	double phi=GetPhi(event[i].px(),event[i].py(),event[i].pz());
	n_particles++;
	PID_v.push_back( pdgId );
	P_X_v.push_back( event[i].px() );
	P_Y_v.push_back( event[i].py() );
	P_Z_v.push_back( event[i].pz() );
	P_T_v.push_back( perp );
	E_v.push_back( event[i].e() );
	M_v.push_back( event[i].m() );
	Eta_v.push_back( eta );
	Phi_v.push_back( phi );
      };
    } // end of loop over particles
    // Cluster jets and sort by pT 
    // Change the jet_def in the head file
    fjcore::ClusterSequence cs(jetParticles, jet_def);
    jets = fjcore::sorted_by_pt(cs.inclusive_jets());
    for ( int j = 0; j < jets.size(); j++ ){
      // Only store jets with pT above JetPtCut 
      // Change the JetPtCut in the head file
      if ( jets.at(j).perp() < JetPtCut ) continue;
      Jet_n++;
      Jet_pt_v.push_back( jets.at(j).perp() );
      Jet_eta_v.push_back( jets.at(j).eta() );
      Jet_phi_v.push_back( jets.at(j).phi() );
      Jet_E_v.push_back( jets.at(j).E() );
    }
    // Fill output MET variables 
    MET_et = MET_p.Et();
    MET_phi = MET_p.phi();
    // Fill tree
    output_tree->Fill();
    // Clean up
  }; // end of loop over events
  // Write out and save
  outFile.cd();
  output_tree->Write();
  outFile.Close();

  //delete output_tree;

  pythia.stat();

  // sigmaGen is meaningful only AFTER stat
  // in unit of mb
  //norm=pythia.info.sigmaGen();
  // the weight of each event is normalized to norm
  //norm=norm/double(iEventshower);
  // change unit from mb to nb
  //norm=norm*1e6;

  return 0;
}

void ResetTreeVars(){
  n_particles = 0;
  MET_et = 0;
  MET_phi = 0;
  Jet_n = 0;
  PID_v.clear();
  P_X_v.clear();
  P_Y_v.clear();
  P_Z_v.clear();
  P_T_v.clear();
  E_v.clear();
  M_v.clear();
  Eta_v.clear();
  Phi_v.clear();
  Jet_E_v.clear();
  Jet_pt_v.clear();
  Jet_eta_v.clear();
  Jet_phi_v.clear();
  jetParticles.clear();
  jets.clear();
  MET_p.reset(0,0,0,0);
}

double GetPseudoRapidity(double Px, double Py, double Pz)
{
  /////////////////////////////////////////////////////////////////////
  // Return the pseudo-rapidity for given momenta
  /////////////////////////////////////////////////////////////////////
  double pAbs = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((pAbs+Pz)/(pAbs-Pz));
}

double GetPt(double Px, double Py)
{
  /////////////////////////////////////////////////////////////////////
  // Return the transverse momentum of the particle:
  /////////////////////////////////////////////////////////////////////
  return sqrt(Px*Px + Py*Py);
}

double GetPhi(double Px, double Py, double Pz)
{
  if ( Px*Px+Py*Py < 0. ){
    if ( Pz >= 0. )return 0.;
    if ( Pz < 0. )return PI;
  };
  double s2= Px*Px/(Px*Px+Py*Py);
  double s=sqrt(s2);
  double res;
  if ( s > 1. ){
    cout<<"GetPhi WARNING S="<<s<<endl;
    res=0.;
    if ( Px < 0. )res=PI;
    return res;
  };
  if ( Px < 0. )s=-s;
  res=acos(s);
  if ( Py < 0. )res=2*PI-res;
  return res;
}

double GetMomentum(double Px, double Py, double Pz)
{
  /////////////////////////////////////////////////////////////////////
  // Return P given Px, Py, Pz
  /////////////////////////////////////////////////////////////////////
  return sqrt(Px*Px + Py*Py + Pz*Pz);
}

double GetInvMass(double E1, double Px1, double Py1, double Pz1,
		  double E2, double Px2, double Py2, double Pz2)
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Return Invariant Mass given E1, Px1, Py1, Pz1, E2, Px2, Py2, Pz2
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  return sqrt( (E1+E2)*(E1+E2) -
               (Px1+Px2)*(Px1+Px2) -
               (Py1+Py2)*(Py1+Py2) -
               (Pz1+Pz2)*(Pz1+Pz2));
}

//_____________ probability (in %) of kaon misidentification with muon
double muonKaonMisidentification(float mom){

  float a = 8.60;
  float b = 0.11;

  return 0.5 + a*exp(-b*mom);


}

//_____________ probability (in %) of pion misidentification with muon
double muonPionMisidentification(float mom){

  float a = 6.63;
  float b = 0.13;

  return 0.5 + a*exp(-b*mom);


}
