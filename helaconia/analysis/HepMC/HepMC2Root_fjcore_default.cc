// -------------------------------------------------------
// Original by Jim Henderson (ATLAS Collaboration) March 2015
// Changed by Hua-Sheng Shao April 2015
// Convert HepMC to simple ROOT tree
// Code finds ALL ".hep" files in current directory and loops through
//--------------------------------------------------------

#include "HepMC2Root_fjcore_default.h"
#include <iostream>
//#include "HepMC/PythiaWrapper.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Units.h"
//#include "HepMC/HEPEVT_Wrapper.h"
#include <TFile.h>

using namespace std;

int main(){

  TFile outFile( "helaconia_hepmc.root", "RECREATE" );
  TTree *output_tree = new TTree("helaconia", "helaconia");
  SetupTree( output_tree );
  getdir(".", files );
  
  for ( int f = 0; f < files.size(); f++ ){
    if ( files.at(f).find(".hep") == std::string::npos ) continue;
    HepMC::IO_GenEvent ascii_io2( files.at(f), std::ios::in );
    cout << "Analysing file: " << files.at(f) << endl;
    int eventCount = 0;
    double weight_tot = 0.; // in some case the weight will be negative
                            // for example in aMC@NLO
    // for calculating the total weight of the events
    cout << "Loop over events for total weight ..."<<endl;
    while ( true ){
      // Read event
      HepMC::GenEvent* genEvt = ascii_io2.read_next_event();
      // If at end of file, break event loop
      if ( !genEvt ) break;
      // Get event weight
      HepMC::WeightContainer eventWeight ( genEvt->weights () );
      weight = *(eventWeight.begin());
      weight_tot = weight_tot+weight;
    };
    cout << "                           ... done  "<<endl;
    cout << "Loop over events for filling ROOT TTree ..."<<endl;
    // reopen the file
    HepMC::IO_GenEvent ascii_io( files.at(f), std::ios::in );
    // Loop over events
    while ( true ){
      // Read event and set the units used
      HepMC::GenEvent* genEvt = ascii_io.read_next_event();
      // If at end of file, break event loop
      if ( !genEvt ) break;
      // Get the units for this file - default is GeV
      // Change it in HepMC2Root_fjcore.h
      if ( !eventCount ){
      	if ( genEvt->momentum_unit () == HepMC::Units::MEV ){
      	  JetPtCut = 1e3*JetPtCut;
      	  ParticlePtCut = 1e3*ParticlePtCut;
      	}
      	else{
      	  JetPtCut =JetPtCut; 
      	  ParticlePtCut = ParticlePtCut;
      	};
	cross_section= genEvt->cross_section()->cross_section();// in unit of pb ;
	cross_section= 1e-3*cross_section; // in unit of nb
      }
      // Reset output variables
      ResetTreeVars();
      // Get event weight
      HepMC::WeightContainer eventWeight ( genEvt->weights () );
      weight = *(eventWeight.begin());
      if ( weight_tot > 0 ){
	evt_weight = cross_section*weight/weight_tot;
      } else {
	evt_weight = 0.;
      };
      if ( eventWeight.size() != 1 ) cerr << "Event Weight vector is greater than one in size, bit weird."  << endl;
      // Loop over all particles in the event
      HepMC::GenEvent::particle_const_iterator pitr;
      for (pitr = genEvt->particles_begin(); pitr != genEvt->particles_end(); ++pitr ) {
	const HepMC::GenParticle* part = (*pitr);
	if ( part->status() != 1 ) continue;
	const HepMC::FourVector partMom = part->momentum();
	int pdgId = part->pdg_id();
	if ( fabs(pdgId) == 12 || fabs(pdgId) == 14 || fabs(pdgId) == 16 ){
	  // Add ALL neutrinos to missing eT
	  MET_p += fjcore::PseudoJet( partMom.px(), partMom.py(), partMom.pz(), partMom.e() );
	}
	else{
	  // Fill all OTHER final state particles into a vector of (pseudo)jets for fastjet to cluster later
	  jetParticles.push_back( fjcore::PseudoJet( partMom.px(), partMom.py(), partMom.pz(), partMom.e() ) );
	}
	// Only store particles with pT above ParticlePtCut
	// Change it in HepMC2Root_fjcore.h
	if ( partMom.perp() > ParticlePtCut ){
	  n_particles++;
	  PID_v.push_back( pdgId );
	  P_X_v.push_back( partMom.px() );
	  P_Y_v.push_back( partMom.py() );
	  P_Z_v.push_back( partMom.pz() );
	  P_T_v.push_back( partMom.perp() );
	  E_v.push_back( partMom.e() );
	  M_v.push_back( partMom.m() );
	  Eta_v.push_back( partMom.eta() );
	  Phi_v.push_back( partMom.phi() );
	}
      }
      // Cluster jets and sort by pT
      // Change jet_def in HepMC2Root_fjcore.h
      fjcore::ClusterSequence cs(jetParticles, jet_def);
      jets = fjcore::sorted_by_pt(cs.inclusive_jets());
      for ( int j = 0; j < jets.size(); j++ ){
	// Only store jets with pT above JetPtCut
	// Change JetPtCut in HepMC2Root_fjcore.h
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
      delete genEvt;
      eventCount++;
    };
    cout<<"... Done. A ROOT TTree is filled in !"<<endl;
  }
  // Write out and save
  outFile.cd();
  output_tree->Write();
  outFile.Close();
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
