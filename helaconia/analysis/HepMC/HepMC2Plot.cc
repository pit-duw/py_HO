// -------------------------------------------------------
// Hua-Sheng Shao April 2015
// 1) Find all of the .hep files in the current directory
// 2) Read information  from HepMC file
// 3) Perform analysis from analysis file
// -------------------------------------------------------

#include <iostream>
//#include "HepMC/PythiaWrapper.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Units.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "fstream"
#include <dirent.h>
#include <errno.h>
#include <string>
#include <sys/types.h>

/*function... might want it in some class?*/
int getdir (std::string dir, std::vector<std::string> &files)
{
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
  }

  while ((dirp = readdir(dp)) != NULL) {
    files.push_back(std::string(dirp->d_name));
  }
  closedir(dp);
  return 0;
}

extern "C" {
  extern struct {
    double EVWGT;
  } cevwgt_;
}
#define cevwgt cevwgt_

// Following are the subroutines defined in the analysis/HepMC
extern "C" {
  void __plot_hepmc_user_MOD_plot_hepmc_begin(int&, char(*)[15]);
  void __plot_hepmc_user_MOD_plot_hepmc_end(double&);
  void __plot_hepmc_user_MOD_plot_hepmc_fill(double(*));
}

int main() {
  
  int cwgtinfo_nn=1;
  char cwgtinfo_weights_info[250][15];
  double cwgt_ww[250];
  std::vector<std::string> files;
  double crosssection;
  
  getdir(".",files);
  std::string cwgtinfo ("central value");
  int stringsize=cwgtinfo.size();
  stringsize=std::min(stringsize,15);
  for (int ii = 0; ii < stringsize; ++ii){
    cwgtinfo_weights_info[0][ii]=cwgtinfo.at(ii);
  };

  // Loop over .hep files
  for ( int f = 0; f < files.size(); f++ ){
    if ( files.at(f).find(".hep") == std::string::npos ) continue;
    HepMC::IO_BaseClass *_hepevtio;
    HepMC::IO_GenEvent ascii_io( files.at(f), std::ios::in );
    std::cout << "Analysing file: "<<files.at(f)<<std::endl;
    __plot_hepmc_user_MOD_plot_hepmc_begin(cwgtinfo_nn,cwgtinfo_weights_info);
    int eventCount = 0;
    
   // Loop over events
    while ( true ){
      // Read event and set the units used
      HepMC::GenEvent* genEvt = ascii_io.read_next_event();
      // If at the end of the file, break event loop
      if ( !genEvt )break;
      // Get the units for this file - default is GeV
      //if ( !eventCount ){
	//  if ( genEvt->momentum_unit () == HepMC::Units::MEV ){
      //	 ptcut = 1e3*ParticlePtCut;
      //  }
      //  else{
      //	 ptcut = ParticlePtCut;
      //  };
      //};
      // Get the total cross section in unit of nb
      if ( !eventCount ){
	crosssection = genEvt->cross_section()->cross_section();// in unit of pb
	crosssection = 1e-3*crosssection; // in unit of nb
      };
      // Get event weight
      HepMC::WeightContainer eventWeight ( genEvt->weights () );
      double weight = *(eventWeight.begin());

      if (eventWeight.size() != 1) std::cerr << "Event Weight vector is greater than one in size !"<< std::endl;
      
      // define the IO_HEPEVT
      // Then we can use HEPMCf90.inc
      _hepevtio = new HepMC::IO_HEPEVT;
      _hepevtio->write_event(genEvt);
      // event weight
      cevwgt.EVWGT=weight;
      cwgt_ww[0]=cevwgt.EVWGT;
      
      // fill in the FORTRAN analysis routine for this event
      __plot_hepmc_user_MOD_plot_hepmc_fill(cwgt_ww);
      
      // Loop over all particles in the event
      //HepMC::GenEvent::particle_const_iterator pitr;
      //for (pitr = genEvt->particles_begin(); pitr != genEvt->particles_end(); ++pitr ){
      //  const HepMC::GenParticle* part = (*pitr);
      // Only include the stable final states ??
      //if ( part->status() != 1) continue;
      //  const HepMC::FourVector partMom = part->momentum();
      //  int pdgId = part->pdg_id();
      // Only store particles with pT above ParticlePtCut
      //  if ( partMom.perp() > ptcut ){
      // 	 n_particles++;
      //  };
      //};
      // Clean up
      delete genEvt;
      eventCount++;

    }; // end of Loop over event
    
    double norm=1.0;
    // in unit of nb
    norm=crosssection/double(eventCount);
    // the weight of each event is normalized to norm
    __plot_hepmc_user_MOD_plot_hepmc_end(norm);
    
  }; // end of Loop over files

  return 0;

}
