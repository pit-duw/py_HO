// Driver for using Pythia8 with lhe file with hepmc

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_BaseClass.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "fstream"
#include "LHEFRead.h"

using namespace Pythia8;

extern "C" {
  extern struct {
    double EVWGT;
  } cevwgt_;
}
#define cevwgt cevwgt_

extern "C" {
  void pyabeg_(int&,char(*)[15]); // pythia8 analysis begin routine
  void pyaend_(int&);             // pythia8 analysis end routine
  void pyanal_(int&,double(*));   // pythia8 analysis routine
}

//==========================================================================

int main() {
  Pythia pythia;
  int cwgtinfo_nn;
  char cwgtinfo_weights_info[250][15];
  double cwgt_ww[250];

  string inputname="Pythia8_lhe.cmnd",outputname="Pythia8_lhe.hep";
  
  pythia.readFile(inputname.c_str());
  pythia.init();
  string filename = pythia.word("Beams:LHEF");
  
  MyReader read(filename);
  read.lhef_read_wgtsinfo_(cwgtinfo_nn,cwgtinfo_weights_info);
  pyabeg_(cwgtinfo_nn,cwgtinfo_weights_info);

  int nAbort=10;

  pythia.start();
  return 0;
}
  
