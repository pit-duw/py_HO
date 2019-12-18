// Driver for Pythia 8. Reads an input file dynamically created on
// the basis of the inputs specified in MCatNLO_MadFKS_PY8.Script 
#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_BaseClass.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "fstream"
#include "LHEFRead.h"
#include <string>
#include <sstream>
//#include "../examples/CombineMatchingInput.h"

using namespace Pythia8;

extern "C" {
  extern struct {
    double EVWGT;
  } cevwgt_;
}
#define cevwgt cevwgt_

// Following are the subroutines defined in the analysis/PYTHIA8
extern "C" { 
  void __plot_mlm_py8_user_MOD_plot_mlm_py8_begin(int&,char(*)[15]);
  void __plot_mlm_py8_user_MOD_plot_mlm_py8_end(double&);
  void __plot_mlm_py8_user_MOD_plot_mlm_py8_fill(double(*));
}

std::string trim(const std::string& str,
		 const std::string& whitespace);

int main() {
  Pythia pythia;

  int cwgtinfo_nn=1;
  char cwgtinfo_weights_info[250][15];
  double cwgt_ww[250];

  string inputname="Pythia8_lhe.cmnd",outputname="Pythia8_lhe.hep";

  pythia.readFile(inputname.c_str());

  pythia.init();
  string filename = pythia.word("Beams:LHEF");

  MyReader read(filename);
  // Read header of event file
  //read.lhef_read_wgtsinfo_(cwgtinfo_nn,cwgtinfo_weights_info);
  std::string cwgtinfo ("central value");
  cwgtinfo=trim(cwgtinfo," \t");
  int stringsize=cwgtinfo.size();
  stringsize=std::min(stringsize,15);
  for (int ii = 0; ii < stringsize; ++ii){
    cwgtinfo_weights_info[0][ii]=cwgtinfo.at(ii);
  };
  if (stringsize < 15 ){
    for (int ii = stringsize; ii < 15; ++ii){
      std::string emptystr (" ");
      cwgtinfo_weights_info[0][ii]=emptystr.at(0); // avoid ^@ in Fortan code
    }
  };
 __plot_mlm_py8_user_MOD_plot_mlm_py8_begin(cwgtinfo_nn,cwgtinfo_weights_info);

  int nAbort=10;
  int nPrintLHA=1;
  int iAbort=0;
  int iPrintLHA=0;
  int nstep=5000;
  int iEventtot=pythia.mode("Main:numberOfEvents");
  int iEventshower=pythia.mode("Main:spareMode1");
  //string evt_norm=pythia.word("Main:spareWord1");
  //int iEventtot_norm=iEventtot;
  //if (evt_norm == "average"){
  //  iEventtot_norm=1;
  //}

  HepMC::IO_BaseClass *_hepevtio;
  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent ascii_io(outputname.c_str(), std::ios::out);
  double nSelected;
  double norm=1.;

  // Cross section
  double sigmaTotal  = 0.;

  for (int iEvent = 0; ; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      break;
    };
    // the number of events read by Pythia so far
    nSelected=double(pythia.info.nSelected());
    // normalisation factor for the default analyses defined in pyanal_
    //norm=iEventtot_norm*iEvent/nSelected;

    if (nSelected >= iEventshower) break;
    if (pythia.info.isLHA() && iPrintLHA < nPrintLHA) {
      pythia.LHAeventList();
      pythia.info.list();
      pythia.process.list();
      pythia.event.list();
      ++iPrintLHA;
    };

    double evtweight = pythia.info.weight();
    double normhepmc;
    // Add the weight of the current event to the cross section.
    normhepmc = 1. / double(iEventshower);
    //if (evt_norm == "average") {
    //  sigmaTotal  += evtweight*normhepmc;
    //} else {
    //  sigmaTotal  += evtweight*normhepmc*iEventtot;
    //}
    sigmaTotal += evtweight*normhepmc;

    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia, hepmcevt ); // fill the pythia event info to hepmcevt

    //define the IO_HEPEVT
    _hepevtio = new HepMC::IO_HEPEVT;
    _hepevtio->write_event(hepmcevt);
    
    //event weight
    cevwgt.EVWGT=hepmcevt->weights()[0];

    //call the FORTRAN analysis for this event
    //read.lhef_read_wgts_(cwgt_ww);
    //plot_py8_fill_(cwgtinfo_nn,cwgt_ww);
    cwgt_ww[0]=cevwgt.EVWGT;
    __plot_mlm_py8_user_MOD_plot_mlm_py8_fill(cwgt_ww);

    //if (iEvent % nstep == 0 && iEvent >= 100){
      //pyaend_(norm);
    //  plot_py8_end_(1.);
    //}
    delete hepmcevt;
  };

  pythia.stat();

  // sigmaGen is meaningful only AFTER stat
  // in unit of mb
  norm=pythia.info.sigmaGen();
  // the weight of each event is normalized to norm
  norm=norm/double(iEventshower);
  // change unit from mb to nb
  norm=norm*1e6;

  __plot_mlm_py8_user_MOD_plot_mlm_py8_end(norm);

  pythia.stat();

  return 0;
}

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
  const std::size_t strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; // no content

  const std::size_t strEnd = str.find_last_not_of(whitespace);
  const std::size_t strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}
