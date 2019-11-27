//
//  Sample test program for running EvtGen
//  

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtStdlibRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"

#ifdef EVTGEN_EXTERNAL
#include "EvtGenExternal/EvtExternalGenList.hh"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <list>

int main(int argc, char** argv) {

  EvtParticle* parent(0);
  int B_NDAUG;
  std::ofstream myfile;
  EvtStdlibRandomEngine eng;
  EvtRandom::setRandomEngine((EvtRandomEngine*)&eng);

  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;

#ifdef EVTGEN_EXTERNAL
  EvtExternalGenList genList;
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();
#endif

  //Initialize the generator - read in the decay table and particle properties
  EvtGen myGenerator("../DECAY_2010.DEC","../evt.pdl", (EvtRandomEngine*)&eng,
  		     radCorrEngine, &extraModels);

  //If I wanted a user decay file, I would read it in now.
  //myGenerator.readUDecay("../user.dec");

  static EvtId UPS4 = EvtPDL::getId(std::string("Upsilon(4S)"));

  int nEvents(100);

  // Loop to create nEvents, starting from an Upsilon(4S)
  int i;
  myfile.open ("myexample.hepmc");
  for (i = 0; i < nEvents; i++) {

    std::cout<<"Event number "<<i<<std::endl;

    // Set up the parent particle
    EvtVector4R pInit(EvtPDL::getMass(UPS4), 0.0, 0.0, 0.0);
    parent = EvtParticleFactory::particleFactory(UPS4, pInit);
    parent->setVectorSpinDensity();      

    // Generate the event
    myGenerator.generateDecay(parent);
    
    // Write out the results
    EvtHepMCEvent theEvent;
    theEvent.constructEvent(parent);
    HepMC::GenEvent* genEvent = theEvent.getEvent();
    // Please look at src/EvtGenBase/EvtHepMCEvent.cpp
    genEvent->print(myfile);

    // Get the momentum of the daughters
    B_NDAUG = parent->getNDaug();
    std::cout<<"Number of daughters:"<<B_NDAUG<<std::endl;

    parent->deleteTree();

  };
  myfile.close();

  return 0;

}
