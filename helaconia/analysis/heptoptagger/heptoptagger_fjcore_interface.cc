#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <exception>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ext/hash_map>

#include "fjcore.hh"

using namespace std;
using namespace fjcore;

#define HEPTOPTAGGER_BEGIN_NAMESPACE namespace HEPTopTaggerInt {
#define HEPTOPTAGGER_END_NAMESPACE   }
#include "HEPTopTagger_fjcore.hh"
#include "output_gnuplot_display.hh"

HEPTOPTAGGER_BEGIN_NAMESPACE

// a namespace for the fortran-wrapper which contains commonly-used
// structures and means to transfer fortran <-> C++

namespace fwrappertoptagger {

  vector<PseudoJet> input_particles, jets;
  JetDefinition jet_def;
  auto_ptr<ClusterSequence> cs;

  // helper routine to transfer fortran input particles into
  void transfer_input_particles(const double * p, const int & npart) {
    input_particles.resize(0);
    input_particles.reserve(npart);
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
	mom[j] = *(p++);
      }
      PseudoJet psjet(mom);
      psjet.set_user_index(i);
      input_particles.push_back(psjet);
    }
  }

  /// helper routine to help transfer jets -> f77jets[4*ijet+0..3]
  void transfer_jets(double * f77jets, int & njets) {
    njets = jets.size();
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
	*f77jets = jets[i][j];
	f77jets++;
      }
    }
  }

  vector<PseudoJet> gran_jets ( vector<PseudoJet> & ori_hadrons,const double & eta_cell, const double & phi_cell, const double & pt_cutoff)
  {
    
    double pi = 3.142592654;
    vector<PseudoJet> granulated_jets;
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
		
		PseudoJet comb_jet( rescaled_px, rescaled_py,rescaled_pz, total_energy);
		
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
    
  };

  void output_vec_pseudojet(  ofstream & fout, vector<PseudoJet> & vec_jets)
  {
    for(unsigned idum=0;idum<vec_jets.size();idum++){
      fout << vec_jets.at(idum).perp() << " ";
      fout << vec_jets.at(idum).eta() << " ";
      fout << vec_jets.at(idum).phi_std() << endl;
    }
    return;
  };

}
HEPTOPTAGGER_END_NAMESPACE

using namespace HEPTopTaggerInt::fwrappertoptagger;

extern "C" {
  
  void heptoptagger_(const double * p, const int & npart,
		     const double & R, const double & ptjetmin,
		     const double & palg, const double & eta_cell,
		     const double & phi_cell, const double * topmass,
		     const double * wmass, double * f77top, 
		     double * f77bsubjet, double * f77Wsubjet1,
		     double * f77Wsubjet2, int & ntop,
		     const bool & print){

  // transfer p[4*ipart+0..3] -> input_particles[i]
    transfer_input_particles(p,npart);
    ofstream fout;
    ofstream fout2;

    if (print) {
      fout.open("sample_event_display.gnu");
      fout2.open("sample_event_display.dat");
      output_vec_pseudojet(fout2,input_particles);
    };

  // jet definition
    if (palg == 1.0) {
      jet_def = JetDefinition(kt_algorithm,R);
    } else if (palg == 0.0) {
      jet_def = JetDefinition(cambridge_algorithm, R);
    } else if (palg == -1.0) {
      jet_def = JetDefinition(antikt_algorithm, R);
    } else {
      jet_def = JetDefinition(genkt_algorithm, R, palg);
    }

  // start of granularization of the hadronic calorimeter to redefine hadrons
    double pt_cutoff=0.5;
    vector<PseudoJet> gran_hadrons
      = gran_jets(input_particles, eta_cell, phi_cell, pt_cutoff);

    // run the jet finding; find the hardest jet
    cs.reset(new ClusterSequence(gran_hadrons, jet_def));

    // extract jets (pt-ordered)
    jets = sorted_by_pt(cs->inclusive_jets(ptjetmin));

    // top quark mass and mass window
    valarray<double> top_mass(3); // top_mass[0..2]
                                  // top_mass[0] central topmass
                                  // top_mass[1] lower limit of top mass range
                                  // top_mass[2] upper limit of top mass range
    for (unsigned i=0;i<=2; i++) {
      top_mass[i]=*(topmass++);
    };

    // w boson mass and mass window
    valarray<double> w_mass(3); // w_mass[0..2]
                                // w_mass[0] central wmass
                                // w_mass[1] lower limit of w mass/top mass range
                                // w_mass[2] upper limit of w mass/top mass range
    for (unsigned i=0;i<=2; i++) {
      w_mass[i]=*(wmass++);
    };

    double mt=top_mass[0];
    double mw=w_mass[0];
    double mtl=top_mass[1];
    double mtu=top_mass[2];
    double rmin=w_mass[1];
    double rmax=w_mass[2];
    
    ntop=0;
    for(unsigned ijet=0; ijet<jets.size(); ijet++)
      {
	HEPTopTagger::HEPTopTagger cm_toptag(*cs,jets[ijet],mt,mw);
	cm_toptag.set_top_range(mtl,mtu);
	cm_toptag.set_mass_ratio_range(rmin,rmax);
	if (print) cout<< "========= Top Tagger ============" << endl;
	cm_toptag.run_tagger();
	if (print) {
	  cout<< "-------- setting  --------" << endl;
	  cm_toptag.get_setting();
	  cout<< "-------- results  --------" << endl;
	  cm_toptag.get_info();
	};

	if(cm_toptag.is_masscut_passed()){
	  if (print) cout << "### masscut_passed ###" << endl;
	  PseudoJet top=cm_toptag.top_candidate();
	  PseudoJet b=cm_toptag.top_subjets().at(0);
	  PseudoJet W1=cm_toptag.top_subjets().at(1); // subjet from W decay
	  PseudoJet W2=cm_toptag.top_subjets().at(2); // subjet from W decay
	  for (unsigned j=0;j<=3;j++){
	    *f77top = top[j];
	    *f77bsubjet = b[j];
	    *f77Wsubjet1 = W1[j];
	    *f77Wsubjet2 = W2[j];
	    f77top++;
	    f77bsubjet++;
	    f77Wsubjet1++;
	    f77Wsubjet2++;
	  };
	  ntop++;
	  if (print) {
	    double rad = max(eta_cell,phi_cell);
	    cout << "top mass: " << top.m() << endl;
	    cout << "bottom mass: "<< b.m() << endl;
	    cout << "W mass: "<< (W1+W2).m() << endl;
	    set_header_for_display(fout);
	    circle(fout,top,rad,1);
	    circle(fout,b,rad,2);
	    circle(fout,W1,rad,3);
	    circle(fout,W2,rad,3);
	    set_footer_for_display(fout);
	  };
	};
      };
  };    
}
