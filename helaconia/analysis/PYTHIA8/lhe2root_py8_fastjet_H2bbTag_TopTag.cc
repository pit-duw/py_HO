#include "Pythia8_lhe2root_Htag_Toptag.h"
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
int ihadr(int id);
double judge_best( double top_mass, double w_mass );

int main() {
  TFile outFile ( "helaconia_hepmc.root", "RECREATE" );
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
    vector<int>  B_INDEX;
    int n_B = 0;
    int n_P = 0;
    vector<int> BjetParticles_index;
    vector<bool> is_b_jet;
    vector<int> n_b_jet;

    // Loop over all particles in the event
    for (int i=1; i < event.size(); ++i){
      // Only final state particles
      if ( !event[i].isFinal() )continue;
      int pdgId = event[i].id();
      int quarkid = ihadr(pdgId);
      if ( fabs(pdgId) == 12 || fabs(pdgId) == 14 || fabs(pdgId) == 16 ){
	// Add ALL neutrinos to missing eT
	MET_p += fastjet::PseudoJet( event[i].px(), event[i].py(), event[i].pz(), event[i].e() );
      } else {
	// Fill all OTHER final state particles into a vector of (pseudo)jets for fastjet to cluster later
	// only hadrons are allowed
	if ( fabs(pdgId) >= 100 ){
	  fastjet::PseudoJet psjet ( event[i].px(), event[i].py(), event[i].pz(), event[i].e() );
	  psjet.set_user_index(n_P); // import for using user_index in the following
	  jetParticles.push_back( psjet );
	  if ( fabs(quarkid) == 5 ){
	    // B hadron for b tagging
	    n_B++;
	    B_INDEX.push_back(i);
	    BjetParticles_index.push_back(n_P);
	  };
	  n_P++;
	};
      }
      double perp=GetPt(event[i].px(),event[i].py());
      // Only store particles with pT above ParticlePtCut 
      // Change ParticlePtCut in the head file
      // and final state hadron or final state neutrino
      if ( perp > ParticlePtCut ){
	//double eta=GetPseudoRapidity(event[i].px(),event[i].py(),event[i].pz());
	//double phi=GetPhi(event[i].px(),event[i].py(),event[i].pz());
	n_particles++;
	//PID_v.push_back( pdgId );
	//P_X_v.push_back( event[i].px() );
	//P_Y_v.push_back( event[i].py() );
	//P_Z_v.push_back( event[i].pz() );
	//P_T_v.push_back( perp );
	//E_v.push_back( event[i].e() );
	//M_v.push_back( event[i].m() );
	//Eta_v.push_back( eta );
	//Phi_v.push_back( phi );
      };
    } // end of loop over particles
    // Cluster jets and sort by pT 
    // Change the jet_def in the head file
    fastjet::ClusterSequence cs(jetParticles, jet_def);
    jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    // Determine which parton/particle ended-up in which jet
    // set all jet entrie to zero first
    vector<int> whichjet;
    for(unsigned int ii=0; ii<n_P; ++ii) whichjet.push_back(0);
    for ( int j = 0; j < jets.size(); j++ ){
      // Only store jets with pT above JetPtCut
      // Change the JetPtCut in the head file
      if ( jets.at(j).perp() < JetPtCut ) continue;
      vector<fastjet::PseudoJet> constit = cs.constituents(jets.at(j));
      for(unsigned int ll=0; ll<constit.size(); ++ll){
	whichjet.at(constit[ll].user_index())=j+1;
      };
      is_b_jet.push_back(false);
      n_b_jet.push_back(0);
      for (unsigned int kk=0; kk<n_B; ++kk){
	if (whichjet.at(BjetParticles_index.at(kk)) == j+1 ){
	  is_b_jet.at(j)=true;
	  n_b_jet.at(j)=n_b_jet.at(j)+1;
	  if ( n_b_jet.at(j) >= 2 )break;
	}
      }
    };
    int nhiggs=0;
    int nhiggs_btag=0;
    int ntop=0;
    int ntop_btag=0;
    int jhntop=0;
    int jhntop_btag=0;
    for ( int j = 0; j < jets.size(); j++ ){
      // Only store jets with pT above JetPtCut 
      // Change the JetPtCut in the head file
      if ( jets.at(j).perp() < JetPtCut ) continue;
      bool htag_done=false;
      bool hbbtag_done=false;
      //bool toptag_done=false;
      //bool toptag_btag_done=false;
      Jet_n++;
      Jet_pt_v.push_back( jets.at(j).perp() );
      Jet_eta_v.push_back( jets.at(j).eta() );
      Jet_phi_v.push_back( jets.at(j).phi() );
      Jet_E_v.push_back( jets.at(j).E() );
      // Higgs tagging via mass drop and filter
      // following the way of arXiv:0802.2470
      if ( (jets.at(j).m() >= hmasslow) && (jets.at(j).m() <= hmassup) ){ 
	fastjet::PseudoJet tagged_jet = htagger(jets.at(j));
	if ( (tagged_jet != 0) ){
	  if ( (htag_mu == 0) && (htag_y == 0) ){
	    htag_mu = tagged_jet.structure_of<fastjet::MassDropTagger>().mu();
	    htag_y =  tagged_jet.structure_of<fastjet::MassDropTagger>().y();
	    htag_mass = tagged_jet.m();
	    htag_pt = tagged_jet.perp();
	    htag_eta = tagged_jet.eta();
            htag_rap = tagged_jet.rapidity();
	    htag_phi = tagged_jet.phi();
	    // using filter
	    vector<fastjet::PseudoJet> tagged_pieces = tagged_jet.pieces();
	    double Rfilt = min(0.3, 0.5*tagged_pieces[0].delta_R(tagged_pieces[1]));
	    fastjet::PseudoJet filtered_tagged_jet = fastjet::Filter(Rfilt, fastjet::SelectorNHardest(3))(tagged_jet);
	    htag_mass_afterfilter = filtered_tagged_jet.m();
	    htag_pt_afterfilter = filtered_tagged_jet.perp();
	    htag_eta_afterfilter = filtered_tagged_jet.eta();
            htag_rap_afterfilter = filtered_tagged_jet.rapidity();
	    htag_phi_afterfilter = filtered_tagged_jet.phi();
	    // using pt as Lambda_hard in 0912.1342
	    fastjet::PseudoJet trimmered_tagged_jet = fastjet::Filter(Rfilt, fastjet::SelectorPtFractionMin(fcut))(tagged_jet);
	    htag_mass_aftertrimmer = trimmered_tagged_jet.m();
	    htag_pt_aftertrimmer = trimmered_tagged_jet.perp();
	    htag_eta_aftertrimmer = trimmered_tagged_jet.eta();
	    htag_rap_aftertrimmer = trimmered_tagged_jet.rapidity();
	    htag_phi_aftertrimmer = trimmered_tagged_jet.phi();
	    htag_done=true;
	    nhiggs++;
	  };
	  // Higgs tagging via mass drop and filter as well as two b tagging
	  if ( ( n_b_jet.at(j) >= 2 ) && (hbbtag_mu == 0) && (hbbtag_y == 0) ){
	    hbbtag_mu = tagged_jet.structure_of<fastjet::MassDropTagger>().mu();
	    hbbtag_y = tagged_jet.structure_of<fastjet::MassDropTagger>().y();
	    hbbtag_mass = tagged_jet.m();
	    hbbtag_pt = tagged_jet.perp();
	    hbbtag_eta = tagged_jet.eta();
            hbbtag_rap = tagged_jet.rapidity();
	    hbbtag_phi = tagged_jet.phi();
	    vector<fastjet::PseudoJet> bbtagged_pieces = tagged_jet.pieces();
            double bbRfilt = min(0.3, 0.5*bbtagged_pieces[0].delta_R(bbtagged_pieces[1]));
	    fastjet::PseudoJet bbfiltered_tagged_jet = fastjet::Filter(bbRfilt, fastjet::SelectorNHardest(3))(tagged_jet);
            hbbtag_mass_afterfilter = bbfiltered_tagged_jet.m();
	    hbbtag_pt_afterfilter = bbfiltered_tagged_jet.perp();
	    hbbtag_eta_afterfilter = bbfiltered_tagged_jet.eta();
            hbbtag_rap_afterfilter = bbfiltered_tagged_jet.rapidity();
	    hbbtag_phi_afterfilter = bbfiltered_tagged_jet.phi();
	    // using pt as Lambda_hard in arXiv:0912.1342
	    fastjet::PseudoJet bbtrimmered_tagged_jet = fastjet::Filter(bbRfilt, fastjet::SelectorPtFractionMin(fcut))(tagged_jet);
            hbbtag_mass_aftertrimmer = bbtrimmered_tagged_jet.m();
            hbbtag_pt_aftertrimmer = bbtrimmered_tagged_jet.perp();
            hbbtag_eta_aftertrimmer = bbtrimmered_tagged_jet.eta();
            hbbtag_rap_aftertrimmer = bbtrimmered_tagged_jet.rapidity();
            hbbtag_phi_aftertrimmer = bbtrimmered_tagged_jet.phi();
	    hbbtag_done=true;
	    nhiggs_btag++;
	  };
	};
	// top tagger via HEPTopTagger (see arXiv:1006.2833)
	bool hcandidate=htag_done && htag_mass_afterfilter<hrecmassup && htag_mass_afterfilter>hrecmasslow;
	bool hbbcandidate=hbbtag_done && hbbtag_mass_afterfilter<hrecmassup && hbbtag_mass_afterfilter>hrecmasslow;
	if (!hcandidate && !hbbcandidate){
	  HEPTopTagger::HEPTopTagger cm_toptag(cs,jets.at(j),topmass,wmass);
	  cm_toptag.set_top_range(topmasslow,topmassup);
	  if ( ntop < 2 || (ntop_btag < 2 && n_b_jet.at(j) >= 1)){
	    cm_toptag.run_tagger();
	    if(cm_toptag.is_masscut_passed()){
	      fastjet::PseudoJet top=cm_toptag.top_candidate();
	      fastjet::PseudoJet b=cm_toptag.top_subjets().at(0);
	      fastjet::PseudoJet W1=cm_toptag.top_subjets().at(1);
	      fastjet::PseudoJet W2=cm_toptag.top_subjets().at(2);
	      fastjet::PseudoJet W=W1+W2;
	      if ( ntop < 2 ){
		if ( ntop == 0 ){
		  heptoptag_tmass1=top.m();
		  heptoptag_wmass1=W.m();
		  heptoptag_bmass1=b.m();
		  heptoptag_t1_pt=top.perp();
		  heptoptag_w1_pt=W.perp();
		  heptoptag_b1_pt=b.perp();
		  heptoptag_t1_eta=top.eta();
                  heptoptag_w1_eta=W.eta();
                  heptoptag_b1_eta=b.eta();
                  heptoptag_t1_rap=top.rapidity();
                  heptoptag_w1_rap=W.rapidity();
                  heptoptag_b1_rap=b.rapidity();
                  heptoptag_t1_phi=top.phi();
                  heptoptag_w1_phi=W.phi();
                  heptoptag_b1_phi=b.phi();
		} else if ( ntop == 1 ){
		  heptoptag_tmass2=top.m();
		  heptoptag_wmass2=W.m();
		  heptoptag_bmass2=b.m();
		  heptoptag_t2_pt=top.perp();
                  heptoptag_w2_pt=W.perp();
                  heptoptag_b2_pt=b.perp();
                  heptoptag_t2_eta=top.eta();
                  heptoptag_w2_eta=W.eta();
                  heptoptag_b2_eta=b.eta();
                  heptoptag_t2_rap=top.rapidity();
                  heptoptag_w2_rap=W.rapidity();
                  heptoptag_b2_rap=b.rapidity();
                  heptoptag_t2_phi=top.phi();
                  heptoptag_w2_phi=W.phi();
                  heptoptag_b2_phi=b.phi();
		};
		ntop++;
		if ( ntop == 2 ){
		  double judge1 = judge_best(heptoptag_tmass1, heptoptag_wmass1);
		  double judge2 = judge_best(heptoptag_tmass2, heptoptag_wmass2);
		  if ( judge1 <= judge2 ){
		    heptoptag_tmass_best = heptoptag_tmass1;
		    heptoptag_wmass_best = heptoptag_wmass1;
		    heptoptag_bmass_best = heptoptag_bmass1;
		    heptoptag_tbest_pt=heptoptag_t1_pt;
		    heptoptag_wbest_pt=heptoptag_w1_pt;
		    heptoptag_bbest_pt=heptoptag_b1_pt;
		    heptoptag_tbest_eta=heptoptag_t1_eta;
		    heptoptag_wbest_eta=heptoptag_w1_eta;
		    heptoptag_bbest_eta=heptoptag_b1_eta;
                    heptoptag_tbest_rap=heptoptag_t1_rap;
                    heptoptag_wbest_rap=heptoptag_w1_rap;
                    heptoptag_bbest_rap=heptoptag_b1_rap;
		    heptoptag_tbest_phi=heptoptag_t1_phi;
		    heptoptag_wbest_phi=heptoptag_w1_phi;
		    heptoptag_bbest_phi=heptoptag_b1_phi;
		  } else {
		    heptoptag_tmass_best = heptoptag_tmass2;
                    heptoptag_wmass_best = heptoptag_wmass2;
                    heptoptag_bmass_best = heptoptag_bmass2;
                    heptoptag_tbest_pt=heptoptag_t2_pt;
                    heptoptag_wbest_pt=heptoptag_w2_pt;
                    heptoptag_bbest_pt=heptoptag_b2_pt;
                    heptoptag_tbest_eta=heptoptag_t2_eta;
                    heptoptag_wbest_eta=heptoptag_w2_eta;
                    heptoptag_bbest_eta=heptoptag_b2_eta;
                    heptoptag_tbest_rap=heptoptag_t2_rap;
                    heptoptag_wbest_rap=heptoptag_w2_rap;
                    heptoptag_bbest_rap=heptoptag_b2_rap;
                    heptoptag_tbest_phi=heptoptag_t2_phi;
                    heptoptag_wbest_phi=heptoptag_w2_phi;
                    heptoptag_bbest_phi=heptoptag_b2_phi;
		  }
		}
	      };
	      if ( ntop_btag < 2 && n_b_jet.at(j) >= 1 ){
		if ( ntop_btag == 0 ){
		  heptoptag_btag_tmass1=top.m();
		  heptoptag_btag_wmass1=W.m();
		  heptoptag_btag_bmass1=b.m();
                  heptoptag_btag_t1_pt=top.perp();
                  heptoptag_btag_w1_pt=W.perp();
                  heptoptag_btag_b1_pt=b.perp();
                  heptoptag_btag_t1_eta=top.eta();
                  heptoptag_btag_w1_eta=W.eta();
                  heptoptag_btag_b1_eta=b.eta();
                  heptoptag_btag_t1_rap=top.rapidity();
                  heptoptag_btag_w1_rap=W.rapidity();
                  heptoptag_btag_b1_rap=b.rapidity();
                  heptoptag_btag_t1_phi=top.phi();
                  heptoptag_btag_w1_phi=W.phi();
                  heptoptag_btag_b1_phi=b.phi();
		} else if ( ntop_btag == 1 ){
		  heptoptag_btag_tmass2=top.m();
		  heptoptag_btag_wmass2=W.m();
		  heptoptag_btag_bmass2=b.m();
                  heptoptag_btag_t2_pt=top.perp();
                  heptoptag_btag_w2_pt=W.perp();
                  heptoptag_btag_b2_pt=b.perp();
                  heptoptag_btag_t2_eta=top.eta();
                  heptoptag_btag_w2_eta=W.eta();
                  heptoptag_btag_b2_eta=b.eta();
                  heptoptag_btag_t2_rap=top.rapidity();
                  heptoptag_btag_w2_rap=W.rapidity();
                  heptoptag_btag_b2_rap=b.rapidity();
                  heptoptag_btag_t2_phi=top.phi();
                  heptoptag_btag_w2_phi=W.phi();
                  heptoptag_btag_b2_phi=b.phi();
		};
		ntop_btag++;
		if ( ntop_btag == 2 ){
		  double bjudge1= judge_best(heptoptag_btag_tmass1, heptoptag_btag_wmass1);
                  double bjudge2= judge_best(heptoptag_btag_tmass2, heptoptag_btag_wmass2);
                  if ( bjudge1 <= bjudge2){
                    heptoptag_btag_tmass_best = heptoptag_btag_tmass1;
                    heptoptag_btag_wmass_best = heptoptag_btag_wmass1;
                    heptoptag_btag_bmass_best = heptoptag_btag_bmass1;
		    heptoptag_btag_tbest_pt=heptoptag_btag_t1_pt;
		    heptoptag_btag_wbest_pt=heptoptag_btag_w1_pt;
		    heptoptag_btag_bbest_pt=heptoptag_btag_b1_pt;
                    heptoptag_btag_tbest_eta=heptoptag_btag_t1_eta;
                    heptoptag_btag_wbest_eta=heptoptag_btag_w1_eta;
                    heptoptag_btag_bbest_eta=heptoptag_btag_b1_eta;
                    heptoptag_btag_tbest_rap=heptoptag_btag_t1_rap;
                    heptoptag_btag_wbest_rap=heptoptag_btag_w1_rap;
                    heptoptag_btag_bbest_rap=heptoptag_btag_b1_rap;
                    heptoptag_btag_tbest_phi=heptoptag_btag_t1_phi;
                    heptoptag_btag_wbest_phi=heptoptag_btag_w1_phi;
                    heptoptag_btag_bbest_phi=heptoptag_btag_b1_phi;
                  } else {
                    heptoptag_btag_tmass_best = heptoptag_btag_tmass2;
                    heptoptag_btag_wmass_best = heptoptag_btag_wmass2;
                    heptoptag_btag_bmass_best = heptoptag_btag_bmass2;
                    heptoptag_btag_tbest_pt=heptoptag_btag_t2_pt;
                    heptoptag_btag_wbest_pt=heptoptag_btag_w2_pt;
                    heptoptag_btag_bbest_pt=heptoptag_btag_b2_pt;
                    heptoptag_btag_tbest_eta=heptoptag_btag_t2_eta;
                    heptoptag_btag_wbest_eta=heptoptag_btag_w2_eta;
                    heptoptag_btag_bbest_eta=heptoptag_btag_b2_eta;
                    heptoptag_btag_tbest_rap=heptoptag_btag_t2_rap;
                    heptoptag_btag_wbest_rap=heptoptag_btag_w2_rap;
                    heptoptag_btag_bbest_rap=heptoptag_btag_b2_rap;
                    heptoptag_btag_tbest_phi=heptoptag_btag_t2_phi;
                    heptoptag_btag_wbest_phi=heptoptag_btag_w2_phi;
                    heptoptag_btag_bbest_phi=heptoptag_btag_b2_phi;
                  }
		}// end of assign b-tagged top for the best one
	      }; // end of b-tagged top tagging
	    };// end of pass mass cut
	  }; // end of deciding HEP top tagger or not
	  // starting JHU top tagger
	  fastjet::JHTopTagger jhtop_tagger(delta_p,delta_r,cos_theta_W_max);
	  // indicate the acceptable range of top, W masses 
	  jhtop_tagger.set_top_selector(fastjet::SelectorMassRange(topmasslow,topmassup));
	  jhtop_tagger.set_W_selector(fastjet::SelectorMassRange(wmasslow, wmassup)); 
	  if ( jhntop < 2 || (jhntop_btag < 2 && n_b_jet.at(j) >= 1)){
	    // try and tag a jet
	    // jet should come from a C/A clustering
	    fastjet::PseudoJet top_candidate = jhtop_tagger(jets.at(j));
	    if ( top_candidate != 0 ){
	      if ( jhntop < 2 ){
                if ( jhntop == 0 ){
                  jhtoptag_tmass1=top_candidate.m();
                  jhtoptag_wmass1=top_candidate.structure_of<fastjet::JHTopTagger>().W().m();
		  jhtoptag_bmass1=top_candidate.structure_of<fastjet::JHTopTagger>().non_W().m();
                } else if ( jhntop == 1 ){
                  jhtoptag_tmass2=top_candidate.m();
                  jhtoptag_wmass2=top_candidate.structure_of<fastjet::JHTopTagger>().W().m();
		  jhtoptag_bmass2=top_candidate.structure_of<fastjet::JHTopTagger>().non_W().m();
                }; // end of assining top and w masses
		jhntop++;
                if ( jhntop == 2 ){
                  double jhjudge1 = judge_best(jhtoptag_tmass1, jhtoptag_wmass1);
                  double jhjudge2 = judge_best(jhtoptag_tmass2, jhtoptag_wmass2);
                  if ( jhjudge1 <= jhjudge2 ){
                    jhtoptag_tmass_best = jhtoptag_tmass1;
                    jhtoptag_wmass_best = jhtoptag_wmass1;
		    jhtoptag_bmass_best = jhtoptag_bmass1;
                  } else {
                    jhtoptag_tmass_best = jhtoptag_tmass2;
                    jhtoptag_wmass_best = jhtoptag_wmass2;
		    jhtoptag_bmass_best = jhtoptag_bmass2;
                  }
                }; // end of assigning best top and w masses when there are 2 tops		
	      }; // end of non-b tagged top tagging
	      if ( jhntop_btag < 2 && n_b_jet.at(j) >= 1 ){
                if ( jhntop_btag == 0 ){
                  jhtoptag_btag_tmass1=top_candidate.m();
                  jhtoptag_btag_wmass1=top_candidate.structure_of<fastjet::JHTopTagger>().W().m();
		  jhtoptag_btag_bmass1=top_candidate.structure_of<fastjet::JHTopTagger>().non_W().m();
                } else if ( jhntop_btag == 1 ){
                  jhtoptag_btag_tmass2=top_candidate.m();
                  jhtoptag_btag_wmass2=top_candidate.structure_of<fastjet::JHTopTagger>().W().m();
		  jhtoptag_btag_bmass2=top_candidate.structure_of<fastjet::JHTopTagger>().non_W().m();
                };// end of assining top and w masses 
                jhntop_btag++;
		if ( jhntop_btag == 2 ){
                  double jhbjudge1= judge_best(jhtoptag_btag_tmass1, jhtoptag_btag_wmass1);
                  double jhbjudge2= judge_best(jhtoptag_btag_tmass2, jhtoptag_btag_wmass2);
                  if ( jhbjudge1 <= jhbjudge2){
                    jhtoptag_btag_tmass_best = jhtoptag_btag_tmass1;
                    jhtoptag_btag_wmass_best = jhtoptag_btag_wmass1;
		    jhtoptag_btag_bmass_best = jhtoptag_btag_bmass1;
                  } else {
                    jhtoptag_btag_tmass_best = jhtoptag_btag_tmass2;
                    jhtoptag_btag_wmass_best = jhtoptag_btag_wmass2;
		    jhtoptag_btag_bmass_best = jhtoptag_btag_bmass2;
                  }
                }// end of assign b-tagged top for the best one
	      }; // end of b-tagged top tagging
	    }; // end of generating a top tagging and passing the mass cut
	  };// end of deciding JH top tagger or not
	} // end of passing jet mass cut
      }; // end of loop over jet
      // in the case of only one top reconstructed
      if ( ntop == 1 ){
	heptoptag_tmass_best = heptoptag_tmass1;
	heptoptag_wmass_best = heptoptag_wmass1;
	heptoptag_bmass_best = heptoptag_bmass1;
	heptoptag_tbest_pt=heptoptag_t1_pt;
	heptoptag_wbest_pt=heptoptag_w1_pt;
	heptoptag_bbest_pt=heptoptag_b1_pt;
	heptoptag_tbest_eta=heptoptag_t1_eta;
	heptoptag_wbest_eta=heptoptag_w1_eta;
	heptoptag_bbest_eta=heptoptag_b1_eta;
        heptoptag_tbest_rap=heptoptag_t1_rap;
        heptoptag_wbest_rap=heptoptag_w1_rap;
        heptoptag_bbest_rap=heptoptag_b1_rap;
	heptoptag_tbest_phi=heptoptag_t1_phi;
	heptoptag_wbest_phi=heptoptag_w1_phi;
	heptoptag_bbest_phi=heptoptag_b1_phi;
      };
      if ( ntop_btag == 1 ){
	heptoptag_btag_tmass_best = heptoptag_btag_tmass1;
	heptoptag_btag_wmass_best = heptoptag_btag_wmass1;
	heptoptag_btag_bmass_best = heptoptag_btag_bmass1;
	heptoptag_btag_tbest_pt=heptoptag_btag_t1_pt;
	heptoptag_btag_wbest_pt=heptoptag_btag_w1_pt;
	heptoptag_btag_bbest_pt=heptoptag_btag_b1_pt;
	heptoptag_btag_tbest_eta=heptoptag_btag_t1_eta;
	heptoptag_btag_wbest_eta=heptoptag_btag_w1_eta;
	heptoptag_btag_bbest_eta=heptoptag_btag_b1_eta;
        heptoptag_btag_tbest_rap=heptoptag_btag_t1_rap;
        heptoptag_btag_wbest_rap=heptoptag_btag_w1_rap;
        heptoptag_btag_bbest_rap=heptoptag_btag_b1_rap;
	heptoptag_btag_tbest_phi=heptoptag_btag_t1_phi;
	heptoptag_btag_wbest_phi=heptoptag_btag_w1_phi;
	heptoptag_btag_bbest_phi=heptoptag_btag_b1_phi;
      };
      if ( jhntop == 1 ){
        jhtoptag_tmass_best = jhtoptag_tmass1;
        jhtoptag_wmass_best = jhtoptag_wmass1;
	jhtoptag_bmass_best = jhtoptag_bmass1;
      };
      if ( jhntop_btag == 1 ){
        jhtoptag_btag_tmass_best = jhtoptag_btag_tmass1;
        jhtoptag_btag_wmass_best = jhtoptag_btag_wmass1;
	jhtoptag_btag_bmass_best = jhtoptag_btag_bmass1;
      };
      // observables defined between top and higgs
      double deltaphi,deltay;
      if ( ntop+nhiggs > 1 ){
	if ( ntop > 1 ){
	  deltaphi=heptoptag_t1_phi-heptoptag_t2_phi;
	  deltay=heptoptag_t1_rap-heptoptag_t2_rap;
	  deltaR_t1heptag_t2heptag=sqrt(deltaphi*deltaphi+deltay*deltay);
	  deltaRmin_httheptag = deltaR_t1heptag_t2heptag;
	};
	if ( nhiggs > 0 && ntop > 0 ){
	  deltaphi=heptoptag_t1_phi-htag_phi_afterfilter;
	  deltay=heptoptag_t1_rap-htag_rap_afterfilter;
	  deltaR_htag_t1heptag = sqrt(deltaphi*deltaphi+deltay*deltay);
	  if ( ntop > 1 ){
	    deltaphi=heptoptag_t2_phi-htag_phi_afterfilter;
	    deltay=heptoptag_t2_rap-htag_rap_afterfilter;
	    deltaR_htag_t2heptag = sqrt(deltaphi*deltaphi+deltay*deltay);
	  };
	  deltaphi=heptoptag_tbest_phi-htag_phi_afterfilter;
	  deltay=heptoptag_tbest_rap-htag_rap_afterfilter;
	  deltaR_htag_tbestheptag = sqrt(deltaphi*deltaphi+deltay*deltay);
	  if ( ntop > 1 ){
	    deltaRmin_httheptag = min(deltaRmin_httheptag, deltaR_htag_t1heptag);
	    deltaRmin_httheptag = min(deltaRmin_httheptag, deltaR_htag_t2heptag);
	  } else {
	    deltaRmin_httheptag = deltaR_htag_t1heptag;
	  };
	};
      };
      if ( ntop_btag+nhiggs_btag > 1 ){
        if ( ntop_btag > 1 ){
          deltaphi=heptoptag_btag_t1_phi-heptoptag_btag_t2_phi;
          deltay=heptoptag_btag_t1_rap-heptoptag_btag_t2_rap;
          deltaR_t1btagheptag_t2btagheptag=sqrt(deltaphi*deltaphi+deltay*deltay);
          deltaRmin_httbtagheptag = deltaR_t1btagheptag_t2btagheptag;
        };
        if ( nhiggs_btag > 0 && ntop_btag > 0 ){
          deltaphi=heptoptag_btag_t1_phi-hbbtag_phi_afterfilter;
          deltay=heptoptag_btag_t1_rap-hbbtag_rap_afterfilter;
          deltaR_hbbtag_t1btagheptag = sqrt(deltaphi*deltaphi+deltay*deltay);
          if ( ntop_btag > 1 ){
            deltaphi=heptoptag_btag_t2_phi-hbbtag_phi_afterfilter;
            deltay=heptoptag_btag_t2_rap-hbbtag_rap_afterfilter;
            deltaR_hbbtag_t2btagheptag = sqrt(deltaphi*deltaphi+deltay*deltay);
          };
          deltaphi=heptoptag_btag_tbest_phi-hbbtag_phi_afterfilter;
          deltay=heptoptag_btag_tbest_rap-hbbtag_rap_afterfilter;
          deltaR_hbbtag_tbestbtagheptag = sqrt(deltaphi*deltaphi+deltay*deltay);
          if ( ntop_btag > 1 ){
            deltaRmin_httbtagheptag = min(deltaRmin_httbtagheptag, deltaR_hbbtag_t1btagheptag);
            deltaRmin_httbtagheptag = min(deltaRmin_httbtagheptag, deltaR_hbbtag_t2btagheptag);
          } else {
            deltaRmin_httbtagheptag = deltaR_hbbtag_t1btagheptag;
          };
        };
      };
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
  //PID_v.clear();
  //P_X_v.clear();
  //P_Y_v.clear();
  //P_Z_v.clear();
  //P_T_v.clear();
  //E_v.clear();
  //M_v.clear();
  //Eta_v.clear();
  //Phi_v.clear();
  Jet_E_v.clear();
  Jet_pt_v.clear();
  Jet_eta_v.clear();
  Jet_phi_v.clear();
  htag_mu = 0;
  htag_y = 0;
  htag_mass = 0;
  htag_mass_afterfilter = 0;
  htag_mass_aftertrimmer = 0;
  htag_pt = 0;
  htag_pt_afterfilter = 0;
  htag_pt_aftertrimmer = 0;
  htag_eta = 0;
  htag_eta_afterfilter = 0;
  htag_eta_aftertrimmer = 0;
  htag_rap = 0;
  htag_rap_afterfilter = 0;
  htag_rap_aftertrimmer = 0;
  htag_phi = 0;
  htag_phi_afterfilter = 0;
  htag_phi_aftertrimmer = 0;
  hbbtag_mu = 0;
  hbbtag_y = 0;
  hbbtag_mass = 0;
  hbbtag_mass_afterfilter = 0;
  hbbtag_mass_aftertrimmer = 0;
  hbbtag_pt = 0;
  hbbtag_pt_afterfilter = 0;
  hbbtag_pt_aftertrimmer = 0;
  hbbtag_eta = 0;
  hbbtag_eta_afterfilter = 0;
  hbbtag_eta_aftertrimmer = 0;
  hbbtag_rap = 0;
  hbbtag_rap_afterfilter = 0;
  hbbtag_rap_aftertrimmer = 0;
  hbbtag_phi = 0;
  hbbtag_phi_afterfilter = 0;
  hbbtag_phi_aftertrimmer = 0;
  heptoptag_tmass1 = 0;
  heptoptag_bmass1 = 0;
  heptoptag_wmass1 = 0;
  heptoptag_t1_pt = 0;
  heptoptag_b1_pt = 0;
  heptoptag_w1_pt = 0;
  heptoptag_t1_eta = 0;
  heptoptag_b1_eta = 0;
  heptoptag_w1_eta = 0;
  heptoptag_t1_rap = 0;
  heptoptag_b1_rap = 0;
  heptoptag_w1_rap = 0;
  heptoptag_t1_phi = 0;
  heptoptag_b1_phi = 0;
  heptoptag_w1_phi = 0;
  heptoptag_tmass2 = 0;
  heptoptag_bmass2 = 0;
  heptoptag_wmass2 = 0;
  heptoptag_t2_pt = 0;
  heptoptag_b2_pt = 0;
  heptoptag_w2_pt = 0;
  heptoptag_t2_eta = 0;
  heptoptag_b2_eta = 0;
  heptoptag_w2_eta = 0;
  heptoptag_t2_rap = 0;
  heptoptag_b2_rap = 0;
  heptoptag_w2_rap = 0;
  heptoptag_t2_phi = 0;
  heptoptag_b2_phi = 0;
  heptoptag_w2_phi = 0;
  heptoptag_tmass_best = 0;
  heptoptag_bmass_best = 0;
  heptoptag_wmass_best = 0;
  heptoptag_tbest_pt = 0;
  heptoptag_bbest_pt = 0;
  heptoptag_wbest_pt = 0;
  heptoptag_tbest_eta = 0;
  heptoptag_bbest_eta = 0;
  heptoptag_wbest_eta = 0;
  heptoptag_tbest_rap = 0;
  heptoptag_bbest_rap = 0;
  heptoptag_wbest_rap = 0;
  heptoptag_tbest_phi = 0;
  heptoptag_bbest_phi = 0;
  heptoptag_wbest_phi = 0;
  heptoptag_btag_tmass1 = 0;
  heptoptag_btag_bmass1 = 0;
  heptoptag_btag_wmass1 = 0;
  heptoptag_btag_t1_pt = 0;
  heptoptag_btag_b1_pt = 0;
  heptoptag_btag_w1_pt = 0;
  heptoptag_btag_t1_eta = 0;
  heptoptag_btag_b1_eta = 0;
  heptoptag_btag_w1_eta = 0;
  heptoptag_btag_t1_rap = 0;
  heptoptag_btag_b1_rap = 0;
  heptoptag_btag_w1_rap = 0;
  heptoptag_btag_t1_phi = 0;
  heptoptag_btag_b1_phi = 0;
  heptoptag_btag_w1_phi = 0;
  heptoptag_btag_tmass2 = 0;
  heptoptag_btag_bmass2 = 0;
  heptoptag_btag_wmass2 = 0;
  heptoptag_btag_t2_pt = 0;
  heptoptag_btag_b2_pt = 0;
  heptoptag_btag_w2_pt = 0;
  heptoptag_btag_t2_eta = 0;
  heptoptag_btag_b2_eta = 0;
  heptoptag_btag_w2_eta = 0;
  heptoptag_btag_t2_rap = 0;
  heptoptag_btag_b2_rap = 0;
  heptoptag_btag_w2_rap = 0;
  heptoptag_btag_t2_phi = 0;
  heptoptag_btag_b2_phi = 0;
  heptoptag_btag_w2_phi = 0;
  heptoptag_btag_tmass_best = 0;
  heptoptag_btag_bmass_best = 0;
  heptoptag_btag_wmass_best = 0;
  heptoptag_btag_tbest_pt = 0;
  heptoptag_btag_bbest_pt = 0;
  heptoptag_btag_wbest_pt = 0;
  heptoptag_btag_tbest_eta = 0;
  heptoptag_btag_bbest_eta = 0;
  heptoptag_btag_wbest_eta = 0;
  heptoptag_btag_tbest_rap = 0;
  heptoptag_btag_bbest_rap = 0;
  heptoptag_btag_wbest_rap = 0;
  heptoptag_btag_tbest_phi = 0;
  heptoptag_btag_bbest_phi = 0;
  heptoptag_btag_wbest_phi = 0;
  jhtoptag_tmass1 = 0;
  jhtoptag_bmass1 = 0;
  jhtoptag_wmass1 = 0;
  jhtoptag_tmass2 = 0;
  jhtoptag_bmass2 = 0;
  jhtoptag_wmass2 = 0;
  jhtoptag_tmass_best = 0;
  jhtoptag_bmass_best = 0;
  jhtoptag_wmass_best = 0;
  jhtoptag_btag_tmass1 = 0;
  jhtoptag_btag_bmass1 = 0;
  jhtoptag_btag_wmass1 = 0;
  jhtoptag_btag_tmass2 = 0;
  jhtoptag_btag_bmass2 = 0;
  jhtoptag_btag_wmass2 = 0;
  jhtoptag_btag_tmass_best = 0;
  jhtoptag_btag_bmass_best = 0;
  jhtoptag_btag_wmass_best = 0;
  deltaR_htag_t1heptag = 0;
  deltaR_htag_t2heptag = 0;
  deltaR_t1heptag_t2heptag = 0;
  deltaR_htag_tbestheptag = 0;
  deltaRmin_httheptag = 0;
  deltaR_hbbtag_t1btagheptag = 0;
  deltaR_hbbtag_t2btagheptag = 0;
  deltaR_t1btagheptag_t2btagheptag = 0;
  deltaR_hbbtag_tbestbtagheptag = 0;
  deltaRmin_httbtagheptag = 0;
  jetParticles.clear();
  jets.clear();
  MET_p.reset(0,0,0,0);
}

// the value of w1*|top_mass-topmass|+w2*|w_mass-wmass|
// where w1,w2,topmass,wmass are defined in the head file
double judge_best( double top_mass, double w_mass )
{
  return w1*fabs(top_mass-topmass)+w2*fabs(w_mass-wmass);
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

//____________ probablity ( in %) of b tagging (mis-)identification
// id1 is the PDG of the heavier quark in the hadron of PDG code ID
double btag_efficiency(int id1){
  if ( id1 == 5 ){
    // it is a beauty jet
    return 60.;
  } else if ( id1 == 4 ){
    // it is a charm jet
    // mis-identification
    return 10.;
  } else {
    // it is a light flavor jet
    // mis-identification
    return 2.;
  }
}

// return the PDG code of the heavier quark in the hadron of PDG code ID
int ihadr(int id){
  if ( id != 0 ){
    int id1=abs(id);
    if ( id1 > 10000 )id1=id1-1000*(id1/1000);
    int nn=int(log10(float(id1)));
    int id2=int(id1/pow(float(10),nn));
    return id2;
  } else {
    return 0;
  };
}
