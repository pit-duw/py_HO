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
#include "fastjet/tools/JHTopTagger.hh"

using namespace std;

//vector<int>     PID_v;
//vector<double>  P_X_v;
//vector<double>  P_Y_v;
//vector<double>  P_Z_v;
//vector<double>  P_T_v;
//vector<double>  E_v;
//vector<double>  M_v;
//vector<double>  Eta_v;
//vector<double>  Phi_v;

vector<double>  Jet_pt_v;
vector<double>  Jet_E_v;
vector<double>  Jet_eta_v;
vector<double>  Jet_phi_v;

fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 1.5 );
vector<fastjet::PseudoJet> jetParticles;
vector<fastjet::PseudoJet> jets;
fastjet::PseudoJet MET_p;
fastjet::MassDropTagger htagger(0.67, 0.15 );
// using trimmer, fcut defined in 0912.1342
double fcut=1e-2;
double hmasslow=0.;
double hmassup=1e30;
double hrecmasslow=115.;
double hrecmassup=135.;
double topmass=173.3;
double wmass=80.4;
double topmasslow=150.;
double topmassup=200.;
// for the toptagging to choose the smallest one of w1*|m_top^rec-m_top|+w2*|m_w^rec-m_w|
double w1=1.;
double w2=1.;
// the following is for JHU Toptagger (see arXiv:0806.0848)
// also see page 49 in FastJet3 manual (arXiv:1111.6097)
double delta_p = 0.10; // subjets must carry at least this fraction of original jet's pt
double delta_r = 0.19; // subjets must be separated by at least this Manhattan distance
double cos_theta_W_max = 0.7; // the maximal allowed value of the W helicity angle
//fastjet::JHTopTagger jhtop_tagger(delta_p,delta_r,cos_theta_W_max);
// indicate the acceptable range of top, W masses
//jhtop_tagger.set_top_selector(fastjet::SelectorMassRange(topmasslow,topmassup));
double wmasslow=65.;
double wmassup=95.;
//jhtop_tagger.set_W_selector(fastjet::SelectorMassRange(wmasslow, wmassup));

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
double htag_mass_aftertrimmer = 0;
double htag_pt = 0;
double htag_pt_afterfilter = 0;
double htag_pt_aftertrimmer = 0;
double htag_eta = 0;
double htag_rap = 0;
double htag_phi = 0;
double htag_eta_afterfilter = 0;
double htag_rap_afterfilter = 0;
double htag_phi_afterfilter = 0;
double htag_eta_aftertrimmer = 0;
double htag_rap_aftertrimmer = 0;
double htag_phi_aftertrimmer = 0;
// using massdrop for higgs tag (with two b tagging)
double hbbtag_mu = 0;
double hbbtag_y = 0;
double hbbtag_mass = 0;
double hbbtag_mass_afterfilter = 0;
double hbbtag_mass_aftertrimmer = 0;
double hbbtag_pt = 0;
double hbbtag_pt_afterfilter = 0;
double hbbtag_pt_aftertrimmer = 0;
double hbbtag_eta = 0;
double hbbtag_rap = 0;
double hbbtag_phi = 0;
double hbbtag_eta_afterfilter = 0;
double hbbtag_rap_afterfilter = 0;
double hbbtag_phi_afterfilter = 0;
double hbbtag_eta_aftertrimmer = 0;
double hbbtag_rap_aftertrimmer = 0;
double hbbtag_phi_aftertrimmer = 0;
// toptagger with HEPTopTagger
// first top in order of jet pt
double heptoptag_tmass1 = 0;
double heptoptag_bmass1 = 0;
double heptoptag_wmass1 = 0;
double heptoptag_t1_pt = 0;
double heptoptag_b1_pt = 0;
double heptoptag_w1_pt = 0;
double heptoptag_t1_eta = 0;
double heptoptag_b1_eta = 0;
double heptoptag_w1_eta = 0;
double heptoptag_t1_rap = 0;
double heptoptag_b1_rap = 0;
double heptoptag_w1_rap = 0;
double heptoptag_t1_phi = 0;
double heptoptag_b1_phi = 0;
double heptoptag_w1_phi = 0;
// second top in order of jet pt
double heptoptag_tmass2 = 0;
double heptoptag_bmass2 = 0;
double heptoptag_wmass2 = 0;
double heptoptag_t2_pt = 0;
double heptoptag_b2_pt = 0;
double heptoptag_w2_pt = 0;
double heptoptag_t2_eta= 0;
double heptoptag_b2_eta= 0;
double heptoptag_w2_eta= 0;
double heptoptag_t2_rap= 0;
double heptoptag_b2_rap= 0;
double heptoptag_w2_rap= 0;
double heptoptag_t2_phi = 0;
double heptoptag_b2_phi = 0;
double heptoptag_w2_phi = 0;
// the closest one with the smallest w1*|m_top^rec-m_top|+w2*|m_w^rec-m_w|
double heptoptag_tmass_best = 0;
double heptoptag_wmass_best = 0;
double heptoptag_bmass_best = 0;
double heptoptag_tbest_pt = 0;
double heptoptag_wbest_pt = 0;
double heptoptag_bbest_pt = 0;
double heptoptag_tbest_eta= 0;
double heptoptag_bbest_eta= 0;
double heptoptag_wbest_eta= 0;
double heptoptag_tbest_rap= 0;
double heptoptag_bbest_rap= 0;
double heptoptag_wbest_rap= 0;
double heptoptag_tbest_phi = 0;
double heptoptag_bbest_phi = 0;
double heptoptag_wbest_phi = 0;
// one toptagger with a b tagging by HEPTopTagger
// first top in order of jet pt
double heptoptag_btag_tmass1 = 0;
double heptoptag_btag_bmass1 = 0;
double heptoptag_btag_wmass1 = 0;
double heptoptag_btag_t1_pt = 0;
double heptoptag_btag_b1_pt = 0;
double heptoptag_btag_w1_pt = 0;
double heptoptag_btag_t1_eta= 0;
double heptoptag_btag_b1_eta= 0;
double heptoptag_btag_w1_eta= 0;
double heptoptag_btag_t1_rap= 0;
double heptoptag_btag_b1_rap= 0;
double heptoptag_btag_w1_rap= 0;
double heptoptag_btag_t1_phi = 0;
double heptoptag_btag_b1_phi = 0;
double heptoptag_btag_w1_phi = 0;
// second top in order of jet pt
double heptoptag_btag_tmass2 = 0;
double heptoptag_btag_bmass2 = 0;
double heptoptag_btag_wmass2 = 0;
double heptoptag_btag_t2_pt = 0;
double heptoptag_btag_b2_pt = 0;
double heptoptag_btag_w2_pt = 0;
double heptoptag_btag_t2_eta= 0;
double heptoptag_btag_b2_eta= 0;
double heptoptag_btag_w2_eta= 0;
double heptoptag_btag_t2_rap= 0;
double heptoptag_btag_b2_rap= 0;
double heptoptag_btag_w2_rap= 0;
double heptoptag_btag_t2_phi = 0;
double heptoptag_btag_b2_phi = 0;
double heptoptag_btag_w2_phi = 0;
// the closest one with the smallest w1*|m_top^rec-m_top|+w2*|m_w^rec-m_w|
double heptoptag_btag_tmass_best = 0;
double heptoptag_btag_bmass_best = 0;
double heptoptag_btag_wmass_best = 0;
double heptoptag_btag_tbest_pt = 0;
double heptoptag_btag_bbest_pt = 0;
double heptoptag_btag_wbest_pt = 0;
double heptoptag_btag_tbest_eta= 0;
double heptoptag_btag_bbest_eta= 0;
double heptoptag_btag_wbest_eta= 0;
double heptoptag_btag_tbest_rap= 0;
double heptoptag_btag_bbest_rap= 0;
double heptoptag_btag_wbest_rap= 0;
double heptoptag_btag_tbest_phi = 0;
double heptoptag_btag_bbest_phi = 0;
double heptoptag_btag_wbest_phi = 0;
// toptagger with JHTopTagger in FastJet3
// first top in order of jet pt
double jhtoptag_tmass1 = 0;
double jhtoptag_bmass1 = 0;
double jhtoptag_wmass1 = 0;
// second top in order of jet pt
double jhtoptag_tmass2 = 0;
double jhtoptag_bmass2 = 0;
double jhtoptag_wmass2 = 0;
// the closest one with the smallest w1*|m_top^rec-m_top|+w2*|m_w^rec-m_w|
double jhtoptag_tmass_best = 0;
double jhtoptag_bmass_best = 0;
double jhtoptag_wmass_best = 0;
// one toptagger with a b tagging by JHTopTagger
// first top in order of jet pt
double jhtoptag_btag_tmass1 = 0;
double jhtoptag_btag_bmass1 = 0;
double jhtoptag_btag_wmass1 = 0;
// second top in order of jet pt
double jhtoptag_btag_tmass2 = 0;
double jhtoptag_btag_bmass2 = 0;
double jhtoptag_btag_wmass2 = 0;
// the closest one with the smallest w1*|m_top^rec-m_top|+w2*|m_w^rec-m_w|
double jhtoptag_btag_tmass_best = 0;
double jhtoptag_btag_bmass_best = 0;
double jhtoptag_btag_wmass_best = 0;
// the ones between top and higgs
double deltaR_htag_t1heptag = 0;
double deltaR_htag_t2heptag = 0;
double deltaR_t1heptag_t2heptag = 0;
double deltaR_htag_tbestheptag = 0;
double deltaRmin_httheptag = 0;
double deltaR_hbbtag_t1btagheptag = 0;
double deltaR_hbbtag_t2btagheptag = 0;
double deltaR_t1btagheptag_t2btagheptag = 0;
double deltaR_hbbtag_tbestbtagheptag = 0;
double deltaRmin_httbtagheptag = 0;

double JetPtCut=200;
double ParticlePtCut=0;

std::vector<std::string> files;

void ResetTreeVars();

void SetupTree( TTree* tree ){
  tree->Branch("n_particles", &(n_particles));
  //tree->Branch("PID", &(PID_v) );
  //tree->Branch("P_X", &(P_X_v) );
  //tree->Branch("P_Y", &(P_Y_v) );
  //tree->Branch("P_Z", &(P_Z_v) );
  //tree->Branch("P_T", &(P_T_v) );
  //tree->Branch("E", &(E_v) );
  //tree->Branch("M", &(M_v) );
  //tree->Branch("Eta", &(Eta_v) );
  //tree->Branch("Phi", &(Phi_v) );
  tree->Branch("Jet_n", &(Jet_n));
  tree->Branch("Jet_E", &(Jet_E_v) );
  tree->Branch("Jet_pt", &(Jet_pt_v) );
  tree->Branch("Jet_eta", &(Jet_eta_v) );
  tree->Branch("Jet_phi", &(Jet_phi_v) );
  tree->Branch("htag_mu", &(htag_mu));
  tree->Branch("htag_y",&(htag_y));
  tree->Branch("htag_mass",&(htag_mass));
  tree->Branch("htag_mass_afterfilter",&(htag_mass_afterfilter));
  tree->Branch("htag_mass_aftertrimmer",&(htag_mass_aftertrimmer));
  tree->Branch("htag_pt",&(htag_pt));
  tree->Branch("htag_pt_afterfilter",&(htag_pt_afterfilter));
  tree->Branch("htag_pt_aftertrimmer",&(htag_pt_aftertrimmer));
  tree->Branch("htag_eta",&(htag_eta));
  tree->Branch("htag_eta_afterfilter",&(htag_eta_afterfilter));
  tree->Branch("htag_eta_aftertrimmer",&(htag_eta_aftertrimmer));
  tree->Branch("htag_rap",&(htag_rap));
  tree->Branch("htag_rap_afterfilter",&(htag_rap_afterfilter));
  tree->Branch("htag_rap_aftertrimmer",&(htag_rap_aftertrimmer));
  tree->Branch("htag_phi",&(htag_phi));
  tree->Branch("htag_phi_afterfilter",&(htag_phi_afterfilter));
  tree->Branch("htag_phi_aftertrimmer",&(htag_phi_aftertrimmer));
  tree->Branch("hbbtag_mu", &(hbbtag_mu));
  tree->Branch("hbbtag_y",&(hbbtag_y));
  tree->Branch("hbbtag_mass",&(hbbtag_mass));
  tree->Branch("hbbtag_mass_afterfilter",&(hbbtag_mass_afterfilter));
  tree->Branch("hbbtag_mass_aftertrimmer",&(hbbtag_mass_aftertrimmer));
  tree->Branch("hbbtag_pt",&(hbbtag_pt));
  tree->Branch("hbbtag_pt_afterfilter",&(hbbtag_pt_afterfilter));
  tree->Branch("hbbtag_pt_aftertrimmer",&(hbbtag_pt_aftertrimmer));
  tree->Branch("hbbtag_eta",&(hbbtag_eta));
  tree->Branch("hbbtag_eta_afterfilter",&(hbbtag_eta_afterfilter));
  tree->Branch("hbbtag_eta_aftertrimmer",&(hbbtag_eta_aftertrimmer));
  tree->Branch("hbbtag_rap",&(hbbtag_rap));
  tree->Branch("hbbtag_rap_afterfilter",&(hbbtag_rap_afterfilter));
  tree->Branch("hbbtag_rap_aftertrimmer",&(hbbtag_rap_aftertrimmer));
  tree->Branch("hbbtag_phi",&(hbbtag_phi));
  tree->Branch("hbbtag_phi_afterfilter",&(hbbtag_phi_afterfilter));
  tree->Branch("hbbtag_phi_aftertrimmer",&(hbbtag_phi_aftertrimmer));
  tree->Branch("heptoptag_tmass1",&(heptoptag_tmass1));
  tree->Branch("heptoptag_bmass1",&(heptoptag_bmass1));
  tree->Branch("heptoptag_wmass1",&(heptoptag_wmass1));
  tree->Branch("heptoptag_t1_pt",&(heptoptag_t1_pt));
  tree->Branch("heptoptag_b1_pt",&(heptoptag_b1_pt));
  tree->Branch("heptoptag_w1_pt",&(heptoptag_w1_pt));
  tree->Branch("heptoptag_t1_eta",&(heptoptag_t1_eta));
  tree->Branch("heptoptag_b1_eta",&(heptoptag_b1_eta));
  tree->Branch("heptoptag_w1_eta",&(heptoptag_w1_eta));
  tree->Branch("heptoptag_t1_rap",&(heptoptag_t1_rap));
  tree->Branch("heptoptag_b1_rap",&(heptoptag_b1_rap));
  tree->Branch("heptoptag_w1_rap",&(heptoptag_w1_rap));
  tree->Branch("heptoptag_t1_phi",&(heptoptag_t1_phi));
  tree->Branch("heptoptag_b1_phi",&(heptoptag_b1_phi));
  tree->Branch("heptoptag_w1_phi",&(heptoptag_w1_phi));
  tree->Branch("heptoptag_tmass2",&(heptoptag_tmass2));
  tree->Branch("heptoptag_bmass2",&(heptoptag_bmass2));
  tree->Branch("heptoptag_wmass2",&(heptoptag_wmass2));
  tree->Branch("heptoptag_t2_pt",&(heptoptag_t2_pt));
  tree->Branch("heptoptag_b2_pt",&(heptoptag_b2_pt));
  tree->Branch("heptoptag_w2_pt",&(heptoptag_w2_pt));
  tree->Branch("heptoptag_t2_eta",&(heptoptag_t2_eta));
  tree->Branch("heptoptag_b2_eta",&(heptoptag_b2_eta));
  tree->Branch("heptoptag_w2_eta",&(heptoptag_w2_eta));
  tree->Branch("heptoptag_t2_rap",&(heptoptag_t2_rap));
  tree->Branch("heptoptag_b2_rap",&(heptoptag_b2_rap));
  tree->Branch("heptoptag_w2_rap",&(heptoptag_w2_rap));
  tree->Branch("heptoptag_t2_phi",&(heptoptag_t2_phi));
  tree->Branch("heptoptag_b2_phi",&(heptoptag_b2_phi));
  tree->Branch("heptoptag_w2_phi",&(heptoptag_w2_phi));
  tree->Branch("heptoptag_tmass_best",&(heptoptag_tmass_best));
  tree->Branch("heptoptag_bmass_best",&(heptoptag_bmass_best));
  tree->Branch("heptoptag_wmass_best",&(heptoptag_wmass_best));
  tree->Branch("heptoptag_tbest_pt",&(heptoptag_tbest_pt));
  tree->Branch("heptoptag_bbest_pt",&(heptoptag_bbest_pt));
  tree->Branch("heptoptag_wbest_pt",&(heptoptag_wbest_pt));
  tree->Branch("heptoptag_tbest_eta",&(heptoptag_tbest_eta));
  tree->Branch("heptoptag_bbest_eta",&(heptoptag_bbest_eta));
  tree->Branch("heptoptag_wbest_eta",&(heptoptag_wbest_eta));
  tree->Branch("heptoptag_tbest_rap",&(heptoptag_tbest_rap));
  tree->Branch("heptoptag_bbest_rap",&(heptoptag_bbest_rap));
  tree->Branch("heptoptag_wbest_rap",&(heptoptag_wbest_rap));
  tree->Branch("heptoptag_tbest_phi",&(heptoptag_tbest_phi));
  tree->Branch("heptoptag_bbest_phi",&(heptoptag_bbest_phi));
  tree->Branch("heptoptag_wbest_phi",&(heptoptag_wbest_phi));
  tree->Branch("heptoptag_btag_tmass1",&(heptoptag_btag_tmass1));
  tree->Branch("heptoptag_btag_bmass1",&(heptoptag_btag_bmass1));
  tree->Branch("heptoptag_btag_wmass1",&(heptoptag_btag_wmass1));
  tree->Branch("heptoptag_btag_t1_pt",&(heptoptag_btag_t1_pt));
  tree->Branch("heptoptag_btag_b1_pt",&(heptoptag_btag_b1_pt));
  tree->Branch("heptoptag_btag_w1_pt",&(heptoptag_btag_w1_pt));
  tree->Branch("heptoptag_btag_t1_eta",&(heptoptag_btag_t1_eta));
  tree->Branch("heptoptag_btag_b1_eta",&(heptoptag_btag_b1_eta));
  tree->Branch("heptoptag_btag_w1_eta",&(heptoptag_btag_w1_eta));
  tree->Branch("heptoptag_btag_t1_rap",&(heptoptag_btag_t1_rap));
  tree->Branch("heptoptag_btag_b1_rap",&(heptoptag_btag_b1_rap));
  tree->Branch("heptoptag_btag_w1_rap",&(heptoptag_btag_w1_rap));
  tree->Branch("heptoptag_btag_t1_phi",&(heptoptag_btag_t1_phi));
  tree->Branch("heptoptag_btag_b1_phi",&(heptoptag_btag_b1_phi));
  tree->Branch("heptoptag_btag_w1_phi",&(heptoptag_btag_w1_phi));
  tree->Branch("heptoptag_btag_tmass2",&(heptoptag_btag_tmass2));
  tree->Branch("heptoptag_btag_bmass2",&(heptoptag_btag_bmass2));
  tree->Branch("heptoptag_btag_wmass2",&(heptoptag_btag_wmass2));
  tree->Branch("heptoptag_btag_t2_pt",&(heptoptag_btag_t2_pt));
  tree->Branch("heptoptag_btag_b2_pt",&(heptoptag_btag_b2_pt));
  tree->Branch("heptoptag_btag_w2_pt",&(heptoptag_btag_w2_pt));
  tree->Branch("heptoptag_btag_t2_eta",&(heptoptag_btag_t2_eta));
  tree->Branch("heptoptag_btag_b2_eta",&(heptoptag_btag_b2_eta));
  tree->Branch("heptoptag_btag_w2_eta",&(heptoptag_btag_w2_eta));
  tree->Branch("heptoptag_btag_t2_rap",&(heptoptag_btag_t2_rap));
  tree->Branch("heptoptag_btag_b2_rap",&(heptoptag_btag_b2_rap));
  tree->Branch("heptoptag_btag_w2_rap",&(heptoptag_btag_w2_rap));
  tree->Branch("heptoptag_btag_t2_phi",&(heptoptag_btag_t2_phi));
  tree->Branch("heptoptag_btag_b2_phi",&(heptoptag_btag_b2_phi));
  tree->Branch("heptoptag_btag_w2_phi",&(heptoptag_btag_w2_phi));
  tree->Branch("heptoptag_btag_tmass_best",&(heptoptag_btag_tmass_best));
  tree->Branch("heptoptag_btag_bmass_best",&(heptoptag_btag_bmass_best));
  tree->Branch("heptoptag_btag_wmass_best",&(heptoptag_btag_wmass_best));
  tree->Branch("heptoptag_btag_tbest_pt",&(heptoptag_btag_tbest_pt));
  tree->Branch("heptoptag_btag_bbest_pt",&(heptoptag_btag_bbest_pt));
  tree->Branch("heptoptag_btag_wbest_pt",&(heptoptag_btag_wbest_pt));
  tree->Branch("heptoptag_btag_tbest_eta",&(heptoptag_btag_tbest_eta));
  tree->Branch("heptoptag_btag_bbest_eta",&(heptoptag_btag_bbest_eta));
  tree->Branch("heptoptag_btag_wbest_eta",&(heptoptag_btag_wbest_eta));
  tree->Branch("heptoptag_btag_tbest_rap",&(heptoptag_btag_tbest_rap));
  tree->Branch("heptoptag_btag_bbest_rap",&(heptoptag_btag_bbest_rap));
  tree->Branch("heptoptag_btag_wbest_rap",&(heptoptag_btag_wbest_rap));
  tree->Branch("heptoptag_btag_tbest_phi",&(heptoptag_btag_tbest_phi));
  tree->Branch("heptoptag_btag_bbest_phi",&(heptoptag_btag_bbest_phi));
  tree->Branch("heptoptag_btag_wbest_phi",&(heptoptag_btag_wbest_phi));
  tree->Branch("jhtoptag_tmass1",&(jhtoptag_tmass1));
  tree->Branch("jhtoptag_bmass1",&(jhtoptag_bmass1));
  tree->Branch("jhtoptag_wmass1",&(jhtoptag_wmass1));
  tree->Branch("jhtoptag_tmass2",&(jhtoptag_tmass2));
  tree->Branch("jhtoptag_bmass2",&(jhtoptag_bmass2));
  tree->Branch("jhtoptag_wmass2",&(jhtoptag_wmass2));
  tree->Branch("jhtoptag_tmass_best",&(jhtoptag_tmass_best));
  tree->Branch("jhtoptag_bmass_best",&(jhtoptag_bmass_best));
  tree->Branch("jhtoptag_wmass_best",&(jhtoptag_wmass_best));
  tree->Branch("jhtoptag_btag_tmass1",&(jhtoptag_btag_tmass1));
  tree->Branch("jhtoptag_btag_bmass1",&(jhtoptag_btag_bmass1));
  tree->Branch("jhtoptag_btag_wmass1",&(jhtoptag_btag_wmass1));
  tree->Branch("jhtoptag_btag_tmass2",&(jhtoptag_btag_tmass2));
  tree->Branch("jhtoptag_btag_bmass2",&(jhtoptag_btag_bmass2));
  tree->Branch("jhtoptag_btag_wmass2",&(jhtoptag_btag_wmass2));
  tree->Branch("jhtoptag_btag_tmass_best",&(jhtoptag_btag_tmass_best));
  tree->Branch("jhtoptag_btag_bmass_best",&(jhtoptag_btag_bmass_best));
  tree->Branch("jhtoptag_btag_wmass_best",&(jhtoptag_btag_wmass_best));
  tree->Branch("deltaR_htag_t1heptag",&(deltaR_htag_t1heptag));
  tree->Branch("deltaR_htag_t2heptag",&(deltaR_htag_t2heptag));
  tree->Branch("deltaR_t1heptag_t2heptag",&(deltaR_t1heptag_t2heptag));
  tree->Branch("deltaR_htag_tbestheptag",&(deltaR_htag_tbestheptag));
  tree->Branch("deltaRmin_httheptag",&(deltaRmin_httheptag));
  tree->Branch("deltaR_hbbtag_t1btagheptag",&(deltaR_hbbtag_t1btagheptag));
  tree->Branch("deltaR_hbbtag_t2btagheptag",&(deltaR_hbbtag_t2btagheptag));
  tree->Branch("deltaR_t1btagheptag_t2btagheptag",&(deltaR_t1btagheptag_t2btagheptag));
  tree->Branch("deltaR_hbbtag_tbestbtagheptag",&(deltaR_hbbtag_tbestbtagheptag));
  tree->Branch("deltaRmin_httbtagheptag",&(deltaRmin_httbtagheptag));
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
