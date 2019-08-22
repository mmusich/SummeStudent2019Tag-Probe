#define myReader_cxx
#include "ttree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <algorithm>
#include <iterator>
#include <TEfficiency.h>
#include <vector>
#include <TGraphAsymmErrors.h>


void myReader::Loop(){
  
	TFile *myfile = new TFile("Teff.root","recreate"); 
	TH1F *invM = new TH1F("invM","stAMuon from the Zinv Mass",50,0,120);
	TH1F *delta_R= new TH1F("delta_R","#stAMuons satisfy the deltaR ",50,0,120);
	TCanvas *c1 = new TCanvas("c1", "Z boson invM",800,600);
  
	//defining the bins
	const int npt_bins = 100; 
	const double pt_lo = 0;
	const double pt_hi = 120.;

  
	// for data				    
	TEfficiency *pEff = new TEfficiency( "eff","Track reconstruction effciency;p_{T};#epsilon", npt_bins, pt_lo, pt_hi);
	TH1D pT_all("pT_all","pT_all", npt_bins, pt_lo, pt_hi);
	TH1D pT_matched("pT_matched","pT_matched", npt_bins, pt_lo, pt_hi);
  
	// for simulation 
	TEfficiency *pEff_sim = new TEfficiency( "eff_sim","Simulated track reconstruction effciency;p_{T};#epsilon", npt_bins, pt_lo, pt_hi);					   
	TH1D pT_all_sim("pT_all_sim","pT_all_sim", npt_bins, pt_lo, pt_hi);
	TH1D pT_matched_sim("pT_matched_sim","pT_matched_sim", npt_bins, pt_lo, pt_hi);
  
	if (fChain == 0) return;
  
	Long64_t nentries = fChain->GetEntriesFast();   
	Long64_t nbytes = 0, nb = 0;
  
	for (Long64_t jentry=0; jentry <nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
    //if (Cut(ientry) < 0) continue;
     
    TLorentzVector gMuon, sMuon, Zboson, Zboson_sim, track, track_sim;
    double deltaR, theInvariantMass, theInvariantMass_sim;
    int nMuon_ZWindows, nMuon_dR;
	int check_passing, check_total;
    
    // choose the highest pT Global Muon from each vector, and loop the big vector 
    
    if (tree_globalMuon_pt->size() == 0) 
      continue;
    
    int MaxElementIndex = std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end()) - tree_globalMuon_pt->begin();
    double MaxElement = *std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end());
    
    gMuon.SetPtEtaPhiM( (*tree_globalMuon_pt)[MaxElementIndex],  (*tree_globalMuon_eta)[MaxElementIndex],  (*tree_globalMuon_phi)[MaxElementIndex] , 0);
    
    if (tree_globalMuon_pt->at(MaxElementIndex) < 25)
      continue;
	check_total++ ;
    pT_all_sim.Fill(gMuon.Pt());
	
    // loop the event with ALL stA Muon
    
    for (unsigned int j = 0 ;  j < tree_staMuon_pt ->size(); j++) {      
      if (tree_staMuon_pt->at(j) < 15.) 
	continue; 
      
      sMuon.SetPtEtaPhiM( tree_staMuon_pt->at(j), tree_staMuon_eta->at(j), tree_staMuon_phi->at(j), 0);
      deltaR = gMuon.DeltaR(sMuon);
      Zboson = gMuon + sMuon;
      double theInvariantMass = Zboson.M();
      
      if (deltaR < 0.5)
	continue;
      
      //select the stAMuons that is within the range of Z inv Mass above, define as total
      if (theInvariantMass >=86 && theInvariantMass  <= 101) {
	
	pT_all.Fill(tree_staMuon_pt->at(j));
	
	int index_dRmin = -1;
	double dR_min = 999;
	
	for (unsigned int k =0 ; k<tree_track_pt->size() ;k++){
	  if((*tree_track_pt)[k] < 15.)  
	    continue;
	  track.SetPtEtaPhiM( (*tree_track_pt)[k], (*tree_track_eta)[k], (*tree_track_phi)[k],0);
	  double dR = sMuon.DeltaR(track);	     
	  if (dR < dR_min){
	    dR_min = dR;
	    index_dRmin = k;
	  }
	}
	
	if (dR_min < 0.1){
	  pT_matched.Fill(tree_staMuon_pt->at(j));	  
	}	
      } // inv mass cut 
    }
	
	
	
    
    // for (unsigned int f = 0 ; f < tree_simtrack_simtrack_pt ->size(); f++) { 
      if (tree_simtrack_simtrack_pt->at(f) < 25.) continue; 
    
	    TLorentzVector track1,track2;
		if (tree_simtrack_simtrack_pt->size() == 2){
		track1.SetPtEtaPhiM( tree_simtrack_simtrack_pt->at(0), tree_simtrack_simtrack_eta->at(0), tree_simtrack_simtrack_phi->at(0), 0);
		track2.SetPtEtaPhiM( tree_simtrack_simtrack_pt->at(1), tree_simtrack_simtrack_eta->at(1), tree_simtrack_simtrack_phi->at(1), 0);
		
		Zboson_sim = track1+track2;
		theInvariantMass_sim =Zboson_sim.M();
		}
		
		if(theInvariantMass_sim >=86 && theInvariantMass_sim  <= 101){
			if(tree_simtrack_isRecoMatched->at(f)== true) {
				check_passing ++ ;
			}
		}
		
	}// end of simulation loop	   
    }//loop over tree simulation
    
		
    
    }//loop over tree
  
	//TEfficiency calculation for data
	if(TEfficiency::CheckConsistency(pT_matched, pT_all)){
    pEff = new TEfficiency(pT_matched, pT_all);    
    pEff->Write();
	}
	
	//TEfficiency calculation for simulation
	if(TEfficiency::CheckConsistency(pT_matched_sim, pT_all_sim)){
    pEff_sim = new TEfficiency(pT_matched_sim, pT_all_sim);
	pEff_sim->Write() ;
    }
	
	// Draw and superimpose
	pEff->Draw();
	pEff->SetLineColor(4);
	pEff_sim->Draw("same");
	
	
	
	
  
 
} // end of void::Loop()
