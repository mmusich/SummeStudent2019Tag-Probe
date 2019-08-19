#define treeReader_cxx

#include "treeReader.h"
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>
#include <iostream>


void treeReader::Loop(){
  
  TFile *myfile = new TFile("Teff.root","recreate"); 
  TH1F *invM = new TH1F("invM","stAMuon from the Zinv Mass",50,0,120);
  TH1F *delta_R= new TH1F("delta_R","#stAMuons satisfy the deltaR ",50,0,120);
  TCanvas *c1 = new TCanvas("c1", "Z boson invM",800,600);
  
  const int npt_bins = 120; 
  const double pt_lo = 0.;
  const double pt_hi = 120;
  
  TEfficiency *pEff = new TEfficiency( "eff","Track reconstruction effciency;p_{T};#epsilon", 
				       npt_bins, pt_lo, pt_hi);
  TH1D pT_all("pT_all","pT_all", npt_bins, pt_lo, pt_hi);
  TH1D pT_matched("pT_matched","pT_matched", npt_bins, pt_lo, pt_hi);
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();   
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry <nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
    //if (Cut(ientry) < 0) continue;
     
    TLorentzVector gMuon, sMuon, Zboson, track;
    double deltaR, theInvariantMass;
    int nMuon_ZWindows, nMuon_dR;
    
    // choose the highest pT Global Muon from each vector, and loop the big vector 
    
    if (tree_globalMuon_pt->size() == 0) 
      continue;
    
    int MaxElementIndex = std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end()) - tree_globalMuon_pt->begin();
    double MaxElement = *std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end());
    
    gMuon.SetPtEtaPhiM( (*tree_globalMuon_pt)[MaxElementIndex],  (*tree_globalMuon_eta)[MaxElementIndex],  (*tree_globalMuon_phi)[MaxElementIndex] , 0);
    
    if (tree_globalMuon_pt->at(MaxElementIndex) < 25)
      continue;
    
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
  } // loop over tree
  
  myfile->cd();    
  if(TEfficiency::CheckConsistency(pT_matched, pT_all)){
    pEff = new TEfficiency(pT_matched, pT_all);    
    pEff->Write();
  }  
}
