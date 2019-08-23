#define myReader_cxx
//#include "ttree.h"
#include "myReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <algorithm>
#include "TLegend.h"
#include <iterator>
#include <TEfficiency.h>
#include <vector>
#include <TGraphAsymmErrors.h>

void myReader::Loop(){
  
  TFile *myfile = new TFile("Teff.root","recreate"); 
  TH1F *invM = new TH1F("invM","stAMuon from the Zinv Mass",50,0,120);
  TH1F *delta_R= new TH1F("delta_R","#stAMuons satisfy the deltaR ",50,0,120);
  //defining the bins
  const int npt_bins = 50; 
  const double pt_lo = 25;
  const double pt_hi = 125.;
  
  const int neta_bins = 50; 
  const double eta_lo = -2.5;
  const double eta_hi = 2.5;

  // for data				    
  TEfficiency *pEff = new TEfficiency( "eff","Track reconstruction effciency;p_{T};#epsilon", npt_bins, pt_lo, pt_hi);
  TH1D pT_all("pT_all","pT_all;p_{T} [GeV];#epsilon", npt_bins, pt_lo, pt_hi);
  TH1D pT_matched("pT_matched","pT_matched;p_{T};#epsilon", npt_bins, pt_lo, pt_hi);

  // for data				    
  TEfficiency *pEff_eta = new TEfficiency( "eff","Track reconstruction effciency;#eta;#epsilon", neta_bins,eta_lo, eta_hi);
  TH1D eta_all("eta_all","eta_all;#eta;#epsilon", neta_bins, eta_lo, eta_hi);
  TH1D eta_matched("eta_matched","eta_matched;#eta;#epsilon", neta_bins, eta_lo, eta_hi);

  // for simulation 
  TEfficiency *pEff_sim = new TEfficiency("eff_sim","Simulated track reconstruction efficiency;p_{T};#epsilon", npt_bins, pt_lo, pt_hi);					   
  TH1D pT_all_sim("pT_all_sim","pT_all_sim;p_{T} [GeV];#epsilon", npt_bins, pt_lo, pt_hi);
  TH1D pT_matched_sim("pT_matched_sim","pT_matched_sim;;p_{T};#epsilon", npt_bins, pt_lo, pt_hi);
  
  // for simulation 
  TEfficiency *pEff_sim_eta = new TEfficiency("eff_eta_sim","Simulated track reconstruction efficiency;#eta;#epsilon", neta_bins, eta_lo, eta_hi);					   
  TH1D eta_all_sim("eta_all_sim","eta_all_sim;#eta;#epsilon", neta_bins, eta_lo, eta_hi);
  TH1D eta_matched_sim("eta_matched_sim","eta_matched_sim;#eta;#epsilon", neta_bins, eta_lo, eta_hi);
  
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
	eta_all.Fill(tree_staMuon_eta->at(j));
	
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
	  eta_matched.Fill(tree_staMuon_eta->at(j));
	}	
      } // inv mass cut 
    }
	  
    for (unsigned int f = 0 ; f < tree_simtrack_simtrack_pt ->size(); f++) { 
      if (tree_simtrack_simtrack_pt->at(f) < 25.) continue;

      TLorentzVector track1,track2;

      /*
      if(tree_simtrack_isRecoMatched->at(f)== true) {
	pT_matched_sim.Fill(tree_simtrack_simtrack_pt->at(f));
	eta_matched_sim.Fill(tree_simtrack_simtrack_eta->at(f));
      }
      */

      for (unsigned int g = 0 ; g < tree_simtrack_simtrack_pt ->size(); g++) { 

      	track1.SetPtEtaPhiM( tree_simtrack_simtrack_pt->at(f), tree_simtrack_simtrack_eta->at(f), tree_simtrack_simtrack_phi->at(f), 0);
      	track2.SetPtEtaPhiM( tree_simtrack_simtrack_pt->at(g), tree_simtrack_simtrack_eta->at(g), tree_simtrack_simtrack_phi->at(g), 0);
	
      	Zboson_sim = track1+track2;
      	theInvariantMass_sim =Zboson_sim.M();
   		
      	if(theInvariantMass_sim >=86 && theInvariantMass_sim  <= 101){

	  pT_all_sim.Fill(tree_simtrack_simtrack_pt->at(f));
	  eta_all_sim.Fill(tree_simtrack_simtrack_eta->at(f));
	  
      	  if(tree_simtrack_isRecoMatched->at(f)== true) {
      	    check_passing ++ ;
      	    pT_matched_sim.Fill(tree_simtrack_simtrack_pt->at(f));
	    eta_matched_sim.Fill(tree_simtrack_simtrack_eta->at(f));
      	    break; // to stop the loop once a match is found
      	  }
      	}
      }
    }// end of simulation loop	   
  }//loop over tree simulation  

  //TEfficiency calculation for data
  if(TEfficiency::CheckConsistency(pT_matched, pT_all)){
    pEff = new TEfficiency(pT_matched, pT_all);    
    //pEff->Write();
  }
  
  //TEfficiency calculation for simulation
  if(TEfficiency::CheckConsistency(pT_matched_sim, pT_all_sim)){
    pEff_sim = new TEfficiency(pT_matched_sim, pT_all_sim);
    //pEff_sim->Write() ;
  }

  //TEfficiency calculation for data
  if(TEfficiency::CheckConsistency(eta_matched, eta_all)){
    pEff_eta = new TEfficiency(eta_matched, eta_all);    
    //pEff->Write();
  }
  
  //TEfficiency calculation for simulation
  if(TEfficiency::CheckConsistency(eta_matched_sim, eta_all_sim)){
    pEff_sim_eta = new TEfficiency(eta_matched_sim, eta_all_sim);
    //pEff_sim->Write() ;
  }

 
  // Draw and superimpose

  TCanvas *c1 = new TCanvas("c1", "c1",800,600);
  TCanvas *c2 = new TCanvas("c2", "c2",800,600);  

  
  c1->cd();
  pEff->SetLineColor(kBlue);
  pEff_sim->SetLineColor(kRed);
  pEff->SetMarkerSize(1.);
  pEff_sim->SetMarkerSize(1.);
  pEff->SetMarkerColor(kBlue);
  pEff_sim->SetMarkerColor(kRed);
  pEff->SetMarkerStyle(20);
  pEff_sim->SetMarkerStyle(21);
  pEff->Draw("AP");
  auto graph = pEff->GetPaintedGraph();
  //graph->SetMinimum(0.75);
  //graph->SetMaximum(1.); 
  //c2->Update();
  
  pEff_sim->Draw("same");

  auto legend = new TLegend(0.1,0.1,0.48,0.3);
  legend->AddEntry(pEff,"Reconstructed efficiency","P");
  legend->AddEntry(pEff_sim,"Simulated efficiency","P");
  legend->Draw("same");


  
  //Draw and superimpose
  c2->cd();
  pEff_eta->SetLineColor(kBlue);
  pEff_sim_eta->SetLineColor(kRed);
  pEff_eta->SetMarkerSize(1.);
  pEff_sim_eta->SetMarkerSize(1.);
  pEff_eta->SetMarkerColor(kBlue);
  pEff_sim_eta->SetMarkerColor(kRed);
  pEff_eta->SetMarkerStyle(20);
  pEff_sim_eta->SetMarkerStyle(21);
  //pEff_eta->GetYAxis()->SetRangeUser(0.75,1.);
  pEff_eta->Draw("AP");
  auto graph1 = pEff_eta->GetPaintedGraph(); 
  //graph1->SetMinimum(0.75);
  //graph1->SetMaximum(1.); 
  c2->Update(); 
  
  pEff_sim_eta->Draw("same");

  //std::cout << eta_all.GetEntries() << std::endl;
  //std::cout << eta_matched.GetEntries() << std::endl;
  auto legend2 = new TLegend(0.1,0.1,0.48,0.3);
  legend2->AddEntry(pEff_eta,"Reconstructed efficiency","P");
  legend2->AddEntry(pEff_sim_eta,"Simulated efficiency","P");
  legend2->Draw("same");

  
} // end of void::Loop()
