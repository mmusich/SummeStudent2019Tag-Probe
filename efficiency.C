#define ttree_cxx
#include "ttree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include<algorithm>
#include<iterator>
#include<TEfficiency.h>


void ttree::Loop()
{
  
   TFile *myfile = new TFile("Teff.root","recreate"); 
   TH1F *invM = new TH1F("invM","stAMuon from the Zinv Mass",50,0,120);
   TH1F *delta_R= new TH1F("delta_R","#stAMuons satisfy the deltaR ",50,0,120);
   TCanvas *c1 = new TCanvas("c1", "Z boson invM",800,600);
   TEfficiency * pEff = new TEfficiency( "eff","Track reconstruction effciency;p_{T};#epsilon", 50, 0 ,120);
   
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast(); 
   Long64_t nbytes = 0, nb = 0;

   //the variables
	int nevents = 0; // initialization
	double dR;
	TLorentzVector gMuon , sMuon, Zboson, track ;
	double deltaR , theInvariantMass , nMuon_ZWindows , nMuon_dR ;
	double pass =0;
	double total = 0;
	bool bpassed;
	double numberofZ=0 ;
	double numberofZ_withMatchedProbe = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry); 
	 
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
     
	  TLorentzVector gMuon , sMuon, Zboson, track ;
	  double deltaR , theInvariantMass  ;
	  int nMuon_ZWindows , nMuon_dR ;
	  
	    // choose the highest pT Global Muon from each vector, and loop the big vector 
	     if ( tree_globalMuon_pt->size() == 0) continue;
			int MaxElementIndex = std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end()) - tree_globalMuon_pt->begin();
			double MaxElement = *std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end());
			gMuon.SetPtEtaPhiM( (*tree_globalMuon_pt)[MaxElementIndex],  (*tree_globalMuon_eta)[MaxElementIndex],  (*tree_globalMuon_phi)[MaxElementIndex] , 0);
			
		// loop the event with ALL stA Muon
		for ( unsigned int j = 0 ;  j < tree_staMuon_pt ->size()  ; j++)
		{
			if ( (*tree_staMuon_pt)[j]< 15.) continue; 
				sMuon.SetPtEtaPhiM( tree_staMuon_pt->at(j), tree_staMuon_eta->at(j), tree_staMuon_phi->at(j), 0);
				deltaR = gMuon.DeltaR(sMuon);
				
			if (deltaR >0.5)
			{
				Zboson = gMuon + sMuon;
				theInvariantMass =Zboson.M();
			}
			
			//select the stAMuons that is within the range of Z inv Mass above, define as total
			if (  theInvariantMass >=86 &&  theInvariantMass  <= 101 ) {total =(*tree_staMuon_pt)[j];};
			break;
		}
		
		//  select the stA from the Z mass windows and associate them with track such that dR < 0.1
		if (theInvariantMass >=86 &&  theInvariantMass  <= 101 )
		{	
			for (unsigned int k =0 ;k< tree_track_pt->size() ;k++)
			{
				if((*tree_track_pt)[k] < 15.)  continue;
				track.SetPtEtaPhiM( (*tree_track_pt)[k], (*tree_track_eta)[k], (*tree_track_phi)[k],0);
				dR= sMuon.DeltaR(track);
				
				// select the staMuon that has delta_R of <0.1 with track from the total stA above
				if(dR< 0.1)  {pass = sMuon.Pt();};
				break;
			}	
		}	
		
		// construct the efficiency
		bpassed = pass < total;
		pEff->Fill(bpassed, sMuon.Pt());
	
	}	
		
		
		
		myfile->WriteTObject(c1);
		myfile->Draw();
		pEff->Draw();
		
}	
	
