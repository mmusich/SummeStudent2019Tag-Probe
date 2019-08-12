#define myReader_cxx
#include "myReader.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <iostream>
#include <TCanvas.h>

// void myReader::Loop()
// {
//   //   In a ROOT session, you can do:
// //      root> .L myReader.C
// //      root> myReader t
// //      root> t.GetEntry(12); // Fill t data members with entry number 12
// //      root> t.Show();       // Show values of entry 12
// //      root> t.Show(16);     // Read and show values of entry 16
// //      root> t.Loop();       // Loop on all entries
// //

// //     This is the loop skeleton where:
// //    jentry is the global entry number in the chain
// //    ientry is the entry number in the current Tree
// //  Note that the argument to GetEntry must be:
// //    jentry for TChain::GetEntry
// //    ientry for TTree::GetEntry and TBranch::GetEntry
// //
// //       To read only selected branches, Insert statements like:
// // METHOD1:
// //    fChain->SetBranchStatus("*",0);  // disable all branches
// //    fChain->SetBranchStatus("branchname",1);  // activate branchname
// // METHOD2: replace line
// //    fChain->GetEntry(jentry);       //read all branches
// //by  b_branchname->GetEntry(ientry); //read only this branch
//    if (fChain == 0) return;

//    Long64_t nentries = fChain->GetEntriesFast();

//    Long64_t nbytes = 0, nb = 0;
//    for (Long64_t jentry=0; jentry<nentries;jentry++) {
//       Long64_t ientry = LoadTree(jentry);
//       if (ientry < 0) break;
//       nb = fChain->GetEntry(jentry);   nbytes += nb;
//       // if (Cut(ientry) < 0) continue;
//    }
// }


void myReader::Loop()
{

   TFile *myfile = new TFile("invM.root","recreate");
   TH1F *invM = new TH1F("invM","Invariant Mass of Z costructed from T&P Muons",100,0,180);
   TH1F *delta_R= new TH1F("delta_R","#Delta R of probe muons and tracks ",100,0,6);
   TCanvas *c1 = new TCanvas("c1", "Z boson invM",800,600);
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast(); 
   Long64_t nbytes = 0, nb = 0;


   double dR;
   TLorentzVector gMuon , sMuon, Zboson, track ;
   double deltaR , theInvariantMass , nMuon_ZWindows , nMuon_dR ;
   
   double numberofZ=0 ;
   double numberofZ_withMatchedProbe = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry); 
     
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     TLorentzVector gMuon , sMuon, Zboson, track ;
     double deltaR , theInvariantMass  ;
     int nMuon_ZWindows , nMuon_dR ;
  
     // choose the highest pT Global Muon from each vector 
 
     if ( tree_globalMuon_pt->size() == 0) continue;
     int MaxElementIndex = std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end()) - tree_globalMuon_pt->begin();
     double MaxElement = *std::max_element(tree_globalMuon_pt->begin(), tree_globalMuon_pt->end());
     gMuon.SetPtEtaPhiM( (*tree_globalMuon_pt)[MaxElementIndex],  (*tree_globalMuon_eta)[MaxElementIndex],  (*tree_globalMuon_phi)[MaxElementIndex] , 0);
     // loop the event with ALL stA Muon
     for ( unsigned int j = 0 ;  j < tree_staMuon_pt ->size()  ; j++){
       if ( (*tree_staMuon_pt)[j]< 15.) continue; 
       sMuon.SetPtEtaPhiM( tree_staMuon_pt->at(j), tree_staMuon_eta->at(j), tree_staMuon_phi->at(j), 0);
       deltaR = gMuon.DeltaR(sMuon);
       if (deltaR >0.5) {
	 Zboson = gMuon + sMuon;
	 theInvariantMass =Zboson.M();
       }
       if (  theInvariantMass >=86 &&  theInvariantMass  <= 101 ) {numberofZ = numberofZ+1;};
       break;
     }
     
     //  select the stA from the Z mass windows and associate them with track such that dR < 0.1
     if (  theInvariantMass >=86 &&  theInvariantMass  <= 101 ) {
       for (unsigned int k =0 ;k< tree_track_pt->size() ;k++) {
	 if((*tree_track_pt)[k] < 15.)  continue;
	 track.SetPtEtaPhiM( (*tree_track_pt)[k], (*tree_track_eta)[k], (*tree_track_phi)[k],0);
	 dR= sMuon.DeltaR(track);
	 if(dR< 0.1)  {numberofZ_withMatchedProbe = numberofZ_withMatchedProbe+1;};
       } 
     }
   } 
   cout << numberofZ << ";" << numberofZ_withMatchedProbe << endl;
   double efficiency = numberofZ_withMatchedProbe /numberofZ;
   cout << "Efficiency : " << efficiency << endl;
} 
