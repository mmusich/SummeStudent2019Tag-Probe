//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 12 11:13:56 2019 by ROOT version 6.08/04
// from TTree ttree/ttree
// found on file: trackingNTuple.root
//////////////////////////////////////////////////////////

#ifndef myReader_h
#define myReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class myReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *tree_track_pt;
   vector<double>  *tree_track_outerPt;
   vector<double>  *tree_track_eta;
   vector<double>  *tree_track_phi;
   vector<int>     *tree_track_nhits;
   vector<double>  *tree_track_NChi2;
   vector<double>  *tree_track_isHighQuality;
   vector<double>  *tree_track_isLoose;
   vector<double>  *tree_track_isTight;
   vector<double>  *tree_track_dxy;
   vector<double>  *tree_track_dxyError;
   vector<double>  *tree_track_dz;
   vector<double>  *tree_track_dzError;
   vector<int>     *tree_track_numberOfLostHits;
   vector<int>     *tree_track_numberOfValidHits;
   vector<unsigned int> *tree_track_originalAlgo;
   vector<unsigned int> *tree_track_algo;
   vector<unsigned short> *tree_track_stopReason;
   vector<bool>    *tree_track_isSimMatched;
   vector<int>     *tree_track_numberOfValidPixelHits;
   vector<int>     *tree_track_numberOfValidStripHits;
   vector<int>     *tree_track_numberOfValidStripTIBHits;
   vector<int>     *tree_track_numberOfValidStripTIDHits;
   vector<int>     *tree_track_numberOfValidStripTOBHits;
   vector<int>     *tree_track_numberOfValidStripTECHits;
   vector<int>     *tree_track_numberOfValidPixelBarrelHits;
   vector<int>     *tree_track_numberOfValidPixelEndcapHits;
   vector<int>     *tree_track_hasValidHitInPixelLayer;
   vector<int>     *tree_track_trackerLayersWithMeasurement;
   vector<int>     *tree_track_pixelLayersWithMeasurement;
   vector<int>     *tree_track_stripTECLayersWithMeasurement;
   vector<int>     *tree_track_stripTIBLayersWithMeasurement;
   vector<int>     *tree_track_stripTIDLayersWithMeasurement;
   vector<int>     *tree_track_stripTOBLayersWithMeasurement;
   vector<int>     *tree_track_nPixel;
   vector<int>     *tree_track_nStrip;
   vector<double>  *tree_track_vx;
   vector<double>  *tree_track_vy;
   vector<double>  *tree_track_vz;
   vector<double>  *tree_track_firsthit_X;
   vector<double>  *tree_track_firsthit_Y;
   vector<double>  *tree_track_firsthit_Z;
   vector<double>  *tree_track_firsthit_phi;
   vector<double>  *tree_track_simtrack_charge;
   vector<double>  *tree_track_simtrack_pt;
   vector<double>  *tree_track_simtrack_eta;
   vector<double>  *tree_track_simtrack_phi;
   vector<bool>    *tree_track_simtrack_longLived;
   vector<int>     *tree_track_simtrack_pdgId;
   vector<int>     *tree_track_simtrack_numberOfTrackerHits;
   vector<int>     *tree_track_simtrack_numberOfTrackerLayers;
   vector<double>  *tree_track_simtrack_mass;
   vector<int>     *tree_track_simtrack_status;
   vector<double>  *tree_track_genVertexPos_X;
   vector<double>  *tree_track_genVertexPos_Y;
   vector<double>  *tree_track_genVertexPos_Z;
   vector<int>     *tree_track_recoVertex_idx;
   vector<int>     *tree_track_recoJet_idx;
   vector<int>     *tree_SiCluster_subDet;
   vector<int>     *tree_SiCluster_PetalSide;
   vector<int>     *tree_SiCluster_LayerNbr;
   vector<int>     *tree_SiCluster_WheelSide;
   vector<int>     *tree_PixCluster_isBPix;
   vector<int>     *tree_PixCluster_LayerNbr;
   vector<bool>    *tree_PixCluster_detID;
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiBlock;
   Float_t         tree_bs_PosX;
   Float_t         tree_bs_PosY;
   Float_t         tree_bs_PosZ;
   vector<float>   *tree_vtx_PosX;
   vector<float>   *tree_vtx_PosY;
   vector<float>   *tree_vtx_PosZ;
   vector<float>   *tree_vtx_PosXError;
   vector<float>   *tree_vtx_PosYError;
   vector<float>   *tree_vtx_PosZError;
   vector<double>  *tree_simtrack_simtrack_charge;
   vector<double>  *tree_simtrack_simtrack_pt;
   vector<double>  *tree_simtrack_simtrack_eta;
   vector<double>  *tree_simtrack_simtrack_phi;
   vector<bool>    *tree_simtrack_simtrack_longLived;
   vector<int>     *tree_simtrack_simtrack_pdgId;
   vector<int>     *tree_simtrack_simtrack_numberOfTrackerHits;
   vector<int>     *tree_simtrack_simtrack_numberOfTrackerLayers;
   vector<double>  *tree_simtrack_simtrack_mass;
   vector<int>     *tree_simtrack_simtrack_status;
   vector<double>  *tree_simtrack_genVertexPos_X;
   vector<double>  *tree_simtrack_genVertexPos_Y;
   vector<double>  *tree_simtrack_genVertexPos_Z;
   vector<bool>    *tree_simtrack_isRecoMatched;
   vector<double>  *tree_simtrack_pca_dxy;
   vector<double>  *tree_simtrack_pca_dz;
   vector<vector<int> > *tree_simtrack_trkIdx;
   vector<float>   *tree_jet_pt;
   vector<float>   *tree_jet_eta;
   vector<float>   *tree_jet_phi;
   vector<double>  *tree_gsftrack_pt;
   vector<double>  *tree_gsftrack_eta;
   vector<double>  *tree_gsftrack_phi;
   vector<int>     *tree_gsftrack_nhits;
   vector<double>  *tree_gsftrack_outerPt;
   vector<double>  *tree_gsftrack_simtrack_pt;
   vector<double>  *tree_gsftrack_simtrack_eta;
   vector<double>  *tree_gsftrack_simtrack_phi;
   vector<double>  *tree_gsftrack_simtrack_nhits;
   vector<double>  *tree_gsftrack_simtrack_outerPt;
   vector<double>  *tree_globalMuon_charge;
   vector<double>  *tree_globalMuon_pt;
   vector<double>  *tree_globalMuon_eta;
   vector<double>  *tree_globalMuon_phi;
   vector<double>  *tree_globalMuon_nChi2;
   vector<double>  *tree_globalMuon_nHit;
   vector<int>     *tree_globalMuon_idXtrack;
   vector<double>  *tree_staMuon_charge;
   vector<double>  *tree_staMuon_pt;
   vector<double>  *tree_staMuon_eta;
   vector<double>  *tree_staMuon_phi;
   vector<double>  *tree_staMuon_nChi2;
   vector<double>  *tree_staMuon_nHit;

   // List of branches
   TBranch        *b_tree_track_pt;   //!
   TBranch        *b_tree_track_outerPt;   //!
   TBranch        *b_tree_track_eta;   //!
   TBranch        *b_tree_track_phi;   //!
   TBranch        *b_tree_track_nhits;   //!
   TBranch        *b_tree_track_NChi2;   //!
   TBranch        *b_tree_track_isHighQuality;   //!
   TBranch        *b_tree_track_isLoose;   //!
   TBranch        *b_tree_track_isTight;   //!
   TBranch        *b_tree_track_dxy;   //!
   TBranch        *b_tree_track_dxyError;   //!
   TBranch        *b_tree_track_dz;   //!
   TBranch        *b_tree_track_dzError;   //!
   TBranch        *b_tree_track_numberOfLostHits;   //!
   TBranch        *b_tree_track_numberOfValidHits;   //!
   TBranch        *b_tree_track_originalAlgo;   //!
   TBranch        *b_tree_track_algo;   //!
   TBranch        *b_tree_track_stopReason;   //!
   TBranch        *b_tree_track_isSimMatched;   //!
   TBranch        *b_tree_track_numberOfValidPixelHits;   //!
   TBranch        *b_tree_track_numberOfValidStripHits;   //!
   TBranch        *b_tree_track_numberOfValidStripTIBHits;   //!
   TBranch        *b_tree_track_numberOfValidStripTIDHits;   //!
   TBranch        *b_tree_track_numberOfValidStripTOBHits;   //!
   TBranch        *b_tree_track_numberOfValidStripTECHits;   //!
   TBranch        *b_tree_track_numberOfValidPixelBarrelHits;   //!
   TBranch        *b_tree_track_numberOfValidPixelEndcapHits;   //!
   TBranch        *b_tree_track_hasValidHitInPixelLayer;   //!
   TBranch        *b_tree_track_trackerLayersWithMeasurement;   //!
   TBranch        *b_tree_track_pixelLayersWithMeasurement;   //!
   TBranch        *b_tree_track_stripTECLayersWithMeasurement;   //!
   TBranch        *b_tree_track_stripTIBLayersWithMeasurement;   //!
   TBranch        *b_tree_track_stripTIDLayersWithMeasurement;   //!
   TBranch        *b_tree_track_stripTOBLayersWithMeasurement;   //!
   TBranch        *b_tree_track_nPixel;   //!
   TBranch        *b_tree_track_nStrip;   //!
   TBranch        *b_tree_track_vx;   //!
   TBranch        *b_tree_track_vy;   //!
   TBranch        *b_tree_track_vz;   //!
   TBranch        *b_tree_track_firsthit_X;   //!
   TBranch        *b_tree_track_firsthit_Y;   //!
   TBranch        *b_tree_track_firsthit_Z;   //!
   TBranch        *b_tree_track_firsthit_phi;   //!
   TBranch        *b_tree_track_simtrack_charge;   //!
   TBranch        *b_tree_track_simtrack_pt;   //!
   TBranch        *b_tree_track_simtrack_eta;   //!
   TBranch        *b_tree_track_simtrack_phi;   //!
   TBranch        *b_tree_track_simtrack_longLived;   //!
   TBranch        *b_tree_track_simtrack_pdgId;   //!
   TBranch        *b_tree_track_simtrack_numberOfTrackerHits;   //!
   TBranch        *b_tree_track_simtrack_numberOfTrackerLayers;   //!
   TBranch        *b_tree_track_simtrack_mass;   //!
   TBranch        *b_tree_track_simtrack_status;   //!
   TBranch        *b_tree_track_genVertexPos_X;   //!
   TBranch        *b_tree_track_genVertexPos_Y;   //!
   TBranch        *b_tree_track_genVertexPos_Z;   //!
   TBranch        *b_tree_track_recoVertex_idx;   //!
   TBranch        *b_tree_track_recoJet_idx;   //!
   TBranch        *b_tree_SiCluster_subDet;   //!
   TBranch        *b_tree_SiCluster_PetalSide;   //!
   TBranch        *b_tree_SiCluster_LayerNbr;   //!
   TBranch        *b_tree_SiCluster_WheelSide;   //!
   TBranch        *b_tree_PixCluster_isBPix;   //!
   TBranch        *b_tree_PixCluster_LayerNbr;   //!
   TBranch        *b_tree_PixCluster_detID;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_tree_bs_PosX;   //!
   TBranch        *b_tree_bs_PosY;   //!
   TBranch        *b_tree_bs_PosZ;   //!
   TBranch        *b_tree_vtx_PosX;   //!
   TBranch        *b_tree_vtx_PosY;   //!
   TBranch        *b_tree_vtx_PosZ;   //!
   TBranch        *b_tree_vtx_PosXError;   //!
   TBranch        *b_tree_vtx_PosYError;   //!
   TBranch        *b_tree_vtx_PosZError;   //!
   TBranch        *b_tree_simtrack_simtrack_charge;   //!
   TBranch        *b_tree_simtrack_simtrack_pt;   //!
   TBranch        *b_tree_simtrack_simtrack_eta;   //!
   TBranch        *b_tree_simtrack_simtrack_phi;   //!
   TBranch        *b_tree_simtrack_simtrack_longLived;   //!
   TBranch        *b_tree_simtrack_simtrack_pdgId;   //!
   TBranch        *b_tree_simtrack_simtrack_numberOfTrackerHits;   //!
   TBranch        *b_tree_simtrack_simtrack_numberOfTrackerLayers;   //!
   TBranch        *b_tree_simtrack_simtrack_mass;   //!
   TBranch        *b_tree_simtrack_simtrack_status;   //!
   TBranch        *b_tree_simtrack_genVertexPos_X;   //!
   TBranch        *b_tree_simtrack_genVertexPos_Y;   //!
   TBranch        *b_tree_simtrack_genVertexPos_Z;   //!
   TBranch        *b_tree_simtrack_isRecoMatched;   //!
   TBranch        *b_tree_simtrack_pca_dxy;   //!
   TBranch        *b_tree_simtrack_pca_dz;   //!
   TBranch        *b_tree_simtrack_trkIdx;   //!
   TBranch        *b_tree_jet_pt;   //!
   TBranch        *b_tree_jet_eta;   //!
   TBranch        *b_tree_jet_phi;   //!
   TBranch        *b_tree_gsftrack_pt;   //!
   TBranch        *b_tree_gsftrack_eta;   //!
   TBranch        *b_tree_gsftrack_phi;   //!
   TBranch        *b_tree_gsftrack_nhits;   //!
   TBranch        *b_tree_gsftrack_outerPt;   //!
   TBranch        *b_tree_gsftrack_simtrack_pt;   //!
   TBranch        *b_tree_gsftrack_simtrack_eta;   //!
   TBranch        *b_tree_gsftrack_simtrack_phi;   //!
   TBranch        *b_tree_gsftrack_simtrack_nhits;   //!
   TBranch        *b_tree_gsftrack_simtrack_outerPt;   //!
   TBranch        *b_tree_globalMuon_charge;   //!
   TBranch        *b_tree_globalMuon_pt;   //!
   TBranch        *b_tree_globalMuon_eta;   //!
   TBranch        *b_tree_globalMuon_phi;   //!
   TBranch        *b_tree_globalMuon_nChi2;   //!
   TBranch        *b_tree_globalMuon_nHit;   //!
   TBranch        *b_tree_globalMuon_idXtrack;   //!
   TBranch        *b_tree_staMuon_charge;   //!
   TBranch        *b_tree_staMuon_pt;   //!
   TBranch        *b_tree_staMuon_eta;   //!
   TBranch        *b_tree_staMuon_phi;   //!
   TBranch        *b_tree_staMuon_nChi2;   //!
   TBranch        *b_tree_staMuon_nHit;   //!

   myReader(TTree *tree=0);
   virtual ~myReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myReader_cxx
myReader::myReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("trackingNTuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("trackingNTuple.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("trackingNTuple.root:/trackingPerf");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

myReader::~myReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   tree_track_pt = 0;
   tree_track_outerPt = 0;
   tree_track_eta = 0;
   tree_track_phi = 0;
   tree_track_nhits = 0;
   tree_track_NChi2 = 0;
   tree_track_isHighQuality = 0;
   tree_track_isLoose = 0;
   tree_track_isTight = 0;
   tree_track_dxy = 0;
   tree_track_dxyError = 0;
   tree_track_dz = 0;
   tree_track_dzError = 0;
   tree_track_numberOfLostHits = 0;
   tree_track_numberOfValidHits = 0;
   tree_track_originalAlgo = 0;
   tree_track_algo = 0;
   tree_track_stopReason = 0;
   tree_track_isSimMatched = 0;
   tree_track_numberOfValidPixelHits = 0;
   tree_track_numberOfValidStripHits = 0;
   tree_track_numberOfValidStripTIBHits = 0;
   tree_track_numberOfValidStripTIDHits = 0;
   tree_track_numberOfValidStripTOBHits = 0;
   tree_track_numberOfValidStripTECHits = 0;
   tree_track_numberOfValidPixelBarrelHits = 0;
   tree_track_numberOfValidPixelEndcapHits = 0;
   tree_track_hasValidHitInPixelLayer = 0;
   tree_track_trackerLayersWithMeasurement = 0;
   tree_track_pixelLayersWithMeasurement = 0;
   tree_track_stripTECLayersWithMeasurement = 0;
   tree_track_stripTIBLayersWithMeasurement = 0;
   tree_track_stripTIDLayersWithMeasurement = 0;
   tree_track_stripTOBLayersWithMeasurement = 0;
   tree_track_nPixel = 0;
   tree_track_nStrip = 0;
   tree_track_vx = 0;
   tree_track_vy = 0;
   tree_track_vz = 0;
   tree_track_firsthit_X = 0;
   tree_track_firsthit_Y = 0;
   tree_track_firsthit_Z = 0;
   tree_track_firsthit_phi = 0;
   tree_track_simtrack_charge = 0;
   tree_track_simtrack_pt = 0;
   tree_track_simtrack_eta = 0;
   tree_track_simtrack_phi = 0;
   tree_track_simtrack_longLived = 0;
   tree_track_simtrack_pdgId = 0;
   tree_track_simtrack_numberOfTrackerHits = 0;
   tree_track_simtrack_numberOfTrackerLayers = 0;
   tree_track_simtrack_mass = 0;
   tree_track_simtrack_status = 0;
   tree_track_genVertexPos_X = 0;
   tree_track_genVertexPos_Y = 0;
   tree_track_genVertexPos_Z = 0;
   tree_track_recoVertex_idx = 0;
   tree_track_recoJet_idx = 0;
   tree_SiCluster_subDet = 0;
   tree_SiCluster_PetalSide = 0;
   tree_SiCluster_LayerNbr = 0;
   tree_SiCluster_WheelSide = 0;
   tree_PixCluster_isBPix = 0;
   tree_PixCluster_LayerNbr = 0;
   tree_PixCluster_detID = 0;
   tree_vtx_PosX = 0;
   tree_vtx_PosY = 0;
   tree_vtx_PosZ = 0;
   tree_vtx_PosXError = 0;
   tree_vtx_PosYError = 0;
   tree_vtx_PosZError = 0;
   tree_simtrack_simtrack_charge = 0;
   tree_simtrack_simtrack_pt = 0;
   tree_simtrack_simtrack_eta = 0;
   tree_simtrack_simtrack_phi = 0;
   tree_simtrack_simtrack_longLived = 0;
   tree_simtrack_simtrack_pdgId = 0;
   tree_simtrack_simtrack_numberOfTrackerHits = 0;
   tree_simtrack_simtrack_numberOfTrackerLayers = 0;
   tree_simtrack_simtrack_mass = 0;
   tree_simtrack_simtrack_status = 0;
   tree_simtrack_genVertexPos_X = 0;
   tree_simtrack_genVertexPos_Y = 0;
   tree_simtrack_genVertexPos_Z = 0;
   tree_simtrack_isRecoMatched = 0;
   tree_simtrack_pca_dxy = 0;
   tree_simtrack_pca_dz = 0;
   tree_simtrack_trkIdx = 0;
   tree_jet_pt = 0;
   tree_jet_eta = 0;
   tree_jet_phi = 0;
   tree_gsftrack_pt = 0;
   tree_gsftrack_eta = 0;
   tree_gsftrack_phi = 0;
   tree_gsftrack_nhits = 0;
   tree_gsftrack_outerPt = 0;
   tree_gsftrack_simtrack_pt = 0;
   tree_gsftrack_simtrack_eta = 0;
   tree_gsftrack_simtrack_phi = 0;
   tree_gsftrack_simtrack_nhits = 0;
   tree_gsftrack_simtrack_outerPt = 0;
   tree_globalMuon_charge = 0;
   tree_globalMuon_pt = 0;
   tree_globalMuon_eta = 0;
   tree_globalMuon_phi = 0;
   tree_globalMuon_nChi2 = 0;
   tree_globalMuon_nHit = 0;
   tree_globalMuon_idXtrack = 0;
   tree_staMuon_charge = 0;
   tree_staMuon_pt = 0;
   tree_staMuon_eta = 0;
   tree_staMuon_phi = 0;
   tree_staMuon_nChi2 = 0;
   tree_staMuon_nHit = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("tree_track_pt", &tree_track_pt, &b_tree_track_pt);
   fChain->SetBranchAddress("tree_track_outerPt", &tree_track_outerPt, &b_tree_track_outerPt);
   fChain->SetBranchAddress("tree_track_eta", &tree_track_eta, &b_tree_track_eta);
   fChain->SetBranchAddress("tree_track_phi", &tree_track_phi, &b_tree_track_phi);
   fChain->SetBranchAddress("tree_track_nhits", &tree_track_nhits, &b_tree_track_nhits);
   fChain->SetBranchAddress("tree_track_NChi2", &tree_track_NChi2, &b_tree_track_NChi2);
   fChain->SetBranchAddress("tree_track_isHighQuality", &tree_track_isHighQuality, &b_tree_track_isHighQuality);
   fChain->SetBranchAddress("tree_track_isLoose", &tree_track_isLoose, &b_tree_track_isLoose);
   fChain->SetBranchAddress("tree_track_isTight", &tree_track_isTight, &b_tree_track_isTight);
   fChain->SetBranchAddress("tree_track_dxy", &tree_track_dxy, &b_tree_track_dxy);
   fChain->SetBranchAddress("tree_track_dxyError", &tree_track_dxyError, &b_tree_track_dxyError);
   fChain->SetBranchAddress("tree_track_dz", &tree_track_dz, &b_tree_track_dz);
   fChain->SetBranchAddress("tree_track_dzError", &tree_track_dzError, &b_tree_track_dzError);
   fChain->SetBranchAddress("tree_track_numberOfLostHits", &tree_track_numberOfLostHits, &b_tree_track_numberOfLostHits);
   fChain->SetBranchAddress("tree_track_numberOfValidHits", &tree_track_numberOfValidHits, &b_tree_track_numberOfValidHits);
   fChain->SetBranchAddress("tree_track_originalAlgo", &tree_track_originalAlgo, &b_tree_track_originalAlgo);
   fChain->SetBranchAddress("tree_track_algo", &tree_track_algo, &b_tree_track_algo);
   fChain->SetBranchAddress("tree_track_stopReason", &tree_track_stopReason, &b_tree_track_stopReason);
   fChain->SetBranchAddress("tree_track_isSimMatched", &tree_track_isSimMatched, &b_tree_track_isSimMatched);
   fChain->SetBranchAddress("tree_track_numberOfValidPixelHits", &tree_track_numberOfValidPixelHits, &b_tree_track_numberOfValidPixelHits);
   fChain->SetBranchAddress("tree_track_numberOfValidStripHits", &tree_track_numberOfValidStripHits, &b_tree_track_numberOfValidStripHits);
   fChain->SetBranchAddress("tree_track_numberOfValidStripTIBHits", &tree_track_numberOfValidStripTIBHits, &b_tree_track_numberOfValidStripTIBHits);
   fChain->SetBranchAddress("tree_track_numberOfValidStripTIDHits", &tree_track_numberOfValidStripTIDHits, &b_tree_track_numberOfValidStripTIDHits);
   fChain->SetBranchAddress("tree_track_numberOfValidStripTOBHits", &tree_track_numberOfValidStripTOBHits, &b_tree_track_numberOfValidStripTOBHits);
   fChain->SetBranchAddress("tree_track_numberOfValidStripTECHits", &tree_track_numberOfValidStripTECHits, &b_tree_track_numberOfValidStripTECHits);
   fChain->SetBranchAddress("tree_track_numberOfValidPixelBarrelHits", &tree_track_numberOfValidPixelBarrelHits, &b_tree_track_numberOfValidPixelBarrelHits);
   fChain->SetBranchAddress("tree_track_numberOfValidPixelEndcapHits", &tree_track_numberOfValidPixelEndcapHits, &b_tree_track_numberOfValidPixelEndcapHits);
   fChain->SetBranchAddress("tree_track_hasValidHitInPixelLayer", &tree_track_hasValidHitInPixelLayer, &b_tree_track_hasValidHitInPixelLayer);
   fChain->SetBranchAddress("tree_track_trackerLayersWithMeasurement", &tree_track_trackerLayersWithMeasurement, &b_tree_track_trackerLayersWithMeasurement);
   fChain->SetBranchAddress("tree_track_pixelLayersWithMeasurement", &tree_track_pixelLayersWithMeasurement, &b_tree_track_pixelLayersWithMeasurement);
   fChain->SetBranchAddress("tree_track_stripTECLayersWithMeasurement", &tree_track_stripTECLayersWithMeasurement, &b_tree_track_stripTECLayersWithMeasurement);
   fChain->SetBranchAddress("tree_track_stripTIBLayersWithMeasurement", &tree_track_stripTIBLayersWithMeasurement, &b_tree_track_stripTIBLayersWithMeasurement);
   fChain->SetBranchAddress("tree_track_stripTIDLayersWithMeasurement", &tree_track_stripTIDLayersWithMeasurement, &b_tree_track_stripTIDLayersWithMeasurement);
   fChain->SetBranchAddress("tree_track_stripTOBLayersWithMeasurement", &tree_track_stripTOBLayersWithMeasurement, &b_tree_track_stripTOBLayersWithMeasurement);
   fChain->SetBranchAddress("tree_track_nPixel", &tree_track_nPixel, &b_tree_track_nPixel);
   fChain->SetBranchAddress("tree_track_nStrip", &tree_track_nStrip, &b_tree_track_nStrip);
   fChain->SetBranchAddress("tree_track_vx", &tree_track_vx, &b_tree_track_vx);
   fChain->SetBranchAddress("tree_track_vy", &tree_track_vy, &b_tree_track_vy);
   fChain->SetBranchAddress("tree_track_vz", &tree_track_vz, &b_tree_track_vz);
   fChain->SetBranchAddress("tree_track_firsthit_X", &tree_track_firsthit_X, &b_tree_track_firsthit_X);
   fChain->SetBranchAddress("tree_track_firsthit_Y", &tree_track_firsthit_Y, &b_tree_track_firsthit_Y);
   fChain->SetBranchAddress("tree_track_firsthit_Z", &tree_track_firsthit_Z, &b_tree_track_firsthit_Z);
   fChain->SetBranchAddress("tree_track_firsthit_phi", &tree_track_firsthit_phi, &b_tree_track_firsthit_phi);
   fChain->SetBranchAddress("tree_track_simtrack_charge", &tree_track_simtrack_charge, &b_tree_track_simtrack_charge);
   fChain->SetBranchAddress("tree_track_simtrack_pt", &tree_track_simtrack_pt, &b_tree_track_simtrack_pt);
   fChain->SetBranchAddress("tree_track_simtrack_eta", &tree_track_simtrack_eta, &b_tree_track_simtrack_eta);
   fChain->SetBranchAddress("tree_track_simtrack_phi", &tree_track_simtrack_phi, &b_tree_track_simtrack_phi);
   fChain->SetBranchAddress("tree_track_simtrack_longLived", &tree_track_simtrack_longLived, &b_tree_track_simtrack_longLived);
   fChain->SetBranchAddress("tree_track_simtrack_pdgId", &tree_track_simtrack_pdgId, &b_tree_track_simtrack_pdgId);
   fChain->SetBranchAddress("tree_track_simtrack_numberOfTrackerHits", &tree_track_simtrack_numberOfTrackerHits, &b_tree_track_simtrack_numberOfTrackerHits);
   fChain->SetBranchAddress("tree_track_simtrack_numberOfTrackerLayers", &tree_track_simtrack_numberOfTrackerLayers, &b_tree_track_simtrack_numberOfTrackerLayers);
   fChain->SetBranchAddress("tree_track_simtrack_mass", &tree_track_simtrack_mass, &b_tree_track_simtrack_mass);
   fChain->SetBranchAddress("tree_track_simtrack_status", &tree_track_simtrack_status, &b_tree_track_simtrack_status);
   fChain->SetBranchAddress("tree_track_genVertexPos_X", &tree_track_genVertexPos_X, &b_tree_track_genVertexPos_X);
   fChain->SetBranchAddress("tree_track_genVertexPos_Y", &tree_track_genVertexPos_Y, &b_tree_track_genVertexPos_Y);
   fChain->SetBranchAddress("tree_track_genVertexPos_Z", &tree_track_genVertexPos_Z, &b_tree_track_genVertexPos_Z);
   fChain->SetBranchAddress("tree_track_recoVertex_idx", &tree_track_recoVertex_idx, &b_tree_track_recoVertex_idx);
   fChain->SetBranchAddress("tree_track_recoJet_idx", &tree_track_recoJet_idx, &b_tree_track_recoJet_idx);
   fChain->SetBranchAddress("tree_SiCluster_subDet", &tree_SiCluster_subDet, &b_tree_SiCluster_subDet);
   fChain->SetBranchAddress("tree_SiCluster_PetalSide", &tree_SiCluster_PetalSide, &b_tree_SiCluster_PetalSide);
   fChain->SetBranchAddress("tree_SiCluster_LayerNbr", &tree_SiCluster_LayerNbr, &b_tree_SiCluster_LayerNbr);
   fChain->SetBranchAddress("tree_SiCluster_WheelSide", &tree_SiCluster_WheelSide, &b_tree_SiCluster_WheelSide);
   fChain->SetBranchAddress("tree_PixCluster_isBPix", &tree_PixCluster_isBPix, &b_tree_PixCluster_isBPix);
   fChain->SetBranchAddress("tree_PixCluster_LayerNbr", &tree_PixCluster_LayerNbr, &b_tree_PixCluster_LayerNbr);
   fChain->SetBranchAddress("tree_PixCluster_detID", &tree_PixCluster_detID, &b_tree_PixCluster_detID);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("tree_bs_PosX", &tree_bs_PosX, &b_tree_bs_PosX);
   fChain->SetBranchAddress("tree_bs_PosY", &tree_bs_PosY, &b_tree_bs_PosY);
   fChain->SetBranchAddress("tree_bs_PosZ", &tree_bs_PosZ, &b_tree_bs_PosZ);
   fChain->SetBranchAddress("tree_vtx_PosX", &tree_vtx_PosX, &b_tree_vtx_PosX);
   fChain->SetBranchAddress("tree_vtx_PosY", &tree_vtx_PosY, &b_tree_vtx_PosY);
   fChain->SetBranchAddress("tree_vtx_PosZ", &tree_vtx_PosZ, &b_tree_vtx_PosZ);
   fChain->SetBranchAddress("tree_vtx_PosXError", &tree_vtx_PosXError, &b_tree_vtx_PosXError);
   fChain->SetBranchAddress("tree_vtx_PosYError", &tree_vtx_PosYError, &b_tree_vtx_PosYError);
   fChain->SetBranchAddress("tree_vtx_PosZError", &tree_vtx_PosZError, &b_tree_vtx_PosZError);
   fChain->SetBranchAddress("tree_simtrack_simtrack_charge", &tree_simtrack_simtrack_charge, &b_tree_simtrack_simtrack_charge);
   fChain->SetBranchAddress("tree_simtrack_simtrack_pt", &tree_simtrack_simtrack_pt, &b_tree_simtrack_simtrack_pt);
   fChain->SetBranchAddress("tree_simtrack_simtrack_eta", &tree_simtrack_simtrack_eta, &b_tree_simtrack_simtrack_eta);
   fChain->SetBranchAddress("tree_simtrack_simtrack_phi", &tree_simtrack_simtrack_phi, &b_tree_simtrack_simtrack_phi);
   fChain->SetBranchAddress("tree_simtrack_simtrack_longLived", &tree_simtrack_simtrack_longLived, &b_tree_simtrack_simtrack_longLived);
   fChain->SetBranchAddress("tree_simtrack_simtrack_pdgId", &tree_simtrack_simtrack_pdgId, &b_tree_simtrack_simtrack_pdgId);
   fChain->SetBranchAddress("tree_simtrack_simtrack_numberOfTrackerHits", &tree_simtrack_simtrack_numberOfTrackerHits, &b_tree_simtrack_simtrack_numberOfTrackerHits);
   fChain->SetBranchAddress("tree_simtrack_simtrack_numberOfTrackerLayers", &tree_simtrack_simtrack_numberOfTrackerLayers, &b_tree_simtrack_simtrack_numberOfTrackerLayers);
   fChain->SetBranchAddress("tree_simtrack_simtrack_mass", &tree_simtrack_simtrack_mass, &b_tree_simtrack_simtrack_mass);
   fChain->SetBranchAddress("tree_simtrack_simtrack_status", &tree_simtrack_simtrack_status, &b_tree_simtrack_simtrack_status);
   fChain->SetBranchAddress("tree_simtrack_genVertexPos_X", &tree_simtrack_genVertexPos_X, &b_tree_simtrack_genVertexPos_X);
   fChain->SetBranchAddress("tree_simtrack_genVertexPos_Y", &tree_simtrack_genVertexPos_Y, &b_tree_simtrack_genVertexPos_Y);
   fChain->SetBranchAddress("tree_simtrack_genVertexPos_Z", &tree_simtrack_genVertexPos_Z, &b_tree_simtrack_genVertexPos_Z);
   fChain->SetBranchAddress("tree_simtrack_isRecoMatched", &tree_simtrack_isRecoMatched, &b_tree_simtrack_isRecoMatched);
   fChain->SetBranchAddress("tree_simtrack_pca_dxy", &tree_simtrack_pca_dxy, &b_tree_simtrack_pca_dxy);
   fChain->SetBranchAddress("tree_simtrack_pca_dz", &tree_simtrack_pca_dz, &b_tree_simtrack_pca_dz);
   fChain->SetBranchAddress("tree_simtrack_trkIdx", &tree_simtrack_trkIdx, &b_tree_simtrack_trkIdx);
   fChain->SetBranchAddress("tree_jet_pt", &tree_jet_pt, &b_tree_jet_pt);
   fChain->SetBranchAddress("tree_jet_eta", &tree_jet_eta, &b_tree_jet_eta);
   fChain->SetBranchAddress("tree_jet_phi", &tree_jet_phi, &b_tree_jet_phi);
   fChain->SetBranchAddress("tree_gsftrack_pt", &tree_gsftrack_pt, &b_tree_gsftrack_pt);
   fChain->SetBranchAddress("tree_gsftrack_eta", &tree_gsftrack_eta, &b_tree_gsftrack_eta);
   fChain->SetBranchAddress("tree_gsftrack_phi", &tree_gsftrack_phi, &b_tree_gsftrack_phi);
   fChain->SetBranchAddress("tree_gsftrack_nhits", &tree_gsftrack_nhits, &b_tree_gsftrack_nhits);
   fChain->SetBranchAddress("tree_gsftrack_outerPt", &tree_gsftrack_outerPt, &b_tree_gsftrack_outerPt);
   fChain->SetBranchAddress("tree_gsftrack_simtrack_pt", &tree_gsftrack_simtrack_pt, &b_tree_gsftrack_simtrack_pt);
   fChain->SetBranchAddress("tree_gsftrack_simtrack_eta", &tree_gsftrack_simtrack_eta, &b_tree_gsftrack_simtrack_eta);
   fChain->SetBranchAddress("tree_gsftrack_simtrack_phi", &tree_gsftrack_simtrack_phi, &b_tree_gsftrack_simtrack_phi);
   fChain->SetBranchAddress("tree_gsftrack_simtrack_nhits", &tree_gsftrack_simtrack_nhits, &b_tree_gsftrack_simtrack_nhits);
   fChain->SetBranchAddress("tree_gsftrack_simtrack_outerPt", &tree_gsftrack_simtrack_outerPt, &b_tree_gsftrack_simtrack_outerPt);
   fChain->SetBranchAddress("tree_globalMuon_charge", &tree_globalMuon_charge, &b_tree_globalMuon_charge);
   fChain->SetBranchAddress("tree_globalMuon_pt", &tree_globalMuon_pt, &b_tree_globalMuon_pt);
   fChain->SetBranchAddress("tree_globalMuon_eta", &tree_globalMuon_eta, &b_tree_globalMuon_eta);
   fChain->SetBranchAddress("tree_globalMuon_phi", &tree_globalMuon_phi, &b_tree_globalMuon_phi);
   fChain->SetBranchAddress("tree_globalMuon_nChi2", &tree_globalMuon_nChi2, &b_tree_globalMuon_nChi2);
   fChain->SetBranchAddress("tree_globalMuon_nHit", &tree_globalMuon_nHit, &b_tree_globalMuon_nHit);
   fChain->SetBranchAddress("tree_globalMuon_idXtrack", &tree_globalMuon_idXtrack, &b_tree_globalMuon_idXtrack);
   fChain->SetBranchAddress("tree_staMuon_charge", &tree_staMuon_charge, &b_tree_staMuon_charge);
   fChain->SetBranchAddress("tree_staMuon_pt", &tree_staMuon_pt, &b_tree_staMuon_pt);
   fChain->SetBranchAddress("tree_staMuon_eta", &tree_staMuon_eta, &b_tree_staMuon_eta);
   fChain->SetBranchAddress("tree_staMuon_phi", &tree_staMuon_phi, &b_tree_staMuon_phi);
   fChain->SetBranchAddress("tree_staMuon_nChi2", &tree_staMuon_nChi2, &b_tree_staMuon_nChi2);
   fChain->SetBranchAddress("tree_staMuon_nHit", &tree_staMuon_nHit, &b_tree_staMuon_nHit);
   Notify();
}

Bool_t myReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myReader_cxx
