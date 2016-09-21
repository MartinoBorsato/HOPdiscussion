//File Edit Options Buffers Tools C++ Help                                                                                                                                                                     
#include<iostream>
using namespace std;
void addHOP()
{

  TString name = "/afs/cern.ch/user/c/cmarinbe/git/TupleToolHOP/DVntuple_newClassPart.root";
  TFile* newFile = new TFile(name, "UPDATE");
  TTree *newTree = (TTree*)newFile->Get("TupleBdKstEE_HOP/DecayTree");
  cout << "Number of events: " << newTree->GetEntries() ;
  
  ULong64_t eventNumber;

  Double_t B0_ENDVERTEX_X ; 
  Double_t B0_ENDVERTEX_Y ; 
  Double_t B0_ENDVERTEX_Z ; 

  Double_t B0_OWNPV_X ; 
  Double_t B0_OWNPV_Y ; 
  Double_t B0_OWNPV_Z ; 

  Double_t B0_FD_OWNPV ; 

  Double_t Kst_PX ; 
  Double_t Kst_PY ; 
  Double_t Kst_PZ ; 
  Double_t Kst_P ; 
  Double_t Kst_PE ; 

  Double_t JPs_PX ; 
  Double_t JPs_PY ; 
  Double_t JPs_PZ ; 
  Double_t JPs_P ; 
  Double_t JPs_MM;

  Double_t E1_PX ; 
  Double_t E1_PY ; 
  Double_t E1_PZ ; 
  Double_t E1_P ; 

  Double_t E2_PX ; 
  Double_t E2_PY ; 
  Double_t E2_PZ ; 
  Double_t E2_P ;
  
  newTree->SetBranchAddress("B0_ENDVERTEX_X" ,&B0_ENDVERTEX_X);
  newTree->SetBranchAddress("B0_ENDVERTEX_Y" ,&B0_ENDVERTEX_Y);
  newTree->SetBranchAddress("B0_ENDVERTEX_Z" ,&B0_ENDVERTEX_Z);

  newTree->SetBranchAddress("B0_OWNPV_X" ,&B0_OWNPV_X);
  newTree->SetBranchAddress("B0_OWNPV_Y" ,&B0_OWNPV_Y);
  newTree->SetBranchAddress("B0_OWNPV_Z" ,&B0_OWNPV_Z);

  newTree->SetBranchAddress("B0_FD_OWNPV" ,&B0_FD_OWNPV);

  newTree->SetBranchAddress("Kst_892_0_PX" ,&Kst_PX) ; 
  newTree->SetBranchAddress("Kst_892_0_PY" ,&Kst_PY) ; 
  newTree->SetBranchAddress("Kst_892_0_PZ" ,&Kst_PZ) ; 
  newTree->SetBranchAddress("Kst_892_0_P" ,&Kst_P) ; 
  newTree->SetBranchAddress("Kst_892_0_PE" ,&Kst_PE) ; 

  newTree->SetBranchAddress("J_psi_1S_PX" ,&JPs_PX) ; 
  newTree->SetBranchAddress("J_psi_1S_PY" ,&JPs_PY) ; 
  newTree->SetBranchAddress("J_psi_1S_PZ" ,&JPs_PZ) ; 
  newTree->SetBranchAddress("J_psi_1S_P"  ,&JPs_P) ;
  newTree->SetBranchAddress("J_psi_1S_MM" ,&JPs_MM) ;

  newTree->SetBranchAddress("eplus_PX" ,&E1_PX) ; 
  newTree->SetBranchAddress("eplus_PY" ,&E1_PY) ; 
  newTree->SetBranchAddress("eplus_PZ" ,&E1_PZ) ; 
  newTree->SetBranchAddress("eplus_P" ,&E1_P) ; 

  newTree->SetBranchAddress("eminus_PX" ,&E2_PX) ; 
  newTree->SetBranchAddress("eminus_PY" ,&E2_PY) ; 
  newTree->SetBranchAddress("eminus_PZ" ,&E2_PZ) ; 
  newTree->SetBranchAddress("eminus_P" ,&E2_P) ; 

  Int_t signumberOfEntries = newTree->GetEntries();

  const Double_t PDG_e_M = 0.510998910 ;
  const Double_t PDG_JPs_M = 3096.916;
  Float_t JPs_SF ;
  TLorentzVector KstarMom ;
  TLorentzVector ScaledL1Mom ;
  TLorentzVector ScaledL2Mom ;
  TLorentzVector ScaledJpsiMom ;
  TLorentzVector ScaledJpsiMom_fixM ;
  TLorentzVector ScaledJpsiMom_measM ;
  TLorentzVector ScaledBMom ;
  TLorentzVector ScaledBMom_fixJPsM ;
  TLorentzVector ScaledBMom_measJPsM ;

  Float_t Hop = 0.0;
  Float_t HopM = 0.0;
  Float_t HopM_fixJPsM = 0.0;
  Float_t HopM_measJPsM = 0.0;
  Float_t HopMee = 0.0;

  TBranch *Hop_branch = newTree->Branch("HOP", &Hop, "HOP") ;
  TBranch *HopM_branch = newTree->Branch("HOP_B0_MM", &HopM, "HOP_B0_MM");
  TBranch *HopM_fixJPsM_branch = newTree->Branch("HOP_B0_MM_fixJPsM", &HopM_fixJPsM, "HOP_B0_MM_fixJPsM");
  TBranch *HopM_measJPsM_branch = newTree->Branch("HOP_B0_MM_measJPsM", &HopM_measJPsM, "HOP_B0_MM_measJPsM");
  TBranch *HopMee_branch = newTree->Branch("HOP_J_psi_1S_MM", &HopMee, "HOP_J_psi_1S_MM");
  for (Int_t loopie=0; loopie < signumberOfEntries; ++loopie)
  {
    newTree->GetEntry(loopie); 
 
    // compute the pT of the K* wrt to the B line of flight 
    
    Double_t cosThetaKstar = ((B0_ENDVERTEX_X-B0_OWNPV_X)*Kst_PX + (B0_ENDVERTEX_Y-B0_OWNPV_Y)*Kst_PY +
                              (B0_ENDVERTEX_Z-B0_OWNPV_Z)*Kst_PZ) / (Kst_P*B0_FD_OWNPV) ;
    Double_t pTKstar = Kst_P*TMath::Sqrt(1.-cosThetaKstar*cosThetaKstar) ;
    Double_t cosThetaY = ((B0_ENDVERTEX_X-B0_OWNPV_X)*JPs_PX + (B0_ENDVERTEX_Y-B0_OWNPV_Y)*JPs_PY + 
                          (B0_ENDVERTEX_Z-B0_OWNPV_Z)*JPs_PZ) / (JPs_P*B0_FD_OWNPV);
    Double_t pTY = JPs_P*TMath::Sqrt(1.-cosThetaY*cosThetaY) ;

    JPs_SF = pTKstar/pTY; 

    Hop = JPs_SF ; 

    ScaledL1Mom.SetXYZM( JPs_SF*E1_PX, JPs_SF*E1_PY, JPs_SF*E1_PZ, PDG_e_M );
    ScaledL2Mom.SetXYZM( JPs_SF*E2_PX, JPs_SF*E2_PY, JPs_SF*E2_PZ, PDG_e_M );
    KstarMom.SetPxPyPzE( Kst_PX, Kst_PY, Kst_PZ, Kst_PE );
    ScaledBMom = ScaledL1Mom + ScaledL2Mom + KstarMom  ;
    ScaledJpsiMom = ScaledL1Mom + ScaledL2Mom ; 
    ScaledJpsiMom_fixM.SetXYZM( JPs_SF*JPs_PX, JPs_SF*JPs_PY, JPs_SF*JPs_PZ, PDG_JPs_M );
    ScaledJpsiMom_measM.SetXYZM( JPs_SF*JPs_PX, JPs_SF*JPs_PY, JPs_SF*JPs_PZ, JPs_MM );
    ScaledBMom_fixJPsM = ScaledJpsiMom_fixM + KstarMom;
    ScaledBMom_measJPsM = ScaledJpsiMom_measM + KstarMom;

    HopM = ScaledBMom.M();
    HopM_fixJPsM = ScaledBMom_fixJPsM.M();
    HopM_measJPsM = ScaledBMom_measJPsM.M();
    HopMee = ScaledJpsiMom.M();

    Hop_branch->Fill();
    HopM_branch->Fill();
    HopM_fixJPsM_branch->Fill();
    HopM_measJPsM_branch->Fill();
    HopMee_branch->Fill();
  }
  newFile->Write();  
  cout << "  new Number of events: " << newTree->GetEntries()<< endl ;
  

}
