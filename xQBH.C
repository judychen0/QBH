#include <vector>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
using namespace std;
#include <iostream>
#include <TProfile.h>

#include "untuplizer.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"


Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}

void xQBH(){
  //***********************Initialization***********************//

  //access EventTree with TreeReader class
  TreeReader data("/home/judy/ntuhep/QBH/MC_sample/ADD_n6_2016/ggtree_mc_*.root");

  //create an output .root file
  TFile *fout_;
  fout_ = new TFile("output_ADD_ggtree_mc.root","RECREATE");

  //create histograms in output .root file
  TH1F *h_npho = new TH1F("h_npho", "n pho", 20, 0., 20);
  TH1F *h_njet = new TH1F("h_njet", "n jet", 20, 0., 20);

  TH1F *h_nMCpho = new TH1F("h_nMCpho", "n MC pho", 40, 0., 40);
  TH1F *h_dr_pho = new TH1F("h_dr_pho", "dr of photon", 100, 0., 0.2);
  TH1F *h_dpt_pho = new TH1F("h_dpt_pho", "dpt of photon", 100, 0., 1);
  TH2F *h_dptdr_pho = new TH2F("h_dptdr_pho", "dptdr of photon", 100, 0., 1, 100, 0., 2);
  TH1F *h_nrealpho = new TH1F("h_nrealpho", "n real photon", 20, 0., 20);
  TH1F *h_ngenpho = new TH1F("h_ngenpho", "n gen-matched pho", 20, 0., 20);
  
  TH1F *h_njetGen = new TH1F("h_njetGen", "n jetGen parton", 20, 0., 20);
  TH1F *h_dr_jet = new TH1F("h_dr_jet", "dr of jet", 100, 0., 1);
  TH1F *h_dpt_jet = new TH1F("h_dpt_jet", "dpt of jet", 100, 0., 2);
  TH2F *h_dptdr_jet = new TH2F("h_dptdr_jet", "dptdr of jet", 100, 0., 2, 100, 0., 4);
  TH1F *h_nrealjet = new TH1F("h_nrealjet", "n real jet", 20, 0., 20);
  TH1F *h_ncjet = new TH1F("h_ncjet", "n c-jet", 20, 0., 20);
  TH1F *h_nLjet = new TH1F("h_nLjet", "n L-jet", 20, 0., 20);
  TH1F *h_nbjet = new TH1F("h_nbjet", "n b-jet", 20, 0., 20);
  
  TH1F *h_MCQBHMass = new TH1F("h_MCQBHMass", "MC QBH mass", 50, 2000, 7000);
  TH1F *h_pjmass_real = new TH1F("h_pjmass_real", "pho jet real inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_HPT = new TH1F("h_pjmass_HPT", "HPT pho jet inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_Sndpho = new TH1F("h_pjmass_Sndpho", "2nd PT pho jet inv mass", 100, 2000, 7000);

  TH1F *h_pjmass_phocjet = new TH1F("h_pjmass_phocjet", "pho cjet inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_phopascjet = new TH1F("h_pjmass_phopascjet", "pho non-cjet inv mass", 100, 2000, 7000);

  TH1F *h_pjmass_phoLjet = new TH1F("h_pjmass_phoLjet", "pho Ljet inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_phobjet = new TH1F("h_pjmass_phobjet", "pho bjet inv mass", 100, 2000, 7000);
  
  TH1F *h_pjmass_HPTcjet = new TH1F("h_pjmass_HPTcjet", "HPT pho cjet inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_HPTLjet = new TH1F("h_pjmass_HPTLjet", "HPT pho Ljet inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_HPTbjet = new TH1F("h_pjmass_HPTbjet", "HPT pho bjet inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_Sndcjet = new TH1F("h_pjmass_Sndcjet", "pho 2nd-cjet inv mass", 100, 2000, 7000);
  TH1F *h_pjmass_SndLjet = new TH1F("h_pjmass_SndLjet", "pho 2nd-Ljet inv mass", 100, 2000, 7000);
  
  
  TH1F *h_realjet_HPT = new TH1F("h_realjet_HPT", "real jet HPT", 100, 0, 4000);
  TH1F *h_realjet_SndPT = new TH1F("h_realjet_SndPT", "real jet 2nd PT", 100, 0, 4000);

  TH1F *h_dr_genphojet = new TH1F("h_dr_genphojet", "dr of gen photon and jet", 100, 0., 1);
  TH1F *h_dr_matchedphojet = new TH1F("h_dr_matchedphojet", "dr of matched photon and jet", 100, 0., 1);
  //TH1F *h_ctageff = new TH1F("h_ctageff", "cjet over Ljet resonance events", 100, 0, 1);
  h_dr_pho->Sumw2();
  h_dpt_pho->Sumw2();
  h_dptdr_pho->Sumw2();
  h_dr_jet->Sumw2();
  h_dpt_jet->Sumw2();
  h_dptdr_jet->Sumw2();
  
  h_pjmass_real->Sumw2();
  h_pjmass_HPT->Sumw2();
  h_pjmass_Sndpho->Sumw2();
  
  h_pjmass_phocjet->Sumw2();
  h_pjmass_phoLjet->Sumw2();
  h_pjmass_phobjet->Sumw2();
  h_realjet_HPT->Sumw2();
  h_realjet_SndPT->Sumw2();

  
  TH1F *h_realphoEB_pt = new TH1F("h_realphoEB_pt", "gen-matched pho pt of EB", 80, 0., 4000);
  TH1F *h_realphoEB_eta = new TH1F("h_realphoEB_eta", "gen-matched pho eta of EB", 100, 0, 5);
  TH1F *h_chIso_rc = new TH1F("h_chIso_rc", "rho corrected chIso", 100, 0, 10);
  TH1F *h_phoIso_rc = new TH1F("h_phoIso_rc", "rho corrected phoIso", 100, 0, 10);
  TH1F *h_nhIso_rc = new TH1F("h_nhIso_rc", "rho corrected nhIso", 100, 0, 10);

  TH1F *h_realphoEB_pt_L = new TH1F("h_realphoEB_pt_L", "gen-matched pho pt of EB with Loose cut", 80, 0., 4000);
  TH1F *h_realphoEB_pt_M = new TH1F("h_realphoEB_pt_M", "gen-matched pho pt of EB with Medium cut", 80, 0., 4000);
  TH1F *h_realphoEB_pt_T = new TH1F("h_realphoEB_pt_T", "gen-matched phopt of EB with Tight cut", 80, 0., 4000);
  //TH1F *h_realphoEB_pt_01 = new TH1F("h_realphoEB_pt_01")
  
  
  //define branch variables
  Bool_t   isData;
  Int_t    run;
  Long64_t event;
  Float_t jetPt_, jetEta_, jetPhi_;
  Int_t   isMatched, isMatchedEle, isConverted;
  Float_t HoverE, sieie, chIso, phoIso, nhIso, eleVeto, rho;
  Float_t chIsoEB_L, chIsoEB_M, chIsoEB_T, phoIsoEB_L, phoIsoEB_M, phoIsoEB_T, nhIsoEB_L, nhIsoEB_M, nhIsoEB_T;
  Float_t sieieFull5x5, sipipFull5x5, sieipFull5x5, r9Full5x5, e2x2Full5x5, e5x5Full5x5;
  Float_t jetCEF_, jetNEF_, jetCHF_, jetNHF_;
  Int_t jetNCH_, jetNNP_;

  TTree *outtree_;
  outtree_ = new TTree("t", "mini tree");

  outtree_->Branch("run", &run, "run/I");
  outtree_->Branch("event", &event, "event/L");
  outtree_->Branch("isData",         &isData,        "isData/O");
  outtree_->Branch("eleVeto",      &eleVeto,      "eleVeto/I");
  outtree_->Branch("HoverE",       &HoverE,       "HoverE/F");
  outtree_->Branch("chIso",        &chIso,        "chIso/F");
  outtree_->Branch("phoIso",       &phoIso,       "phoIso/F");
  outtree_->Branch("nhIso",        &nhIso,        "nhIso/F");

  outtree_->Branch("sieieFull5x5",        &sieieFull5x5,        "sieieFull5x5/F");
  outtree_->Branch("sieipFull5x5",        &sieipFull5x5,        "sieipFull5x5/F");
  outtree_->Branch("sipipFull5x5",        &sipipFull5x5,        "sipipFull5x5/F");
  outtree_->Branch("r9Full5x5",        &r9Full5x5,        "r9Full5x5/F");
  outtree_->Branch("e2x2Full5x5",        &e2x2Full5x5,        "e2x2Full5x5/F");
  outtree_->Branch("e5x5Full5x5",        &e5x5Full5x5,        "e5x5Full5x5/F");
  outtree_->Branch("jetCEF", &jetCEF_, "jetCEF/F");
  outtree_->Branch("jetNEF", &jetNEF_, "jetNEF/F");
  outtree_->Branch("jetCHF", &jetCHF_, "jetCHF/F");
  outtree_->Branch("jetNHF", &jetNHF_, "jetNHF/F");
  outtree_->Branch("jetNCH", &jetNCH_, "jetNCH/I");
  outtree_->Branch("jetNNP", &jetNNP_, "jetNNP/I");
  

  //***********************Loop***********************//

  // printf("processing entries %lli \n", data.GetEntriesFast());
  
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
     //for (Long64_t ev = 0; ev <10000; ev++) {
   //for (Long64_t ev = 0; ev < data.GetEntriesFast()/2.; ev++) {
      
    TLorentzVector phoP4, jetP4, pjP4;

    if (ev % 100000 == 0){
      fprintf(stderr, "Processing event %lli of %lli (%.3f \%)\n", ev+1, data.GetEntriesFast(), (ev+1)*100./data.GetEntriesFast());
    }
    data.GetEntry(ev);
    run     = data.GetInt("run");
    event   = data.GetLong64("event"); 
    isData = data.GetBool("isData");
    
    Int_t nPho     = data.GetInt("nPho");
    Int_t nJet     = data.GetInt("nJet");
    h_npho->Fill(nPho);
    h_njet->Fill(nJet);
    if(nJet <1) continue;

    //reco QBH
    Float_t* phoE = data.GetPtrFloat("phoE");
    Float_t* phoEt = data.GetPtrFloat("phoEt");
    Float_t* phoEta = data.GetPtrFloat("phoEta");
    Float_t* phoPhi = data.GetPtrFloat("phoPhi");
    Float_t* phoR9  = data.GetPtrFloat("phoR9");
    Float_t* phoHoverE = data.GetPtrFloat("phoHoverE");
    Float_t* phoPFChIso          = data.GetPtrFloat("phoPFChIso");
    Float_t* phoPFNeuIso         = data.GetPtrFloat("phoPFNeuIso");
    Float_t* phoPFPhoIso         = data.GetPtrFloat("phoPFPhoIso");
    Int_t*   phoEleVeto          = data.GetPtrInt("phoEleVeto");
    Float_t* phoSigmaIEtaIEtaFull5x5  = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
    Float_t* phoSigmaIEtaIPhiFull5x5  = data.GetPtrFloat("phoSigmaIEtaIPhiFull5x5");
    Float_t* phoSigmaIPhiIPhiFull5x5  = data.GetPtrFloat("phoSigmaIPhiIPhiFull5x5");
    Float_t* phoR9Full5x5           = data.GetPtrFloat("phoR9Full5x5");
    Float_t  rho                    = data.GetFloat("rho");
        
    Float_t* jetPt = data.GetPtrFloat("jetPt");
    Float_t* jetEta = data.GetPtrFloat("jetEta");
    Float_t* jetPhi = data.GetPtrFloat("jetPhi");
    Float_t* jetEn = data.GetPtrFloat("jetEn");
    Float_t* jetCEF = data.GetPtrFloat("jetCEF");
    Float_t* jetNEF = data.GetPtrFloat("jetNEF");
    Float_t* jetCHF = data.GetPtrFloat("jetCHF");
    Float_t* jetNHF = data.GetPtrFloat("jetNHF");
    Int_t*   jetNCH = data.GetPtrInt("jetNCH");
    Int_t*   jetNNP = data.GetPtrInt("jetNNP");

    //MC QBH
    Int_t    nMC   =0;     
    Int_t*   mcPID =0;
    Int_t*   mcMomPID =0;
    Float_t* mcPt      =0;
    Float_t* mcEta     =0;
    Float_t* mcPhi     =0;
    Float_t* mcE       =0;
    Float_t* mcEt      =0;
    Float_t* mcMass    =0;
    Int_t*   mcStatus  =0;

    Float_t* mcMomMass =0;
    Float_t* mcMomPt   =0;
    Float_t* mcMomEta  =0;
    Float_t* mcMomPhi  =0;
    Float_t* mcMomE    =0;
    Float_t* mcCalIsoDR04 =0;
    
    Int_t*   jetGenID    = 0;
    Int_t*   jetGenMomID = 0;
    Float_t* jetGenJetPt    = 0;
    Float_t* jetGenJetEn    = 0;
    Float_t* jetGenJetEta   = 0;
    Float_t* jetGenJetPhi   = 0;

    TLorentzVector mc_phoP4, mc_pQBHP4;

    if(!isData){
      nMC       = data.GetInt("nMC");
      mcPID     = data.GetPtrInt("mcPID");
      
      mcPt      = data.GetPtrFloat("mcPt");
      mcEta     = data.GetPtrFloat("mcEta");
      mcPhi     = data.GetPtrFloat("mcPhi");
      mcE       = data.GetPtrFloat("mcE");
      mcMass    = data.GetPtrFloat("mcMass");
      mcStatus  = data.GetPtrInt("mcStatus");

      mcMomPID  = data.GetPtrInt("mcMomPID");
      mcMomMass = data.GetPtrFloat("mcMomMass");
      mcMomPt   = data.GetPtrFloat("mcMomPt");
      mcMomEta  = data.GetPtrFloat("mcMomEta");
      mcMomPhi  = data.GetPtrFloat("mcMomPhi");
      mcCalIsoDR04 = data.GetPtrFloat("mcCalIsoDR04");

      vector <Int_t> match;
      vector <Int_t> converted;
      vector <Float_t> mmcPt;
      vector <Float_t> mmcEta;
      vector <Float_t> mmcPhi;

      //TLorentzVector genphoP4, genjetP4;
      
      //get mc photon id
      vector <Int_t> mc_phoid;
      Int_t nMCpho = 0;
      for(Int_t k=0; k < nMC; k++){
  	if(mcPID[k] == 22 && mcPt[k] >= 15 && abs(mcMomPID[k]) <= 22 ){
  	  mc_phoid.push_back(k);
  	  nMCpho++;
  	}
      }

      //create real pho list and gen pho list
      vector <Int_t> realpho_list;
      Int_t nrealpho=0;
      for(Int_t ipho=0; ipho < nPho; ipho++){
  	//if(phoEt[ipho] < 165.) continue;
  	isMatched     = -1;
  	isMatchedEle  = -1;
  	isConverted   = -1;

  	for(Int_t nn=0; nn < nMCpho; nn++){
  	  Int_t k = mc_phoid[nn];
  	  Float_t dr = deltaR(phoEta[ipho], phoPhi[ipho], mcEta[k], mcPhi[k]);
  	  Float_t dpt = fabs((phoEt[ipho] - mcPt[k])/mcPt[k]);
  	  h_dptdr_pho->Fill(dr, dpt);
  	  if(dr < 0.1) h_dpt_pho->Fill(dpt);
  	  if(dpt < 0.2) h_dr_pho->Fill(dr);
  	  if(dr < 0.1 && dpt < 0.2){
  	    isMatched = 1;
  	    //printf("MC phomatched !");
  	    realpho_list.push_back(ipho);
	    nrealpho++;
  	    break;
  	  }
  	}
      }
      //if(nPho == 0) continue;

      h_nMCpho->Fill(nMCpho);
      h_nrealpho->Fill(nrealpho);
      //printf("h_nrealpho histo saved \n");
      //h_ngenpho->Fill(ngenpho);

      
      //get rho corrected Iso
      vector <Float_t> chIso_rc;
      vector <Float_t> phoIso_rc;
      vector <Float_t> nhIso_rc;
      Float_t EAch[7] = {0.0112, 0.0108, 0.0106, 0.01002, 0.0098, 0.0089, 0.000087};
      Float_t EAnh[7] = {0.0668, 0.1054, 0.0786, 0.0223, 0.0078, 0.0028, 0.0137};
      Float_t EApho[7] = {0.1113, 0.0953, 0.0619, 0.0837, 0.1070, 0.1212, 0.1446};

      for(Int_t ii=0; ii < nrealpho; ii++){
	Int_t ipho = realpho_list[ii];
	if(fabs(phoEta[ipho]) < 1.0){
	  chIso_rc.push_back(phoPFChIso[ipho] - rho*EAch[0]);
	  phoIso_rc.push_back(phoPFPhoIso[ipho] - rho*EApho[0]);
	  nhIso_rc.push_back(phoPFNeuIso[ipho] - rho*EAnh[0]);
	}
	
	else if(fabs(phoEta[ipho]) > 1.0 && fabs(phoEta[ipho]) < 1.479){
	  chIso_rc.push_back(phoPFChIso[ipho] - rho*EAch[1]);
	  phoIso_rc.push_back(phoPFPhoIso[ipho] - rho*EApho[1]);
	  nhIso_rc.push_back(phoPFNeuIso[ipho] - rho*EAnh[1]);
	}
	else if(fabs(phoEta[ipho]) > 1.479 && fabs(phoEta[ipho]) < 2.0){
	  chIso_rc.push_back(phoPFChIso[ipho] - rho*EAch[2]);
	  phoIso_rc.push_back(phoPFPhoIso[ipho] - rho*EApho[2]);
	  nhIso_rc.push_back(phoPFNeuIso[ipho] - rho*EAnh[2]);
	}
	else if(fabs(phoEta[ipho]) > 2.0 && fabs(phoEta[ipho]) < 2.2){
	  chIso_rc.push_back(phoPFChIso[ipho] - rho*EAch[3]);
	  phoIso_rc.push_back(phoPFPhoIso[ipho] - rho*EApho[3]);
	  nhIso_rc.push_back(phoPFNeuIso[ipho] - rho*EAnh[3]);
	}
	else if(fabs(phoEta[ipho]) > 2.2 && fabs(phoEta[ipho]) < 2.3){
	  chIso_rc.push_back(phoPFChIso[ipho] - rho*EAch[4]);
	  phoIso_rc.push_back(phoPFPhoIso[ipho] - rho*EApho[4]);
	  nhIso_rc.push_back(phoPFNeuIso[ipho] - rho*EAnh[4]);
	}
	else if(fabs(phoEta[ipho]) > 2.3 && fabs(phoEta[ipho]) < 2.4){
	  chIso_rc.push_back(phoPFChIso[ipho] - rho*EAch[5]);
	  phoIso_rc.push_back(phoPFPhoIso[ipho] - rho*EApho[5]);
	  nhIso_rc.push_back(phoPFNeuIso[ipho] - rho*EAnh[5]);
	}
	else if(fabs(phoEta[ipho]) > 2.4){
	  chIso_rc.push_back(phoPFChIso[ipho] - rho*EAch[6]);
	  phoIso_rc.push_back(phoPFPhoIso[ipho] - rho*EApho[6]);
	  nhIso_rc.push_back(phoPFNeuIso[ipho] - rho*EAnh[6]);
	}
	
      }
      
      
      //get genjet id
      jetGenID = data.GetPtrInt("jetGenPartonID");
      jetGenMomID = data.GetPtrInt("jetGenPartonMomID");
      jetGenJetPt = data.GetPtrFloat("jetGenJetPt");
      jetGenJetEn = data.GetPtrFloat("jetGenJetEn");
      jetGenJetEta = data.GetPtrFloat("jetGenJetEta");
      jetGenJetPhi = data.GetPtrFloat("jetGenJetPhi");

      vector <Int_t> jetGenid;
      Int_t njetGen=0;
      for(Int_t j=0; j < nJet; j++){
  	if(fabs(jetGenID[j]) <= 6 || jetGenID[j] == 21){
  	  jetGenid.push_back(j);
  	  njetGen++;
  	}
      }
      h_njetGen->Fill(njetGen);

      //create realjet list
      vector <Int_t> realjet_list;
      Int_t nrealjet=0;
      for(Int_t ijet=0; ijet < nJet; ijet++){
  	//if(fabs(jetPt[ijet]) < 1000) continue;
  	for(Int_t jj=0; jj < njetGen; jj++){
  	  Int_t j = jetGenid[jj];
  	  Float_t dr = deltaR(jetEta[ijet], jetPhi[ijet], jetGenJetEta[j], jetGenJetPhi[j]);
  	  Float_t dpt = fabs((jetPt[ijet]-jetGenJetPt[j])/jetGenJetPt[j]);
  	  h_dptdr_jet->Fill(dpt, dr);
  	  if(dpt < 0.5) h_dr_jet->Fill(dr);
  	  if(dr < 0.3) h_dpt_jet->Fill(dpt);
  	  if(dr < 0.3 && dpt < 0.5){
  	    realjet_list.push_back(ijet);
  	    nrealjet++;
  	    break;
  	  }
  	}
      }
      h_nrealjet->Fill(nrealjet);

      /*
      //deltaR of genpho and genjet
      TLorentzVector genphoP4, genjetP4;
      for(Int_t nn=0; nn<nMCpho; nn++){
  	genphoP4.SetPtEtaPhiM(mcPt[mc_phoid[nn]], mcEta[mc_phoid[nn]], mcPhi[mc_phoid[nn]], 0.);
  	for(Int_t jj=0; jj<njetGen; jj++){
  	  genjetP4.SetPtEtaPhiM(jetGenJetPt[jetGenid[jj]], jetGenJetEta[jetGenid[jj]], jetGenJetPhi[jetGenid[jj]], jetGenJetEn[jetGenid[jj]]);
	  
  	  h_dr_genphojet->Fill(genphoP4.DeltaR(genjetP4));
  	}
      }
      */
      
      //create c-jet list
      vector <Int_t> cjetGenid;
      Int_t ncjetGen=0;
      for(Int_t j=0; j < nJet; j++){
  	if(fabs(jetGenID[j]) == 4){
  	  cjetGenid.push_back(j);
  	  ncjetGen++;
  	}
      }
      
      vector <Int_t> cjet_list;
      Int_t ncjet=0;
      for(Int_t ijet=0; ijet < nJet; ijet++){
  	for(Int_t jj=0; jj < ncjetGen; jj++){
  	  Int_t j = cjetGenid[jj];
  	  Float_t dr = deltaR(jetEta[ijet], jetPhi[ijet], jetGenJetEta[j], jetGenJetPhi[j]);
  	  Float_t dpt = fabs((jetPt[ijet]-jetGenJetPt[j])/jetGenJetPt[j]);
  	  if(dr < 0.2 && dpt < 0.5){
  	    cjet_list.push_back(ijet);
  	    ncjet++;
  	    break;
  	  }
  	}
      }
      //if(ncjet ==0) continue;
      h_ncjet->Fill(ncjet);

      //create L-jet list
      vector <Int_t> LjetGenid;
      Int_t nLjetGen=0;
      for(Int_t j=0; j < nJet; j++){
  	if(fabs(jetGenID[j]) == 1 || fabs(jetGenID[j] == 2)){
  	  LjetGenid.push_back(j);
  	  nLjetGen++;
  	}
      }

      // vector <Int_t> Ljet_list;  
      vector <Int_t> Ljet_list;
      Int_t nLjet=0;
      for(Int_t ijet=0; ijet < nJet; ijet++){
  	for(Int_t jj=0; jj < nLjetGen; jj++){
  	  Int_t j = LjetGenid[jj];
  	  Float_t dr = deltaR(jetEta[ijet], jetPhi[ijet], jetGenJetEta[j], jetGenJetPhi[j]);
  	  Float_t dpt = fabs((jetPt[ijet]-jetGenJetPt[j])/jetGenJetPt[j]);
  	  if(dr < 0.2 && dpt < 0.5){
  	    Ljet_list.push_back(ijet);
  	    nLjet++;
  	    break;
  	  }
  	}
      }
      //if(nLjet ==0) continue;
      h_nLjet->Fill(nLjet);
      
      //create c-jet list
      vector <Int_t> bjetGenid;
      Int_t nbjetGen=0;
      for(Int_t j=0; j < nJet; j++){
  	if(fabs(jetGenID[j]) == 5){
  	  bjetGenid.push_back(j);
  	  nbjetGen++;
  	}
      }
      
      vector <Int_t> bjet_list;
      Int_t nbjet=0;
      for(Int_t ijet=0; ijet < nJet; ijet++){
  	for(Int_t jj=0; jj < nbjetGen; jj++){
  	  Int_t j = bjetGenid[jj];
  	  Float_t dr = deltaR(jetEta[ijet], jetPhi[ijet], jetGenJetEta[j], jetGenJetPhi[j]);
  	  Float_t dpt = fabs((jetPt[ijet]-jetGenJetPt[j])/jetGenJetPt[j]);
  	  if(dr < 0.2 && dpt < 0.5){
  	    bjet_list.push_back(ijet);
  	    nbjet++;
  	    break;
  	  }
  	}
      }
      //if(nbjet ==0) continue;
      h_nbjet->Fill(nbjet);
      
      //Int_t ncjet=0;
      //Int_t nLjet=0;
      
      if(isMatched ==1){
  	TLorentzVector pjP4_tmp, jetP4, cjetP4, LjetP4, bjetP4;
  	//Int_t ncjet=0;
  	//Int_t nLjet=0;
  	//Int_t nbjet=0;
	
  	for(Int_t iipho=0; iipho < nrealpho; iipho++){
  	  if(fabs(phoEta[realpho_list[iipho]]) > 1.4442) continue;
  	  if(phoEt[realpho_list[iipho]] < 1000) continue;
  	  phoP4.SetPtEtaPhiM(phoEt[realpho_list[iipho]], phoEta[realpho_list[iipho]], phoPhi[realpho_list[iipho]], 0.);

  	  for(Int_t iijet=0; iijet < nrealjet; iijet++){
  	    if(fabs(jetPt[realjet_list[iijet]]) < 1000) continue;
  	    if(iijet == 0) {
  	      h_realjet_HPT->Fill(jetPt[realjet_list[iijet]]);
  	    }
  	    else{
  	      h_realjet_SndPT->Fill(jetPt[realjet_list[iijet]]);
  	    }
	    
	    
  	    jetP4.SetPtEtaPhiE(jetPt[realjet_list[iijet]], jetEta[realjet_list[iijet]], jetPhi[realjet_list[iijet]], jetEn[realjet_list[iijet]]);
  	    //if(fabs(mcMomPID[iipho]) != fabs(jetGenMomID[iijet])) continue;
  	    pjP4_tmp = phoP4 + jetP4;

	    
  	    if(phoP4.DeltaR(jetP4) < 0.4) continue;
  	    h_pjmass_real->Fill(pjP4_tmp.M());
	    
  	    if(iipho ==0 && iijet == 0){
  	      h_pjmass_HPT->Fill(pjP4_tmp.M());
  	    }
  	    else{
  	      h_pjmass_Sndpho->Fill(pjP4_tmp.M());
  	    }
	    
  	    if(fabs(jetGenID[realjet_list[iijet]]) == 1 || fabs(jetGenID[realjet_list[iijet]]) ==2){
  	      h_pjmass_phoLjet->Fill(pjP4_tmp.M());
  	      h_nLjet->Fill(jetPt[realjet_list[iijet]]);
  	    }
  	    else if(fabs(jetGenID[realjet_list[iijet]]) == 4){
  	      h_pjmass_phocjet->Fill(pjP4_tmp.M());
  	      h_ncjet->Fill(jetPt[realjet_list[iijet]]);
  	    }
  	    else if(fabs(jetGenID[realjet_list[iijet]]) == 5){
  	      h_pjmass_phobjet->Fill(pjP4_tmp.M());
  	      h_nbjet->Fill(jetPt[realjet_list[iijet]]);
  	    }
  	    else{
  	      h_pjmass_phopascjet->Fill(pjP4_tmp.M());
  	    }

	    

  	    if(iipho ==0 && iijet ==0 && jetGenID[realjet_list[iijet]] == 4) h_pjmass_HPTcjet->Fill(pjP4_tmp.M());
  	    if(iipho ==0 && iijet ==0 && (jetGenID[realjet_list[iijet]] ==2 || jetGenID[realjet_list[iijet]] ==1)) h_pjmass_HPTLjet->Fill(pjP4_tmp.M());
  	    if(iipho ==0 && iijet ==1 && jetGenID[realjet_list[iijet]] ==4){
  	      h_pjmass_Sndcjet->Fill(pjP4_tmp.M());
  	    }
	    
  	    //else{
  	    //h_pjmass_SndLjet->Fill(pjP4_tmp.M());
  	    //}
  	    //if(iipho ==0 && jetGenID[realjet_list[iijet]] > 6) h_pjmass_SndLjet->Fill(pjP4_tmp.M());
	    
	    
  	  }
  	}
	

  	/*
  	for(Int_t iipho=0; iipho < nrealpho; iipho++){
  	  if(fabs(phoEta[realpho_list[iipho]]) > 1.4442) continue;
  	  if(phoEt[realpho_list[iipho]] < 1000) continue;
  	  phoP4.SetPtEtaPhiM(phoEt[realpho_list[iipho]], phoEta[realpho_list[iipho]], phoPhi[realpho_list[iipho]], 0.);

  	  for(Int_t iijet=0; iijet < ncjet; iijet++){
  	    if(jetPt[cjet_list[iijet]] < 1000) continue;
  	    cjetP4.SetPtEtaPhiM(jetPt[cjet_list[iijet]], jetEta[cjet_list[iijet]], jetPhi[cjet_list[iijet]], jetEn[cjet_list[iijet]]);

  	    pjP4_tmp = phoP4 + cjetP4;
  	    h_pjmass_phocjet->Fill(pjP4_tmp.M());

  	    if(iipho ==0 && iijet ==0) h_pjmass_HPTcjet->Fill(pjP4_tmp.M()); 
  	  }
  	}

  	for(Int_t iipho=0; iipho < nrealpho; iipho++){
  	  if(fabs(phoEta[realpho_list[iipho]]) > 1.4442) continue;
  	  if(phoEt[realpho_list[iipho]] < 1000) continue;
  	  phoP4.SetPtEtaPhiM(phoEt[realpho_list[iipho]], phoEta[realpho_list[iipho]], phoPhi[realpho_list[iipho]], 0.);

  	  for(Int_t iijet=0; iijet < nLjet; iijet++){
  	    if(jetPt[Ljet_list[iijet]] < 1000) continue;
  	    LjetP4.SetPtEtaPhiM(jetPt[Ljet_list[iijet]], jetEta[Ljet_list[iijet]], jetPhi[Ljet_list[iijet]], jetEn[Ljet_list[iijet]]);

  	    pjP4_tmp = phoP4 + LjetP4;
  	    h_pjmass_phoLjet->Fill(pjP4_tmp.M());

  	    if(iipho ==0 && iijet ==0) h_pjmass_HPTLjet->Fill(pjP4_tmp.M()); 
  	  }
  	}

  	for(Int_t iipho=0; iipho < nrealpho; iipho++){
  	  if(fabs(phoEta[realpho_list[iipho]]) > 1.4442) continue;
  	  if(phoEt[realpho_list[iipho]] < 1000) continue;
  	  phoP4.SetPtEtaPhiM(phoEt[realpho_list[iipho]], phoEta[realpho_list[iipho]], phoPhi[realpho_list[iipho]], 0.);

  	  for(Int_t iijet=0; iijet < nbjet; iijet++){
  	    if(jetPt[bjet_list[iijet]] < 1000) continue;
  	    bjetP4.SetPtEtaPhiM(jetPt[bjet_list[iijet]], jetEta[bjet_list[iijet]], jetPhi[bjet_list[iijet]], jetEn[bjet_list[iijet]]);

  	    pjP4_tmp = phoP4 + bjetP4;
  	    h_pjmass_phobjet->Fill(pjP4_tmp.M());

  	    if(iipho ==0 && iijet ==0) h_pjmass_HPTbjet->Fill(pjP4_tmp.M()); 
  	  }
  	}
  	*/
	

  	//Float_t c_eff=ncjet/nLjet;
   	//h_ctageff->Fill();
   	//printf("Efficiency : %i", ncjet)	
      }

      
      

      //fill matched phoId & jetId par
      for(Int_t ii=0; ii < nrealpho; ii++){
	Int_t ipho = realpho_list[ii];
	eleVeto = phoEleVeto[ipho];
	HoverE = phoHoverE[ipho];
	sieieFull5x5 = phoSigmaIEtaIEtaFull5x5[ipho];
	sieipFull5x5 = phoSigmaIEtaIPhiFull5x5[ipho];
	sipipFull5x5 = phoSigmaIPhiIPhiFull5x5[ipho];
	r9Full5x5 = phoR9Full5x5[ipho];
	
	/*
	  chIso = phoPFChIso[ipho];
	  phoIso = phoPFPhoIso[ipho];
	  nhIso = phoPFNeuIso[ipho];
	*/
	
	
	chIso = chIso_rc[ipho];
	phoIso = phoIso_rc[ipho];
	nhIso = nhIso_rc[ipho];
	
	
	//phoR9 = phoR9[ipho];
	
	
	for(Int_t jj=0; jj < nrealjet; jj++){
	  Int_t ijet = realjet_list[jj];
	  if(jj != 0) continue;
	  jetCEF_ = jetCEF[ijet];
	  jetNEF_ = jetNEF[ijet];
	  jetCHF_ = jetCHF[ijet];
	  jetNHF_ = jetNHF[ijet];
	  jetNCH_ = jetNCH[ijet];
	  jetNNP_ = jetNNP[ijet];
	}
	outtree_->Fill();
	
	
	h_realphoEB_pt->Fill(phoEt[ipho]);
	h_chIso_rc->Fill(chIso);
	h_phoIso_rc->Fill(phoIso);
	h_nhIso_rc->Fill(nhIso);

	//Loose EB cut (Fall2017)
	chIsoEB_L = 1.694;
	phoIsoEB_L = 2.876 + 0.004017*phoEt[ipho];
	nhIsoEB_L = 24.032 + 0.01512*phoEt[ipho] + (2.259*phoEt[ipho]*phoEt[ipho])/100000;
	if( HoverE < 0.04596 && sieieFull5x5 < 0.0106 && chIso < chIsoEB_L && phoIso < phoIsoEB_L && nhIso < nhIsoEB_L ){
	  h_realphoEB_pt_L->Fill(phoEt[ipho]);
	}

	//Medium EB cut (Fall2017)
	chIsoEB_M = 1.141;
	phoIsoEB_M = 2.08 + 0.004017*phoEt[ipho];
	nhIsoEB_M = 1.189 + 0.01512*phoEt[ipho] + (2.259*phoEt[ipho]*phoEt[ipho])/100000;
	if( HoverE < 0.02197 && sieieFull5x5 < 0.01015 && chIso < chIsoEB_M && phoIso < phoIsoEB_M && nhIso < nhIsoEB_M ){
	  h_realphoEB_pt_M->Fill(phoEt[ipho]);
       	}

	//Tight EB cut (Fall2017)
	chIsoEB_T = 0.65;
	phoIsoEB_T = 2.044 + 0.004017*phoEt[ipho];
	nhIsoEB_T = 0.317 + 0.01512*phoEt[ipho] + 2.259*phoEt[ipho]*phoEt[ipho]/100000;
	if( HoverE < 0.02148 && sieieFull5x5 < 0.0996 && chIso < chIsoEB_T && phoIso < phoIsoEB_T && nhIso < nhIsoEB_T ){
	  h_realphoEB_pt_T->Fill(phoEt[ipho]);
	}
      }
      
    }
  }
  //****************END LOOP**********************//

  //****************Terminate*********************//
  fout_->cd();
  outtree_->Write();

  h_npho->Write();
  h_njet->Write();
  //h_phoEt->Write();

  h_nMCpho->Write();
  h_dr_pho->Write();
  h_dpt_pho->Write();
  h_dptdr_pho->Write();
  h_nrealpho->Write();
  h_ngenpho->Write();
  
  h_njetGen->Write();
  h_dr_jet->Write();
  h_dpt_jet->Write();
  h_nrealjet->Write();
  h_dptdr_jet->Write();
  h_ncjet->Write();
  h_nLjet->Write();
  h_nbjet->Write();
  
  h_MCQBHMass->Write();
  h_pjmass_real->Write();
  h_pjmass_phocjet->Write();
  h_pjmass_phopascjet->Write();
  h_pjmass_phoLjet->Write();
  h_pjmass_phobjet->Write();
  
  h_pjmass_HPT->Write();
  h_pjmass_Sndpho ->Write();
  h_pjmass_HPTcjet->Write();
  h_pjmass_HPTLjet->Write();
  //h_pjmass_HPTbjet->Write();
  h_pjmass_Sndcjet->Write();
  h_pjmass_SndLjet->Write();

  h_realjet_HPT->Write();
  h_realjet_SndPT->Write();
  h_dr_genphojet->Write();
  //h_ctageff->Write();

  h_realphoEB_pt->Write();
  h_chIso_rc->Write();
  h_phoIso_rc->Write();
  h_nhIso_rc->Write();
  h_realphoEB_pt_L->Write();
  h_realphoEB_pt_M->Write();
  h_realphoEB_pt_T->Write();
  
  fout_->Close();
  fprintf(stderr, "Processed all events\n");
  
}
