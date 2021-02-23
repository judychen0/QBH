 void xDraw(){
  TFile *fopen = new TFile("output_ADD_ggtree_mc.root");
  
  TCanvas *c1 = new TCanvas("c1");
  TTree *t =(TTree*)fopen->Get("t");

  const char *title;
  //draw TH1F non-logY histograms
  TH1F *H_npho = (TH1F*)fopen->Get("h_npho");
  H_npho->GetXaxis()->SetTitle("photon number");
  H_npho->GetYaxis()->SetTitle("Events");
  H_npho->Draw();
  title = H_npho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  TH1F *H_njet = (TH1F*)fopen->Get("h_njet");
  H_njet->Draw();
  title = H_njet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  /*
  TH1F *H_pjmass = (TH1F*)fopen->Get("h_pjmass");
  H_pjmass->Draw();
  title = H_pjmass->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_nMCpho = (TH1F*)fopen->Get("h_nMCpho");
  H_nMCpho->Draw();
  title = H_nMCpho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_nMCQstar = (TH1F*)fopen->Get("h_nMCQstar");
  H_nMCQstar->Draw();
  title = H_nMCQstar->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_MCQstarMass = (TH1F*)fopen->Get("h_MCQstarMass");
  H_MCQstarMass->Draw();
  title = H_MCQstarMass->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  */

  TH1F *H_nrealpho = (TH1F*)fopen->Get("h_nrealpho");
  H_nrealpho->Draw();
  title = H_nrealpho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  TH1F *H_njetGen = (TH1F*)fopen->Get("h_njetGen");
  H_njetGen->Draw();
  title = H_njetGen->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_nrealjet = (TH1F*)fopen->Get("h_nrealjet");
  H_nrealjet->Draw();
  title = H_nrealjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  /*
  TH1F *H_nrealjet_HPT = (TH1F*)fopen->Get("h_nrealjet_HPT");
  H_nrealjet_HPT->Draw("HIST");
  title = H_nrealjet_HPT->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  */
  
  TH1F *H_pjmass_real = (TH1F*)fopen->Get("h_pjmass_real");
  H_pjmass_real->Draw("HIST");
  title = H_pjmass_real->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  TH1F *H_pjmass_phocjet = (TH1F*)fopen->Get("h_pjmass_phocjet");
  H_pjmass_phocjet->Draw("HIST");
  title = H_pjmass_phocjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_pjmass_phoLjet = (TH1F*)fopen->Get("h_pjmass_phoLjet");
  H_pjmass_phoLjet->Draw("HIST");
  title = H_pjmass_phoLjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_pjmass_phobjet = (TH1F*)fopen->Get("h_pjmass_phobjet");
  H_pjmass_phobjet->Draw("HIST");
  title = H_pjmass_phobjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  
  TH1F *H_pjmass_HPT = (TH1F*)fopen->Get("h_pjmass_HPT");
  H_pjmass_HPT->Draw("HIST");
  title = H_pjmass_HPT->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  TH1F *H_pjmass_Sndpho = (TH1F*)fopen->Get("h_pjmass_Sndpho");
  H_pjmass_Sndpho->Draw("HIST");
  title = H_pjmass_Sndpho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  

  
  /*
  TH1F *H_pjmass_HPTcjet = (TH1F*)fopen->Get("h_pjmass_HPTcjet");
  H_pjmass_HPTcjet->Draw();
  title = H_pjmass_HPTcjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_pjmass_HPTLjet = (TH1F*)fopen->Get("h_pjmass_HPTLjet");
  H_pjmass_HPTLjet->Draw();
  title = H_pjmass_HPTLjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_pjmass_Sndcjet = (TH1F*)fopen->Get("h_pjmass_Sndcjet");
  H_pjmass_Sndcjet->Draw();
  title = H_pjmass_Sndcjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_pjmass_SndLjet = (TH1F*)fopen->Get("h_pjmass_SndLjet");
  H_pjmass_SndLjet->Draw();
  title = H_pjmass_SndLjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  */
  
  //draw multi histograms
  H_pjmass_real->SetFillColorAlpha(kBlack, 0.3);
  H_pjmass_real->Draw("HIST");
  H_pjmass_phoLjet->SetFillColorAlpha(kGreen+1, 0.3);
  H_pjmass_phoLjet->Draw("HISTSAME");
  H_pjmass_phocjet->SetFillColorAlpha(kAzure+1, 0.3);
  H_pjmass_phocjet->Draw("HISTSAME");
  H_pjmass_phobjet->SetFillColorAlpha(kYellow-7, 0.3);
  H_pjmass_phobjet->Draw("HISTSAME");
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", "h_pjmass_jet_relation"));

  H_pjmass_real->SetFillColorAlpha(kBlack, 0.3);
  H_pjmass_real->Draw("HIST");
  H_pjmass_phocjet->SetFillColorAlpha(kAzure+1, 0.3);
  H_pjmass_phocjet->Draw("HISTSAME");
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", "h_pjmass_cjet_relation"));

  H_pjmass_real->Draw();
  H_pjmass_HPT->Draw("HISTSAME");
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", "h_pjmass_pt_relation"));
  
  /*
  H_pjmass_Sndcjet->SetFillColorAlpha(kBlack);
  H_pjmass_Sndcjet->Draw("HIST");
  H_pjmass_HPTcjet->SetFillColorAlpha(kAzure);
  H_pjmass_HPTcjet->Draw("HISTSAME");
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", "h_pjmass_HPTandSnd_cjet_relation"));
  
  
  H_pjmass_SndLjet->SetFillColorAlpha(kBlack);
  H_pjmass_SndLjet->Draw("HIST");
  H_pjmass_HPTLjet->SetFillColorAlpha(kAzure);
  H_pjmass_HPTLjet->Draw("HISTSAME");
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", "h_pjmass_HPTandSnd_Ljet_relation"));

  H_pjmass_HPTLjet->SetFillColorAlpha(kBlack);
  H_pjmass_HPTLjet->Draw("HIST");
  H_pjmass_HPTcjet->SetFillColorAlpha(kAzure);
  H_pjmass_HPTcjet->Draw("HISTSAME");
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", "h_pjmass_HPT_jet_relation"));
  
  H_pjmass_SndLjet->SetFillColorAlpha(kBlack);
  H_pjmass_SndLjet->Draw("HIST");
  H_pjmass_Sndcjet->SetFillColorAlpha(kAzure);
  H_pjmass_Sndcjet->Draw("HISTSAME");
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", "h_pjmass_Snd_jet_relation"));
  */
  
  //draw TH1F logY histograms
  c1->SetLogy();
  c1->Update();
  
  
  TH1F *H_dr_pho = (TH1F*)fopen->Get("h_dr_pho");
  H_dr_pho->Draw();
  title = H_dr_pho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_dpt_pho = (TH1F*)fopen->Get("h_dpt_pho");
  H_dpt_pho->Draw();
  title = H_dpt_pho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_dr_jet = (TH1F*)fopen->Get("h_dr_jet");
  H_dr_jet->Draw();
  title = H_dr_jet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_dpt_jet = (TH1F*)fopen->Get("h_dpt_jet");
  H_dpt_jet->Draw();
  title = H_dpt_jet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_sieieFull5x5 = new TH1F("H_sieieFull5x5", "H_sieieFull5x5", 18, 0.003, 0.012);
  H_sieieFull5x5->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("sieieFull5x5>>H_sieieFull5x5");
  title = H_sieieFull5x5->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  TH1F *H_r9Full5x5 = new TH1F("H_r9Full5x5", "H_r9Full5x5", 5, 0.5, 1);
  H_r9Full5x5->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("r9Full5x5>>H_r9Full5x5");
  title = H_r9Full5x5->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  TH1F *H_HoverE = new TH1F("H_HoverE", "H_HoverE", 13, 0., 0.026);
  H_HoverE->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("HoverE>>H_HoverE");
  title = H_HoverE->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_chIso = new TH1F("H_chIso", "H_chIso", 15, 0, 0.6);
  H_chIso->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("chIso>>H_chIso");
  title = H_chIso->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_phoIso = new TH1F("H_phoIso", "H_phoIso", 50, 0, 10);
  H_phoIso->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("phoIso>>H_phoIso");
  title = H_phoIso->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_nhIso = new TH1F("H_nhIso", "H_nhIso", 29, 0, 81);
  H_nhIso->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("nhIso>>H_nhIso");
  title = H_nhIso->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_jetCEF = new TH1F("H_jetCEF", "H_jetCEF", 40, 0, 1);
  H_jetCEF->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("jetCEF>>H_jetCEF");
  title = H_jetCEF->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_jetNEF = new TH1F("H_jetNEF", "H_jetNEF", 40, 0, 1);
  H_jetNEF->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("jetNEF>>H_jetCE");
  title = H_jetNEF->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_jetCHF = new TH1F("H_jetCHF", "H_jetCHF", 40, 0, 1);
  H_jetCHF->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("jetCHF>>H_jetCHF");
  title = H_jetCHF->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH1F *H_jetNHF = new TH1F("H_jetNHF", "H_jetNHF", 40, 0, 1);
  H_jetNHF->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("jetNHF>>H_jetNHF");
  title = H_jetNHF->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  TH1F *H_jetNCH = new TH1F("H_jetNCH", "H_jetNCH", 100, 0, 100);
  H_jetNCH->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("jetNCH>>H_jetNCH");
  title = H_jetNCH->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  TH1F *H_jetNNP = new TH1F("H_jetNNP", "H_jetNNP", 100, 0, 100);
  H_jetNNP->GetYaxis()->SetRangeUser(.0001, 1000000);
  t->Draw("jetNNP>>H_jetNNP");
  title = H_jetNNP->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
  c1->SetLogy(0);
  c1->SetLogz();
  c1->Update();

  TH2F *H_dptdr_pho = (TH2F*)fopen->Get("h_dptdr_pho");
  H_dptdr_pho->Draw("colz");
  title = H_dptdr_pho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));

  TH2F *H_dptdr_jet = (TH2F*)fopen->Get("h_dptdr_jet");
  H_dptdr_jet->Draw("colz");
  title = H_dptdr_jet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/QBH/graph/%s.png", title));
  
}
