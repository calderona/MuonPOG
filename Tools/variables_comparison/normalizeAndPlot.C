#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

void normalizeAndPlot() {

  TString refproc("DY_25ns_weights"); 
  //TString refproc("Data_50ns_Run2015C_weights"); 

  // plottext:  plot name, min |eta|, max |eta|, x-axis title, y-axis title. 
  // plotnumb:  min range for norm., max range for norm., log y, rebin, under/overflow.   

  //TString plottext[5] = {"STAmuonTime"         ,  ""  , ""   , "Muon time [ns]", "Entries"}; Float_t plotnumb[5] = {  -2.5,    2.5, 1, 1, 0}; 
  //TString plottext[5] = {"STAmuonTimeBarrel"   ,  ""  , ""   , "Muon time [ns]", "Entries"}; Float_t plotnumb[5] = {  -2.5,    2.5, 1, 1, 0}; 
  //TString plottext[5] = {"STAmuonTimeEndcap"   ,  ""  , ""   , "Muon time [ns]", "Entries"}; Float_t plotnumb[5] = {  -2.5,    2.5, 1, 1, 0}; 
  //TString plottext[5] = {"UnbSTAmuonTime"      ,  ""  , ""   , "Muon time [ns]", "Entries"}; Float_t plotnumb[5] = {  -2.5,    2.5, 1, 1, 0}; 
  //TString plottext[5] = {"UnbSTAmuonTimeBarrel",  ""  , ""   , "Muon time [ns]", "Entries"}; Float_t plotnumb[5] = {  -2.5,    2.5, 1, 1, 0}; 
  //TString plottext[5] = {"UnbSTAmuonTimeEndcap",  ""  , ""   , "Muon time [ns]", "Entries"}; Float_t plotnumb[5] = {  -2.5,    2.5, 1, 1, 0}; 

  //TString plottext[5] = {"control/invMass"       , ""   , "", "Inv. mass [GeV]", "Entries"}; Float_t plotnumb[5] = {  85,  105, 1, 1, 0}; 
  //TString plottext[5] = {"control/invMassInRange", ""   , "", "Inv. mass [GeV]", "Entries"}; Float_t plotnumb[5] = {-101, -101, 1, 1, 0}; 

  //TString plottext[5] = {"probePt"             , "0.0", "2.4", "Muon p_{T} [GeV]", "Entries"}; Float_t plotnumb[5] = {-101, -101, 1, 1, 1}; 
  //TString plottext[5] = {"probeEta"            , "0.0", "2.4", "Muon #eta"       , "Entries"}; Float_t plotnumb[5] = {-101, -101, 0, 1, 0}; 
  //TString plottext[5] = {"probePhi"            , "0.0", "2.4", "Muon #phi [rad]" , "Entries"}; Float_t plotnumb[5] = {-101, -101, 0, 2, 0}; 
  TString plottext[5] = {"probeDxy"            , "0.0", "2.4", "Muon d_{xy} [cm]", "Entries"}; Float_t plotnumb[5] = {-101, -101, 1, 1, 1}; 
  //TString plottext[5] = {"probeDz"             , "0.0", "2.4", "Muon d_{z} [cm]" , "Entries"}; Float_t plotnumb[5] = {-101, -101, 1, 1, 1}; 

  //TString plottext[5] = {"matchedStation"      , "0.0", "2.4", "Matched muon stations", "Entries"}; Float_t plotnumb[5] = {-101, -101, 1, 1, 0}; 
  //TString plottext[5] = {"muonValidHitsGLB"    , "0.0", "2.4", "Valid muon hits"      , "Entries"}; Float_t plotnumb[5] = {-101, -101, 0, 1, 0}; 
  //TString plottext[5] = {"PixelHitsTRK"        , "0.0", "2.4", "Pixel hits"           , "Entries"}; Float_t plotnumb[5] = {-101, -101, 1, 1, 0}; 
  //TString plottext[5] = {"TrackerLayersTRK"    , "0.0", "2.4", "Tracker layers"       , "Entries"}; Float_t plotnumb[5] = {-101, -101, 0, 1, 0}; 

  //TString plottext[5] = {"dBetaRelIso"         , "0.0", "2.4", "PF comb. rel. iso"    , "Entries"}; Float_t plotnumb[5] = {-101, -101, 1, 1, 1}; 

  TString etatag(""); 
  if(plottext[1].Length()>0 && plottext[2].Length()>0) 
    etatag = "_fEtaMin" + plottext[1] + "_fEtaMax" + plottext[2]; 

  TFile *ref = TFile::Open((refproc+"/results.root").Data()); 
  if(ref==0) {std::cout << "No ref file!" << std::endl; return;} 

  TFile *trg = TFile::Open("Data_25ns_Run2015C/results.root"); 
  if(trg==0) {std::cout << "No trg file!" << std::endl; return;} 

  TH1F *hr; 
  if(refproc.BeginsWith("Data")) 
    hr = (TH1F*)ref->Get(("Data/" + plottext[0] + "_Data" + etatag).Data()); 
  else 
    hr = (TH1F*)ref->Get(("DY/"+plottext[0]+"_DY" + etatag).Data()); 
  if(hr==0) {std::cout << "No hr histo!" << std::endl; return;} 

  TH1F *ht = (TH1F*)trg->Get(("Data/"+plottext[0]+"_Data" + etatag).Data()); 
  if(ht==0) {std::cout << "No ht histo!" << std::endl; return;} 

  TCanvas *ccc = new TCanvas("ccc", "ccc", 600, 600); 
  ccc->cd(); 
  ccc->SetLogy(Int_t(plotnumb[2])); 

  hr->SetLineWidth(1); 
  hr->SetLineColor(kBlack); 
  if(refproc.BeginsWith("Data")) 
    hr->SetFillColor(kGray); 
  else 
    hr->SetFillColor(kAzure+7); 
  hr->GetXaxis()->SetTitle(plottext[3].Data()); 
  hr->GetYaxis()->SetTitle(plottext[4].Data()); 
  hr->GetXaxis()->SetTitleSize(0.06); 
  hr->GetYaxis()->SetTitleSize(0.06); 
  hr->GetXaxis()->SetTitleOffset(1.05); 
  hr->GetYaxis()->SetTitleOffset(1.30); 
  if(Int_t(plotnumb[4])==1) {
    hr->SetBinContent(1, hr->GetBinContent(1)+hr->GetBinContent(0)) ; 
    hr->SetBinError(1, sqrt(hr->GetBinError(1)*hr->GetBinError(1) + hr->GetBinError(0)*hr->GetBinError(0))); 
    hr->SetBinContent(0, 0) ; 
    hr->SetBinError(0, 0) ; 

    unsigned int lastb = hr->GetNbinsX(); 
    hr->SetBinContent(lastb, hr->GetBinContent(lastb)+hr->GetBinContent(lastb+1)) ; 
    hr->SetBinError(lastb, sqrt(hr->GetBinError(lastb)*hr->GetBinError(lastb) + hr->GetBinError(lastb+1)*hr->GetBinError(lastb+1))); 
    hr->SetBinContent(lastb+1, 0) ; 
    hr->SetBinError(lastb+1, 0) ; 
  } 

  ht->SetMarkerStyle(20); 
  ht->SetMarkerColor(kBlack); 
  //ht->SetLineWidth(1); 
  //ht->SetLineColor(kBlack); 
  //ht->SetLineStyle(7); 
  if(Int_t(plotnumb[4])==1) {
    ht->SetBinContent(1, ht->GetBinContent(1)+ht->GetBinContent(0)) ; 
    ht->SetBinError(1, sqrt(ht->GetBinError(1)*ht->GetBinError(1) + ht->GetBinError(0)*ht->GetBinError(0))); 
    ht->SetBinContent(0, 0) ; 
    ht->SetBinError(0, 0) ; 

    unsigned int lastb = ht->GetNbinsX(); 
    ht->SetBinContent(lastb, ht->GetBinContent(lastb)+ht->GetBinContent(lastb+1)) ; 
    ht->SetBinError(lastb, sqrt(ht->GetBinError(lastb)*ht->GetBinError(lastb) + ht->GetBinError(lastb+1)*ht->GetBinError(lastb+1))); 
    ht->SetBinContent(lastb+1, 0) ; 
    ht->SetBinError(lastb+1, 0) ; 
  } 

  Int_t firstbin = plotnumb[0]<-100. ? 1               : hr->FindBin(plotnumb[0]); 
  Int_t lastbin  = plotnumb[1]<-100. ? hr->GetNbinsX() : hr->FindBin(plotnumb[1]); 
  Float_t hrintegral = hr->Integral(firstbin, lastbin); 
  Float_t htintegral = ht->Integral(firstbin, lastbin); 

  if(hrintegral>0. && htintegral>0.) 
    hr->Scale(htintegral/hrintegral); 

  if(Int_t(plotnumb[3])>1) { 
    hr->Rebin(Int_t(plotnumb[3])); 
    ht->Rebin(Int_t(plotnumb[3])); 
  } 

  //hr->GetYaxis()->SetRangeUser(0.05,20000.); 
  hr->Draw("hist"); 
  ht->Draw("e1same"); 

  TLegend *ll = new TLegend(0.73,0.74,0.93,0.90); 
  ll->SetBorderSize(0); 
  ll->SetLineWidth(0); 
  ll->SetFillColor(0); 
  ll->SetFillStyle(0); 
  if(refproc.BeginsWith("Data")) 
    ll->AddEntry(hr, "#splitline{SingleMuon}{2015C, 50 ns}", "LF"); 
  else 
    ll->AddEntry(hr, "#splitline{Z/#gamma* #rightarrow #mu#mu}{25 ns}", "LF"); 
  ll->AddEntry(hr, "#splitline{SingleMuon}{2015C, 25 ns}", "PE"); 
  ll->Draw(); 
  if(plottext[0].Contains("/")==0) ccc->SaveAs((plottext[0]+".png").Data()); 
  
  return; 
}

