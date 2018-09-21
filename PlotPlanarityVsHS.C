#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TPad.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TLatex.h>
#include <TDirectory.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TH2F.h>

#endif

const int xaxistitle[] = {4,5,6,7,8,9};
//const int xaxistitle[] = {2,3,4,5,6,7,8,9};

int PlotPlanarityVsHS(const char* _path="");
void ComputePlanarityRMSandMean(TGraph* graph, double &planarity, double &zRMS, double &zmean);

int PlotPlanarityVsHS(const char* _path) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const int nFiles = sizeof(xaxistitle)/sizeof(int);
  TGraph* gResToNomPlaneHSleft_Z[nFiles];
  TGraph* gResToNomPlaneHSright_Z[nFiles];
  TH1F* hResToNominalHSleft_X[nFiles];
  TH1F* hResToNominalHSright_X[nFiles];
  TH1F* hResToNominalHSleft_Y[nFiles];
  TH1F* hResToNominalHSright_Y[nFiles];
  TGraph* gResToNomPlaneStaveHSleft_Z[nFiles];
  TGraph* gResToNomPlaneStaveHSright_Z[nFiles];
  TH1F* hResToNominalStave_HSleft_X[nFiles];
  TH1F* hResToNominalStave_HSright_X[nFiles];
  TH1F* hResToNominalStave_HSleft_Y[nFiles];
  TH1F* hResToNominalStave_HSright_Y[nFiles];
  double planarityZ_HSleft[nFiles];
  double RMSZ_HSleft[nFiles];
  double meanZ_HSleft[nFiles];
  double planarityZ_HSright[nFiles];
  double RMSZ_HSright[nFiles];
  double meanZ_HSright[nFiles];
  double planarityZ_StaveHSleft[nFiles];
  double RMSZ_StaveHSleft[nFiles];
  double meanZ_StaveHSleft[nFiles];
  double planarityZ_StaveHSright[nFiles];
  double RMSZ_StaveHSright[nFiles];
  double meanZ_StaveHSright[nFiles];

  const TString sBasePath(_path);
  for(int iFile=0; iFile<nFiles; iFile++) {
    TString inFileName;
    inFileName.Form("%s/T-OL-Stave-%03d_MeasAndExtrap.root", sBasePath.Data(), xaxistitle[iFile]);
    Printf("%s", inFileName.Data());
    TFile* infile = TFile::Open(inFileName.Data());
    if(!infile) return -1;
    gResToNomPlaneHSleft_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneHSleft_Z");
    if(!gResToNomPlaneHSleft_Z[iFile]) return -2;
    gResToNomPlaneHSright_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneHSright_Z");
    if(!gResToNomPlaneHSright_Z[iFile]) return -3;
    hResToNominalHSleft_X[iFile] = (TH1F*)infile->Get("hResToNominalHSleft_X");
    if(!hResToNominalHSleft_X[iFile]) return -4;
    hResToNominalHSleft_X[iFile]->SetDirectory(0);
    hResToNominalHSright_X[iFile] = (TH1F*)infile->Get("hResToNominalHSright_X");
    if(!hResToNominalHSright_X[iFile]) return -5;
    hResToNominalHSright_X[iFile]->SetDirectory(0);
    hResToNominalHSleft_Y[iFile] = (TH1F*)infile->Get("hResToNominalHSleft_Y");
    if(!hResToNominalHSleft_Y[iFile]) return -6;
    hResToNominalHSleft_Y[iFile]->SetDirectory(0);
    hResToNominalHSright_Y[iFile] = (TH1F*)infile->Get("hResToNominalHSright_Y");
    if(!hResToNominalHSright_Y[iFile]) return -7;
    hResToNominalHSright_Y[iFile]->SetDirectory(0);
    ComputePlanarityRMSandMean(gResToNomPlaneHSleft_Z[iFile],planarityZ_HSleft[iFile],RMSZ_HSleft[iFile],meanZ_HSleft[iFile]);
    ComputePlanarityRMSandMean(gResToNomPlaneHSright_Z[iFile],planarityZ_HSright[iFile],RMSZ_HSright[iFile],meanZ_HSright[iFile]);
    gResToNomPlaneStaveHSleft_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneStaveHSleft_Z");
    if(!gResToNomPlaneStaveHSleft_Z[iFile]) return -8;
    gResToNomPlaneStaveHSright_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneStaveHSright_Z");
    if(!gResToNomPlaneStaveHSright_Z[iFile]) return -9;
    hResToNominalStave_HSleft_X[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSleft_X");
    if(!hResToNominalStave_HSleft_X[iFile]) return -10;
    hResToNominalStave_HSleft_X[iFile]->SetDirectory(0);
    hResToNominalStave_HSright_X[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSright_X");
    if(!hResToNominalStave_HSright_X[iFile]) return -11;
    hResToNominalStave_HSright_X[iFile]->SetDirectory(0);
    hResToNominalStave_HSleft_Y[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSleft_Y");
    if(!hResToNominalStave_HSleft_Y[iFile]) return -12;
    hResToNominalStave_HSleft_Y[iFile]->SetDirectory(0);
    hResToNominalStave_HSright_Y[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSright_Y");
    if(!hResToNominalStave_HSright_Y[iFile]) return -13;
    hResToNominalStave_HSright_Y[iFile]->SetDirectory(0);
    ComputePlanarityRMSandMean(gResToNomPlaneStaveHSleft_Z[iFile],planarityZ_StaveHSleft[iFile],RMSZ_StaveHSleft[iFile],meanZ_StaveHSleft[iFile]);
    ComputePlanarityRMSandMean(gResToNomPlaneStaveHSright_Z[iFile],planarityZ_StaveHSright[iFile],RMSZ_StaveHSright[iFile],meanZ_StaveHSright[iFile]);

    infile->Close();

  }
//****************************************|||||||||RMS||||||||*****************************************************
  TH1F* hRMSZVsHSleft = new TH1F("hRMSZVsHSleft","HSleft vs. RMS(z);;RMS(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  TH1F* hRMSXVsHSleft = new TH1F("hRMSXVsHSleft","HSleft vs. RMS(x);;RMS(x) (#mum)",nFiles,-0.5,nFiles-0.5);
  TH1F* hRMSYVsHSleft = new TH1F("hRMSYVsHSleft","HSleft vs. RMS(y);;RMS(y) (#mum)",nFiles,-0.5,nFiles-0.5);
  TH1F* hRMSZVsHSright = new TH1F("hRMSZVsHSright","HSright vs. RMS(z);;RMS(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  TH1F* hRMSXVsHSright = new TH1F("hRMSXVsHSright","HSright vs. RMS(x);;RMS(x) (#mum)",nFiles,-0.5,nFiles-0.5);
  TH1F* hRMSYVsHSright = new TH1F("hRMSYVsHSright","HSright vs. RMS(y);;RMS(y) (#mum)",nFiles,-0.5,nFiles-0.5);
  TH1F* hRMSZVsStaveHSleft = new TH1F("hRMSZVsStaveHSleft","StaveHSleft vs. RMS(z);;RMS(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hRMSZVsStaveHSleft->SetLineColor(kRed);
  TH1F* hRMSXVsStaveHSleft = new TH1F("hRMSXVsStaveHSleft","StaveHSleft vs. RMS(x);StaveHSL;RMS(x) (#mum)",nFiles,-0.5,nFiles-0.5);
  hRMSXVsStaveHSleft->SetLineColor(kRed);
  TH1F* hRMSYVsStaveHSleft = new TH1F("hRMSYVsStaveHSleft","StaveHSleft vs. RMS(y);StaveHSL;RMS(y) (#mum)",nFiles,-0.5,nFiles-0.5);
  hRMSYVsStaveHSleft->SetLineColor(kRed);
  TH1F* hRMSZVsStaveHSright = new TH1F("hRMSZVsStaveHSright","StaveHSright vs. RMS(z);StaveHSR;RMS(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hRMSZVsStaveHSright->SetLineColor(kRed);
  TH1F* hRMSXVsStaveHSright = new TH1F("hRMSXVsStaveHSright","StaveHSright vs. RMS(x);StaveHSR;RMS(x) (#mum)",nFiles,-0.5,nFiles-0.5);
  hRMSXVsStaveHSright->SetLineColor(kRed);
  TH1F* hRMSYVsStaveHSright = new TH1F("hRMSYVsStaveHSright","StaveHSright vs. RMS(y);StaveHSR;RMS(y) (#mum)",nFiles,-0.5,nFiles-0.5);
  hRMSYVsStaveHSright->SetLineColor(kRed);
//************************************************PLANARITY******************************************************
  TH1F* hPlanZVsHSleft = new TH1F("hPlanZVsHSleft","HSleft vs. Planarity(z);;Planarity(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hPlanZVsHSleft->SetLineColor(kGreen+2);
  TH1F* hPlanZVsHSright = new TH1F("hPlanZVsHSright","HSright vs. Planarity(z);;Planarity(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hPlanZVsHSright->SetLineColor(kGreen+2);
  TH1F* hPlanZVsStaveHSleft = new TH1F("hPlanZVsStaveHSleft","StaveHSleft vs. Planarity(z);;Planarity(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hPlanZVsStaveHSleft->SetLineColor(kMagenta+1);
  TH1F* hPlanZVsStaveHSright = new TH1F("hPlanZVsStaveHSright","StaveHSright vs. Planarity(z);StaveHSR;Planarity(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hPlanZVsStaveHSright->SetLineColor(kMagenta+1);
  //************************************************MEAN******************************************************
  TH1F* hMeanZVsHSleft = new TH1F("hMeanZVsHSleft","HSleft vs. Mean(z);;Mean(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hMeanZVsHSleft->SetLineColor(kBlack);
  TH1F* hMeanZVsHSright = new TH1F("hMeanZVsHSright","HSright vs. Mean(z);;Mean(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hMeanZVsHSright->SetLineColor(kBlack);
  TH1F* hMeanZVsStaveHSleft = new TH1F("hMeanZVsStaveHSleft","StaveHSleft vs. Mean(z);;Mean(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hMeanZVsStaveHSleft->SetLineColor(kRed);
  TH1F* hMeanZVsStaveHSright = new TH1F("hMeanZVsStaveHSright","StaveHSright vs. Mean(z);StaveHSR;Mean(z) (#mum)",nFiles,-0.5,nFiles-0.5);
  hMeanZVsStaveHSright->SetLineColor(kRed);

  for(int iFile=0; iFile<nFiles; iFile++) {
    //****************************************|||||||||RMS||||||||*****************************************************
    hRMSZVsHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("HSL%d",xaxistitle[iFile]));
    hRMSXVsHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("HSL%d",xaxistitle[iFile]));
    hRMSYVsHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("HSL%d",xaxistitle[iFile]));
    hRMSZVsHSleft->SetBinContent(iFile+1,RMSZ_HSleft[iFile]);
    hRMSXVsHSleft->SetBinContent(iFile+1,hResToNominalHSleft_X[iFile]->GetRMS());
    hRMSYVsHSleft->SetBinContent(iFile+1,hResToNominalHSleft_Y[iFile]->GetRMS());
    hRMSZVsHSright->GetXaxis()->SetBinLabel(iFile+1,Form("HSR%d",xaxistitle[iFile]));
    hRMSXVsHSright->GetXaxis()->SetBinLabel(iFile+1,Form("HSR%d",xaxistitle[iFile]));
    hRMSYVsHSright->GetXaxis()->SetBinLabel(iFile+1,Form("HSR%d",xaxistitle[iFile]));
    hRMSZVsHSright->SetBinContent(iFile+1,RMSZ_HSright[iFile]);
    hRMSXVsHSright->SetBinContent(iFile+1,hResToNominalHSright_X[iFile]->GetRMS());
    hRMSYVsHSright->SetBinContent(iFile+1,hResToNominalHSright_Y[iFile]->GetRMS());
    hRMSZVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSL%d",xaxistitle[iFile]));
    hRMSXVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSL%d",xaxistitle[iFile]));
    hRMSYVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSL%d",xaxistitle[iFile]));
    hRMSZVsStaveHSleft->SetBinContent(iFile+1,RMSZ_StaveHSleft[iFile]);
    hRMSXVsStaveHSleft->SetBinContent(iFile+1,hResToNominalStave_HSleft_X[iFile]->GetRMS());
    hRMSYVsStaveHSleft->SetBinContent(iFile+1,hResToNominalStave_HSleft_Y[iFile]->GetRMS());
    hRMSZVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSR%d",xaxistitle[iFile]));
    hRMSXVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSR%d",xaxistitle[iFile]));
    hRMSYVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSR%d",xaxistitle[iFile]));
    hRMSZVsStaveHSright->SetBinContent(iFile+1,RMSZ_StaveHSright[iFile]);
    hRMSXVsStaveHSright->SetBinContent(iFile+1,hResToNominalStave_HSright_X[iFile]->GetRMS());
    hRMSYVsStaveHSright->SetBinContent(iFile+1,hResToNominalStave_HSright_Y[iFile]->GetRMS());
    //************************************************PLANARITY******************************************************
    hPlanZVsHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("HSL%d",xaxistitle[iFile]));
    hPlanZVsHSleft->SetBinContent(iFile+1,planarityZ_HSleft[iFile]);
    hPlanZVsHSright->GetXaxis()->SetBinLabel(iFile+1,Form("HSR%d",xaxistitle[iFile]));
    hPlanZVsHSright->SetBinContent(iFile+1,planarityZ_HSright[iFile]);
    hPlanZVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSL%d",xaxistitle[iFile]));
    hPlanZVsStaveHSleft->SetBinContent(iFile+1,planarityZ_StaveHSleft[iFile]);
    hPlanZVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSR%d",xaxistitle[iFile]));
    hPlanZVsStaveHSright->SetBinContent(iFile+1,planarityZ_StaveHSright[iFile]);
    //************************************************MEAN******************************************************
    hMeanZVsHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("HSL%d",xaxistitle[iFile]));
    hMeanZVsHSleft->SetBinContent(iFile+1,meanZ_HSleft[iFile]);
    hMeanZVsHSright->GetXaxis()->SetBinLabel(iFile+1,Form("HSR%d",xaxistitle[iFile]));
    hMeanZVsHSright->SetBinContent(iFile+1,meanZ_HSright[iFile]);
    hMeanZVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSL%d",xaxistitle[iFile]));
    hMeanZVsStaveHSleft->SetBinContent(iFile+1,meanZ_StaveHSleft[iFile]);
    hMeanZVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,Form("StaveHSR%d",xaxistitle[iFile]));
    hMeanZVsStaveHSright->SetBinContent(iFile+1,meanZ_StaveHSright[iFile]);

}

//************************************************MEAN******************************************************
  TLegend* legMeanZleft = new TLegend(0.885,0.7,0.7,0.825);
  legMeanZleft->SetTextSize(0.04);
  legMeanZleft->AddEntry(hMeanZVsHSleft,"HSL","l");
  legMeanZleft->AddEntry(hMeanZVsStaveHSleft,"StaveHSL","l");
  TLegend* legMeanZright = new TLegend(0.885,0.7,0.7,0.825);
  legMeanZright->SetTextSize(0.04);
  legMeanZright->AddEntry(hMeanZVsHSright,"HSR","l");
  legMeanZright->AddEntry(hMeanZVsStaveHSright,"StaveHSR","l");
  TCanvas* cMeanZ_HSLeftRight_StaveHSLeftRight = new TCanvas("cMeanZ_HSLeftRight_StaveHSLeftRight","",1150,450);
  cMeanZ_HSLeftRight_StaveHSLeftRight->Divide(2,1);
  cMeanZ_HSLeftRight_StaveHSLeftRight->cd(1);
  hMeanZVsHSleft->GetYaxis()->SetRangeUser(-150.,100.);
  hMeanZVsHSleft->Draw();
  hMeanZVsStaveHSleft->Draw("same");
  legMeanZleft->Draw("same");
  cMeanZ_HSLeftRight_StaveHSLeftRight->cd(2);
  hMeanZVsHSright->GetYaxis()->SetRangeUser(-150.,250.);
  hMeanZVsHSright->Draw();
  hMeanZVsStaveHSright->Draw("same");
  legMeanZright->Draw("same");
  //************************************************PLANARITY******************************************************
  TLegend* legPlanZleft = new TLegend(0.885,0.56,0.7,0.45);
  legPlanZleft->SetTextSize(0.04);
  legPlanZleft->AddEntry(hPlanZVsHSleft,"HSL","l");
  legPlanZleft->AddEntry(hPlanZVsStaveHSleft,"StaveHSL","l");
  TLegend* legPlanZright = new TLegend(0.885,0.7,0.7,0.825);
  legPlanZright->SetTextSize(0.04);
  legPlanZright->AddEntry(hPlanZVsHSright,"HSR","l");
  legPlanZright->AddEntry(hPlanZVsStaveHSright,"StaveHSR","l");
  TCanvas* cPlanZ_HSLeftRight_StaveHSLeftRight = new TCanvas("cPlanZ_HSLeftRight_StaveHSLeftRight","",1150,450);
  cPlanZ_HSLeftRight_StaveHSLeftRight->Divide(2,1);
  cPlanZ_HSLeftRight_StaveHSLeftRight->cd(1);
  hPlanZVsHSleft->GetYaxis()->SetRangeUser(0.,450.);
  hPlanZVsHSleft->Draw();
  hPlanZVsStaveHSleft->Draw("same");
  legPlanZleft->Draw("same");
  cPlanZ_HSLeftRight_StaveHSLeftRight->cd(2);
  hPlanZVsHSright->GetYaxis()->SetRangeUser(0.,500.);
  hPlanZVsHSright->Draw();
  hPlanZVsStaveHSright->Draw("same");
  legPlanZright->Draw("same");
  //****************************************|||||||||RMS||||||||*****************************************************
  TLegend* legRMSZright = new TLegend(0.885,0.7,0.7,0.825);
  legRMSZright->SetTextSize(0.04);
  legRMSZright->AddEntry(hRMSZVsHSright,"HSR","l");
  legRMSZright->AddEntry(hRMSZVsStaveHSright,"StaveHSR","l");
  TLegend* legRMSXright = new TLegend(0.7,0.7,0.5,0.825);
  legRMSXright->SetTextSize(0.04);
  legRMSXright->AddEntry(hRMSXVsHSright,"HSR","l");
  legRMSXright->AddEntry(hRMSXVsStaveHSright,"StaveHSR","l");
  TLegend* legRMSYright = new TLegend(0.885,0.7,0.7,0.825);
  legRMSYright->SetTextSize(0.04);
  legRMSYright->AddEntry(hRMSYVsHSright,"HSR","l");
  legRMSYright->AddEntry(hRMSYVsStaveHSright,"StaveHSR","l");
  TCanvas* cRMSZXY_HSright_StaveHSright = new TCanvas("cRMSZXY_HSright_StaveHSright","",900,900);
  cRMSZXY_HSright_StaveHSright->Divide(2,2);
  cRMSZXY_HSright_StaveHSright->cd(1);
  hRMSZVsHSright->GetYaxis()->SetRangeUser(0.,100.);
  hRMSZVsHSright->Draw();
  hRMSZVsStaveHSright->Draw("same");
  legRMSZright->Draw("same");
  cRMSZXY_HSright_StaveHSright->cd(2);
  hRMSXVsHSright->GetYaxis()->SetRangeUser(0.,120.);
  hRMSXVsHSright->Draw();
  hRMSXVsStaveHSright->Draw("same");
  legRMSXright->Draw("same");
  cRMSZXY_HSright_StaveHSright->cd(3);
  hRMSYVsHSright->GetYaxis()->SetRangeUser(0.,35.);
  hRMSYVsHSright->Draw();
  hRMSYVsStaveHSright->Draw("same");
  legRMSYright->Draw("same");

  TLegend* legRMSZleft = new TLegend(0.885,0.56,0.7,0.45);
  legRMSZleft->SetTextSize(0.04);
  legRMSZleft->AddEntry(hRMSZVsHSleft,"HSL","l");
  legRMSZleft->AddEntry(hRMSZVsStaveHSleft,"StaveHSL","l");
  TLegend* legRMSXleft = new TLegend(0.885,0.7,0.7,0.825);
  legRMSXleft->SetTextSize(0.04);
  legRMSXleft->AddEntry(hRMSXVsHSleft,"HSL","l");
  legRMSXleft->AddEntry(hRMSXVsStaveHSleft,"StaveHSL","l");
  TLegend* legRMSYleft = new TLegend(0.885,0.7,0.7,0.825);
  legRMSYleft->SetTextSize(0.04);
  legRMSYleft->AddEntry(hRMSYVsHSleft,"HSL","l");
  legRMSYleft->AddEntry(hRMSYVsStaveHSleft,"StaveHSL","l");
  TCanvas* cRMSZXY_HSleft_StaveHSleft = new TCanvas("cRMSZXY_HSleft_StaveHSleft","",900,900);
  cRMSZXY_HSleft_StaveHSleft->Divide(2,2);
  cRMSZXY_HSleft_StaveHSleft->cd(1);
  hRMSZVsHSleft->GetYaxis()->SetRangeUser(0.,100.);
  hRMSZVsHSleft->Draw();
  hRMSZVsStaveHSleft->Draw("same");
  legRMSZleft->Draw("same");
  cRMSZXY_HSleft_StaveHSleft->cd(2);
  hRMSXVsHSleft->GetYaxis()->SetRangeUser(0.,100.);
  hRMSXVsHSleft->Draw();
  hRMSXVsStaveHSleft->Draw("same");
  legRMSXleft->Draw("same");
  cRMSZXY_HSleft_StaveHSleft->cd(3);
  hRMSYVsHSleft->GetYaxis()->SetRangeUser(0.,25.);
  hRMSYVsHSleft->Draw();
  hRMSYVsStaveHSleft->Draw("same");
  legRMSYleft->Draw("same");

  return 0;
}

//______________________________________________________________________________________________
void ComputePlanarityRMSandMean(TGraph* graph, double &planarity, double &zRMS, double &zmean)
{

  double zmax = -1.e20;
  double zmin = 1.e20;
  zmean = 0;
  zRMS = 0;

  for(int iEntry=0; iEntry<graph->GetN(); iEntry++) {
    double y,z; //planarity measurements are z vs. y
    graph->GetPoint(iEntry,y,z);
    if(z>zmax) zmax=z;
    if(z<zmin) zmin=z;
    zmean+=z;
  }
  zmean /= graph->GetN();

  for(int iEntry=0; iEntry<graph->GetN(); iEntry++) {
    double y,z; //planarity measurements are z vs. y
    graph->GetPoint(iEntry,y,z);
    zRMS += (z-zmean)*(z-zmean);
  }

  zRMS /= (graph->GetN()-1);
  zRMS = TMath::Sqrt(zRMS);

  planarity = zmax-zmin;
}
