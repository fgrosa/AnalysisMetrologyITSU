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

const int xaxistitle[] = {2,3,4,5,6,7,8,9,10};

int PlotPlanarityVsHS(const char* _path="");
void ComputePlanarityRMSandMean(TGraph* graph, double &planarity, double &zRMS, double &zmean);

int PlotPlanarityVsHS(const char* _path) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.3);

  const int nFiles = sizeof(xaxistitle)/sizeof(int);
  TGraph* gResToNomPlaneHSleft_Z[nFiles];
  TGraph* gResToNomPlaneHSright_Z[nFiles];
  TGraph* gResToAvPlaneHSleft_Z[nFiles];
  TGraph* gResToAvPlaneHSright_Z[nFiles];
  TH1F* hResToNominalHSleft_X[nFiles];
  TH1F* hResToNominalHSright_X[nFiles];
  TH1F* hResToNominalHSleft_Y[nFiles];
  TH1F* hResToNominalHSright_Y[nFiles];
  TGraph* gResToNomPlaneStaveHSleft_Z[nFiles];
  TGraph* gResToNomPlaneStaveHSright_Z[nFiles];
  TGraph* gResToAvPlaneStaveHSleft_Z[nFiles];
  TGraph* gResToAvPlaneStaveHSright_Z[nFiles];
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
  double dummy;

  const TString sBasePath(_path);
  for(int iFile=0; iFile<nFiles; iFile++) {
    TString inFileName;
    inFileName.Form("%s/T-OL-Stave-%03d_MeasAndExtrap.root", sBasePath.Data(), xaxistitle[iFile]);
    Printf("%s", inFileName.Data());
    TFile* infile = TFile::Open(inFileName.Data());
    int _err = 0;
    if(!infile) return --_err;
    gResToNomPlaneHSleft_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneHSleft_Z");
    if(!gResToNomPlaneHSleft_Z[iFile]) return --_err;
    gResToNomPlaneHSright_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneHSright_Z");
    if(!gResToNomPlaneHSright_Z[iFile]) return --_err;
    gResToAvPlaneHSleft_Z[iFile] = (TGraph*)infile->Get("gResToAvPlaneHSleft_Z");
    if(!gResToAvPlaneHSleft_Z[iFile]) return --_err;
    gResToAvPlaneHSright_Z[iFile] = (TGraph*)infile->Get("gResToAvPlaneHSright_Z");
    if(!gResToAvPlaneHSright_Z[iFile]) return --_err;
    hResToNominalHSleft_X[iFile] = (TH1F*)infile->Get("hResToNominalHSleft_X");
    if(!hResToNominalHSleft_X[iFile]) return --_err;
    hResToNominalHSleft_X[iFile]->SetDirectory(0);
    hResToNominalHSright_X[iFile] = (TH1F*)infile->Get("hResToNominalHSright_X");
    if(!hResToNominalHSright_X[iFile]) return --_err;
    hResToNominalHSright_X[iFile]->SetDirectory(0);
    hResToNominalHSleft_Y[iFile] = (TH1F*)infile->Get("hResToNominalHSleft_Y");
    if(!hResToNominalHSleft_Y[iFile]) return --_err;
    hResToNominalHSleft_Y[iFile]->SetDirectory(0);
    hResToNominalHSright_Y[iFile] = (TH1F*)infile->Get("hResToNominalHSright_Y");
    if(!hResToNominalHSright_Y[iFile]) return --_err;
    hResToNominalHSright_Y[iFile]->SetDirectory(0);
    ComputePlanarityRMSandMean(gResToNomPlaneHSleft_Z[iFile],  dummy, dummy,meanZ_HSleft[iFile]);
    ComputePlanarityRMSandMean(gResToNomPlaneHSright_Z[iFile], dummy, dummy,meanZ_HSright[iFile]);
    ComputePlanarityRMSandMean(gResToAvPlaneHSleft_Z[iFile],  planarityZ_HSleft[iFile],  RMSZ_HSleft[iFile],  dummy);
    ComputePlanarityRMSandMean(gResToAvPlaneHSright_Z[iFile], planarityZ_HSright[iFile], RMSZ_HSright[iFile], dummy);
    gResToNomPlaneStaveHSleft_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneStaveHSleft_Z");
    if(!gResToNomPlaneStaveHSleft_Z[iFile]) return --_err;
    gResToNomPlaneStaveHSright_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneStaveHSright_Z");
    if(!gResToNomPlaneStaveHSright_Z[iFile]) return --_err;
    gResToAvPlaneStaveHSleft_Z[iFile] = (TGraph*)infile->Get("gResToAvPlaneStaveHSleft_Z");
    if(!gResToAvPlaneStaveHSleft_Z[iFile]) return --_err;
    gResToAvPlaneStaveHSright_Z[iFile] = (TGraph*)infile->Get("gResToAvPlaneStaveHSright_Z");
    if(!gResToAvPlaneStaveHSright_Z[iFile]) return --_err;;
    hResToNominalStave_HSleft_X[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSleft_X");
    if(!hResToNominalStave_HSleft_X[iFile]) return --_err;
    hResToNominalStave_HSleft_X[iFile]->SetDirectory(0);
    hResToNominalStave_HSright_X[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSright_X");
    if(!hResToNominalStave_HSright_X[iFile]) return --_err;
    hResToNominalStave_HSright_X[iFile]->SetDirectory(0);
    hResToNominalStave_HSleft_Y[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSleft_Y");
    if(!hResToNominalStave_HSleft_Y[iFile]) return --_err;
    hResToNominalStave_HSleft_Y[iFile]->SetDirectory(0);
    hResToNominalStave_HSright_Y[iFile] = (TH1F*)infile->Get("hResToNominalStave_HSright_Y");
    if(!hResToNominalStave_HSright_Y[iFile]) return --_err;
    hResToNominalStave_HSright_Y[iFile]->SetDirectory(0);
    ComputePlanarityRMSandMean(gResToNomPlaneStaveHSleft_Z[iFile],  dummy, dummy, meanZ_StaveHSleft[iFile]);
    ComputePlanarityRMSandMean(gResToNomPlaneStaveHSright_Z[iFile], dummy, dummy, meanZ_StaveHSright[iFile]);
    ComputePlanarityRMSandMean(gResToAvPlaneStaveHSleft_Z[iFile],  planarityZ_StaveHSleft[iFile],  RMSZ_StaveHSleft[iFile],  dummy);
    ComputePlanarityRMSandMean(gResToAvPlaneStaveHSright_Z[iFile], planarityZ_StaveHSright[iFile], RMSZ_StaveHSright[iFile], dummy);

    infile->Close();
  }
  TString xaxis_name("T-OL-Stave"), yaxis_name, hist_title;
//****************************************|||||||||RMS||||||||*****************************************************
  hist_title="HSleft vs. RMS(z)";
  yaxis_name="RMS(z) (#mum)";
  TH1F* hRMSZVsHSleft = new TH1F("hRMSZVsHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSZVsHSleft->SetLineColor(kRed);
  hist_title="HSleft vs. RMS(x)";
  yaxis_name="RMS(x) (#mum)";
  TH1F* hRMSXVsHSleft = new TH1F("hRMSXVsHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSXVsHSleft->SetLineColor(kRed);
  hist_title="HSleft vs. RMS(y)";
  yaxis_name="RMS(y) (#mum)";
  TH1F* hRMSYVsHSleft = new TH1F("hRMSYVsHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSYVsHSleft->SetLineColor(kRed);
  hist_title="HSright vs. RMS(z)";
  yaxis_name="RMS(z) (#mum)";
  TH1F* hRMSZVsHSright = new TH1F("hRMSZVsHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSZVsHSright->SetLineColor(kBlue);
  hist_title="HSright vs. RMS(x)";
  yaxis_name="RMS(x) (#mum)";
  TH1F* hRMSXVsHSright = new TH1F("hRMSXVsHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSZVsHSright->SetLineColor(kBlue);
  hist_title="HSright vs. RMS(y)";
  yaxis_name="RMS(y) (#mum)";
  TH1F* hRMSYVsHSright = new TH1F("hRMSYVsHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSZVsHSright->SetLineColor(kBlue);

  hist_title="StaveHSleft vs. RMS(z)";
  yaxis_name="RMS(z) (#mum)";
  TH1F* hRMSZVsStaveHSleft = new TH1F("hRMSZVsStaveHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSZVsStaveHSleft->SetLineColor(kRed);
  hist_title="StaveHSleft vs. RMS(x)";
  yaxis_name="RMS(x) (#mum)";
  TH1F* hRMSXVsStaveHSleft = new TH1F("hRMSXVsStaveHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSXVsStaveHSleft->SetLineColor(kRed);
  hist_title="StaveHSleft vs. RMS(y)";
  yaxis_name="RMS(y) (#mum)";
  TH1F* hRMSYVsStaveHSleft = new TH1F("hRMSYVsStaveHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSYVsStaveHSleft->SetLineColor(kRed);
  hist_title="StaveHSright vs. RMS(z)";
  yaxis_name="RMS(z) (#mum)";
  TH1F* hRMSZVsStaveHSright = new TH1F("hRMSZVsStaveHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSZVsStaveHSright->SetLineColor(kBlue);
  hist_title="StaveHSright vs. RMS(x)";
  yaxis_name="RMS(x) (#mum)";
  TH1F* hRMSXVsStaveHSright = new TH1F("hRMSXVsStaveHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSXVsStaveHSright->SetLineColor(kBlue);
  hist_title="StaveHSright vs. RMS(y)";
  yaxis_name="RMS(y) (#mum)";
  TH1F* hRMSYVsStaveHSright = new TH1F("hRMSYVsStaveHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hRMSYVsStaveHSright->SetLineColor(kBlue);

//************************************************PLANARITY******************************************************
  hist_title="HSleft vs. Planarity(z)";
  yaxis_name="Planarity(z) (#mum)";
  TH1F* hPlanZVsHSleft = new TH1F("hPlanZVsHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hPlanZVsHSleft->SetLineColor(kGreen+2);
  hist_title="HSright vs. Planarity(z)";
  yaxis_name="Planarity(z) (#mum)";
  TH1F* hPlanZVsHSright = new TH1F("hPlanZVsHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hPlanZVsHSright->SetLineColor(kMagenta+1);
  hist_title="StaveHSleft vs. Planarity(z)";
  yaxis_name="Planarity(z) (#mum)";
  TH1F* hPlanZVsStaveHSleft = new TH1F("hPlanZVsStaveHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hPlanZVsStaveHSleft->SetLineColor(kGreen+2);
  hist_title="StaveHSright vs. Planarity(z)";
  yaxis_name="Planarity(z) (#mum)";
  TH1F* hPlanZVsStaveHSright = new TH1F("hPlanZVsStaveHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hPlanZVsStaveHSright->SetLineColor(kMagenta+1);
  //************************************************MEAN******************************************************
  hist_title="HSleft vs. Mean(z)";
  yaxis_name="Mean(z) (#mum)";
  TH1F* hMeanZVsHSleft = new TH1F("hMeanZVsHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hMeanZVsHSleft->SetLineColor(kBlack);
  hist_title="HSright vs. Mean(z)";
  yaxis_name="Mean(z) (#mum)";
  TH1F* hMeanZVsHSright = new TH1F("hMeanZVsHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hMeanZVsHSright->SetLineColor(kRed);
  hist_title="StaveHSleft vs. Mean(z)";
  yaxis_name="Mean(z) (#mum)";
  TH1F* hMeanZVsStaveHSleft = new TH1F("hMeanZVsStaveHSleft", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hMeanZVsStaveHSleft->SetLineColor(kBlack);
  hist_title="StaveHSright vs. Mean(z)";
  yaxis_name="Mean(z) (#mum)";
  TH1F* hMeanZVsStaveHSright = new TH1F("hMeanZVsStaveHSright", Form("%s;%s;%s", hist_title.Data(), xaxis_name.Data(), yaxis_name.Data()), nFiles, -.5, nFiles-.5);
  hMeanZVsStaveHSright->SetLineColor(kRed);

  for(int iFile=0; iFile<nFiles; iFile++) {
    //****************************************|||||||||RMS||||||||*****************************************************
    TString xaxis_label(Form("%03d",xaxistitle[iFile]));
    hRMSZVsHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSXVsHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSYVsHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSZVsHSleft->SetBinContent(iFile+1,RMSZ_HSleft[iFile]);
    hRMSXVsHSleft->SetBinContent(iFile+1,hResToNominalHSleft_X[iFile]->GetRMS());
    hRMSYVsHSleft->SetBinContent(iFile+1,hResToNominalHSleft_Y[iFile]->GetRMS());
    hRMSZVsHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSXVsHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSYVsHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSZVsHSright->SetBinContent(iFile+1,RMSZ_HSright[iFile]);
    hRMSXVsHSright->SetBinContent(iFile+1,hResToNominalHSright_X[iFile]->GetRMS());
    hRMSYVsHSright->SetBinContent(iFile+1,hResToNominalHSright_Y[iFile]->GetRMS());
    hRMSZVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSXVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSYVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSZVsStaveHSleft->SetBinContent(iFile+1,RMSZ_StaveHSleft[iFile]);
    hRMSXVsStaveHSleft->SetBinContent(iFile+1,hResToNominalStave_HSleft_X[iFile]->GetRMS());
    hRMSYVsStaveHSleft->SetBinContent(iFile+1,hResToNominalStave_HSleft_Y[iFile]->GetRMS());
    hRMSZVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSXVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSYVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hRMSZVsStaveHSright->SetBinContent(iFile+1,RMSZ_StaveHSright[iFile]);
    hRMSXVsStaveHSright->SetBinContent(iFile+1,hResToNominalStave_HSright_X[iFile]->GetRMS());
    hRMSYVsStaveHSright->SetBinContent(iFile+1,hResToNominalStave_HSright_Y[iFile]->GetRMS());
    //************************************************PLANARITY******************************************************
    hPlanZVsHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hPlanZVsHSleft->SetBinContent(iFile+1,planarityZ_HSleft[iFile]);
    hPlanZVsHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hPlanZVsHSright->SetBinContent(iFile+1,planarityZ_HSright[iFile]);
    hPlanZVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hPlanZVsStaveHSleft->SetBinContent(iFile+1,planarityZ_StaveHSleft[iFile]);
    hPlanZVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hPlanZVsStaveHSright->SetBinContent(iFile+1,planarityZ_StaveHSright[iFile]);
    //************************************************MEAN******************************************************
    hMeanZVsHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hMeanZVsHSleft->SetBinContent(iFile+1,meanZ_HSleft[iFile]);
    hMeanZVsHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hMeanZVsHSright->SetBinContent(iFile+1,meanZ_HSright[iFile]);
    hMeanZVsStaveHSleft->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hMeanZVsStaveHSleft->SetBinContent(iFile+1,meanZ_StaveHSleft[iFile]);
    hMeanZVsStaveHSright->GetXaxis()->SetBinLabel(iFile+1,xaxis_label.Data());
    hMeanZVsStaveHSright->SetBinContent(iFile+1,meanZ_StaveHSright[iFile]);

}

//************************************************MEAN******************************************************
  TLegend* legMeanZHS = new TLegend(.4,.7,0.6,.85);
  legMeanZHS->SetTextSize(0.04);
  legMeanZHS->AddEntry(hMeanZVsHSleft,  "HSL","l");
  legMeanZHS->AddEntry(hMeanZVsHSright, "HSR","l");
  TLegend* legPlanZHS = new TLegend(.4,.7,0.6,.85);
  legPlanZHS->SetTextSize(0.04);
  legPlanZHS->AddEntry(hPlanZVsHSleft,  "HSL","l");
  legPlanZHS->AddEntry(hPlanZVsHSright, "HSR","l");
  TCanvas* cMeanZ_PlanZ_HSLeftRight = new TCanvas("cMeanZ_PlanZ_HSLeftRight","",1150,450);
  cMeanZ_PlanZ_HSLeftRight->Divide(2,1);
  cMeanZ_PlanZ_HSLeftRight->cd(1);
  hMeanZVsHSleft->GetYaxis()->SetRangeUser(-250.,350.);
  hMeanZVsHSleft->Draw();
  hMeanZVsHSright->Draw("same");
  legMeanZHS->Draw("same");
  cMeanZ_PlanZ_HSLeftRight->cd(2);
  hPlanZVsHSleft->GetYaxis()->SetRangeUser(0.,600.);
  hPlanZVsHSleft->Draw();
  hPlanZVsHSright->Draw("same");
  legPlanZHS->Draw("same");
  //************************************************PLANARITY******************************************************
  TLegend* legMeanZStave = new TLegend(.3,.7,.5,.82);
  legMeanZStave->SetTextSize(0.04);
  legMeanZStave->AddEntry(hMeanZVsStaveHSleft, "StaveHSL","l");
  legMeanZStave->AddEntry(hMeanZVsStaveHSright,"StaveHSR","l");
  TLegend* legPlanZStave = new TLegend(.3,.7,.5,.82);
  legPlanZStave->SetTextSize(0.04);
  legPlanZStave->AddEntry(hPlanZVsStaveHSleft, "StaveHSL", "l");
  legPlanZStave->AddEntry(hPlanZVsStaveHSright,"StaveHSR", "l");
  TCanvas* cMeanZ_PlanZ_StaveHSLeftRight = new TCanvas("cMeanZ_PlanZ_StaveHSLeftRight","",1150,450);
  cMeanZ_PlanZ_StaveHSLeftRight->Divide(2,1);
  cMeanZ_PlanZ_StaveHSLeftRight->cd(1);
  hMeanZVsStaveHSleft->GetYaxis()->SetRangeUser(-250.,350.);
  hMeanZVsStaveHSleft->Draw();
  hMeanZVsStaveHSright->Draw("same");
  legMeanZStave->Draw("same");
  cMeanZ_PlanZ_StaveHSLeftRight->cd(2);
  hPlanZVsStaveHSleft->GetYaxis()->SetRangeUser(0.,600.);
  hPlanZVsStaveHSleft->Draw();
  hPlanZVsStaveHSright->Draw("same");
  legPlanZStave->Draw("same");
  //****************************************|||||||||RMS||||||||*****************************************************
  TLegend* legRMSZHS = new TLegend(.65,.7,.83,.83);
  legRMSZHS->SetTextSize(0.04);
  legRMSZHS->AddEntry(hRMSZVsHSleft,  "HSL","l");
  legRMSZHS->AddEntry(hRMSZVsHSright, "HSR","l");
  TLegend* legRMSXHS = new TLegend(.65,.7,.83,.83);
  legRMSXHS->SetTextSize(0.04);
  legRMSXHS->AddEntry(hRMSXVsHSleft,  "HSL","l");
  legRMSXHS->AddEntry(hRMSXVsHSright, "HSR","l");
  TLegend* legRMSYHS = new TLegend(.65,.7,.83,.83);
  legRMSYHS->SetTextSize(0.04);
  legRMSYHS->AddEntry(hRMSYVsHSleft,  "HSL","l");
  legRMSYHS->AddEntry(hRMSYVsHSright, "HSR","l");
  TCanvas* cRMSZXY_HS = new TCanvas("cRMSZXY_HS","",900,900);
  cRMSZXY_HS->Divide(2,2);
  cRMSZXY_HS->cd(1);
  hRMSZVsHSleft->GetYaxis()->SetRangeUser(0.,160.);
  hRMSZVsHSleft->Draw();
  hRMSZVsHSright->Draw("same");
  legRMSZHS->Draw("same");
  cRMSZXY_HS->cd(2);
  hRMSXVsHSleft->GetYaxis()->SetRangeUser(0.,160.);
  hRMSXVsHSleft->Draw();
  hRMSXVsHSright->Draw("same");
  legRMSXHS->Draw("same");
  cRMSZXY_HS->cd(3);
  hRMSYVsHSleft->GetYaxis()->SetRangeUser(0.,35.);
  hRMSYVsHSleft->Draw();
  hRMSYVsHSright->Draw("same");
  legRMSYHS->Draw("same");

  TLegend* legRMSZStave = new TLegend(.62,.72,0.86,.86);
  legRMSZStave->SetTextSize(0.04);
  legRMSZStave->AddEntry(hRMSZVsStaveHSleft , "StaveHSL","l");
  legRMSZStave->AddEntry(hRMSZVsStaveHSright, "StaveHSR","l");
  TLegend* legRMSXStave = new TLegend(.62,.72,0.86,.86);
  legRMSXStave->SetTextSize(0.04);
  legRMSXStave->AddEntry(hRMSXVsStaveHSleft,  "StaveHSL","l");
  legRMSXStave->AddEntry(hRMSXVsStaveHSright, "StaveHSR","l");
  TLegend* legRMSYStave = new TLegend(.62,.72,0.86,.86);
  legRMSYStave->SetTextSize(0.04);
  legRMSYStave->AddEntry(hRMSYVsStaveHSleft,  "StaveHSL","l");
  legRMSYStave->AddEntry(hRMSYVsStaveHSright, "StaveHSR","l");
  TCanvas* cRMSZXY_Stave = new TCanvas("cRMSZXY_Stave","",900,900);
  cRMSZXY_Stave->Divide(2,2);
  cRMSZXY_Stave->cd(1);
  hRMSZVsStaveHSleft->GetYaxis()->SetRangeUser(0.,160.);
  hRMSZVsStaveHSleft->Draw();
  hRMSZVsStaveHSright->Draw("same");
  legRMSZStave->Draw("same");
  cRMSZXY_Stave->cd(2);
  hRMSXVsStaveHSleft->GetYaxis()->SetRangeUser(0.,160.);
  hRMSXVsStaveHSleft->Draw();
  hRMSXVsStaveHSright->Draw("same");
  legRMSXStave->Draw("same");
  cRMSZXY_Stave->cd(3);
  hRMSYVsStaveHSleft->GetYaxis()->SetRangeUser(0.,35.);
  hRMSYVsStaveHSleft->Draw();
  hRMSYVsStaveHSright->Draw("same");
  legRMSYStave->Draw("same");

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
