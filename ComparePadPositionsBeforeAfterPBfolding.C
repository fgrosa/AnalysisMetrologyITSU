#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TPad.h>
#include <TString.h>
#include <TLine.h>
#include <TLatex.h>
#include <TBox.h>

#endif

//***************************************************************************************//
//                                                                                       //
//    Main Functions: ComparePadPositionsBeforeAfterPBfolding(), DrawMeanAndRMSvsHS()    //
//                                                                                       //
//***************************************************************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIBALES
const bool fMitutoyoFile=false;

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
int ComparePadPositionsBeforeAfterPBfolding(TString infileAfter = "PadPos/BeforeAfterPB/Stave7/PADPOS_after_PB_T-OL-Stave-007.txt",
                                            TString infileBefore = "PadPos/BeforeAfterPB/Stave7/PADPOS_before_PB_T-OL-Stave-007.txt",
                                            TString outfilename = "PadPos/BeforeAfterPB/Stave7/PADPOS_before_after_PB_T-OL-Stave-007.root");
int DrawMeanAndRMSvsHS();
bool ReadFile(const TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset=0.);
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset=0.);
bool MatchMeas(double xmeas1, double ymeas1, double zmeas1, double xmeas2, double ymeas2, double zmeas2, double accradius);
void SetStyle();

//_____________________________________________________________________________________________
int ComparePadPositionsBeforeAfterPBfolding(TString infileAfter, TString infileBefore, TString outfilename) {

  SetStyle();

  std::vector<double> x_1;
  std::vector<double> y_1;
  std::vector<double> z_1;
  std::vector<double> x_2;
  std::vector<double> y_2;
  std::vector<double> z_2;

  bool readMeas1=false;
  bool readMeas2=false;
  if(!fMitutoyoFile) {
    readMeas1=ReadFile(infileAfter,x_1,y_1,z_1);
    readMeas2=ReadFile(infileBefore,x_2,y_2,z_2);
  }
  else {
    readMeas1=ReadDatFileMitutoyo(infileAfter,x_1,y_1,z_1);
    readMeas2=ReadDatFileMitutoyo(infileBefore,x_2,y_2,z_2);
  }
  if(!readMeas1) {
    cerr << "Impossibile to find the input fileAfter for measured positions. Exit" << endl;
    return 1;
  }
  if(!readMeas2) {
    cerr << "Impossibile to find the input fileBefore for measured positions. Exit" << endl;
    return 2;
  }

  std::vector<double> resx_HSR;
  std::vector<double> resx_HSL;
  std::vector<double> resy_HSR;
  std::vector<double> resy_HSL;
  std::vector<double> resz_HSR;
  std::vector<double> resz_HSL;
  std::vector<double> x1_HSR;
  std::vector<double> x1_HSL;
  std::vector<double> y1_HSR;
  std::vector<double> y1_HSL;
  std::vector<double> z1_HSR;
  std::vector<double> z1_HSL;
  std::vector<double> x2_HSR;
  std::vector<double> x2_HSL;
  std::vector<double> y2_HSR;
  std::vector<double> y2_HSL;
  std::vector<double> z2_HSR;
  std::vector<double> z2_HSL;
  std::vector<double> ycommon_HSL;
  std::vector<double> ycommon_HSR;
  TH1F* hResX_HSR = new TH1F("hResX_HSR","HSR;#DeltaX = X_{(after PB)}-X_{(before PB)} (#mum);Counts",100,-100,100);
  TH1F* hResY_HSR = new TH1F("hResY_HSR","HSR;#DeltaY = Y_{(after PB)}-Y_{(before PB)} (#mum);Counts",100,-100,100);
  TH1F* hResZ_HSR = new TH1F("hResZ_HSR","HSR;#DeltaZ = Z_{(after PB)}-Z_{(before PB)} (#mum);Counts",100,-200,200);
  TH1F* hResX_HSL = new TH1F("hResX_HSL","HSL;#DeltaX = X_{(after PB)}-X_{(before PB)} (#mum);Counts",100,-100,100);
  TH1F* hResY_HSL = new TH1F("hResY_HSL","HSL;#DeltaY = Y_{(after PB)}-Y_{(before PB)} (#mum);Counts",100,-100,100);
  TH1F* hResZ_HSL = new TH1F("hResZ_HSL","HSL;#DeltaZ = Z_{(after PB)}-Z_{(before PB)} (#mum);Counts",100,-200,200);

  //match measurements
  int meassinglepoint=0;
  for(unsigned int iEntry1=0; iEntry1<x_1.size(); iEntry1++) {
    for(unsigned int iEntry2=0; iEntry2<x_2.size(); iEntry2++) {
      if(y_2[iEntry2]<-735 && y_2[iEntry2]>-740 && z_2[iEntry2]<10 && meassinglepoint==0) {
        x2_HSR.push_back(x_2[iEntry2]);
        y2_HSR.push_back(y_2[iEntry2]);
        z2_HSR.push_back(z_2[iEntry2]);
        meassinglepoint++;
      }

      bool match = MatchMeas(x_1[iEntry1],y_1[iEntry1],z_1[iEntry1],x_2[iEntry2],y_2[iEntry2],z_2[iEntry2],0.5);
      if(match) {
        if(z_2[iEntry2]>12) {
          resx_HSL.push_back(x_1[iEntry1]-x_2[iEntry2]);
          resy_HSL.push_back(y_1[iEntry1]-y_2[iEntry2]);
          resz_HSL.push_back(z_1[iEntry1]-z_2[iEntry2]);
          x1_HSL.push_back(x_1[iEntry1]);
          y1_HSL.push_back(y_1[iEntry1]);
          z1_HSL.push_back(z_1[iEntry1]);
          x2_HSL.push_back(x_2[iEntry2]);
          y2_HSL.push_back(y_2[iEntry2]);
          z2_HSL.push_back(z_2[iEntry2]);
          ycommon_HSL.push_back((y_1[iEntry1]+y_2[iEntry2])/2);
          hResX_HSL->Fill((x_1[iEntry1]-x_2[iEntry2])*1000);
          hResY_HSL->Fill((y_1[iEntry1]-y_2[iEntry2])*1000);
          hResZ_HSL->Fill((z_1[iEntry1]-z_2[iEntry2])*1000);
        }
        else{
          if(y_1[iEntry1]<-735 && y_1[iEntry1]>-740) continue;
          resx_HSR.push_back(x_1[iEntry1]-x_2[iEntry2]);
          resy_HSR.push_back(y_1[iEntry1]-y_2[iEntry2]);
          resz_HSR.push_back(z_1[iEntry1]-z_2[iEntry2]);
          x1_HSR.push_back(x_1[iEntry1]);
          y1_HSR.push_back(y_1[iEntry1]);
          z1_HSR.push_back(z_1[iEntry1]);
          x2_HSR.push_back(x_2[iEntry2]);
          y2_HSR.push_back(y_2[iEntry2]);
          z2_HSR.push_back(z_2[iEntry2]);
          ycommon_HSR.push_back((y_1[iEntry1]+y_2[iEntry2])/2);
          hResX_HSR->Fill((x_1[iEntry1]-x_2[iEntry2])*1000);
          hResY_HSR->Fill((y_1[iEntry1]-y_2[iEntry2])*1000);
          hResZ_HSR->Fill((z_1[iEntry1]-z_2[iEntry2])*1000);
        }
      }
    }
  }

  TGraph *gZ1_HSL = new TGraph(0);
  gZ1_HSL->SetMarkerStyle(kOpenCircle);
  gZ1_HSL->SetMarkerSize(1.5);
  gZ1_HSL->SetMarkerColor(kRed);
  gZ1_HSL->SetName("gZafter_HSL");
  gZ1_HSL->SetTitle("HSL (after PB)");

  TGraph *gZ1_HSR = new TGraph(0);
  gZ1_HSR->SetMarkerStyle(kFullCircle);
  gZ1_HSR->SetMarkerSize(1.5);
  gZ1_HSR->SetMarkerColor(kRed);
  gZ1_HSR->SetName("gZafter_HSR");
  gZ1_HSR->SetTitle("HSR (after PB)");

  TGraph *gZ2_HSL = new TGraph(0);
  gZ2_HSL->SetMarkerStyle(kOpenCircle);
  gZ2_HSL->SetMarkerSize(1.5);
  gZ2_HSL->SetMarkerColor(kBlack);
  gZ1_HSL->SetName("gZbefore_HSL");
  gZ1_HSL->SetTitle("HSL (before PB)");

  TGraph *gZ2_HSR = new TGraph(0);
  gZ2_HSR->SetMarkerStyle(kFullCircle);
  gZ2_HSR->SetMarkerSize(1.5);
  gZ2_HSR->SetMarkerColor(kBlack);
  gZ2_HSR->SetName("gZbefore_HSR");
  gZ2_HSR->SetTitle("HSR (before PB)");

  TGraph *gDeltaZ_HSL = new TGraph(0);
  gDeltaZ_HSL->SetMarkerStyle(kFullCircle);
  gDeltaZ_HSL->SetMarkerColor(kBlack);
  gDeltaZ_HSL->SetName("gDeltaZ_HSL");
  gDeltaZ_HSL->SetTitle("HSL (after-before PB)");

  TGraph *gDeltaZ_HSR = new TGraph(0);
  gDeltaZ_HSR->SetMarkerStyle(kFullCircle);
  gDeltaZ_HSR->SetMarkerColor(kBlack);
  gDeltaZ_HSR->SetName("gDeltaZ_HSR");
  gDeltaZ_HSR->SetTitle("HSR (after-before PB)");

  for(unsigned int iEntry1=0; iEntry1<x1_HSL.size(); iEntry1++) {
    gZ1_HSL->SetPoint(iEntry1,y1_HSL[iEntry1],z1_HSL[iEntry1]);
  }
  for(unsigned int iEntry1=0; iEntry1<x1_HSR.size(); iEntry1++) {
    gZ1_HSR->SetPoint(iEntry1,y1_HSR[iEntry1],z1_HSR[iEntry1]);
  }

  for(unsigned int iEntry2=0; iEntry2<x2_HSL.size(); iEntry2++) {
    gZ2_HSL->SetPoint(iEntry2,y2_HSL[iEntry2],z2_HSL[iEntry2]);
  }
  for(unsigned int iEntry2=0; iEntry2<x2_HSR.size(); iEntry2++) {
    gZ2_HSR->SetPoint(iEntry2,y2_HSR[iEntry2],z2_HSR[iEntry2]);
  }

  for(unsigned int iEntry=0; iEntry<resz_HSL.size(); iEntry++) {
    gDeltaZ_HSL->SetPoint(iEntry,ycommon_HSL[iEntry],resz_HSL[iEntry]*1000);
  }
  for(unsigned int iEntry=0; iEntry<resz_HSR.size(); iEntry++) {
    gDeltaZ_HSR->SetPoint(iEntry,ycommon_HSR[iEntry],resz_HSR[iEntry]*1000);
  }

  TCanvas* cRes_HSL = new TCanvas("cRes_HSL","",1200,600);
  cRes_HSL->Divide(3,1);
  cRes_HSL->cd(1);
  hResX_HSL->SetLineWidth(2);
  hResX_HSL->Draw();
  cRes_HSL->cd(2);
  hResY_HSL->SetLineWidth(2);
  hResY_HSL->Draw();
  cRes_HSL->cd(3);
  hResZ_HSL->SetLineWidth(2);
  hResZ_HSL->Draw();

  TCanvas* cRes_HSR = new TCanvas("cRes_HSR","",1200,600);
  cRes_HSR->Divide(3,1);
  cRes_HSR->cd(1);
  hResX_HSR->SetLineWidth(2);
  hResX_HSR->Draw();
  cRes_HSR->cd(2);
  hResY_HSR->SetLineWidth(2);
  hResY_HSR->Draw();
  cRes_HSR->cd(3);
  hResZ_HSR->SetLineWidth(2);
  hResZ_HSR->Draw();

  TF1* fZ1_HSL = new TF1("fZ1_HSL","pol2",-999,999);
  TF1* fZ1_HSR = new TF1("fZ1_HSR","pol2",-999,999);
  TF1* fZ2_HSL = new TF1("fZ2_HSL","pol2",-999,999);
  TF1* fZ2_HSR = new TF1("fZ2_HSR","pol2",-999,999);
  fZ2_HSL->SetLineColor(kBlack);
  fZ2_HSR->SetLineColor(kBlack);

  TLegend* leg = new TLegend(0.4,0.325,0.7,0.575);
  leg->SetTextSize(0.05);
  leg->AddEntry(fZ1_HSL,"Before PB folding","l");
  leg->AddEntry(fZ2_HSL,"After PB folding","l");
  leg->AddEntry(gZ2_HSL,"HSL","p");
  leg->AddEntry(gZ2_HSR,"HSR","p");

  TCanvas* cZVsY = new TCanvas("cZVsY","",1920,1080);
  cZVsY->cd(1)->DrawFrame(-999,9,999,15,";Y (mm);Z (mm)");
  gZ2_HSL->Draw("P");
  gZ2_HSL->Fit("fZ2_HSL");
  gZ2_HSR->Draw("P");
  gZ2_HSR->Fit("fZ2_HSR");
  gZ1_HSL->Draw("P");
  gZ1_HSL->Fit("fZ1_HSL");
  gZ1_HSR->Draw("P");
  gZ1_HSR->Fit("fZ1_HSR");
  leg->Draw("same");

  TCanvas* cDeltaZVsY = new TCanvas("cDeltaZVsY","",1920,1080);
  cDeltaZVsY->Divide(2,1);
  cDeltaZVsY->cd(1)->DrawFrame(-999,-500,999,500,"HSL;Y (mm);#DeltaZ = Z_{(after PB)}-Z_{(before PB)} (#mum)");
  gDeltaZ_HSL->Draw("P");
  cDeltaZVsY->cd(2)->DrawFrame(-999,-500,999,500,"HSR;Y (mm);#DeltaZ = Z_{(after PB)}-Z_{(before PB)} (#mum)");
  gDeltaZ_HSR->Draw("P");

  TFile outfile(outfilename.Data(),"RECREATE");
  hResX_HSL->Write();
  hResY_HSL->Write();
  hResZ_HSL->Write();
  hResX_HSR->Write();
  hResY_HSR->Write();
  hResZ_HSR->Write();
  gZ1_HSL->Write();
  gZ2_HSL->Write();
  gZ1_HSR->Write();
  gZ2_HSR->Write();
  gDeltaZ_HSR->Write();
  gDeltaZ_HSL->Write();
  cDeltaZVsY->Write();
  cZVsY->Write();
  cRes_HSL->Write();
  cRes_HSR->Write();
  outfile.Close();
  
  outfilename.ReplaceAll(".root",".pdf");
  cZVsY->Print(Form("%s[",outfilename.Data()));
  cZVsY->Print(outfilename.Data());
  cDeltaZVsY->Print(outfilename.Data());
  cRes_HSL->Print(outfilename.Data());
  cRes_HSR->Print(outfilename.Data());
  cRes_HSR->Print(Form("%s]",outfilename.Data()));
  
  return 0;
}

//_____________________________________________________________________________________________
int DrawMeanAndRMSvsHS() {
  
  SetStyle();
  
  TString infilenames[] = {"PadPos/BeforeAfterPB/Stave3/PADPOS_before_after_PB_T-OL-Stave-003.root","PadPos/BeforeAfterPB/Stave2/PADPOS_before_after_PB_T-OL-Stave-002.root","PadPos/BeforeAfterPB/Stave4/PADPOS_before_after_PB_T-OL-Stave-004.root","PadPos/BeforeAfterPB/Stave5/PADPOS_before_after_PB_T-OL-Stave-005.root","PadPos/BeforeAfterPB/Stave6/PADPOS_before_after_PB_T-OL-Stave-006.root","PadPos/BeforeAfterPB/Stave7/PADPOS_before_after_PB_T-OL-Stave-007.root","PadPos/BeforeAfterPB/Stave8/PADPOS_before_after_PB_T-OL-Stave-008.root","PadPos/BeforeAfterPB/Stave9/PADPOS_before_after_PB_T-OL-Stave-009.root"};

  TString HStitles[] = {"HS2","HS3","HS4","HS5","HS6","HS7","HS8","HS9"};

  int nFiles = sizeof(infilenames)/sizeof(infilenames[0]);
  TH1F* hResX_HSL[nFiles];
  TH1F* hResY_HSL[nFiles];
  TH1F* hResZ_HSL[nFiles];
  TH1F* hResX_HSR[nFiles];
  TH1F* hResY_HSR[nFiles];
  TH1F* hResZ_HSR[nFiles];
  for(int iFile=0; iFile<nFiles; iFile++) {
    TFile* infile = TFile::Open(infilenames[iFile].Data());
    if(!infile) {cerr << "Error: input file "<<iFile << " not found! Check the path." << endl; return iFile;}
    hResX_HSL[iFile] = (TH1F*)infile->Get("hResX_HSL");
    hResY_HSL[iFile] = (TH1F*)infile->Get("hResY_HSL");
    hResZ_HSL[iFile] = (TH1F*)infile->Get("hResZ_HSL");
    hResX_HSR[iFile] = (TH1F*)infile->Get("hResX_HSR");
    hResY_HSR[iFile] = (TH1F*)infile->Get("hResY_HSR");
    hResZ_HSR[iFile] = (TH1F*)infile->Get("hResZ_HSR");
    hResX_HSL[iFile]->SetDirectory(0);
    hResY_HSL[iFile]->SetDirectory(0);
    hResZ_HSL[iFile]->SetDirectory(0);
    hResX_HSR[iFile]->SetDirectory(0);
    hResY_HSR[iFile]->SetDirectory(0);
    hResZ_HSR[iFile]->SetDirectory(0);
    infile->Close();
  }
  
  TH1F* hRMSX_HSR_VsHS = new TH1F("hResX_HSR_VsHS",";;RMS(#DeltaX) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hRMSY_HSR_VsHS = new TH1F("hResY_HSR_VsHS",";;RMS(#DeltaY) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hRMSZ_HSR_VsHS = new TH1F("hResZ_HSR_VsHS",";;RMS(#DeltaZ) (#mum)",nFiles,0.5,nFiles+0.5);
  hRMSX_HSR_VsHS->SetStats(0);
  hRMSY_HSR_VsHS->SetStats(0);
  hRMSZ_HSR_VsHS->SetStats(0);

  TH1F* hRMSX_HSL_VsHS = new TH1F("hResX_HSL_VsHS",";;RMS(#DeltaX) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hRMSY_HSL_VsHS = new TH1F("hResY_HSL_VsHS",";;RMS(#DeltaY) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hRMSZ_HSL_VsHS = new TH1F("hResZ_HSL_VsHS",";;RMS(#DeltaZ) (#mum)",nFiles,0.5,nFiles+0.5);
  hRMSX_HSL_VsHS->SetLineColor(kRed);
  hRMSY_HSL_VsHS->SetLineColor(kRed);
  hRMSZ_HSL_VsHS->SetLineColor(kRed);
  hRMSX_HSL_VsHS->SetStats(0);
  hRMSY_HSL_VsHS->SetStats(0);
  hRMSZ_HSL_VsHS->SetStats(0);

  TH1F* hMeanX_HSR_VsHS = new TH1F("hMeanX_HSR_VsHS",";;Mean(#DeltaX) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hMeanY_HSR_VsHS = new TH1F("hMeanY_HSR_VsHS",";;Mean(#DeltaY) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hMeanZ_HSR_VsHS = new TH1F("hMeanZ_HSR_VsHS",";;Mean(#DeltaZ) (#mum)",nFiles,0.5,nFiles+0.5);
  hMeanX_HSR_VsHS->SetStats(0);
  hMeanY_HSR_VsHS->SetStats(0);
  hMeanZ_HSR_VsHS->SetStats(0);

  TH1F* hMeanX_HSL_VsHS = new TH1F("hMeanX_HSL_VsHS",";;Mean(#DeltaX) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hMeanY_HSL_VsHS = new TH1F("hMeanY_HSL_VsHS",";;Mean(#DeltaY) (#mum)",nFiles,0.5,nFiles+0.5);
  TH1F* hMeanZ_HSL_VsHS = new TH1F("hMeanZ_HSL_VsHS",";;Mean(#DeltaZ) (#mum)",nFiles,0.5,nFiles+0.5);
  hMeanX_HSL_VsHS->SetLineColor(kRed);
  hMeanY_HSL_VsHS->SetLineColor(kRed);
  hMeanZ_HSL_VsHS->SetLineColor(kRed);
  hMeanX_HSL_VsHS->SetStats(0);
  hMeanY_HSL_VsHS->SetStats(0);
  hMeanZ_HSL_VsHS->SetStats(0);

  for(int iFile=0; iFile<nFiles; iFile++) {
    hRMSX_HSR_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hRMSY_HSR_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hRMSZ_HSR_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hRMSX_HSL_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hRMSY_HSL_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hRMSZ_HSL_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hMeanX_HSR_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hMeanY_HSR_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hMeanZ_HSR_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hMeanX_HSL_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hMeanY_HSL_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    hMeanZ_HSL_VsHS->GetXaxis()->SetBinLabel(iFile+1,HStitles[iFile].Data());
    
    hRMSX_HSR_VsHS->SetBinContent(iFile+1,hResX_HSR[iFile]->GetRMS());
    hRMSY_HSR_VsHS->SetBinContent(iFile+1,hResY_HSR[iFile]->GetRMS());
    hRMSZ_HSR_VsHS->SetBinContent(iFile+1,hResZ_HSR[iFile]->GetRMS());
    hRMSX_HSL_VsHS->SetBinContent(iFile+1,hResX_HSL[iFile]->GetRMS());
    hRMSY_HSL_VsHS->SetBinContent(iFile+1,hResY_HSL[iFile]->GetRMS());
    hRMSZ_HSL_VsHS->SetBinContent(iFile+1,hResZ_HSL[iFile]->GetRMS());
    hMeanX_HSR_VsHS->SetBinContent(iFile+1,hResX_HSR[iFile]->GetMean());
    hMeanY_HSR_VsHS->SetBinContent(iFile+1,hResY_HSR[iFile]->GetMean());
    hMeanZ_HSR_VsHS->SetBinContent(iFile+1,hResZ_HSR[iFile]->GetMean());
    hMeanX_HSL_VsHS->SetBinContent(iFile+1,hResX_HSL[iFile]->GetMean());
    hMeanY_HSL_VsHS->SetBinContent(iFile+1,hResY_HSL[iFile]->GetMean());
    hMeanZ_HSL_VsHS->SetBinContent(iFile+1,hResZ_HSL[iFile]->GetMean());
  }
  
  TLegend* leg = new TLegend(0.55,0.75,0.85,0.85);
  leg->SetTextSize(0.05);
  leg->AddEntry(hRMSX_HSL_VsHS,"HS Left","l");
  leg->AddEntry(hRMSX_HSR_VsHS,"HS Right","l");

  TCanvas* cRMS = new TCanvas("cRMS","",1200,500);
  cRMS->Divide(3,1);
  cRMS->cd(1);
  hRMSX_HSL_VsHS->GetYaxis()->SetRangeUser(0,50);
  hRMSX_HSL_VsHS->Draw();
  hRMSX_HSR_VsHS->Draw("same");
  leg->Draw("same");
  cRMS->cd(2);
  hRMSY_HSL_VsHS->GetYaxis()->SetRangeUser(0,50);
  hRMSY_HSL_VsHS->Draw();
  hRMSY_HSR_VsHS->Draw("same");
  leg->Draw("same");
  cRMS->cd(3);
  hRMSZ_HSL_VsHS->GetYaxis()->SetRangeUser(0,50);
  hRMSZ_HSL_VsHS->Draw();
  hRMSZ_HSR_VsHS->Draw("same");
  leg->Draw("same");
  
  TCanvas* cMean = new TCanvas("cMean","",1200,500);
  cMean->Divide(3,1);
  cMean->cd(1);
  hMeanX_HSL_VsHS->GetYaxis()->SetRangeUser(-150,150);
  hMeanX_HSL_VsHS->Draw();
  hMeanX_HSR_VsHS->Draw("same");
  leg->Draw("same");
  cMean->cd(2);
  hMeanY_HSL_VsHS->GetYaxis()->SetRangeUser(-150,150);
  hMeanY_HSL_VsHS->Draw();
  hMeanY_HSR_VsHS->Draw("same");
  leg->Draw("same");
  cMean->cd(3);
  hMeanZ_HSL_VsHS->GetYaxis()->SetRangeUser(-150,150);
  hMeanZ_HSL_VsHS->Draw();
  hMeanZ_HSR_VsHS->Draw("same");
  leg->Draw("same");

  cMean->SaveAs("PadPos/BeforeAfterPB/Mean_PADPOS_before_after_PB_allstaves.pdf");
  cRMS->SaveAs("PadPos/BeforeAfterPB/RMS_PADPOS_before_after_PB_allstaves.pdf");

  return 0;
}

//_____________________________________________________________________________________________
bool ReadFile(const TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset)
{
  ifstream inSet(FileName.Data());
  if(!inSet) {
    cerr<<"File "<<FileName.Data() <<" does not exists. "<<endl;
    return false;
  }

  double xcoord;
  double ycoord;
  double zcoord;
  while(inSet>>xcoord>>ycoord>>zcoord) {
    x.push_back(xcoord+xoffset);
    y.push_back(ycoord);
    z.push_back(zcoord);
  }

  inSet.close();

  return true;
}

//_____________________________________________________________________________________________
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset)
{
  if(!FileName.Contains(".txt") && !FileName.Contains(".dat") && !FileName.Contains(".csv")) {
    cerr << "Wrong file format. Exit." << endl;
    return false;
  }

  string valueSeparator = ";";
  ifstream inSet(FileName.Data());
  if(!inSet) {
    cerr<<"Please check if "<<FileName.Data() <<" is the right path. Exit."<<endl;
    return false;
  }

  bool islineok=true;
  int linecounter=0;
  if (inSet.is_open()) {
    while(!inSet.eof())
    {
      islineok=true;
      std::vector<string> values;
      string line;
      getline(inSet, line);
      TString test = line;
      if(test.Contains("Mean") || linecounter>0) {
        islineok=false;
        linecounter++;
      }
      if(linecounter==5) {
        linecounter=0;
        islineok=true;
      }
      if(test.Contains("MarkerCenter")) {islineok=false;}

      if(islineok && linecounter<5) {
        size_t pos = 0;
        while((pos = line.find(valueSeparator)) != string::npos)
        {
          double xtmp;
          double ytmp;
          double ztmp;
          values.push_back(line.substr(0,pos));
          if(values.size()==5) {
            stringstream convert(values[4]);
            if( !(convert >> xtmp) ) xtmp = -1;
            x.push_back(xtmp+xoffset);
          }
          else if(values.size()==6) {
            stringstream convert(values[5]);
            if( !(convert >> ytmp) ) ytmp = -1;
            y.push_back(ytmp);
          }
          else if(values.size()==7) {
            stringstream convert(values[6]);
            if( !(convert >> ztmp) ) ztmp = -1;
            z.push_back(ztmp);
          }
          line = line.substr(pos + 1);
        }
      }
    }
    inSet.close();
  }

  return true;
}

//______________________________________________________________________________________________
bool MatchMeas(double xmeas1, double ymeas1, double zmeas1, double xmeas2, double ymeas2, double zmeas2, double accradius) {

  double radius = TMath::Sqrt((xmeas1-xmeas2)*(xmeas1-xmeas2)+(ymeas1-ymeas2)*(ymeas1-ymeas2)+(zmeas1-zmeas2)*(zmeas1-zmeas2));
  if(radius<accradius) return true;

  return false;
}

//______________________________________________________________________________________________
void SetStyle() {
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.035);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleFont(42);
  gStyle->SetLabelFont(42);
  gStyle->SetTitleSize(0.05,"t");
  gStyle->SetLabelSize(0.045,"xyzt");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetNdivisions(505);
  gStyle->SetOptStat(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadGridX(1);

}
