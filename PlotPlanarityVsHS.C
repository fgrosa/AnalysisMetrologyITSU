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

#endif

const TString infilenames[] = {"/home/cecilia/Desktop/STAGE/StavesMetrology/Stave5/Stave/STAVE_MARKERPOS_2018_6_6_T-OL-Stave-5_T-OL-HS-U-005_ALC-0312-01_109_T-OL-HS-L-005_ALC-0312-01_248_MeasAndExtrap.root"};
const int xaxistitle[] = {5};

int PlotPlanarityVsHS();
void ComputePlanarityRMSandMean(TGraph* graph, double &planarity, double &zRMS, double &zmean);

int PlotPlanarityVsHS() {

  gStyle->SetOptStat(0);

  const int nFiles = sizeof(infilenames)/sizeof(infilenames[0]);
  TGraph* gResToNomPlaneHSleft_Z[nFiles];
  TH1F* hResToNominalHSleft_X[nFiles];
  double planarityZ_HSleft[nFiles];
  double RMSZ_HSleft[nFiles];
  double meanZ_HSleft[nFiles];

  for(int iFile=0; iFile<nFiles; iFile++) {
    TFile* infile = TFile::Open(infilenames[iFile].Data());
    if(!infile) return -1;
    gResToNomPlaneHSleft_Z[iFile] = (TGraph*)infile->Get("gResToNomPlaneHSleft_Z");
    if(!gResToNomPlaneHSleft_Z[iFile]) return -2;
    hResToNominalHSleft_X[iFile] = (TH1F*)infile->Get("hResToNominalHSleft_X");
    if(!hResToNominalHSleft_X[iFile]) return -3;
    hResToNominalHSleft_X[iFile]->SetDirectory(0);
    ComputePlanarityRMSandMean(gResToNomPlaneHSleft_Z[iFile],planarityZ_HSleft[iFile],RMSZ_HSleft[iFile],meanZ_HSleft[iFile]);

    infile->Close();
  }

  TH1F* hRMSZVsHSleft = new TH1F("hRMSZVsHSleft","RMS(z) vs. HSleft;;RMS (#mum)",nFiles,-0.5,nFiles-0.5);
  TH1F* hRMSXVsHSleft = new TH1F("hRMSXVsHSleft","RMS(x) vs. HSleft;;RMS (#mum)",nFiles,-0.5,nFiles-0.5);
  for(int iFile=0; iFile<nFiles; iFile++) {
    hRMSZVsHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("HSL%d",xaxistitle[iFile]));
    hRMSXVsHSleft->GetXaxis()->SetBinLabel(iFile+1,Form("HSL%d",xaxistitle[iFile]));
    hRMSZVsHSleft->SetBinContent(iFile+1,RMSZ_HSleft[iFile]);
    hRMSXVsHSleft->SetBinContent(iFile+1,hResToNominalHSleft_X[iFile]->GetRMS());
  }

  TCanvas* cRMSZVsHSleft = new TCanvas("cRMSZVsHSleft","",800,800);
  hRMSZVsHSleft->Draw();

  TCanvas* cRMSXVsHSleft = new TCanvas("cRMSXVsHSleft","",800,800);
  hRMSXVsHSleft->Draw();

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
