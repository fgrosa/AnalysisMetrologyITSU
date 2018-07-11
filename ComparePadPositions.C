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

//*********************************************//
//                                             //
//    Main Functions: ComparePadPositions()    //
//                                             //
//*********************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIBALES
const TString infile1 = "STAVE4/STAVE_PADPOS_2018_5_7_T-OL-Stave-4_T-OL-HS-L-4_ALC-0312-01_246_T-OL-HS-R-4_ALC-0312-01_210_RS212.dat";
const TString infile2 = "STAVE4/STAVE_PADPOS_2018_5_7_T-OL-Stave-4_T-OL-HS-L-4_ALC-0312-01_246_T-OL-HS-R-4_ALC-0312-01_210_RS222.dat";
bool fMitutoyoFile=true;
bool fChangeRS1=true;

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
int ComparePadPositions();
bool ReadFile(const TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset=0.);
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset=0.);
bool MatchMeas(double xmeas1, double ymeas1, double xmeas2, double ymeas2, double accradius);
void SetStyle();

//_____________________________________________________________________________________________
int ComparePadPositions() {

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
    readMeas1=ReadFile(infile1,x_1,y_1,z_1);
    readMeas2=ReadFile(infile2,x_2,y_2,z_2);
  }
  else {
    readMeas1=ReadDatFileMitutoyo(infile1,x_1,y_1,z_1);//read from 212
    readMeas2=ReadDatFileMitutoyo(infile2,x_2,y_2,z_2);//read from 222
  }
  if(!readMeas1) {
    cerr << "Impossibile to find the input file1 for measured positions. Exit" << endl;
    return 1;
  }
  if(!readMeas2) {
    cerr << "Impossibile to find the input file2 for measured positions. Exit" << endl;
    return 2;
  }

  //transform RS212 into RS222
  std::vector<double> x_1_tr;
  std::vector<double> z_1_tr;
  for(unsigned int iEntry=0; iEntry<x_1.size(); iEntry++) {
    x_1_tr.push_back(-x_1[iEntry]);
    z_1_tr.push_back(-z_1[iEntry]);
  }

  std::vector<double> resx;
  std::vector<double> resy;
  TH1F* hResX = new TH1F("hResX",";#DeltaX (#mum);Counts",100,-100,100);
  TH1F* hResY = new TH1F("hResY",";#DeltaY (#mum);Counts",100,-100,100);
  TH1F* hResZ = new TH1F("hResZ",";#DeltaZ (#mum);Counts",100,-200,200);

  //match measurements
  for(unsigned int iEntry1=0; iEntry1<x_1.size(); iEntry1++) {
    for(unsigned int iEntry2=0; iEntry2<x_2.size(); iEntry2++) {
      bool match = false;
      if(fChangeRS1) match = MatchMeas(x_1_tr[iEntry1],y_1[iEntry1],x_2[iEntry2],y_2[iEntry2],1.);
      else match = MatchMeas(x_1[iEntry1],y_1[iEntry1],x_2[iEntry2],y_2[iEntry2],1.);
      if(match) {
        if(fChangeRS1) {
          resx.push_back(x_1_tr[iEntry1]-x_2[iEntry2]);
          hResX->Fill((x_1_tr[iEntry1]-x_2[iEntry2])*1000);
            hResZ->Fill((z_1_tr[iEntry1]-z_2[iEntry2])*1000);
        }
        else {
          resx.push_back(x_1[iEntry1]-x_2[iEntry2]);
          hResX->Fill((x_1[iEntry1]-x_2[iEntry2])*1000);
            hResZ->Fill((z_1[iEntry1]-z_2[iEntry2])*1000);
        }
        resy.push_back(y_1[iEntry1]-y_2[iEntry2]);
        hResY->Fill((y_1[iEntry1]-y_2[iEntry2])*1000);

      }
    }
  }

  TString title1 = "RS222";
  TString title2 = "RS212";
  TString title3 = "RS212 reversed";
  TString title4 = "RS222 + R212 reversed";
  if(!fChangeRS1) {
    title2 = "RS222 (2)";
    title3 = "RS222 (2)";
    title3 = "RS222 comparison";
  }

  TGraph *gZ1_Up = new TGraph(0);
  gZ1_Up->SetTitle(title2.Data());
  gZ1_Up->SetMarkerStyle(kFullCircle);
  gZ1_Up->SetMarkerColor(kRed);

  TGraph *gZ1_Down = new TGraph(0);
  gZ1_Down->SetTitle(title2.Data());
  gZ1_Down->SetMarkerStyle(kFullCircle);
  gZ1_Down->SetMarkerColor(kRed);

  TGraph *gZ1rev_Up = new TGraph(0);
  gZ1rev_Up->SetTitle(title3.Data());
  gZ1rev_Up->SetMarkerStyle(kFullCircle);
  gZ1rev_Up->SetMarkerColor(kRed);

  TGraph *gZ1rev_Down = new TGraph(0);
  gZ1rev_Down->SetTitle(title3.Data());
  gZ1rev_Down->SetMarkerStyle(kFullCircle);
  gZ1rev_Down->SetMarkerColor(kRed);

  TGraph *gZ2_Up = new TGraph(0);
  gZ2_Up->SetTitle(title1.Data());
  gZ2_Up->SetMarkerStyle(kFullCircle);
  gZ2_Up->SetMarkerColor(kBlack);

  TGraph *gZ2_Down = new TGraph(0);
  gZ2_Down->SetTitle(title1.Data());
  gZ2_Down->SetMarkerStyle(kFullCircle);
  gZ2_Down->SetMarkerColor(kBlack);

  int point_up = 0;
  int point_down = 0;
  for(unsigned int iEntry1=0; iEntry1<x_1.size(); iEntry1++) {
    if(z_1[iEntry1]>-10.5 && fChangeRS1) {
      gZ1_Up->SetPoint(point_up,y_1[iEntry1],z_1[iEntry1]);
      gZ1rev_Up->SetPoint(point_up,y_1[iEntry1],z_1_tr[iEntry1]);
      point_up++;
    }
    else if(z_1[iEntry1]<-10.5 && fChangeRS1){
      gZ1_Down->SetPoint(point_down,y_1[iEntry1],z_1[iEntry1]);
      gZ1rev_Down->SetPoint(point_down,y_1[iEntry1],z_1_tr[iEntry1]);
      point_down++;
    }
    else if(z_1[iEntry1]>12.5 && !fChangeRS1) {
      gZ1_Up->SetPoint(point_up,y_1[iEntry1],z_1[iEntry1]);
      gZ1rev_Up->SetPoint(point_up,y_1[iEntry1],z_1[iEntry1]);
      point_up++;
    }
    else if(z_1[iEntry1]<12.5 && !fChangeRS1) {
      gZ1_Down->SetPoint(point_down,y_1[iEntry1],z_1[iEntry1]);
      gZ1rev_Down->SetPoint(point_down,y_1[iEntry1],z_1[iEntry1]);
      point_down++;
    }
  }
  point_up = 0;
  point_down = 0;
  for(unsigned int iEntry2=0; iEntry2<x_2.size(); iEntry2++) {
    if(z_2[iEntry2]>12.5) {
      gZ2_Up->SetPoint(point_up,y_2[iEntry2],z_2[iEntry2]);
      point_up++;
    }
    else {
      gZ2_Down->SetPoint(point_down,y_2[iEntry2],z_2[iEntry2]);
      point_down++;
    }
  }

  TCanvas* cRes = new TCanvas("cRes","",1200,600);
  cRes->Divide(3,1);
  cRes->cd(1);
  hResX->SetLineWidth(2);
  hResX->Draw();
  cRes->cd(2);
  hResY->SetLineWidth(2);
  hResY->Draw();
  cRes->cd(3);
  hResZ->SetLineWidth(2);
  hResZ->Draw();

  TH2F* h2FramePos = new TH2F("h2FramePos",Form("%s;y (mm);z (mm)",title1.Data()),1000,-1000,1000,100,9,15);
  h2FramePos->SetStats(0);
  h2FramePos->SetNdivisions(505);
  TH2F* h2FramePos2 = new TH2F("h2FramePos",Form("%s;y (mm);z (mm)",title4.Data()),1000,-1000,1000,100,9,15);
  h2FramePos2->SetStats(0);
  h2FramePos2->SetNdivisions(505);
  TH2F* h2FrameNeg = new TH2F("h2FramePos",Form("%s;y (mm);z (mm)",title2.Data()),1000,-1000,1000,100,-15,-9);
  h2FrameNeg->SetStats(0);
  h2FrameNeg->SetNdivisions(505);

  TF1* fZ1_Up = new TF1("fZ1_Up","pol2",-1000,1000);
  TF1* fZ1_Down = new TF1("fZ1_Down","pol2",-1000,1000);
  TF1* fZ1rev_Up = new TF1("fZ1rev_Up","pol2",-1000,1000);
  TF1* fZ1rev_Down = new TF1("fZ1rev_Down","pol2",-1000,1000);
  TF1* fZ2_Up = new TF1("fZ2_Up","pol2",-1000,1000);
  fZ2_Up->SetLineColor(kBlack);
  TF1* fZ2_Down = new TF1("fZ2_Down","pol2",-1000,1000);
  fZ2_Down->SetLineColor(kBlack);

  TCanvas* cZVsY = new TCanvas("cZVsY","",1500,500);
  cZVsY->Divide(3,1);
  cZVsY->cd(1);
  if(fChangeRS1) h2FrameNeg->Draw();
  else h2FramePos->Draw();
  gZ1_Up->Draw("P");
  gZ1_Up->Fit("fZ1_Up");
  gZ1_Down->Draw("P");
  gZ1_Down->Fit("fZ1_Down");
  cZVsY->cd(2);
  h2FramePos->Draw();
  gZ2_Up->Draw("P");
  gZ2_Down->Draw("P");
  cZVsY->cd(3);

  h2FramePos2->Draw();
  gZ2_Up->Draw("P");
  gZ2_Up->Fit("fZ2_Up");
  gZ2_Down->Draw("P");
  gZ2_Down->Fit("fZ2_Down");
  gZ1rev_Up->Draw("P");
  gZ1rev_Up->Fit("fZ1rev_Up");
  gZ1rev_Down->Draw("P");
  gZ1rev_Down->Fit("fZ1rev_Down");

  TString outfilename = infile1;
  outfilename.ReplaceAll("RS202","ResXYZ");
  outfilename.ReplaceAll("RS203","ResXYZ");
  outfilename.ReplaceAll("RS212","ResXYZ");
  outfilename.ReplaceAll("RS222","ResXYZ");
  outfilename.ReplaceAll("txt","pdf");
  outfilename.ReplaceAll("dat","pdf");
  cRes->SaveAs(outfilename.Data());
  outfilename.ReplaceAll("ResXYZ","Zprofile");
  cZVsY->SaveAs(outfilename.Data());

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
bool MatchMeas(double xmeas1, double ymeas1, double xmeas2, double ymeas2, double accradius) {

  if(TMath::Sqrt((xmeas1-xmeas2)*(xmeas1-xmeas2)+(ymeas1-ymeas2)*(ymeas1-ymeas2))<accradius) return true;

  return false;
}

//______________________________________________________________________________________________
void SetStyle() {
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleFont(42);
  gStyle->SetLabelFont(42);
  gStyle->SetTitleSize(0.05,"t");
  gStyle->SetLabelSize(0.045,"xyzt");
  gStyle->SetTitleOffset(1.,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetOptStat(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadGridX(1);

}
