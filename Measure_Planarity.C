#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <TInterpreter.h>
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
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TLatex.h>
#include <TDirectory.h>
#include <TBox.h>
#include <TVirtualFitter.h>
#include <TPaveText.h>

#endif

//*******************************************//
//                                           //
//    Main Functions: Measure_Planarity()    //
//                                           //
//*******************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIBALES

//all values are in mm

//input file
//const TString FileName="/home/bianca/Downloads/StaveMetrology/StaveMetrology/STAVE4/HS_left/ALC-0312-01_246_PLANARITY_ALLMODULES_2018_3_29_beforeUarms.dat";
TString FileName="./";
//flag to set the type of file (simple text or from Mitutoyo)
const bool fMitutoyoFile=true;

//number of series of measurements with a given x coordinate
const int nXcoord=2;
const double xnom[nXcoord]={-15,15};
const double xtoll = 3.5;

//flag to set type of planarity (wrt average plane or wrt to nominal plane)
const bool fPlanarityWRTnominal=true;
const double Znominal = 0.41;

//suction planes positions
const int nplanes = 12;
const double yleft[nplanes] = {-745.07,-685.07,-585.07,-435.07,-235.07,-75.07,34.93,194.93,394.93,534.93,644.93,744.93};
const double yright[nplanes] = {-744.07,-655.07,-555.07,-395.07,-205.07,-35.07,74.93,224.93,434.93,584.93,684.93,745.93};
const bool drawplanes=false; //draw suction planes to visualize where they are
const bool onlywithinplanes=false; //draw only measurements within suction planes

//modules positions
const int nmodules = 8;
const double ymodlims[nmodules] = {-738.5,-527.5,-316.5,-105.5,105.5,316.5,527.5,738.5};
const bool drawmodulespositions=false; //draw module positions

//scale factor for z coordinate (to change measure unit)
const int scalefactor = 1000; //from mm to mum (only z)

//coordinates for canvas dimension
const double LY = 1500;//total length along y
const double LX = 30;//total length along x
const double xmin=-20;//minimum x for plotting
const double xmax=20;//maximum x for plotting
const double ymin=-1000;//minimum y for plotting
const double ymax=1000;//maximum y for plotting
const double zmin=-1000;//minimum z for plotting (residuals)
const double zmax=1000;//maximum z for plotting (residuals)
const int colors[4]={kRed,kBlue,kGreen+2,kBlack};
const int markers [4]={kFullCircle,kFullSquare,kFullDiamond,kFullTriangleUp};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
int Measure_Planarity(const char* flname="./");
bool ReadFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
void CreateTxtFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
bool IsInside(double y);
void GetMeanSigmaAndPlanarity(std::vector<double> z, double& mean, double& sigma, double& planarity);
void GetDirCosines(TF2* fplane, double dircos[3]);
void SetStyle();

//_____________________________________________________________________________________________
//FUNCTION FOR PLANARITY MEASUREMENT
int Measure_Planarity(const char* flname)
{
  FileName=flname;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  bool read=false;
  if(!fMitutoyoFile) read=ReadFile(FileName,x,y,z);
  else read=ReadDatFileMitutoyo(FileName,x,y,z);

  if(!read) {
    cerr << "Impossibile to find the input file. Exit" << endl;
    return 1;
  }
  if(x.size()<=0 || y.size()<=0 || z.size()<=0) {
    cerr << "The number of points is <=0. Exit" << endl;
    return 2;
  }
  if(y.size()!=z.size() || x.size()!=z.size()) {
    cerr << "There is a discrepacy between the number of entries for the different coordinates. Exit." << endl;
    return 3;
  }

  TGraph2D* g = new TGraph2D(y.size());
  for(unsigned int iEntry=0; iEntry<y.size(); iEntry++) {
    g->SetPoint(iEntry,x[iEntry],y[iEntry],z[iEntry]);
  }

  TF2* fplane = new TF2("fplane","[0]+[1]*x+[2]*y",xmin,xmax,ymin,ymax);

  TCanvas* c = new TCanvas("c","",1200,900);
  TH3F* hFrame3D = new TH3F("hFrame3D",";x (mm);y (mm);z (#mum)",100,-20,20,100,-1000,1000,100,0,2000);
  c->cd()->SetGrid();
  g->SetMarkerStyle(20);
  fplane->GetXaxis()->SetTitle("x (mm)");
  fplane->GetYaxis()->SetTitle("y (mm)");
  fplane->GetZaxis()->SetTitle("z (#mum)");
  g->SetTitle("");
  g->Fit("fplane");
  fplane->SetTitle("");
  hFrame3D->Draw();
  fplane->Draw("surf same");
  g->Draw("P0 same");
  c->Update();

  double para = fplane->GetParameter(0);
  double parb = fplane->GetParameter(1);
  double parc = fplane->GetParameter(2);
  if(fPlanarityWRTnominal) {
    para=Znominal*scalefactor;
    parb=0.;
    parc=0.;
  }

  int nPoints=0;

  for(unsigned int iEntry=0; iEntry<y.size(); iEntry++) {
    if(onlywithinplanes && !IsInside(y[iEntry]))
      continue;
    nPoints++;
  }

  TGraph2D* gcorr = new TGraph2D(nPoints);
  gcorr->SetName("gcorr");
  TGraph** gcorr_x = new TGraph*[nXcoord];
  for(int iX=0; iX<nXcoord; iX++) {
    gcorr_x[iX] = new TGraph();
    gcorr_x[iX]->SetName(Form("gcorr_%d",iX));
  }
  TGraph* gcorr_xint = new TGraph(0);
  gcorr_xint->SetName("gcorr_xint");

  std::vector<double> zcorr;
  std::vector<double> zinside;
  unsigned int poscounter=0;
  unsigned int iPoint[nXcoord];
  for(int iX=0; iX<nXcoord; iX++) {iPoint[iX]=0;}
  unsigned int counter=0;
  for(unsigned int iEntry=0; iEntry<nPoints; iEntry++) {

    zcorr.push_back(z[iEntry]-(para+parb*x[iEntry]+parc*y[iEntry]));

    if(onlywithinplanes && !IsInside(y[iEntry]))
      continue;

    zinside.push_back(zcorr[iEntry]);
    gcorr->SetPoint(iEntry,x[iEntry],y[iEntry],zcorr[iEntry]);
    gcorr_xint->SetPoint(iEntry,y[iEntry],zcorr[iEntry]);
    for(int iX=0; iX<nXcoord; iX++) {
      if(x[iEntry]>xnom[iX]-xtoll && x[iEntry]<xnom[iX]+xtoll) {
        gcorr_x[iX]->SetPoint(iPoint[iX],y[iEntry],zcorr[iEntry]);
        iPoint[iX]++;
      }
    }
  }

  double mean=0;
  double sigma=0;
  double planarity=0;
  GetMeanSigmaAndPlanarity(zinside,mean,sigma,planarity);
  SetStyle();

  TBox* box[nplanes];
  for(int iPlane=0; iPlane<nplanes; iPlane++) {
    box[iPlane] = new TBox(yleft[iPlane],-100,yright[iPlane],100);
    box[iPlane]->SetLineColor(kRed);
    box[iPlane]->SetLineWidth(1);
    box[iPlane]->SetFillStyle(0);
  }

  TLine* line[nmodules];
  for(int iMod=0; iMod<nmodules; iMod++) {
    line[iMod] = new TLine(ymodlims[iMod],-100,ymodlims[iMod],100);
    line[iMod]->SetLineColor(kRed);
    line[iMod]->SetLineWidth(1);
  }

  TCanvas* ccorr2D = new TCanvas("ccorr2D","",1200,900);
  ccorr2D->cd()->SetGrid();
  gcorr->SetTitle("");
  gcorr->GetXaxis()->SetTitle("x (mm)");
  gcorr->GetYaxis()->SetTitle("y (mm)");
  gcorr->GetZaxis()->SetTitle("y (#mum)");
  gcorr->GetZaxis()->SetRangeUser(zmin,zmax);
  gcorr->SetMarkerStyle(20);
  gcorr->SetTitle("");
  gcorr->Draw("AP");

  TLegend* l = new TLegend(0.15,0.15,0.4,0.30);
  l->SetFillColor(kWhite);
  l->SetTextSize(0.045);
  for(int iX=0; iX<nXcoord; iX++) {
    l->AddEntry(gcorr_x[iX],Form("x = %0.f mm",xnom[iX]),"p");
  }

  TPaveText* info = new TPaveText(0.6,0.7,0.89,0.9,"NDC");
  info->SetTextSize(0.045);
  info->SetTextFont(42);
  info->SetFillColor(kWhite);
  info->SetBorderSize(1);
  info->AddText(Form("mean = %0.f #mum",mean));
  info->AddText(Form("RMS = %0.f #mum",sigma));
  info->AddText(Form("planarity = %0.f #mum",planarity));

  TCanvas* ccorr = new TCanvas("ccorr","",1200,900);
  ccorr->cd()->SetGrid();
  TString drawopt="AP";
  for(int iX=0; iX<nXcoord; iX++) {
    gcorr_x[iX]->SetTitle("");
    gcorr_x[iX]->GetXaxis()->SetTitle("y (mm)");
    gcorr_x[iX]->GetYaxis()->SetTitle("z_{meas} - z_{plane} (#mum)");
    gcorr_x[iX]->GetYaxis()->SetRangeUser(zmin,zmax);
    gcorr_x[iX]->SetMarkerStyle(markers[iX]);
    gcorr_x[iX]->SetMarkerColor(colors[iX]);
    gcorr_x[iX]->SetTitle("");
    gcorr_x[iX]->Draw(drawopt.Data());
    drawopt="P";
  }

  if(drawplanes) {
    for(int iPlane=0; iPlane<nplanes; iPlane++) {
      box[iPlane]->Draw("same");
    }
  }
  if(drawmodulespositions) {
    TLatex* lat = new TLatex();
    lat->SetTextFont(42);
    lat->SetTextSize(0.045);
    lat->SetTextColor(kRed);
    for(int iMod=0; iMod<nmodules; iMod++) {
      line[iMod]->Draw("same");
      cout << line[iMod]->GetX1() << endl;
      if(iMod<nmodules-1) lat->DrawLatex(line[iMod]->GetX1()+18,line[iMod]->GetY1()*1.2,Form("MOD%d",iMod+1));
    }
  }
  l->Draw("same");
  info->Draw("same");

  TString OutFileNamewoext=FileName;
  OutFileNamewoext.ReplaceAll(".txt","_corr");

  if(!onlywithinplanes) {
    ccorr2D->SaveAs(Form("%s2D.pdf",OutFileNamewoext.Data()));
    c->SaveAs(Form("%s2Dfit.pdf",OutFileNamewoext.Data()));
    ccorr->SaveAs(Form("%s.pdf",OutFileNamewoext.Data()));
    TString outfiletxt = Form("%s.txt",OutFileNamewoext.Data());
    CreateTxtFile(outfiletxt,x,y,zcorr);
  } else {
    ccorr2D->SaveAs(Form("%s_insideplanes_2D.pdf",OutFileNamewoext.Data()));
    ccorr->SaveAs(Form("%s_insideplanes_.pdf",OutFileNamewoext.Data()));
    TString outfiletxt = Form("%s_insideplanes_.txt",OutFileNamewoext.Data());
    CreateTxtFile(outfiletxt,x,y,zcorr);
  }
  TFile outfile(Form("%s.root",OutFileNamewoext.Data()),"RECREATE");
  gcorr->Write();
  for(int iX=0; iX<nXcoord; iX++) {
    gcorr_x[iX]->Write();
  }
  gcorr_xint->Write();
  fplane->Write();
  outfile.Close();

  double InclY = TMath::ATan(fplane->GetParameter(2));
  double InclYgrad = InclY/TMath::Pi()*180;
  double InclX = TMath::ATan(fplane->GetParameter(1));
  double InclXgrad = InclX/TMath::Pi()*180;
  cout << "\nInclination:" << endl;
  cout << "Inclination in x = " << InclXgrad << " deg"<< endl;
  cout << "Inclination in y = " << InclYgrad << " deg"<< endl;

  double DeltaX = TMath::Abs(LX/2*TMath::Sin(InclX));
  double DeltaY = TMath::Abs(LY/2*TMath::Sin(InclY));
  cout << "\nDelta:" << endl;
  cout << "DeltaZ (x-axis) = " << DeltaX << " mum"<< endl;
  cout << "DeltaZ (y-axis) = " << DeltaY << " mum"<< endl;
  cout << "\n" << endl;

  return 0;
}

//_____________________________________________________________________________________________
bool ReadFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
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
    x.push_back(xcoord);
    y.push_back(ycoord);
    z.push_back(zcoord*scalefactor);
  }

  inSet.close();

  return true;
}

//_____________________________________________________________________________________________
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
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
      vector<string> values;
      string line;
      getline(inSet, line);
      TString test = line;
      if(test.Contains("Gauss") || test.Contains("Mean") || linecounter>0) {
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
            x.push_back(xtmp);
          }
          else if(values.size()==6) {
            stringstream convert(values[5]);
            if( !(convert >> ytmp) ) ytmp = -1;
            y.push_back(ytmp);
          }
          else if(values.size()==7) {
            stringstream convert(values[6]);
            if( !(convert >> ztmp) ) ztmp = -1;
            z.push_back(ztmp*scalefactor);
          }
          line = line.substr(pos + 1);
        }
      }
    }
    inSet.close();
  }

  return true;
}

//_____________________________________________________________________________________________
void CreateTxtFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
{
  ofstream outSet;
  outSet.open(FileName.Data());

  for(unsigned int iEntry=0; iEntry<z.size(); iEntry++) {
    outSet << x[iEntry] << " " << y[iEntry] << " " << z[iEntry] << endl;
  }

  outSet.close();
}

bool IsInside(double y)
{

  bool isinside=true;
  int iPlane=0;

  while(isinside==true && iPlane<nplanes) {
    if(y<yleft[iPlane] || y>yright[iPlane])
      isinside=false;
    iPlane++;
  }

  return isinside;

}

//_____________________________________________________________________________________________
void GetMeanSigmaAndPlanarity(std::vector<double> z, double& mean, double& sigma, double& planarity)
{

  mean=0;
  sigma=0;
  double max=z[0];
  double min=z[0];

  for(unsigned int iEntry=0; iEntry<z.size(); iEntry++) {
    mean += z[iEntry];
    if(z[iEntry]>max)
      max=z[iEntry];
    if(z[iEntry]<min)
      min=z[iEntry];
  }
  cout << min << "  " << max << endl;
  mean /= z.size();

  for(unsigned int iEntry=0; iEntry<z.size(); iEntry++) {
    sigma += (z[iEntry]-mean)*(z[iEntry]-mean);
  }
  sigma /= (z.size()-1);
  sigma = TMath::Sqrt(sigma);

  planarity = max-min;

}

//_____________________________________________________________________________________________
void GetDirCosines(TF2* fplane, double dircos[3])
{

  double c = fplane->GetParameter(0);
  double a = fplane->GetParameter(1);
  double b = fplane->GetParameter(2);
  double norm = TMath::Sqrt(a*a+b*b+1);

  if(c>0) norm -= norm;

  // double norm = TMath::Sqrt();
  dircos[0] = a/norm;
  dircos[1] = b/norm;
  dircos[2] = -1/norm;

}

//______________________________________________________________________________________________
void SetStyle() {
  cout << "Setting drawing style!" << endl;
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleFont(42);
  gStyle->SetLabelFont(42);
  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetLabelSize(0.045,"xyzt");
  gStyle->SetTitleOffset(1.,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
}
