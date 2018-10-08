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
//    Main Functions: ShopwCP_Planarity()    //
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

const int colors[4]={kRed,kBlue,kGreen+2,kBlack};
const int markers [4]={kFullCircle,kFullSquare,kFullDiamond,kFullTriangleUp};

enum EType {
  kPlaneCP,
  kPlaneHSR,
  kPlaneHSL
};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
int Show_CPPlanarity_HSPlane(const char* flname="./", EType _type=kPlaneCP);
bool ReadFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
void SetStyle();

//_____________________________________________________________________________________________
//FUNCTION FOR PLANARITY MEASUREMENT
int Show_CPPlanarity_HSPlane(const char* flname, EType _type)
{
  SetStyle();
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

  unsigned int nPlot=0;
  float z_min = 0.;
  float z_max = 0.;
  //scale factor for z coordinate (to change measure unit)
  unsigned int scalefactor = 1000; //from mm to mum (only z)

  if (_type==kPlaneCP){
    nPlot=3;
    z_min = 0.;
    z_max = 500.;
  } else if (_type==kPlaneHSR){
    nPlot=3;
    scalefactor = 1;
    z_min = -12.;
    z_max = -8.;
  } else if (_type==kPlaneHSL){
    nPlot=2;
    scalefactor = 1;
    z_min = -15.;
    z_max = -11.;
  }

  TLegend* l = new TLegend(0.15,0.15,0.4,0.30);
  l->SetFillColor(kWhite);
  l->SetTextSize(0.045);
  TGraph* graph[static_cast<const unsigned int>(nPlot)];
  for(unsigned int i=0; i<nPlot; ++i){
    graph[i]=new TGraph(y.size()/nPlot);
    graph[i]->SetMarkerStyle(20);
    graph[i]->SetMarkerColor(colors[i]);
    graph[i]->SetName(Form("gr_%.0fmm", x[i]));
    graph[i]->SetTitle(Form("gr_%.0fmm", x[i]));
    l->AddEntry(graph[i], Form("x = %.0f mm", x[i]), "p");
  }
  for(unsigned int iEntry=0; iEntry<y.size(); iEntry++) {
    Printf("%f, %f, %f", x[iEntry], y[iEntry], z[iEntry]);
    for(unsigned int iPlot=0; iPlot<nPlot; ++iPlot){
      if(iEntry%nPlot==iPlot){
        graph[iPlot]->SetPoint(iEntry/nPlot,y[iEntry],z[iEntry]*scalefactor);
      }
    }
  }
//
//  TF2* fplane = new TF2("fplane","[0]+[1]*x+[2]*y",xmin,xmax,ymin,ymax);
//
  TCanvas* c = new TCanvas("c","",1200,900);
  TH2F* hFrame2D = new TH2F("hFrame2D", Form(";y (mm);z (%s)",(scalefactor==1?"mm":"#mum")),100,-1000,1000,100,z_min,z_max);
  c->cd()->SetGrid();

  hFrame2D->Draw();
  for(unsigned int i=0; i<nPlot; ++i)
    graph[i]->Draw("P0 same");
  l->Draw("same");
  c->Update();

  TString prefix_fl_name(FileName);
  prefix_fl_name.ReplaceAll(".csv","");
  c->SaveAs(Form("%s.pdf", prefix_fl_name.Data()));

  TFile* root_fl = TFile::Open(Form("%s.root", prefix_fl_name.Data()), "recreate");
  for(unsigned int iPlot=0; iPlot<nPlot; ++iPlot)  {
    graph[iPlot]->Write();
  }

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
    z.push_back(zcoord);
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
      vector<string> values;
      string line;
      getline(inSet, line);
      TString test = line;
      if(test.Contains("Piano") || test.Contains("Plane"))
        islineok=false;
      if(islineok) {
        size_t pos = 0;
        values.clear();
        while((pos = line.find(valueSeparator)) != string::npos)
        {
          double xtmp=-1.;
          double ytmp=-1.;
          double ztmp=-1.;
          values.push_back(line.substr(0,pos));
          if(values.size()==5) {
            stringstream convert(values[4]);
            if( !(convert >> xtmp) ) xtmp = -2;
            x.push_back(xtmp);
          }
          else if(values.size()==6) {
            stringstream convert(values[5]);
            if( !(convert >> ytmp) ) ytmp = -2;
            y.push_back(ytmp);
          }
          else if(values.size()==7) {
            stringstream convert(values[6]);
            if( !(convert >> ztmp) ) ztmp = -2;
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
void SetStyle() {
  cout << "Setting drawing style!" << endl;
  gStyle->SetOptStat(0000);
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
