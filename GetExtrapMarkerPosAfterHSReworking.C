#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TSystem.h>
#include <Riostream.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TPad.h>
#include <TString.h>
#include <TLine.h>
#include <TPaveText.h>

#endif

//**********************************************************************//
//                                                                      //
//        Main Function: GetExtrapMarkerPosAfterHSReworking()           //
//                                                                      //
//**********************************************************************//

//author: Fabrizio Grosa, INFN Torino
//e-mail: fabrizio.grosa@to.infn.it

using namespace std;
enum {kHSL, kHSR};

//_____________________________________________________________________________________________
//METHOD PROTOTYPES
int GetExtrapMarkerPosAfterHSReworking(TString file_before_rew = "", TString fileAfterRew = "", int LeftOrRight = kHSL);
void DoInvisibleMarkerExtrapolation(int LeftOrRight, vector<double> x_nominal, vector<double> y_nominal, vector<double> x_after_rew_filled, vector<double> y_after_rew_filled, vector<double> z_after_rew_filled, vector<double> x_before_rew_filled, vector<double> y_before_rew_filled, vector<double> z_before_rew_filled, vector<double>& x_extrap, vector<double>& y_extrap, vector<double>& z_extrap);
bool MatchMeasToNominal(double xmeas, double ymeas, double xnom, double ynom, double accdeltax, double accdeltay, double yshift);
void FillVectorsOfMatchedPositions(vector<double> xmeas, vector<double> ymeas, vector<double> zmeas, vector<double> xnom, vector<double> ynom, vector<double>& xmeasfilled, vector<double>& ymeasfilled, vector<double>& zmeasfilled);
int FindOppositeNominalMarkerPosition(double x, double y, vector<double> xvec, vector<double> yvec);
bool ReadDatFile(TString FileName, vector<double>& x, vector<double>& y, vector<double>& z);
bool ReadDatFileMitutoyo(TString FileName, vector<double>& x, vector<double>& y, vector<double>& z);
void CreateDatFileMitutoyo(TString FileName, vector<double> x_after_rew_filled, vector<double> y_after_rew_filled, vector<double> z_after_rew_filled, vector<double> x_extrap, vector<double> y_extrap, vector<double> z_extrap);
void SetStyle();

//_____________________________________________________________________________________________
int GetExtrapMarkerPosAfterHSReworking(TString file_before_rew, TString file_after_rew, int LeftOrRight) {

  vector<double> x_before_rew, y_before_rew, z_before_rew;
  vector<double> x_after_rew, y_after_rew, z_after_rew;
  vector<double> x_nominal, y_nominal, z_nominal;

  bool read_before_rew = ReadDatFileMitutoyo(file_before_rew,x_before_rew,y_before_rew,z_before_rew);
  if(!read_before_rew) return 1;
  bool read_after_rew = ReadDatFileMitutoyo(file_after_rew,x_after_rew,y_after_rew,z_after_rew);
  if(!read_after_rew) return 2;
  bool read_nominal_pos = ReadDatFile("OL_marker_nominal_positions_HS.dat",x_nominal,y_nominal,z_nominal);

  vector<double> x_after_rew_filled, y_after_rew_filled, z_after_rew_filled;
  vector<double> x_before_rew_filled, y_before_rew_filled, z_before_rew_filled;
  FillVectorsOfMatchedPositions( x_after_rew, y_after_rew, z_after_rew, x_nominal, y_nominal, x_after_rew_filled, y_after_rew_filled, z_after_rew_filled);
  FillVectorsOfMatchedPositions( x_before_rew, y_before_rew, z_before_rew, x_nominal, y_nominal, x_before_rew_filled, y_before_rew_filled, z_before_rew_filled);
  
  vector<double> x_extrap, y_extrap, z_extrap;
  DoInvisibleMarkerExtrapolation(LeftOrRight,x_nominal,y_nominal,x_after_rew_filled,y_after_rew_filled,z_after_rew_filled,x_before_rew_filled,y_before_rew,z_before_rew_filled,x_extrap,y_extrap,z_extrap);
  
  //define graphs and histos
  TGraph2D* g_before_rew = new TGraph2D(x_before_rew.size());
  g_before_rew->SetMarkerColor(kRed);
  g_before_rew->SetMarkerStyle(kFullCircle);
  for(unsigned int iMarker=0; iMarker<x_before_rew.size(); iMarker++) {
    g_before_rew->SetPoint(iMarker,x_before_rew[iMarker],y_before_rew[iMarker],z_before_rew[iMarker]);
  }

  TGraph2D* g_after_rew = new TGraph2D(x_after_rew.size());
  g_after_rew->SetMarkerColor(kBlue);
  g_after_rew->SetMarkerStyle(kFullSquare);
  for(unsigned int iMarker=0; iMarker<x_after_rew.size(); iMarker++) {
    g_after_rew->SetPoint(iMarker,x_after_rew[iMarker],y_after_rew[iMarker],z_after_rew[iMarker]);
  }

  TGraph* g_before_rew_posX = new TGraph(0);
  TGraph* g_before_rew_negX = new TGraph(0);
  g_before_rew_posX->SetMarkerColor(kRed);
  g_before_rew_posX->SetMarkerStyle(kFullCircle);
  g_before_rew_negX->SetMarkerColor(kBlue);
  g_before_rew_negX->SetMarkerStyle(kFullSquare);
  int iPosMarker = 0;
  int iNegMarker = 0;
  for(unsigned int iMarker=0; iMarker<x_before_rew.size(); iMarker++) {
    if(x_before_rew[iMarker]>0) {
      g_before_rew_posX->SetPoint(iPosMarker,y_before_rew[iMarker],z_before_rew[iMarker]*1000);
      iPosMarker++;
    }
    else {
      g_before_rew_negX->SetPoint(iNegMarker,y_before_rew[iMarker],z_before_rew[iMarker]*1000);
      iNegMarker++;
    }
  }

  TGraph* g_after_rew_posX = new TGraph(0);
  TGraph* g_after_rew_negX = new TGraph(0);
  g_after_rew_posX->SetMarkerColor(kRed);
  g_after_rew_posX->SetMarkerStyle(kFullCircle);
  g_after_rew_negX->SetMarkerColor(kBlue);
  g_after_rew_negX->SetMarkerStyle(kFullSquare);
  iPosMarker = 0;
  iNegMarker = 0;
  for(unsigned int iMarker=0; iMarker<x_after_rew.size(); iMarker++) {
    if(x_after_rew[iMarker]>0) {
      g_after_rew_posX->SetPoint(iPosMarker,y_after_rew[iMarker],z_after_rew[iMarker]*1000);
      iPosMarker++;
    }
    else {
      g_after_rew_negX->SetPoint(iNegMarker,y_after_rew[iMarker],z_after_rew[iMarker]*1000);
      iNegMarker++;
    }
  }

  TGraph* g_after_rew_withextr_posX = new TGraph(0);
  TGraph* g_after_rew_withextr_negX = new TGraph(0);
  g_after_rew_withextr_posX->SetMarkerColor(kRed);
  g_after_rew_withextr_posX->SetMarkerStyle(kFullCircle);
  g_after_rew_withextr_negX->SetMarkerColor(kBlue);
  g_after_rew_withextr_negX->SetMarkerStyle(kFullSquare);
  iPosMarker = 0;
  iNegMarker = 0;
  for(unsigned int iMarker=0; iMarker<x_nominal.size(); iMarker++) {
    if(x_after_rew_filled[iMarker]>0 && x_after_rew_filled[iMarker]!=-10000.) {
      g_after_rew_withextr_posX->SetPoint(iPosMarker,y_after_rew_filled[iMarker],z_after_rew_filled[iMarker]*1000);
      iPosMarker++;
    }
    else if(x_after_rew_filled[iMarker]<0 && x_after_rew_filled[iMarker]!=-10000.) {
      g_after_rew_withextr_negX->SetPoint(iNegMarker,y_after_rew_filled[iMarker],z_after_rew_filled[iMarker]*1000);
      iNegMarker++;
    }
    else if(x_extrap[iMarker]>0 && x_extrap[iMarker]!=-10000.) {
      g_after_rew_withextr_posX->SetPoint(iPosMarker,y_extrap[iMarker],z_extrap[iMarker]*1000);
      iPosMarker++;
    }
    else if(x_extrap[iMarker]<0 && x_extrap[iMarker]!=-10000.) {
      g_after_rew_withextr_negX->SetPoint(iNegMarker,y_extrap[iMarker],z_extrap[iMarker]*1000);
      iNegMarker++;
    }
  }

  TH1F* h_extrap_res_X = new TH1F("h_extrap_res_X","Extrapolated - measured marker position X;x_{extr}-x_{meas} (#mum);Entries",100,-50,50);
  TH1F* h_extrap_res_Y = new TH1F("h_extrap_res_Y","Extrapolated - measured marker position Y;y_{extr}-y_{meas} (#mum);Entries",100,-50,50);
  TH1F* h_extrap_res_Z = new TH1F("h_extrap_res_Z","Extrapolated - measured marker position Z;z_{extr}-z_{meas} (#mum);Entries",100,-50,50);
  h_extrap_res_X->SetLineWidth(2);
  h_extrap_res_Y->SetLineWidth(2);
  h_extrap_res_Z->SetLineWidth(2);
  for(unsigned int iMarker=0; iMarker<x_nominal.size(); iMarker++) {
    if(x_after_rew_filled[iMarker]!=-10000. && x_extrap[iMarker]!=-10000.) h_extrap_res_X->Fill((x_extrap[iMarker]-x_after_rew_filled[iMarker])*1000);
    if(y_after_rew_filled[iMarker]!=-10000. && y_extrap[iMarker]!=-10000.) h_extrap_res_Y->Fill((y_extrap[iMarker]-y_after_rew_filled[iMarker])*1000);
    if(z_after_rew_filled[iMarker]!=-10000. && z_extrap[iMarker]!=-10000.) h_extrap_res_Z->Fill((z_extrap[iMarker]-z_after_rew_filled[iMarker])*1000);
  }

  //plots
  SetStyle();

  TH3F* hFrameHS = new TH3F("hFrame",";x (mm);y (mm);z (mm)",100,-50.,50.,100,-1000,1000,100,0.,1.);
  hFrameHS->SetStats(0);

  TLegend* leg_before_after = new TLegend(0.6,0.75,0.9,0.9);
  leg_before_after->SetTextSize(0.045);
  leg_before_after->AddEntry(g_before_rew,"Before reworking","p");
  leg_before_after->AddEntry(g_after_rew,"After reworking","p");

  TLegend* leg_pos_neg = new TLegend(0.8,0.75,0.9,0.9);
  leg_pos_neg->SetTextSize(0.05);
  leg_pos_neg->AddEntry(g_before_rew_posX,"X = 15 mm","p");
  leg_pos_neg->AddEntry(g_before_rew_negX,"X = -15 mm","p");

  TCanvas* cHS3D = new TCanvas("cHS3D","",1920,1080);
  hFrameHS->Draw();
  g_before_rew->Draw("Psame");
  g_after_rew->Draw("Psame");
  leg_before_after->Draw();

  TCanvas* cHS = new TCanvas("cHS","",1920,1080);
  cHS->Divide(1,3);
  cHS->cd(1)->DrawFrame(-800,200,800,800,"Before reworking;y (mm);z (#mum)");
  g_before_rew_posX->Draw("P");
  g_before_rew_negX->Draw("P");
  leg_pos_neg->Draw();
  cHS->cd(2)->DrawFrame(-800,200,800,800,"After reworking;y (mm);z (#mum)");
  g_after_rew_posX->Draw("P");
  g_after_rew_negX->Draw("P");
  leg_pos_neg->Draw();
  cHS->cd(3)->DrawFrame(-800,200,800,800,"After reworking with extr;y (mm);z (#mum)");
  g_after_rew_withextr_posX->Draw("P");
  g_after_rew_withextr_negX->Draw("P");
  leg_pos_neg->Draw();

  TCanvas* cExtrapQuality = new TCanvas("cExtrapQuality","",1920,1080);
  cExtrapQuality->Divide(3,1);
  cExtrapQuality->cd(1);
  h_extrap_res_X->Draw();
  cExtrapQuality->cd(2);
  h_extrap_res_Y->Draw();
  cExtrapQuality->cd(3);
  h_extrap_res_Z->Draw();

  //output files
  TString outfile_root = file_after_rew;
  TString outfile_pdf = file_after_rew;
  TString outfile_dat = file_after_rew;
  outfile_root.ReplaceAll(".dat","_AfterRew_WithExtrap.root");
  outfile_pdf.ReplaceAll(".dat","_AfterRew_WithExtrap.pdf");
  outfile_dat.ReplaceAll(".dat","_AfterRew_WithExtrap.dat");

  TFile outfile(outfile_root.Data(),"recreate");
  g_before_rew->Write("gHSBeforeRework");
  g_after_rew->Write("gHSAfterRework");
  g_before_rew_posX->Write("gHSBeforeRework_PosX");
  g_before_rew_negX->Write("gHSBeforeRework_NegX");
  g_after_rew_posX->Write("gHSAfterRework_PosX");
  g_after_rew_negX->Write("gHSAfterRework_NegX");
  g_after_rew_withextr_posX->Write("gHSAfterReworkWithExtrap_PosX");
  g_after_rew_withextr_negX->Write("gHSAfterReworkWithExtrap_NegX");
  h_extrap_res_X->Write("hExtrapResX");
  h_extrap_res_Y->Write("hExtrapResY");
  h_extrap_res_Z->Write("hExtrapResZ");
  cHS3D->Write();
  cHS->Write();
  cExtrapQuality->Write();
  outfile.Close();

  cHS3D->Print(Form("%s[",outfile_pdf.Data()));
  cHS3D->Print(Form("%s",outfile_pdf.Data()));
  cHS->Print(Form("%s",outfile_pdf.Data()));
  cExtrapQuality->Print(Form("%s",outfile_pdf.Data()));
  cExtrapQuality->Print(Form("%s]",outfile_pdf.Data()));

  CreateDatFileMitutoyo(outfile_dat,x_after_rew_filled,y_after_rew_filled,z_after_rew_filled,x_extrap,y_extrap,z_extrap);

  return 0;
}

//______________________________________________________________________________________________
void DoInvisibleMarkerExtrapolation(int LeftOrRight, vector<double> x_nominal, vector<double> y_nominal, vector<double> x_after_rew_filled, vector<double> y_after_rew_filled, vector<double> z_after_rew_filled, vector<double> x_before_rew_filled, vector<double> y_before_rew_filled, vector<double> z_before_rew_filled, vector<double>& x_extrap, vector<double>& y_extrap, vector<double>& z_extrap) {

  //positions of initial markers in vector
  vector<unsigned int>  modposcornermarkers = {14,27,42,55,70,83,98,111,126,139,154,167,182,195};
  vector<unsigned int>  modnegcornermarkers = {0,13,28,41,56,69,84,97,112,125,140,153,168,181};

  vector<unsigned int>::iterator it;
  bool isposcorner = false;
  bool isnegcorner = false;

  for(unsigned int iMarker=0; iMarker<x_nominal.size(); iMarker++) {
    if((LeftOrRight==kHSL && x_nominal[iMarker]>0) || (LeftOrRight==kHSR && x_nominal[iMarker]<0)) {
      int posoppmarker = FindOppositeNominalMarkerPosition(x_nominal[iMarker],y_nominal[iMarker],x_nominal,y_nominal);

      it = find(modposcornermarkers.begin(), modposcornermarkers.end(), iMarker);
      if(it!=modposcornermarkers.end()) isposcorner=true;
      if(!isposcorner) {
        it = find(modnegcornermarkers.begin(), modnegcornermarkers.end(), iMarker);
        if(it!=modnegcornermarkers.end()) isnegcorner=true;
      }

      const int noppmarkers=3;
      int oppmarkarray[noppmarkers] = {posoppmarker-1,posoppmarker,posoppmarker+1};
      if(isposcorner) {
        oppmarkarray[0] = posoppmarker;
        oppmarkarray[1] = posoppmarker-1;
        oppmarkarray[2] = posoppmarker-2;
      }
      else if(isnegcorner) {
        oppmarkarray[0] = posoppmarker;
        oppmarkarray[1] = posoppmarker+1;
        oppmarkarray[2] = posoppmarker+2;
      }

  	  double deltax[noppmarkers] = {0,0,0};
  	  double deltay[noppmarkers]= {0,0,0};
  	  double deltaz[noppmarkers]= {0,0,0};

      double rawexrap_x[noppmarkers] = {-10000,-10000,-10000};
      double rawexrap_y[noppmarkers] = {-10000,-10000,-10000};
      double rawexrap_z[noppmarkers] = {-10000,-10000,-10000};

  	  for (int j=0; j<noppmarkers; j++) {
        if(x_before_rew_filled[oppmarkarray[j]]!=-10000 && x_before_rew_filled[iMarker]!=-10000 && x_after_rew_filled[oppmarkarray[j]]!=-10000) {
          deltax[j]= x_before_rew_filled[iMarker]-x_before_rew_filled[oppmarkarray[j]];
          deltay[j] = y_before_rew_filled[iMarker]-y_before_rew_filled[oppmarkarray[j]];
     			deltaz[j] = z_before_rew_filled[iMarker]-z_before_rew_filled[oppmarkarray[j]];
		      rawexrap_x[j]=x_after_rew_filled[oppmarkarray[j]]+deltax[j];
			    rawexrap_y[j]=y_after_rew_filled[oppmarkarray[j]]+deltay[j];
			    rawexrap_z[j]=z_after_rew_filled[oppmarkarray[j]]+deltaz[j];
		    }
	    }

      double sumextrapolation_x=0;
      double sumextrapolation_y=0;
      double sumextrapolation_z=0;
      int countpoint=0;

      if(!((TMath::Abs(rawexrap_x[0]-rawexrap_x[1])>0.2 && TMath::Abs(rawexrap_x[0]-rawexrap_x[2])>0.2) || (TMath::Abs(rawexrap_y[0]-rawexrap_y[1])>0.2 && TMath::Abs(rawexrap_y[0]-rawexrap_y[2])>0.2) || (TMath::Abs(rawexrap_z[0]-rawexrap_z[1])>0.2 && TMath::Abs(rawexrap_z[2]-rawexrap_z[2])>0.2))) {
				sumextrapolation_x+=rawexrap_x[0];
				sumextrapolation_y+=rawexrap_y[0];
				sumextrapolation_z+=rawexrap_z[0];
				countpoint++;
			}
      if(!((TMath::Abs(rawexrap_x[1]-rawexrap_x[0])>0.2 && TMath::Abs(rawexrap_x[1]-rawexrap_x[2])>0.2) || (TMath::Abs(rawexrap_y[1]-rawexrap_y[0])>0.2 && TMath::Abs(rawexrap_y[1]-rawexrap_y[2])>0.2) || (TMath::Abs(rawexrap_z[1]-rawexrap_z[0])>0.2 && TMath::Abs(rawexrap_z[1]-rawexrap_z[2])>0.2))) {
				sumextrapolation_x+=rawexrap_x[1];
				sumextrapolation_y+=rawexrap_y[1];
				sumextrapolation_z+=rawexrap_z[1];
				countpoint++;
			}
      if(!((TMath::Abs(rawexrap_x[2]-rawexrap_x[0])>0.2 && TMath::Abs(rawexrap_x[1]-rawexrap_x[2])>0.2) || (TMath::Abs(rawexrap_y[2]-rawexrap_y[0])>0.2 && TMath::Abs(rawexrap_y[1]-rawexrap_y[2])>0.2) || (TMath::Abs(rawexrap_z[2]-rawexrap_z[0])>0.2 && TMath::Abs(rawexrap_z[1]-rawexrap_z[2])>0.2))) {
        sumextrapolation_x+=rawexrap_x[2];
        sumextrapolation_y+=rawexrap_y[2];
        sumextrapolation_z+=rawexrap_z[2];
        countpoint++;
      }
			if(countpoint==0){ //if all the three distant from each other, use the first
				sumextrapolation_x+=rawexrap_x[0];
				sumextrapolation_y+=rawexrap_y[0];
				sumextrapolation_z+=rawexrap_z[0];
				countpoint=1;
			}

      x_extrap.push_back(sumextrapolation_x/countpoint);
      y_extrap.push_back(sumextrapolation_y/countpoint);
      z_extrap.push_back(sumextrapolation_z/countpoint);
    }
    else {
      x_extrap.push_back(-10000.);
      y_extrap.push_back(-10000.);
      z_extrap.push_back(-10000.);
    }
    isposcorner=false;
    isnegcorner=false;
  }
}

//______________________________________________________________________________________________
bool MatchMeasToNominal(double xmeas, double ymeas, double xnom, double ynom, double accdeltax, double accdeltay, double yshift) {

  if(TMath::Abs(xmeas-xnom)<accdeltax && TMath::Abs((ymeas-yshift)-ynom)<accdeltay) return true;

  return false;
}

//______________________________________________________________________________________________
void FillVectorsOfMatchedPositions(vector<double> xmeas, vector<double> ymeas, vector<double> zmeas, vector<double> xnom, vector<double> ynom, vector<double>& xmeasfilled, vector<double>& ymeasfilled, vector<double>& zmeasfilled) {

  //positions of initial markers in vector
  vector<unsigned int>  modinitmarkers = {0,28,56,84,112,140,168};

  double yshift=0.;
  vector<unsigned int>::iterator it;
  int countermeasmarkers=0;
  int matchedmarkerpos=-1;

  for(unsigned int iNomEntry=0; iNomEntry<xnom.size(); iNomEntry++) {
    bool measfound=false;
    it = find(modinitmarkers.begin(), modinitmarkers.end(), iNomEntry);
    if(it!=modinitmarkers.end()) {
      yshift=0.; //reset y shift before starting each module
      matchedmarkerpos=countermeasmarkers;
      measfound=true;
      countermeasmarkers++;
    }
    else {
      for(unsigned int iMeasEntry=0; iMeasEntry<xmeas.size(); iMeasEntry++) {
        if(MatchMeasToNominal(xmeas[iMeasEntry],ymeas[iMeasEntry],xnom[iNomEntry],ynom[iNomEntry],0.800,0.350,yshift)) {
          measfound=true;
          countermeasmarkers++;
          matchedmarkerpos=iMeasEntry;
          break;
        }
      }
    }
    if(measfound) {
      xmeasfilled.push_back(xmeas[matchedmarkerpos]);
      ymeasfilled.push_back(ymeas[matchedmarkerpos]);
      zmeasfilled.push_back(zmeas[matchedmarkerpos]);

      if(it!=modinitmarkers.end()) {yshift = ymeas[matchedmarkerpos]-ynom[iNomEntry];} //if first marker of a module is shifted all the others will be shifted
    }
    else {
      xmeasfilled.push_back(-10000.);
      ymeasfilled.push_back(-10000.);
      zmeasfilled.push_back(-10000.);
    }
  }
}

//______________________________________________________________________________________________
int FindOppositeNominalMarkerPosition(double x, double y, vector<double> xvec, vector<double> yvec) {
  int posmarker=-1;
  for(unsigned int iMarker=0; iMarker<xvec.size(); iMarker++) {
    if(yvec[iMarker]==y && xvec[iMarker]!=x) {
      posmarker=iMarker;
      break;
    }
  }

  return posmarker;
}

//_____________________________________________________________________________________________
bool ReadDatFile(TString FileName, vector<double>& x, vector<double>& y, vector<double>& z)
{
  if(!FileName.Contains("txt") && !FileName.Contains("dat") && !FileName.Contains(".csv")) {
    cerr << "Wrong file format. Exit." << endl;
    return false;
  }

  ifstream inSet(FileName.Data());
  if(!inSet) {
    cerr<<"Please check if "<<FileName.Data() <<" is the right path. Exit."<<endl;
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
bool ReadDatFileMitutoyo(TString FileName, vector<double>& x, vector<double>& y, vector<double>& z)
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
void CreateDatFileMitutoyo(TString FileName, vector<double> x_after_rew_filled, vector<double> y_after_rew_filled, vector<double> z_after_rew_filled, vector<double> x_extrap, vector<double> y_extrap, vector<double> z_extrap) {

  ofstream outFile;
  outFile.open(FileName.Data());

  int iPoint = 0;
  int Num = 0;
  for(unsigned int iMarker=0; iMarker<x_after_rew_filled.size(); iMarker++) {
    iPoint++;
    if(iPoint==1)
      Num = 15;
    else if(iPoint>=2 && iPoint<=13)
      Num = 39;
    else if(iPoint==14)
      Num = 59;
    else if(iPoint==15)
      Num = 74;
    else if(iPoint>=16 && iPoint<=27)
      Num = 98;
    else if(iPoint==28)
      Num = 117;
      
    if(x_after_rew_filled[iMarker]!=-10000.) outFile << Form("%d;%d;point ;1;%.4f;%.4f;%.4f;;%.4f",iPoint,Num,x_after_rew_filled[iMarker],y_after_rew_filled[iMarker],z_after_rew_filled[iMarker],0.) << endl;
    else if(x_extrap[iMarker]!=-10000.) outFile << Form("%d;%d;point ;1;%.4f;%.4f;%.4f;;%.4f",iPoint,Num,x_extrap[iMarker],y_extrap[iMarker],z_extrap[iMarker],0.) << endl;
    if((iMarker+1)%28==0) {
      iPoint = 0;
      outFile << "100;121;Plane_Module Gauss;28;0.000;0.0000;0.0000;;0.0000" << endl;
      outFile << ";;;;0.0000;0.0000;1.0000;;" << endl;
      outFile << "100;122;Plane_Module;;;0.1000;0.0000;0.0000;" << endl;
      outFile << ";;PlanaritÖ;;;;;;     ****--" << endl;
    }
  }

  cout << "File " << FileName.Data() << " saved." << endl;
  outFile.close();
}

//______________________________________________________________________________________________
void SetStyle() {
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleSize(0.045,"xyz");
  gStyle->SetTitleSize(0.045,"t");
  gStyle->SetLabelSize(0.04,"xyzt");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetOptStat(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetHistLineWidth(2);
}
