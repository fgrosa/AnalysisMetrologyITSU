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

//************************************************************//
//                                                            //
//    Main Functions: ComputeResidualsToNominalPositions()    //
//                                                            //
//************************************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIBALES

//all values are in mm
const double ymodlims[8] = {-738.75,-527.75,-316.75,-105.55,105.55,316.75,527.75,738.75};
const bool fMitutoyoFile = true;

const TString FileNameMeas = "PatternRecoTest/T-OB-HS-L-006_ALC-0312-01_233_PLANARITY_ALLMODULES_2018_7_6_beforeUarms_AUT.dat";
const TString FileNameNom = "OL_marker_nominal_positions_HS.dat";

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
int ComputeResidualsToNominalPositions();
bool ReadFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
bool MatchMeasToNominal(double xmeas, double ymeas, double xnom, double ynom, double accdeltax, double accdeltay, double yshift);
void SetStyle();

//_____________________________________________________________________________________________
int ComputeResidualsToNominalPositions()
{
  SetStyle();

  std::vector<double> x_nom;
  std::vector<double> y_nom;
  std::vector<double> z_nom;
  std::vector<double> x_meas;
  std::vector<double> y_meas;
  std::vector<double> z_meas;

  bool readNom=ReadFile(FileNameNom,x_nom,y_nom,z_nom);
  if(!readNom) {
    cerr << "Impossibile to find the input file for nominal positions. Exit" << endl;
    return 1;
  }
  if(x_nom.size()<=0 || y_nom.size()<=0 || z_nom.size()<=0) {
    cerr << "The number of points is <=0. Exit" << endl;
    return 2;
  }
  if(y_nom.size()!=z_nom.size() || x_nom.size()!=z_nom.size()) {
    cerr << "There is a discrepacy between the number of entries for the different coordinates. Exit." << endl;
    return 3;
  }

  bool readMeas=false;
  if(!fMitutoyoFile) readMeas=ReadFile(FileNameMeas,x_meas,y_meas,z_meas);
  else readMeas=ReadDatFileMitutoyo(FileNameMeas,x_meas,y_meas,z_meas);
  
  if(!readMeas) {
    cerr << "Impossibile to find the input file for measured positions. Exit" << endl;
    return 1;
  }
  if(x_meas.size()<=0 || y_meas.size()<=0 || z_meas.size()<=0) {
    cerr << "The number of points is <=0. Exit" << endl;
    return 2;
  }
  if(y_meas.size()!=z_meas.size() || x_meas.size()!=z_meas.size()) {
    cerr << "There is a discrepacy between the number of entries for the different coordinates. Exit." << endl;
    return 3;
  }

  TGraph* gNomPos = new TGraph(y_nom.size());
  for(unsigned int iEntryNom=0; iEntryNom<y_nom.size(); iEntryNom++) {
    gNomPos->SetPoint(iEntryNom,x_nom[iEntryNom],y_nom[iEntryNom]);
  }
  TGraph* gMeasPos = new TGraph(y_meas.size());
  for(unsigned int iEntryMeas=0; iEntryMeas<y_meas.size(); iEntryMeas++) {
    gMeasPos->SetPoint(iEntryMeas,x_meas[iEntryMeas],y_meas[iEntryMeas]);
  }

  TBox* ModBox[7];
  for(int iMod=0; iMod<7; iMod++) {
    ModBox[iMod] = new TBox(-15.1,ymodlims[iMod],15.1,ymodlims[iMod+1]);
    ModBox[iMod]->SetLineColor(kRed);
    ModBox[iMod]->SetLineWidth(1);
    ModBox[iMod]->SetFillColor(kWhite);
    ModBox[iMod]->SetFillStyle(0);
  }

  TLegend* leg = new TLegend(0.25,0.83,0.75,0.92);
  leg->SetTextSize(0.045);
  leg->SetFillColor(kWhite);
  leg->AddEntry(gNomPos,"Nominal positions","p");
  leg->AddEntry(gMeasPos,"Measured positions","p");

  TCanvas* cPos = new TCanvas("cPos","",500,1000);
  gNomPos->SetTitle("");
  gMeasPos->SetTitle("");
  gNomPos->GetYaxis()->SetTitle("y (mm)");
  gNomPos->GetXaxis()->SetTitle("x (mm)");
  gNomPos->GetYaxis()->SetRangeUser(-800.,1100);
  gMeasPos->GetYaxis()->SetTitle("y (mm)");
  gMeasPos->GetXaxis()->SetTitle("x (mm)");
  gMeasPos->GetYaxis()->SetRangeUser(-800.,1100);
  gNomPos->SetMarkerStyle(kFullCircle);
  gNomPos->SetMarkerSize(1);
  gMeasPos->SetMarkerStyle(kFullDiamond);
  gMeasPos->SetMarkerColor(kBlue);
  gNomPos->SetMarkerSize(1);
  gNomPos->Draw("AP");
  gMeasPos->Draw("P");
  TLatex* lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.045);
  lat->SetTextColor(kRed);
  for(int iMod=0; iMod<7; iMod++) {
    ModBox[iMod]->Draw("same");
    lat->DrawLatex(-3,(ymodlims[iMod+1]+ymodlims[iMod])/2,Form("MOD%d",iMod+1));
  }
  leg->Draw("same");

  std::vector<double> x_res;
  std::vector<double> y_res;
  std::vector<double> xy_res;
  std::vector<int> nom_index;
  
  //positions of initial markers in vector
  std::vector<unsigned int>  modinitmarkers;
  modinitmarkers.push_back(0);
  modinitmarkers.push_back(28);
  modinitmarkers.push_back(56);
  modinitmarkers.push_back(84);
  modinitmarkers.push_back(112);
  modinitmarkers.push_back(140);
  modinitmarkers.push_back(168);
  
  double yshift=0.;
  std::vector<unsigned int>::iterator it;
  int countermeasmarkers=0;
  int matchedmarkerpos=-1;
  
  for(unsigned int iNomEntry=0; iNomEntry<x_nom.size(); iNomEntry++) {
    bool measfound=false;
    it = find(modinitmarkers.begin(), modinitmarkers.end(), iNomEntry);
    if(it!=modinitmarkers.end()) { // it found the index in the first-marker verctor
      yshift=0.; //reset y shift before starting each module
      matchedmarkerpos=countermeasmarkers;
      measfound=true;
      countermeasmarkers++;
    }
    else {
      for(unsigned int iMeasEntry=0; iMeasEntry<x_meas.size(); iMeasEntry++) {
        if(MatchMeasToNominal(x_meas[iMeasEntry],y_meas[iMeasEntry],x_nom[iNomEntry],y_nom[iNomEntry],0.800,0.350,yshift)) {
          measfound=true;
          countermeasmarkers++;
          matchedmarkerpos=iMeasEntry;
          break;
        }
      }
    }
    if(measfound) {
      x_res.push_back((x_meas[matchedmarkerpos]-x_nom[iNomEntry])*1000);
      y_res.push_back((y_meas[matchedmarkerpos]-y_nom[iNomEntry])*1000);
      double dist = TMath::Sqrt((x_meas[matchedmarkerpos]-x_nom[iNomEntry])*(x_meas[matchedmarkerpos]-x_nom[iNomEntry])+(y_meas[matchedmarkerpos]-y_nom[iNomEntry])*(y_meas[matchedmarkerpos]-y_nom[iNomEntry]));
      xy_res.push_back(dist*1000);
      nom_index.push_back(iNomEntry);

      if(it!=modinitmarkers.end()) {yshift = y_meas[matchedmarkerpos]-y_nom[iNomEntry];} //if first marker of a module is shifted all the others will be shifted
    }
  }

  TH1F* hResX = new TH1F("hResX","",100,-100,100);
  hResX->GetYaxis()->SetTitle("Entries");
  hResX->GetXaxis()->SetTitle("X_{meas}-X_{nom} (#mum)");
  TH1F* hResY = new TH1F("hResY","",100,-100,100);
  hResY->GetYaxis()->SetTitle("Entries");
  hResY->GetXaxis()->SetTitle("Y_{meas}-Y_{nom} (#mum)");
  TH1F* hResXY = new TH1F("hResXY","",100,0,200);
  hResXY->GetYaxis()->SetTitle("Entries");
  hResXY->GetXaxis()->SetTitle("dist_{XY} (#mum)");

  TGraph* gResXvsYnom_Xpos = new TGraph();
  TGraph* gResXvsYnom_Xneg = new TGraph();
  TGraph* gResXvsYnom_Xint = new TGraph();
  TGraph* gResYvsYnom_Xpos = new TGraph();
  TGraph* gResYvsYnom_Xneg = new TGraph();
  TGraph* gResYvsYnom_Xint = new TGraph();

  int iPointPos=0;
  int iPointNeg=0;
  for(unsigned int iRes=0; iRes<y_res.size(); iRes++) {
    hResX->Fill(x_res[iRes]);
    hResY->Fill(y_res[iRes]);
    hResXY->Fill(xy_res[iRes]);

    if(x_nom[nom_index[iRes]]>0) {
      gResXvsYnom_Xpos->SetPoint(iPointPos,y_nom[nom_index[iRes]],x_res[iRes]);
      gResYvsYnom_Xpos->SetPoint(iPointPos,y_nom[nom_index[iRes]],y_res[iRes]);
      iPointPos++;
    }
    else {
      gResXvsYnom_Xneg->SetPoint(iPointNeg,y_nom[nom_index[iRes]],x_res[iRes]);
      gResYvsYnom_Xneg->SetPoint(iPointNeg,y_nom[nom_index[iRes]],y_res[iRes]);
      iPointNeg++;
    }
    gResXvsYnom_Xint->SetPoint(iPointNeg,y_nom[nom_index[iRes]],x_res[iRes]);
    gResYvsYnom_Xint->SetPoint(iPointPos,y_nom[nom_index[iRes]],y_res[iRes]);

  }

  TPaveText* infox = new TPaveText(0.6,0.7,0.9,0.9,"NDC");
  infox->SetTextFont(42);
  infox->SetTextSize(0.045);
  infox->SetTextColor(kBlack);
  infox->SetBorderSize(1);
  infox->SetFillColor(kWhite);
  infox->AddText(Form("mean = %0.1f (#mum)",hResX->GetMean()));
  infox->AddText(Form("RMS = %0.1f (#mum)",hResX->GetRMS()));

  TPaveText* infoy = new TPaveText(0.6,0.7,0.9,0.9,"NDC");
  infoy->SetTextFont(42);
  infoy->SetTextSize(0.045);
  infoy->SetTextColor(kBlack);
  infoy->SetBorderSize(1);
  infoy->SetFillColor(kWhite);
  infoy->AddText(Form("mean = %0.1f (#mum)",hResY->GetMean()));
  infoy->AddText(Form("RMS = %0.1f (#mum)",hResY->GetRMS()));

  TCanvas* cResX = new TCanvas("cResX","",800,800);
  cResX->SetGrid();
  hResX->Draw();
  infox->Draw("same");
  TCanvas* cResY = new TCanvas("cResY","",800,800);
  cResY->SetGrid();
  hResY->Draw();
  infoy->Draw("same");
  TCanvas* cResXY = new TCanvas("cResXY","",800,800);
  cResXY->SetGrid();
  hResXY->Draw();

  TLegend* leg2 = new TLegend(0.6,0.75,0.89,0.87);
  leg2->SetTextSize(0.045);
  leg2->SetFillColor(kWhite);
  leg2->AddEntry(gResXvsYnom_Xpos,"x = 15 mm","p");
  leg2->AddEntry(gResXvsYnom_Xneg,"x = -15 mm","p");

  TCanvas* cResXvsYnom = new TCanvas("cResXvsYnom","",1200,900);
  cResXvsYnom->SetGrid();
  gResXvsYnom_Xpos->GetXaxis()->SetTitle("Y_{nom} (mm)");
  gResXvsYnom_Xpos->GetYaxis()->SetTitle("X_{meas}-X_{nom} (#mum)");
  gResXvsYnom_Xneg->GetXaxis()->SetTitle("Y_{nom} (mm)");
  gResXvsYnom_Xneg->GetYaxis()->SetTitle("X_{meas}-X_{nom} (#mum)");
  gResXvsYnom_Xneg->GetYaxis()->SetRangeUser(-150,150);
  gResXvsYnom_Xpos->GetYaxis()->SetRangeUser(-150,150);
  gResXvsYnom_Xpos->SetMarkerStyle(kFullCircle);
  gResXvsYnom_Xpos->SetMarkerSize(1);
  gResXvsYnom_Xpos->SetMarkerColor(kBlue+1);
  gResXvsYnom_Xneg->SetMarkerStyle(kFullSquare);
  gResXvsYnom_Xneg->SetMarkerSize(1);
  gResXvsYnom_Xneg->SetMarkerColor(kRed+1);
  gResXvsYnom_Xpos->Draw("AP");
  gResXvsYnom_Xneg->Draw("P");
  leg2->Draw("same");

  TCanvas* cResYvsYnom = new TCanvas("cResYvsYnom","",1200,900);
  cResYvsYnom->SetGrid();
  gResYvsYnom_Xpos->GetXaxis()->SetTitle("Y_{nom} (mm)");
  gResYvsYnom_Xpos->GetYaxis()->SetTitle("Y_{meas}-Y_{nom} (#mum)");
  gResYvsYnom_Xneg->GetXaxis()->SetTitle("Y_{nom} (mm)");
  gResYvsYnom_Xneg->GetYaxis()->SetTitle("Y_{meas}-Y_{nom} (#mum)");
  gResYvsYnom_Xneg->GetYaxis()->SetRangeUser(-150,150);
  gResYvsYnom_Xpos->GetYaxis()->SetRangeUser(-150,150);
  gResYvsYnom_Xpos->SetMarkerStyle(kFullCircle);
  gResYvsYnom_Xpos->SetMarkerSize(1);
  gResYvsYnom_Xpos->SetMarkerColor(kBlue+1);
  gResYvsYnom_Xneg->SetMarkerStyle(kFullSquare);
  gResYvsYnom_Xneg->SetMarkerSize(1);
  gResYvsYnom_Xneg->SetMarkerColor(kRed+1);
  gResYvsYnom_Xpos->Draw("AP");
  gResYvsYnom_Xneg->Draw("P");
  leg2->Draw("same");

  TString OutFileName = FileNameMeas;
  OutFileName.ReplaceAll(".txt","_MeasPos.pdf");
  OutFileName.ReplaceAll(".dat","_MeasPos.pdf");
  cPos->SaveAs(OutFileName.Data());

  OutFileName.ReplaceAll("_MeasPos.pdf","_ResX.pdf");
  cResX->SaveAs(OutFileName.Data());
  OutFileName.ReplaceAll("_ResX.pdf","_ResY.pdf");
  cResY->SaveAs(OutFileName.Data());
  OutFileName.ReplaceAll("_ResY.pdf","_DistXY.pdf");
  cResXY->SaveAs(OutFileName.Data());
  OutFileName.ReplaceAll("_DistXY.pdf","_ResX_vs_Y.pdf");
  cResXvsYnom->SaveAs(OutFileName.Data());
  OutFileName.ReplaceAll("_ResX_vs_Y.pdf","_ResY_vs_Y.pdf");
  cResYvsYnom->SaveAs(OutFileName.Data());

  OutFileName.ReplaceAll("_ResY_vs_Y.pdf","_Res.root");
  TFile outfile(OutFileName.Data(),"RECREATE");
  gResYvsYnom_Xint->SetName("gResYvsYnom_Xint");
  gResXvsYnom_Xint->SetName("gResXvsYnom_Xint");
  gResYvsYnom_Xint->Write();
  gResXvsYnom_Xint->Write();
  hResX->Write();
  hResY->Write();
  outfile.Close();
  
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
      islineok=true;
      std::vector<string> values;
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
bool MatchMeasToNominal(double xmeas, double ymeas, double xnom, double ynom, double accdeltax, double accdeltay, double yshift) {
  
  if(TMath::Abs(xmeas-xnom)<accdeltax && TMath::Abs((ymeas-yshift)-ynom)<accdeltay) return true;
  
  return false;
}

//______________________________________________________________________________________________
void SetStyle() {
  cout << "Setting drawing style!" << endl;
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetTitleFont(42);
  gStyle->SetLabelFont(42);
  gStyle->SetTitleOffset(1.,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetOptStat(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
}
