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
#include <TPaveText.h>
#include <TBox.h>
#include <RQ_OBJECT.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGTab.h>
#include <TGFileDialog.h>
#include <TRootEmbeddedCanvas.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>

#endif

//**********************************************************************//
//                                                                      //
//        Main Function: FinalMarkerPositionExtrapolationAndQA()        //
//                                                                      //
//**********************************************************************//

//author: Fabrizio Grosa, INFN Torino
//e-mail: fabrizio.grosa@to.infn.it

//_____________________________________________________________________________________________
//METHOD PROTOTYPES

void FinalMarkerPositionExtrapolationAndQA(const bool _fileFromCMM=true);
int MetrologyAndExtrapolation(TString infilename_HSleft, TString infilename_HSright, TString infilename_FinalStave, TString initdir, bool _fileFromCMM, TCanvas** fCanvas);
bool ReadDatFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset=0.);
bool ReadDatFileMitutoyo(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset=0.);
void CreateDatFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);
bool MatchMeasToNominal(double xmeas, double ymeas, double xnom, double ynom, double accdeltax, double accdeltay, double yshift);
void FillVectorsOfMatchedPositions(std::vector<double> xmeas, std::vector<double> ymeas, std::vector<double> zmeas, std::vector<double> xnom, std::vector<double> ynom, std::vector<double>& xmeasfilled, std::vector<double>& ymeasfilled, std::vector<double>& zmeasfilled, std::vector<double>& xres, std::vector<double>& yres);
int FindOppositeNominalMarkerPosition(double x, double y, std::vector<double> xvec, std::vector<double> yvec);
double ComputeResidualToPlane(double x, double y, double z, double *pars);
void ComputePlanarityRMSandMean(TGraph* graph, double &planarity, double &zRMS, double &zmean);
void DoInvisibleMarkerExtrapolation(int LeftOrRight, int &nextrap, std::vector<double> x_nominal, std::vector<double> y_nominal, std::vector<double> x_filled_HS, std::vector<double> y_filled_HS, std::vector<double> z_filled_HS, std::vector<double> x_filled_Stave, std::vector<double> y_filled_Stave, std::vector<double> z_filled_Stave, std::vector<double>& x_extrap, std::vector<double>& y_extrap, std::vector<double>& z_extrap, bool correctfortilt, double* parsHSplane, double* parsStaveplane);
bool PrintMatchedCoordinates(std::vector<double> xmeas, std::vector<double> ymeas, std::vector<double> xnom, std::vector<double> ynom, double shiftx=0.);
void SetStyle();

const char *fileformats[] = { "All files",     "*",
                              "ROOT files",    "*.root",
                              "ROOT macros",   "*.C",
                              "Text files",    "*.[tT][xX][tT]",
                              0,               0 };

class GroupBox : public TGGroupFrame {
  private:
    TGTextButton *fButton;
    TGTextEntry  *fEntry;

  public:
    GroupBox(const TGWindow *p, const char *name);
    TGTextEntry*  GetEntry()  { return fEntry; }
    TGTextButton* GetButton() { return fButton; }

    ClassDef(GroupBox, 1);
};

//______________________________________________________________________________
GroupBox::GroupBox(const TGWindow *p, const char *name) :
   TGGroupFrame(p, name)
{
   // Group frame containing text entry and text button.
   fEntry = new TGTextEntry(this);
   fEntry->SetDefaultSize(fEntry->GetWidth()*8, fEntry->GetHeight());
   AddFrame(fEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   TGHorizontalFrame *horz = new TGHorizontalFrame(this);
   AddFrame(horz, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

   fButton = new TGTextButton(horz, "&Open");
   horz->AddFrame(fButton, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 3,3,3,3));
}

//______________________________________________________________________________________________
//GUI IMPLEMENTATION
class MainFrame {
  RQ_OBJECT("MainFrame")

 private:
  TGMainFrame*          fMain;
  TGTab*                fTab;
  GroupBox*             fGboxHSL;
  GroupBox*             fGboxHSR;
  GroupBox*             fGboxStave;
  TString               fFileNameHSL;
  TString               fFileNameHSR;
  TString               fFileNameStave;
  TString               fInitDir;
  TRootEmbeddedCanvas*  fEcanvas;
  bool                  fIsFileFromCMM;

 public:
   enum HStype{kHSL, kHSR, kStave};

  MainFrame(const TGWindow* p, unsigned int w, unsigned int h, const bool _isFromCMM);
  virtual ~MainFrame();
  TString GetFileName(int ftype) {
    if(ftype==kHSL)        return fFileNameHSL;
    else if(ftype==kHSR)   return fFileNameHSR;
    else if(ftype==kStave) return fFileNameStave;
    else return ".";
  }
  void DoOpenHSLFile();
  void DoOpenHSRFile();
  void DoOpenStaveFile();
  int DoMetrologyAndExtrapolation();
};

//_______________________________________________________________________________
MainFrame::MainFrame(const TGWindow* p, unsigned int w, unsigned int h, const bool _isFromCMM):
  fIsFileFromCMM(_isFromCMM),
  fFileNameHSL("."),
  fFileNameHSR("."),
  fFileNameStave(".")
{
  fInitDir = gSystem->WorkingDirectory();

  fMain = new TGMainFrame(p,w,h);
  fMain->DontCallClose(); // to avoid double deletions.
  // use hierarchical cleaning
  fMain->SetCleanup(kDeepCleanup);
  //
  TGLayoutHints* fL1=new TGLayoutHints(kLHintsTop|kLHintsLeft,
                                       2, 2, 2, 2);
  TGLayoutHints* fL2=new TGLayoutHints(kLHintsTop|kLHintsRight,
                                       2, 2, 2, 2);
  TGLayoutHints* fL3=new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX,
                                       2, 2, 2, 2);
  TGLayoutHints* fL4=new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX,
                                       5, 5, 5, 5);
  TGLayoutHints* fL5=new TGLayoutHints(kLHintsTop|kLHintsRight,
                                       2, 2, 5, 1);
  TGLayoutHints* fL6=new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,
                                       5, 5, 5, 5);
  // Create a horiz frame widget button and combobox
  TGVerticalFrame* hFrameTop= new TGVerticalFrame(fMain,40,80);
  fGboxHSL = new GroupBox(fMain, "HS-L");
  hFrameTop->AddFrame(fGboxHSL, fL1);
  fGboxHSL->GetButton()->Connect("Clicked()","MainFrame",this,"DoOpenHSLFile()");
  fGboxHSR = new GroupBox(fMain, "HS-R");
  hFrameTop->AddFrame(fGboxHSR, fL1);
  fGboxHSR->GetButton()->Connect("Clicked()","MainFrame",this,"DoOpenHSRFile()");
  fGboxStave = new GroupBox(fMain, "HS-Stave");
  hFrameTop->AddFrame(fGboxStave, fL1);
  fGboxStave->GetButton()->Connect("Clicked()","MainFrame",this,"DoOpenStaveFile()");

  TGHorizontalFrame* hFrameBot = new TGHorizontalFrame(hFrameTop);
  hFrameTop->AddFrame(hFrameBot, fL5);

  TGTextButton* dometrology=new TGTextButton(hFrameBot,"&DoMetrology");
  dometrology->Connect("Clicked()","MainFrame",this,"DoMetrologyAndExtrapolation()");
  hFrameBot->AddFrame(dometrology, fL1);

  TGTextButton* exit = new TGTextButton(hFrameBot,"&Exit","gApplication->Terminate(0)");
  hFrameBot->AddFrame(exit, fL1);


  fMain->AddFrame(hFrameTop, fL4);

  //--------- create Tab widget and some composite frames
  fTab = new TGTab(fMain, 300, 300);
  const int kNtab=11;
  TRootEmbeddedCanvas* fEcanvas[kNtab];
  TString tabnames[kNtab] = {"HS","HSL planarity","HSR planarity","HSL residuals","HSR residuals","Stave planarity (HSL)","Stave planarity (HSR)","Stave residuals (HSL)","Stave residuals (HSR)","Stave extrapolation","Stave extrapolation goodness"};
  for(Int_t i=0; i<kNtab; ++i) {
    TGCompositeFrame* tf = fTab->AddTab(tabnames[i].Data());
   }
  fMain->AddFrame(fTab, fL6);

  // Set a name to the main frame
  fMain->SetWindowName("Analysis of metrological surveys and extrapolation of invisible positions");

  // Map all subwindows of main frame
  fMain->MapSubwindows();

  // Initialize the layout algorithm
  fMain->Resize(fMain->GetDefaultSize());

  // Map main frame
  fMain->MapWindow();
}

//_______________________________________________________________________________
MainFrame::~MainFrame()
{
  // Clean up used widgets: frames, buttons, layout hints
  fMain->Cleanup();
  delete fMain;
}

//_______________________________________________________________________________
void MainFrame::DoOpenHSLFile()
{
  TGFileInfo fi;
  fi.fFileTypes = fileformats;
  fi.fIniDir = StrDup(fFileNameHSL);
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  fFileNameHSL=fi.fFilename;
  TGTextEntry* entry = fGboxHSL->GetEntry();
  entry->SetText(fFileNameHSL.Data());
  std::cout << "\n\n\nOpen HSL file: "<<  fFileNameHSL << std::endl;

  return;
}

//_______________________________________________________________________________
void MainFrame::DoOpenHSRFile()
{
  TGFileInfo fi;
  fi.fFileTypes = fileformats;
  fi.fIniDir = StrDup(fFileNameHSR);
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  fFileNameHSR=fi.fFilename;
  TGTextEntry* entry = fGboxHSR->GetEntry();
  entry->SetText(fFileNameHSR.Data());
  std::cout << "Open HSR file: "<<  fFileNameHSR << std::endl;

  return;
}

//_______________________________________________________________________________
void MainFrame::DoOpenStaveFile()
{
  TGFileInfo fi;
  fi.fFileTypes = fileformats;
  fi.fIniDir = StrDup(fFileNameStave);
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  fFileNameStave=fi.fFilename;
  TGTextEntry* entry = fGboxStave->GetEntry();
  entry->SetText(fFileNameStave.Data());
  std::cout << "Open Stave file: "<<  fFileNameStave << "\n\n" <<std::endl;

  return;
}

//_______________________________________________________________________________
int MainFrame::DoMetrologyAndExtrapolation()
{
  SetStyle();

  TGLayoutHints* fL3=new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY,
                                       1, 1, 1, 1);
  const int kNtab=11;
  TRootEmbeddedCanvas* fEcanvas[kNtab];
  TCanvas* fCanvas[kNtab];
  for(Int_t i=0; i<kNtab; ++i) {
    TGCompositeFrame* tb=fTab->GetTabContainer(i);
    tb->RemoveAll();
    fEcanvas[i] = new TRootEmbeddedCanvas(Form("efc%d",i), tb, 1200, 800);
    tb->AddFrame(fEcanvas[i], fL3);
    fEcanvas[i]->GetCanvas()->SetBorderMode(0);
    fCanvas[i]=fEcanvas[i]->GetCanvas();
  }

  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();

  return MetrologyAndExtrapolation(fFileNameHSL,fFileNameHSR,fFileNameStave,fInitDir,fIsFileFromCMM,fCanvas);
}

//_______________________________________________________________________________
//METHOD IMPLEMENTATION
void FinalMarkerPositionExtrapolationAndQA(const bool _fileFromCMM)
{
  new MainFrame(gClient->GetRoot(), 200, 200, _fileFromCMM);
}

//_____________________________________________________________________________________________
int MetrologyAndExtrapolation(TString infilename_HSleft, TString infilename_HSright, TString infilename_FinalStave, TString initdir, bool _fileFromCMM, TCanvas** fCanvas) {

  //load measured marker positions HS_left BEFORE U-ARMS
  std::vector<double> x_HSleft;
  std::vector<double> y_HSleft;
  std::vector<double> z_HSleft;
  bool loadfile = _fileFromCMM ? ReadDatFileMitutoyo(infilename_HSleft,x_HSleft,y_HSleft,z_HSleft,-12.9) : ReadDatFile(infilename_HSleft,x_HSleft,y_HSleft,z_HSleft,-12.9);
  if(!loadfile) {return 1;}

  //load measured marker positions HS_right BEFORE U-ARMS
  std::vector<double> x_HSright;
  std::vector<double> y_HSright;
  std::vector<double> z_HSright;
  loadfile = _fileFromCMM ? ReadDatFileMitutoyo(infilename_HSright,x_HSright,y_HSright,z_HSright,12.9) : ReadDatFile(infilename_HSright,x_HSright,y_HSright,z_HSright,12.9);
  if(!loadfile) {return 2;}

  //load measured marker positions Final Stave Metrology
  std::vector<double> x_FinalStave;
  std::vector<double> y_FinalStave;
  std::vector<double> z_FinalStave;
  loadfile = _fileFromCMM ? ReadDatFileMitutoyo(infilename_FinalStave,x_FinalStave,y_FinalStave,z_FinalStave) : ReadDatFile(infilename_FinalStave,x_FinalStave,y_FinalStave,z_FinalStave);
  if(!loadfile) {return 3;}

  //divide final positions in HSleft and HSright
  std::vector<double> x_FinalStave_HSleft;
  std::vector<double> y_FinalStave_HSleft;
  std::vector<double> z_FinalStave_HSleft;
  std::vector<double> x_FinalStave_HSright;
  std::vector<double> y_FinalStave_HSright;
  std::vector<double> z_FinalStave_HSright;
  for(unsigned int iEntry=0; iEntry<x_FinalStave.size(); iEntry++) {
    if((x_FinalStave[iEntry]>-32 && x_FinalStave[iEntry]<-24) || (x_FinalStave[iEntry]>0 && x_FinalStave[iEntry]<6)) {
      x_FinalStave_HSleft.push_back(x_FinalStave[iEntry]);
      y_FinalStave_HSleft.push_back(y_FinalStave[iEntry]);
      z_FinalStave_HSleft.push_back(z_FinalStave[iEntry]);
    }
    else {
      x_FinalStave_HSright.push_back(x_FinalStave[iEntry]);
      y_FinalStave_HSright.push_back(y_FinalStave[iEntry]);
      z_FinalStave_HSright.push_back(z_FinalStave[iEntry]);
    }
  }

  //load nominal marker positions
  const TString nominalposition_HS = Form("%s/OL_marker_nominal_positions_HS.dat",initdir.Data());

  //HSL shifted of -12.9 cm in x, as it is in the satve
  std::vector<double> x_HSleft_nominal;
  std::vector<double> y_HSleft_nominal;
  std::vector<double> z_HSleft_nominal;
  loadfile=ReadDatFile(nominalposition_HS,x_HSleft_nominal,y_HSleft_nominal,z_HSleft_nominal,-12.9);
  if(!loadfile) {return 4;}

  //HSR shifted of +12.9 cm in x, as it is in the satve
  std::vector<double> x_HSright_nominal;
  std::vector<double> y_HSright_nominal;
  std::vector<double> z_HSright_nominal;
  loadfile=ReadDatFile(nominalposition_HS,x_HSright_nominal,y_HSright_nominal,z_HSright_nominal,12.9);
  if(!loadfile) {return 5;}

  //create vectors with nominal number of markers; fill with -10000 if marker not measured
  //compute nominal if != -10000
  //HS_left
  std::vector<double> x_HSleft_filled;
  std::vector<double> y_HSleft_filled;
  std::vector<double> z_HSleft_filled;

  std::vector<double> x_HSleft_restonominal;
  std::vector<double> y_HSleft_restonominal;

  FillVectorsOfMatchedPositions(x_HSleft,y_HSleft,z_HSleft,x_HSleft_nominal,y_HSleft_nominal,x_HSleft_filled,y_HSleft_filled,z_HSleft_filled,x_HSleft_restonominal,y_HSleft_restonominal);

  std::cout << "\n\n\n\n**************************************** HSL first metrology on bss ******************************************" << std::endl;
  bool okmatch = PrintMatchedCoordinates(x_HSleft_filled,y_HSleft_filled,x_HSleft_nominal,y_HSleft_nominal,12.9);
  if(!okmatch) {return -1;}
  std::cout << "\nNumber of measured markers: " << x_HSleft.size() << ", Number of not measured markers: " << x_HSleft_filled.size() - x_HSleft.size() << std::endl;
  std::cout << "\n***************************************************************************************************************\n\n\n" << std::endl;

  //HS_right
  std::vector<double> x_HSright_filled;
  std::vector<double> y_HSright_filled;
  std::vector<double> z_HSright_filled;

  std::vector<double> x_HSright_restonominal;
  std::vector<double> y_HSright_restonominal;

  FillVectorsOfMatchedPositions(x_HSright,y_HSright,z_HSright,x_HSright_nominal,y_HSright_nominal,x_HSright_filled,y_HSright_filled,z_HSright_filled,x_HSright_restonominal,y_HSright_restonominal);

  std::cout << "\n\n\n\n**************************************** HSR first metrology on bss ******************************************" << std::endl;
  okmatch = PrintMatchedCoordinates(x_HSright_filled,y_HSright_filled,x_HSright_nominal,y_HSright_nominal,-12.9);
  if(!okmatch) {return -1;}
  std::cout << "\nNumber of measured markers: " << x_HSright.size() << ", Number of not measured markers: " << x_HSright_filled.size() - x_HSright.size() << std::endl;
  std::cout << "\n***************************************************************************************************************\n\n\n" << std::endl;

  //FinalStave_HS_left
  std::vector<double> x_FinalStave_HSleft_filled;
  std::vector<double> y_FinalStave_HSleft_filled;
  std::vector<double> z_FinalStave_HSleft_filled;

  std::vector<double> x_FinalStave_HSleft_restonominal;
  std::vector<double> y_FinalStave_HSleft_restonominal;

  FillVectorsOfMatchedPositions(x_FinalStave_HSleft,y_FinalStave_HSleft,z_FinalStave_HSleft,x_HSleft_nominal,y_HSleft_nominal,x_FinalStave_HSleft_filled,y_FinalStave_HSleft_filled,z_FinalStave_HSleft_filled,x_FinalStave_HSleft_restonominal,y_FinalStave_HSleft_restonominal);

  std::cout << "\n\n\n\n**************************************** HSL final metrology on Meas ******************************************" << std::endl;
  okmatch = PrintMatchedCoordinates(x_FinalStave_HSleft_filled,y_FinalStave_HSleft_filled,x_HSleft_nominal,y_HSleft_nominal,0.);
  if(!okmatch) {return -1;}
  std::cout << "\nNumber of measured markers: " << x_FinalStave_HSleft.size() << ", Number of not measured markers: " << x_FinalStave_HSleft_filled.size() - x_FinalStave_HSleft.size() << std::endl;
  std::cout << "\n***************************************************************************************************************\n\n\n" << std::endl;

  //FinalStave_HS_right
  std::vector<double> x_FinalStave_HSright_filled;
  std::vector<double> y_FinalStave_HSright_filled;
  std::vector<double> z_FinalStave_HSright_filled;

  std::vector<double> x_FinalStave_HSright_restonominal;
  std::vector<double> y_FinalStave_HSright_restonominal;

  FillVectorsOfMatchedPositions(x_FinalStave_HSright,y_FinalStave_HSright,z_FinalStave_HSright,x_HSright_nominal,y_HSright_nominal,x_FinalStave_HSright_filled,y_FinalStave_HSright_filled,z_FinalStave_HSright_filled,x_FinalStave_HSright_restonominal,y_FinalStave_HSright_restonominal);

  std::cout << "\n\n\n\n**************************************** HSR final metrology on Meas ******************************************" << std::endl;
  okmatch = PrintMatchedCoordinates(x_FinalStave_HSright_filled,y_FinalStave_HSright_filled,x_HSright_nominal,y_HSright_nominal,0.);
  if(!okmatch) {return -1;}
  std::cout << "\nNumber of measured markers: " << x_FinalStave_HSright.size() << ", Number of not measured markers: " << x_FinalStave_HSright_filled.size() - x_FinalStave_HSright.size() << std::endl;
  std::cout << "\n***************************************************************************************************************\n\n\n" << std::endl;

  //Fill graphs / histograms
  TF2* fPlaneHSleft = new TF2("fPlaneHSleft","[0]+[1]*x+[2]*y",-15,15,-1000,1000);
  TF2* fPlaneHSright = new TF2("fPlaneHSright","[0]+[1]*x+[2]*y",-15,15,-1000,1000);
  TF2* fPlaneStaveHSleft = new TF2("fPlaneStaveHSleft","[0]+[1]*x+[2]*y",-30,5,-1000,1000);
  TF2* fPlaneStaveHSright = new TF2("fPlaneStaveHSright","[0]+[1]*x+[2]*y",-5,30,-1000,1000);
  fPlaneStaveHSright->FixParameter(1,0); //only one set of points for one x -> no degree of freedom in x

  //HS left
  TGraph2D* gHSleftMeas = new TGraph2D(x_HSleft.size());
  gHSleftMeas->SetName("gHSleftMeas");
  gHSleftMeas->SetMarkerStyle(kFullCircle);
  gHSleftMeas->SetMarkerColor(kRed);
  for(unsigned int iEntry=0; iEntry<x_HSleft.size(); iEntry++) {
    gHSleftMeas->SetPoint(iEntry,x_HSleft[iEntry]+12.9,y_HSleft[iEntry],z_HSleft[iEntry]);
  }
  gHSleftMeas->Fit("fPlaneHSleft");

  //HS right
  TGraph2D* gHSrightMeas = new TGraph2D(x_HSright.size());
  gHSrightMeas->SetName("gHSrightMeas");
  gHSrightMeas->SetMarkerStyle(kFullSquare);
  gHSrightMeas->SetMarkerColor(kBlue);
  for(unsigned int iEntry=0; iEntry<x_HSright.size(); iEntry++) {
    gHSrightMeas->SetPoint(iEntry,x_HSright[iEntry]-12.9,y_HSright[iEntry],z_HSright[iEntry]);
  }
  gHSrightMeas->Fit("fPlaneHSright");

  //Final stave (all markers)
  TGraph2D* gStaveMeas = new TGraph2D(x_FinalStave.size());
  gStaveMeas->SetName("gStaveMeas");
  gStaveMeas->SetMarkerStyle(kFullCircle);
  gStaveMeas->SetMarkerColor(kBlack);
  for(unsigned int iEntry=0; iEntry<x_FinalStave.size(); iEntry++) {
    gStaveMeas->SetPoint(iEntry,x_FinalStave[iEntry],y_FinalStave[iEntry],z_FinalStave[iEntry]);
  }

  //HS left (final stave)
  TGraph2D* gStaveMeasHSleft = new TGraph2D(x_FinalStave_HSleft.size());
  gStaveMeasHSleft->SetName("gStaveMeasHSleft");
  gStaveMeasHSleft->SetMarkerStyle(kFullCircle);
  gStaveMeasHSleft->SetMarkerColor(kRed);
  for(unsigned int iEntry=0; iEntry<x_FinalStave_HSleft.size(); iEntry++) {
    gStaveMeasHSleft->SetPoint(iEntry,x_FinalStave_HSleft[iEntry],y_FinalStave_HSleft[iEntry],z_FinalStave_HSleft[iEntry]);
  }
  gStaveMeasHSleft->Fit("fPlaneStaveHSleft");

  //HS right (final stave)
  TGraph2D* gStaveHSrightMeas = new TGraph2D(x_FinalStave_HSright.size());
  gStaveHSrightMeas->SetName("gStaveHSrightMeas");
  gStaveHSrightMeas->SetMarkerStyle(kFullSquare);
  gStaveHSrightMeas->SetMarkerColor(kBlue);
  for(unsigned int iEntry=0; iEntry<x_FinalStave_HSright.size(); iEntry++) {
    gStaveHSrightMeas->SetPoint(iEntry,x_FinalStave_HSright[iEntry],y_FinalStave_HSright[iEntry],z_FinalStave_HSright[iEntry]);
  }
  gStaveHSrightMeas->Fit("fPlaneStaveHSright");

  //compute X and Y residuals w.r.t. nominal positions
  //HS left
  TH1F* hResToNominalHSleft_X = new TH1F("hResToNominalHSleft_X","Residuals to X nominal positions on HSL;x_{meas}-x_{nom} (#mum); Entries",500,-1000,1000);
  TH1F* hResToNominalHSleft_Y = new TH1F("hResToNominalHSleft_Y","Residuals to Y nominal positions on HSL;y_{meas}-y_{nom} (#mum); Entries",500,-1000,1000);
  for(unsigned int iEntry=0; iEntry<x_HSleft_restonominal.size(); iEntry++) {
    hResToNominalHSleft_X->Fill(x_HSleft_restonominal[iEntry]);
    hResToNominalHSleft_Y->Fill(y_HSleft_restonominal[iEntry]);
  }

  //HS right
  TH1F* hResToNominalHSright_X = new TH1F("hResToNominalHSright_X","Residuals to X nominal positions on HSR;x_{meas}-x_{nom} (#mum); Entries",500,-1000,1000);
  TH1F* hResToNominalHSright_Y = new TH1F("hResToNominalHSright_Y","Residuals to Y nominal positions on HSR;y_{meas}-y_{nom} (#mum); Entries",500,-1000,1000);
  for(unsigned int iEntry=0; iEntry<x_HSright_restonominal.size(); iEntry++) {
    hResToNominalHSright_X->Fill(x_HSright_restonominal[iEntry]);
    hResToNominalHSright_Y->Fill(y_HSright_restonominal[iEntry]);
  }

  //Final Stave - HS left
  TH1F* hResToNominalStave_HSleft_X = new TH1F("hResToNominalStave_HSleft_X","Final Stave - Residuals to X nominal positions on HSL;x_{meas}-x_{nom} (#mum); Entries",500,-1000,1000);
  TH1F* hResToNominalStave_HSleft_Y = new TH1F("hResToNominalStave_HSleft_Y","Final Stave - Residuals to Y nominal positions on HSL;y_{meas}-y_{nom} (#mum); Entries",500,-1000,1000);
  for(unsigned int iEntry=0; iEntry<x_FinalStave_HSleft_restonominal.size(); iEntry++) {
    hResToNominalStave_HSleft_X->Fill(x_FinalStave_HSleft_restonominal[iEntry]);
    hResToNominalStave_HSleft_Y->Fill(y_FinalStave_HSleft_restonominal[iEntry]);
  }

  //Final Stave - HS right
  TH1F* hResToNominalStave_HSright_X = new TH1F("hResToNominalStave_HSright_X","Final Stave - Residuals to X nominal positions on HSR;x_{meas}-x_{nom} (#mum); Entries",500,-1000,1000);
  TH1F* hResToNominalStave_HSright_Y = new TH1F("hResToNominalStave_HSright_Y","Final Stave - Residuals to Y nominal positions on HSR;y_{meas}-y_{nom} (#mum); Entries",500,-1000,1000);
  for(unsigned int iEntry=0; iEntry<x_FinalStave_HSright_restonominal.size(); iEntry++) {
    hResToNominalStave_HSright_X->Fill(x_FinalStave_HSright_restonominal[iEntry]);
    hResToNominalStave_HSright_Y->Fill(y_FinalStave_HSright_restonominal[iEntry]);
  }

  //compute Z residuals w.r.t. average and nominal planes and planarity
  //HS left
  TGraph* gResToNomPlaneHSleft_Z = new TGraph(x_HSleft.size());
  gResToNomPlaneHSleft_Z->SetName("gResToNomPlaneHSleft_Z");
  gResToNomPlaneHSleft_Z->SetTitle("HSL planarity (bss metrology)");
  gResToNomPlaneHSleft_Z->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneHSleft_Z->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneHSleft_Z->SetMarkerStyle(kFullCircle);
  gResToNomPlaneHSleft_Z->SetMarkerSize(1.);
  gResToNomPlaneHSleft_Z->SetMarkerColor(kBlack);
  TGraph* gResToNomPlaneHSleft_Z_posX = new TGraph(0);
  gResToNomPlaneHSleft_Z_posX->SetName("gResToNomPlaneHSleft_Z_posX");
  gResToNomPlaneHSleft_Z_posX->SetTitle("HSL planarity (bss metrology) - pos X");
  gResToNomPlaneHSleft_Z_posX->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneHSleft_Z_posX->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneHSleft_Z_posX->SetMarkerStyle(kFullCircle);
  gResToNomPlaneHSleft_Z_posX->SetMarkerSize(1.);
  gResToNomPlaneHSleft_Z_posX->SetMarkerColor(kRed);
  TGraph* gResToNomPlaneHSleft_Z_negX = new TGraph(0);
  gResToNomPlaneHSleft_Z_negX->SetName("gResToNomPlaneHSleft_Z_negX");
  gResToNomPlaneHSleft_Z_negX->SetTitle("HSL planarity (bss metrology) - neg X");
  gResToNomPlaneHSleft_Z_negX->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneHSleft_Z_negX->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneHSleft_Z_negX->SetMarkerStyle(kFullCircle);
  gResToNomPlaneHSleft_Z_negX->SetMarkerSize(1.);
  gResToNomPlaneHSleft_Z_negX->SetMarkerColor(kBlue);
  TGraph* gResToAvPlaneHSleft_Z = new TGraph(x_HSleft.size());
  gResToAvPlaneHSleft_Z->SetName("gResToAvPlaneHSleft_Z");
  gResToAvPlaneHSleft_Z->SetTitle("HSL planarity (bss metrology)");
  gResToAvPlaneHSleft_Z->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneHSleft_Z->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneHSleft_Z->SetMarkerStyle(kFullCircle);
  gResToAvPlaneHSleft_Z->SetMarkerSize(1.);
  gResToAvPlaneHSleft_Z->SetMarkerColor(kBlack);
  TGraph* gResToAvPlaneHSleft_Z_posX = new TGraph(0);
  gResToAvPlaneHSleft_Z_posX->SetName("gResToAvPlaneHSleft_Z_posX");
  gResToAvPlaneHSleft_Z_posX->SetTitle("HSL planarity (bss metrology) - pos X");
  gResToAvPlaneHSleft_Z_posX->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneHSleft_Z_posX->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneHSleft_Z_posX->SetMarkerStyle(kFullCircle);
  gResToAvPlaneHSleft_Z_posX->SetMarkerSize(1.);
  gResToAvPlaneHSleft_Z_posX->SetMarkerColor(kRed);
  TGraph* gResToAvPlaneHSleft_Z_negX = new TGraph(0);
  gResToAvPlaneHSleft_Z_negX->SetName("gResToAvPlaneHSleft_Z_negX");
  gResToAvPlaneHSleft_Z_negX->SetTitle("HSL planarity (bss metrology) - neg X");
  gResToAvPlaneHSleft_Z_negX->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneHSleft_Z_negX->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneHSleft_Z_negX->SetMarkerStyle(kFullCircle);
  gResToAvPlaneHSleft_Z_negX->SetMarkerSize(1.);
  gResToAvPlaneHSleft_Z_negX->SetMarkerColor(kBlue);

  double planepars_base[3] = {0.41,0.,0.};
  int xposcounter=0;
  int xnegcounter=0;
  double deltazmax_HSL=0;
  for(unsigned int iEntry=0; iEntry<x_HSleft.size(); iEntry++) {
    double zrestonominal = ComputeResidualToPlane(x_HSleft[iEntry]+12.9,y_HSleft[iEntry],z_HSleft[iEntry],planepars_base);
    double zrestoaverage = ComputeResidualToPlane(x_HSleft[iEntry]+12.9,y_HSleft[iEntry],z_HSleft[iEntry],fPlaneHSleft->GetParameters());
    if(TMath::Abs(zrestonominal)>deltazmax_HSL) deltazmax_HSL=TMath::Abs(zrestonominal);
    gResToNomPlaneHSleft_Z->SetPoint(iEntry,y_HSleft[iEntry],zrestonominal*1000);
    gResToAvPlaneHSleft_Z->SetPoint(iEntry,y_HSleft[iEntry],zrestoaverage*1000);
    if(x_HSleft[iEntry]+12.9>0) {
      gResToNomPlaneHSleft_Z_posX->SetPoint(xposcounter,y_HSleft[iEntry],zrestonominal*1000);
      gResToAvPlaneHSleft_Z_posX->SetPoint(xposcounter,y_HSleft[iEntry],zrestoaverage*1000);
      xposcounter++;
    }
    else {
      gResToNomPlaneHSleft_Z_negX->SetPoint(xnegcounter,y_HSleft[iEntry],zrestonominal*1000);
      gResToAvPlaneHSleft_Z_negX->SetPoint(xnegcounter,y_HSleft[iEntry],zrestoaverage*1000);
      xnegcounter++;
    }
  }
  TF1* fAvPlaneHSleft_posX = new TF1("fAvPlaneHSleft_posX","pol1",-800,800);
  fAvPlaneHSleft_posX->SetParameters((fPlaneHSleft->GetParameter(0)+fPlaneHSleft->GetParameter(1)*(14.999)-planepars_base[0])*1000,fPlaneHSleft->GetParameter(2)*1000);
  fAvPlaneHSleft_posX->SetLineColor(kRed);
  TF1* fAvPlaneHSleft_negX = new TF1("fAvPlaneHSleft_negX","pol1",-800,800);
  fAvPlaneHSleft_negX->SetParameters((fPlaneHSleft->GetParameter(0)+fPlaneHSleft->GetParameter(1)*(-14.999)-planepars_base[0])*1000,fPlaneHSleft->GetParameter(2)*1000);
  fAvPlaneHSleft_negX->SetLineColor(kBlue);

  //HS right
  TGraph* gResToNomPlaneHSright_Z = new TGraph(x_HSright.size());
  gResToNomPlaneHSright_Z->SetName("gResToNomPlaneHSright_Z");
  gResToNomPlaneHSright_Z->SetTitle("HSR planarity (bss metrology)");
  gResToNomPlaneHSright_Z->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneHSright_Z->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneHSright_Z->SetMarkerStyle(kFullCircle);
  gResToNomPlaneHSright_Z->SetMarkerSize(1.);
  gResToNomPlaneHSright_Z->SetMarkerColor(kBlack);
  TGraph* gResToNomPlaneHSright_Z_posX = new TGraph(0);
  gResToNomPlaneHSright_Z_posX->SetName("gResToNomPlaneHSright_Z_posX");
  gResToNomPlaneHSright_Z_posX->SetTitle("HSR planarity (bss metrology) - pos X");
  gResToNomPlaneHSright_Z_posX->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneHSright_Z_posX->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneHSright_Z_posX->SetMarkerStyle(kFullCircle);
  gResToNomPlaneHSright_Z_posX->SetMarkerSize(1.);
  gResToNomPlaneHSright_Z_posX->SetMarkerColor(kRed);
  TGraph* gResToNomPlaneHSright_Z_negX = new TGraph(0);
  gResToNomPlaneHSright_Z_negX->SetName("gResToNomPlaneHSright_Z_negX");
  gResToNomPlaneHSright_Z_negX->SetTitle("HSR planarity (bss metrology) - neg X");
  gResToNomPlaneHSright_Z_negX->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneHSright_Z_negX->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneHSright_Z_negX->SetMarkerStyle(kFullCircle);
  gResToNomPlaneHSright_Z_negX->SetMarkerSize(1.);
  gResToNomPlaneHSright_Z_negX->SetMarkerColor(kBlue);
  TGraph* gResToAvPlaneHSright_Z = new TGraph(x_HSright.size());
  gResToAvPlaneHSright_Z->SetName("gResToAvPlaneHSright_Z");
  gResToAvPlaneHSright_Z->SetTitle("HSR planarity (bss metrology)");
  gResToAvPlaneHSright_Z->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneHSright_Z->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneHSright_Z->SetMarkerStyle(kFullCircle);
  gResToAvPlaneHSright_Z->SetMarkerSize(1.);
  gResToAvPlaneHSright_Z->SetMarkerColor(kBlack);
  TGraph* gResToAvPlaneHSright_Z_posX = new TGraph(0);
  gResToAvPlaneHSright_Z_posX->SetName("gResToAvPlaneHSright_Z_posX");
  gResToAvPlaneHSright_Z_posX->SetTitle("HSR planarity (bss metrology) - pos X");
  gResToAvPlaneHSright_Z_posX->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneHSright_Z_posX->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneHSright_Z_posX->SetMarkerStyle(kFullCircle);
  gResToAvPlaneHSright_Z_posX->SetMarkerSize(1.);
  gResToAvPlaneHSright_Z_posX->SetMarkerColor(kRed);
  TGraph* gResToAvPlaneHSright_Z_negX = new TGraph(0);
  gResToAvPlaneHSright_Z_negX->SetName("gResToAvPlaneHSright_Z_negX");
  gResToAvPlaneHSright_Z_negX->SetTitle("HSR planarity (bss metrology) - neg X");
  gResToAvPlaneHSright_Z_negX->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneHSright_Z_negX->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneHSright_Z_negX->SetMarkerStyle(kFullCircle);
  gResToAvPlaneHSright_Z_negX->SetMarkerSize(1.);
  gResToAvPlaneHSright_Z_negX->SetMarkerColor(kBlue);
  xposcounter=0;
  xnegcounter=0;
  double deltazmax_HSR=0;
  for(unsigned int iEntry=0; iEntry<x_HSright.size(); iEntry++) {
    double zrestonominal = ComputeResidualToPlane(x_HSright[iEntry]-12.9,y_HSright[iEntry],z_HSright[iEntry],planepars_base);
    double zresaverage = ComputeResidualToPlane(x_HSright[iEntry]-12.9,y_HSright[iEntry],z_HSright[iEntry],fPlaneHSright->GetParameters());
    if(TMath::Abs(zrestonominal)>deltazmax_HSR) deltazmax_HSR=TMath::Abs(zrestonominal);
    gResToNomPlaneHSright_Z->SetPoint(iEntry,y_HSright[iEntry],zrestonominal*1000);
    gResToAvPlaneHSright_Z->SetPoint(iEntry,y_HSright[iEntry],zresaverage*1000);
    if(x_HSright[iEntry]-12.9>0) {
      gResToNomPlaneHSright_Z_posX->SetPoint(xposcounter,y_HSright[iEntry],zrestonominal*1000);
      gResToAvPlaneHSright_Z_posX->SetPoint(xposcounter,y_HSright[iEntry],zresaverage*1000);
      xposcounter++;
    }
    else {
      gResToNomPlaneHSright_Z_negX->SetPoint(xnegcounter,y_HSright[iEntry],zrestonominal*1000);
      gResToAvPlaneHSright_Z_negX->SetPoint(xnegcounter,y_HSright[iEntry],zresaverage*1000);
      xnegcounter++;
    }
  }
  TF1* fAvPlaneHSright_posX = new TF1("fAvPlaneHSright_posX","pol1",-800,800);
  fAvPlaneHSright_posX->SetParameters((fPlaneHSright->GetParameter(0)+fPlaneHSright->GetParameter(1)*(14.999)-planepars_base[0])*1000,fPlaneHSright->GetParameter(2)*1000);
  fAvPlaneHSright_posX->SetLineColor(kRed);
  TF1* fAvPlaneHSright_negX = new TF1("fAvPlaneHSright_negX","pol1",-800,800);
  fAvPlaneHSright_negX->SetParameters((fPlaneHSright->GetParameter(0)+fPlaneHSright->GetParameter(1)*(-14.999)-planepars_base[0])*1000,fPlaneHSright->GetParameter(2)*1000);
  fAvPlaneHSright_negX->SetLineColor(kBlue);

  //HS left (final stave)
  TGraph* gResToNomPlaneStaveHSleft_Z = new TGraph(x_FinalStave_HSleft.size());
  gResToNomPlaneStaveHSleft_Z->SetName("gResToNomPlaneStaveHSleft_Z");
  gResToNomPlaneStaveHSleft_Z->SetTitle("HSL planarity (final stave metrology)");
  gResToNomPlaneStaveHSleft_Z->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneStaveHSleft_Z->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneStaveHSleft_Z->SetMarkerStyle(kFullCircle);
  gResToNomPlaneStaveHSleft_Z->SetMarkerSize(1.);
  gResToNomPlaneStaveHSleft_Z->SetMarkerColor(kBlack);
  TGraph* gResToNomPlaneStaveHSleft_Z_posX = new TGraph(0);
  gResToNomPlaneStaveHSleft_Z_posX->SetName("gResToNomPlaneStaveHSleft_Z_posX");
  gResToNomPlaneStaveHSleft_Z_posX->SetTitle("HSL planarity (final stave metrology) - pos X");
  gResToNomPlaneStaveHSleft_Z_posX->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneStaveHSleft_Z_posX->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneStaveHSleft_Z_posX->SetMarkerStyle(kFullCircle);
  gResToNomPlaneStaveHSleft_Z_posX->SetMarkerSize(1.);
  gResToNomPlaneStaveHSleft_Z_posX->SetMarkerColor(kRed);
  TGraph* gResToNomPlaneStaveHSleft_Z_negX = new TGraph(0);
  gResToNomPlaneStaveHSleft_Z_negX->SetName("gResToNomPlaneStaveHSleft_Z_negX");
  gResToNomPlaneStaveHSleft_Z_negX->SetTitle("HSL planarity (final stave metrology) - neg X");
  gResToNomPlaneStaveHSleft_Z_negX->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneStaveHSleft_Z_negX->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneStaveHSleft_Z_negX->SetMarkerStyle(kFullCircle);
  gResToNomPlaneStaveHSleft_Z_negX->SetMarkerSize(1.);
  gResToNomPlaneStaveHSleft_Z_negX->SetMarkerColor(kBlue);
  TGraph* gResToAvPlaneStaveHSleft_Z = new TGraph(x_FinalStave_HSleft.size());
  gResToAvPlaneStaveHSleft_Z->SetName("gResToAvPlaneStaveHSleft_Z");
  gResToAvPlaneStaveHSleft_Z->SetTitle("HSL planarity (final stave metrology)");
  gResToAvPlaneStaveHSleft_Z->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneStaveHSleft_Z->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneStaveHSleft_Z->SetMarkerStyle(kFullCircle);
  gResToAvPlaneStaveHSleft_Z->SetMarkerSize(1.);
  gResToAvPlaneStaveHSleft_Z->SetMarkerColor(kBlack);
  TGraph* gResToAvPlaneStaveHSleft_Z_posX = new TGraph(0);
  gResToAvPlaneStaveHSleft_Z_posX->SetName("gResToAvPlaneStaveHSleft_Z_posX");
  gResToAvPlaneStaveHSleft_Z_posX->SetTitle("HSL planarity (final stave metrology) - pos X");
  gResToAvPlaneStaveHSleft_Z_posX->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneStaveHSleft_Z_posX->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneStaveHSleft_Z_posX->SetMarkerStyle(kFullCircle);
  gResToAvPlaneStaveHSleft_Z_posX->SetMarkerSize(1.);
  gResToAvPlaneStaveHSleft_Z_posX->SetMarkerColor(kRed);
  TGraph* gResToAvPlaneStaveHSleft_Z_negX = new TGraph(0);
  gResToAvPlaneStaveHSleft_Z_negX->SetName("gResToAvPlaneStaveHSleft_Z_negX");
  gResToAvPlaneStaveHSleft_Z_negX->SetTitle("HSL planarity (final stave metrology) - neg X");
  gResToAvPlaneStaveHSleft_Z_negX->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneStaveHSleft_Z_negX->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneStaveHSleft_Z_negX->SetMarkerStyle(kFullCircle);
  gResToAvPlaneStaveHSleft_Z_negX->SetMarkerSize(1.);
  gResToAvPlaneStaveHSleft_Z_negX->SetMarkerColor(kBlue);
  double planepars_HSleft[3] = {13.27,0.,0.};
  xposcounter=0;
  xnegcounter=0;
  double deltazmax_stave_HSL=0;
  for(unsigned int iEntry=0; iEntry<x_FinalStave_HSleft.size(); iEntry++) {
    double zrestonominal = ComputeResidualToPlane(x_FinalStave_HSleft[iEntry],y_FinalStave_HSleft[iEntry],z_FinalStave_HSleft[iEntry],planepars_HSleft);
    double zrestoaverage = ComputeResidualToPlane(x_FinalStave_HSleft[iEntry],y_FinalStave_HSleft[iEntry],z_FinalStave_HSleft[iEntry],fPlaneStaveHSleft->GetParameters());
    if(TMath::Abs(zrestonominal)>deltazmax_stave_HSL) deltazmax_stave_HSL=TMath::Abs(zrestonominal);
    gResToNomPlaneStaveHSleft_Z->SetPoint(iEntry,y_FinalStave_HSleft[iEntry],zrestonominal*1000);
    gResToAvPlaneStaveHSleft_Z->SetPoint(iEntry,y_FinalStave_HSleft[iEntry],zrestoaverage*1000);
    if(x_FinalStave_HSleft[iEntry]>0) {
      gResToNomPlaneStaveHSleft_Z_posX->SetPoint(xposcounter,y_FinalStave_HSleft[iEntry],zrestonominal*1000);
      gResToAvPlaneStaveHSleft_Z_posX->SetPoint(xposcounter,y_FinalStave_HSleft[iEntry],zrestoaverage*1000);
      xposcounter++;
    }
    else {
      gResToNomPlaneStaveHSleft_Z_negX->SetPoint(xnegcounter,y_FinalStave_HSleft[iEntry],zrestonominal*1000);
      gResToAvPlaneStaveHSleft_Z_negX->SetPoint(xnegcounter,y_FinalStave_HSleft[iEntry],zrestoaverage*1000);
      xnegcounter++;
    }
  }
  TF1* fAvPlaneStaveHSleft_posX = new TF1("fAvPlaneStaveHSleft_posX","pol1",-800,800);
  fAvPlaneStaveHSleft_posX->SetParameters((fPlaneStaveHSleft->GetParameter(0)+fPlaneStaveHSleft->GetParameter(1)*(2.099)-planepars_HSleft[0])*1000,fPlaneStaveHSleft->GetParameter(2)*1000);
  fAvPlaneStaveHSleft_posX->SetLineColor(kRed);
  TF1* fAvPlaneStaveHSleft_negX = new TF1("fAvPlaneStaveHSleft_negX","pol1",-800,800);
  fAvPlaneStaveHSleft_negX->SetParameters((fPlaneStaveHSleft->GetParameter(0)+fPlaneStaveHSleft->GetParameter(1)*(-27.899)-planepars_HSleft[0])*1000,fPlaneStaveHSleft->GetParameter(2)*1000);
  fAvPlaneStaveHSleft_negX->SetLineColor(kBlue);

  //HS right (final stave)
  TGraph* gResToNomPlaneStaveHSright_Z = new TGraph(x_FinalStave_HSright.size());
  gResToNomPlaneStaveHSright_Z->SetName("gResToNomPlaneStaveHSright_Z");
  gResToNomPlaneStaveHSright_Z->SetTitle("HSR planarity (final stave metrology)");
  gResToNomPlaneStaveHSright_Z->GetXaxis()->SetTitle("y (mm)");
  gResToNomPlaneStaveHSright_Z->GetYaxis()->SetTitle("z_{meas} - z_{nom plane} (#mum)");
  gResToNomPlaneStaveHSright_Z->SetMarkerStyle(kFullCircle);
  gResToNomPlaneStaveHSright_Z->SetMarkerSize(1.);
  gResToNomPlaneStaveHSright_Z->SetMarkerColor(kBlack);
  TGraph* gResToAvPlaneStaveHSright_Z = new TGraph(x_FinalStave_HSright.size());
  gResToAvPlaneStaveHSright_Z->SetName("gResToAvPlaneStaveHSright_Z");
  gResToAvPlaneStaveHSright_Z->SetTitle("HSR planarity (final stave metrology)");
  gResToAvPlaneStaveHSright_Z->GetXaxis()->SetTitle("y (mm)");
  gResToAvPlaneStaveHSright_Z->GetYaxis()->SetTitle("z_{meas} - z_{average plane} (#mum)");
  gResToAvPlaneStaveHSright_Z->SetMarkerStyle(kFullCircle);
  gResToAvPlaneStaveHSright_Z->SetMarkerSize(1.);
  gResToAvPlaneStaveHSright_Z->SetMarkerColor(kBlack);
  double planepars_HSright[3] = {9.67,0.,0.};
  double deltazmax_stave_HSR=0;
  for(unsigned int iEntry=0; iEntry<x_FinalStave_HSright.size(); iEntry++) {
    double zrestonominal = ComputeResidualToPlane(x_FinalStave_HSright[iEntry],y_FinalStave_HSright[iEntry],z_FinalStave_HSright[iEntry],planepars_HSright);
    double zrestoaverage = ComputeResidualToPlane(x_FinalStave_HSright[iEntry],y_FinalStave_HSright[iEntry],z_FinalStave_HSright[iEntry],fPlaneStaveHSright->GetParameters());
    if(TMath::Abs(zrestonominal)>deltazmax_stave_HSR) deltazmax_stave_HSR=TMath::Abs(zrestonominal);
    gResToNomPlaneStaveHSright_Z->SetPoint(iEntry,y_FinalStave_HSright[iEntry],zrestonominal*1000);
    gResToAvPlaneStaveHSright_Z->SetPoint(iEntry,y_FinalStave_HSright[iEntry],zrestoaverage*1000);
  }
  TF1* fAvPlaneStaveHSright = new TF1("fAvPlaneStaveHSright","pol1",-800,800);
  fAvPlaneStaveHSright->SetParameters((fPlaneStaveHSright->GetParameter(0)+fPlaneStaveHSright->GetParameter(1)*(27.899)-planepars_HSright[0])*1000,fPlaneStaveHSright->GetParameter(2)*1000);
  fAvPlaneStaveHSright->SetLineColor(kBlack);

  //extrapolate position of markers that are not visible in the last metrology survey from the measurements of the two HS (before U-arms)
  std::vector<double> x_extrap_HSleft;
  std::vector<double> y_extrap_HSleft;
  std::vector<double> z_extrap_HSleft;

  std::vector<double> x_extrap_HSright;
  std::vector<double> y_extrap_HSright;
  std::vector<double> z_extrap_HSright;

  int nextrap_left=0;
  int nextrap_right=0;
  DoInvisibleMarkerExtrapolation(MainFrame::kHSL,nextrap_left,x_HSleft_nominal,y_HSleft_nominal,x_HSleft_filled,y_HSleft_filled,z_HSleft_filled,x_FinalStave_HSleft_filled,y_FinalStave_HSleft_filled,z_FinalStave_HSleft_filled,x_extrap_HSleft,y_extrap_HSleft,z_extrap_HSleft,true,fPlaneHSleft->GetParameters(),fPlaneStaveHSleft->GetParameters());
  DoInvisibleMarkerExtrapolation(MainFrame::kHSR,nextrap_right,x_HSright_nominal,y_HSright_nominal,x_HSright_filled,y_HSright_filled,z_HSright_filled,x_FinalStave_HSright_filled,y_FinalStave_HSright_filled,z_FinalStave_HSright_filled,x_extrap_HSright,y_extrap_HSright,z_extrap_HSright,true,fPlaneHSleft->GetParameters(),fPlaneStaveHSleft->GetParameters());

  //vectors of final positions (measured when available, otherwise extrapolated)
  std::vector<double> x_FinalStave_MeasPlusExtr;
  std::vector<double> y_FinalStave_MeasPlusExtr;
  std::vector<double> z_FinalStave_MeasPlusExtr;
  //first HSL
  for(unsigned int iNomEntry=0; iNomEntry<x_HSleft_nominal.size(); iNomEntry++) {
    if(x_FinalStave_HSleft_filled[iNomEntry]!=-10000.) {
      x_FinalStave_MeasPlusExtr.push_back(x_FinalStave_HSleft_filled[iNomEntry]);
      y_FinalStave_MeasPlusExtr.push_back(y_FinalStave_HSleft_filled[iNomEntry]);
      z_FinalStave_MeasPlusExtr.push_back(z_FinalStave_HSleft_filled[iNomEntry]);
    }
    else if(x_extrap_HSright[iNomEntry]!=-10000.) {
      x_FinalStave_MeasPlusExtr.push_back(x_extrap_HSright[iNomEntry]);
      y_FinalStave_MeasPlusExtr.push_back(y_extrap_HSright[iNomEntry]);
      z_FinalStave_MeasPlusExtr.push_back(z_extrap_HSright[iNomEntry]);
    }
  }
  //and then HSR
  for(unsigned int iNomEntry=0; iNomEntry<x_HSright_nominal.size(); iNomEntry++) {
    if(x_FinalStave_HSright_filled[iNomEntry]!=-10000.) {
      x_FinalStave_MeasPlusExtr.push_back(x_FinalStave_HSright_filled[iNomEntry]);
      y_FinalStave_MeasPlusExtr.push_back(y_FinalStave_HSright_filled[iNomEntry]);
      z_FinalStave_MeasPlusExtr.push_back(z_FinalStave_HSright_filled[iNomEntry]);
    }
    else if(x_extrap_HSright[iNomEntry]!=-10000.) {
      x_FinalStave_MeasPlusExtr.push_back(x_extrap_HSright[iNomEntry]);
      y_FinalStave_MeasPlusExtr.push_back(y_extrap_HSright[iNomEntry]);
      z_FinalStave_MeasPlusExtr.push_back(z_extrap_HSright[iNomEntry]);
    }
  }

  //Fill graphs / histograms
  //extrapolated measurements
  TGraph2D* gStaveExtrap = new TGraph2D(nextrap_right+nextrap_left);
  gStaveExtrap->SetName("gStaveExtrap");
  gStaveExtrap->SetMarkerStyle(kFullSquare);
  gStaveExtrap->SetMarkerColor(kRed);
  int nentryleft=0;
  for(unsigned int iEntryLeft=0; iEntryLeft<x_extrap_HSleft.size(); iEntryLeft++) {
    if(x_extrap_HSleft[iEntryLeft]!=-10000.) {
      gStaveExtrap->SetPoint(nentryleft,x_extrap_HSleft[iEntryLeft],y_extrap_HSleft[iEntryLeft],z_extrap_HSleft[iEntryLeft]);
      nentryleft++;
    }
  }
  int nentryright=0;
  for(unsigned int iEntryRight=0; iEntryRight<x_extrap_HSright.size(); iEntryRight++) {
    if(x_extrap_HSright[iEntryRight]!=-10000.) {
      gStaveExtrap->SetPoint(nentryright+nentryleft,x_extrap_HSright[iEntryRight],y_extrap_HSright[iEntryRight],z_extrap_HSright[iEntryRight]);
      nentryright++;
    }
  }

  TH1F* hStaveMeasExtrapResX = new TH1F("hStaveMeasExtrapResX","Extrapolated - measured marker position X;x_{extr}-x_{meas} (#mum);Entries",100,-50.,50.);
  TH1F* hStaveMeasExtrapResY = new TH1F("hStaveMeasExtrapResY","Extrapolated - measured marker position Y;y_{extr}-y_{meas} (#mum);Entries",100,-50.,50.);
  TH1F* hStaveMeasExtrapResZ = new TH1F("hStaveMeasExtrapResZ","Extrapolated - measured marker position Z;z_{extr}-z_{meas} (#mum);Entries",100,-500.,500.);
  for(unsigned int iEntry=0; iEntry<x_extrap_HSleft.size(); iEntry++) {
    if(x_extrap_HSleft[iEntry]!=-10000. && x_FinalStave_HSleft_filled[iEntry]!=-10000.) {
      hStaveMeasExtrapResX->Fill((x_extrap_HSleft[iEntry]-x_FinalStave_HSleft_filled[iEntry])*1000);
      hStaveMeasExtrapResY->Fill((y_extrap_HSleft[iEntry]-y_FinalStave_HSleft_filled[iEntry])*1000);
      hStaveMeasExtrapResZ->Fill((z_extrap_HSleft[iEntry]-z_FinalStave_HSleft_filled[iEntry])*1000);
    }
  }

  //Draw results
  TH3F* hFrameHS = new TH3F("hFrame",";x (mm);y (mm);z (mm)",100,-50.,50.,100,-1000,1000,100,0.,1.);
  hFrameHS->SetStats(0);
  TH3F* hFrameStave = new TH3F("hFrame",";x (mm);y (mm);z (mm)",100,-50.,50.,100,-1000,1000,100,8.,15.);
  hFrameStave->SetStats(0);

  TLegend* legHS = new TLegend(0.6,0.725,0.7,0.825);
  legHS->SetTextSize(0.045);
  legHS->AddEntry(gHSleftMeas,"HSL","p");
  legHS->AddEntry(gHSrightMeas,"HSR","p");

  TLegend* legStave = new TLegend(0.6,0.725,0.75,0.825);
  legStave->SetTextSize(0.045);
  legStave->AddEntry(gStaveMeas,"Measured","p");
  legStave->AddEntry(gStaveExtrap,"Extrapolated","p");

  Double_t x1 = .57, y1 = .65, x2 = .93, y2 = .85;
  TLegend* legHS_X = new TLegend(x1,y1,x2,y2);
  legHS_X->SetTextSize(0.045);
  legHS_X->AddEntry(gResToNomPlaneHSright_Z_posX,"x = 15 mm","p");
  legHS_X->AddEntry(gResToNomPlaneHSright_Z_negX,"x = -15 mm","p");
  legHS_X->AddEntry(fAvPlaneHSright_posX,"average plane x = 15 mm","l");
  legHS_X->AddEntry(fAvPlaneHSright_negX,"average plane x = -15 mm","l");

  TLegend* legStaveHSL_X = new TLegend(x1,y1,x2,y2);
  legStaveHSL_X->SetTextSize(0.045);
  legStaveHSL_X->AddEntry(gResToNomPlaneStaveHSleft_Z_posX,"x = 2 mm","p");
  legStaveHSL_X->AddEntry(gResToNomPlaneStaveHSleft_Z_negX,"x = -28 mm","p");
  legStaveHSL_X->AddEntry(fAvPlaneStaveHSleft_posX,"average plane x = 2 mm","l");
  legStaveHSL_X->AddEntry(fAvPlaneStaveHSleft_negX,"average plane x = -28 mm","l");

  TLegend* legStaveHSR_X = new TLegend(x1,y1,x2,y2);
  legStaveHSR_X->SetTextSize(0.045);
  legStaveHSR_X->AddEntry(gResToNomPlaneStaveHSright_Z,"x = 28 mm","p");
  legStaveHSR_X->AddEntry(fAvPlaneStaveHSright,"average plane x = 28 mm","l");

  TPaveText* text[4];
  for(int iText=0; iText<4; iText++) {
    text[iText] = new TPaveText(0.15,0.65,0.53,0.85,"NDC");
    text[iText]->SetTextSize(0.045);
    text[iText]->SetTextColor(kBlack);
    text[iText]->SetTextFont(42);
    text[iText]->SetFillColor(kWhite);
    text[iText]->SetBorderSize(1);
  }

  TCanvas*& cHS = fCanvas[0];
  cHS->cd();
  hFrameHS->Draw();
  gHSleftMeas->Draw("Psame");
  gHSrightMeas->Draw("Psame");
  legHS->Draw("same");
  cHS->Update();

  TCanvas*& cHSleft_plan = fCanvas[1];
  cHSleft_plan->cd();
  gResToNomPlaneHSleft_Z->GetYaxis()->SetRangeUser(-500.,800);
  gResToNomPlaneHSleft_Z->Draw("AP");
  gResToNomPlaneHSleft_Z_posX->Draw("P");
  gResToNomPlaneHSleft_Z_negX->Draw("P");
  fAvPlaneHSleft_posX->Draw("same");
  fAvPlaneHSleft_negX->Draw("same");
  double planaritytonominal=0, meantonominal=0, RMStonominal=0, planaritytoaverage=0, meantoaverage=0, RMStoaverage=0;
  ComputePlanarityRMSandMean(gResToNomPlaneHSleft_Z,planaritytonominal,RMStonominal,meantonominal);
  ComputePlanarityRMSandMean(gResToAvPlaneHSleft_Z,planaritytoaverage,RMStoaverage,meantoaverage);
  text[0]->AddText(Form("planarity (to average) = %0.1f #mum",planaritytoaverage));
  text[0]->AddText(Form("RMS (to average) = %0.1f #mum",RMStoaverage));
  text[0]->AddText(Form("shift (to nominal)  = %0.1f #mum",meantonominal));
  text[0]->AddText(Form("Max #DeltaZ (to nominal) = %0.1f #mum",deltazmax_HSL*1000));
  text[0]->Draw("same");
  legHS_X->Draw("same");
  cHSleft_plan->Update();

  TCanvas*& cHSright_plan = fCanvas[2];
  cHSright_plan->cd();
  gResToNomPlaneHSright_Z->GetYaxis()->SetRangeUser(-500.,800);
  gResToNomPlaneHSright_Z->Draw("AP");
  gResToNomPlaneHSright_Z_posX->Draw("P");
  gResToNomPlaneHSright_Z_negX->Draw("P");
  fAvPlaneHSright_posX->Draw("same");
  fAvPlaneHSright_negX->Draw("same");
  ComputePlanarityRMSandMean(gResToNomPlaneHSright_Z,planaritytonominal,RMStonominal,meantonominal);
  ComputePlanarityRMSandMean(gResToAvPlaneHSright_Z,planaritytoaverage,RMStoaverage,meantoaverage);
  text[1]->AddText(Form("planarity (to average) = %0.1f #mum",planaritytoaverage));
  text[1]->AddText(Form("RMS (to average) = %0.1f #mum",RMStoaverage));
  text[1]->AddText(Form("mean (to nominal) = %0.1f #mum",meantonominal));
  text[1]->AddText(Form("Max #DeltaZ (to nominal) = %0.1f #mum",deltazmax_HSR*1000));
  text[1]->Draw("same");
  legHS_X->Draw("same");
  cHSright_plan->Update();

  TCanvas*& cHSleft_res = fCanvas[3];
  cHSleft_res->Divide(2,1);
  cHSleft_res->cd(1);
  hResToNominalHSleft_X->Draw();
  cHSleft_res->cd(2);
  hResToNominalHSleft_Y->Draw();
  cHSleft_res->Update();

  TCanvas*& cHSright_res = fCanvas[4];
  cHSright_res->Divide(2,1);
  cHSright_res->cd(1);
  hResToNominalHSright_X->Draw();
  cHSright_res->cd(2);
  hResToNominalHSright_Y->Draw();
  cHSright_res->Update();

  TCanvas*& cStave = fCanvas[9];
  cStave->cd();
  hFrameStave->Draw();
  gStaveMeas->Draw("Psame");
  gStaveExtrap->Draw("Psame");
  legStave->Draw("same");
  cStave->Update();

  TCanvas*& cStaveHSleft_plan = fCanvas[5];
  cStaveHSleft_plan->cd();
  gResToNomPlaneStaveHSleft_Z->GetYaxis()->SetRangeUser(-800,800);
  gResToNomPlaneStaveHSleft_Z->Draw("AP");
  gResToNomPlaneStaveHSleft_Z_posX->Draw("P");
  gResToNomPlaneStaveHSleft_Z_negX->Draw("P");
  fAvPlaneStaveHSleft_posX->Draw("same");
  fAvPlaneStaveHSleft_negX->Draw("same");
  ComputePlanarityRMSandMean(gResToNomPlaneStaveHSleft_Z,planaritytonominal,RMStonominal,meantonominal);
  ComputePlanarityRMSandMean(gResToAvPlaneStaveHSleft_Z,planaritytoaverage,RMStoaverage,meantoaverage);
  text[2]->AddText(Form("planarity (to average) = %0.1f #mum",planaritytoaverage));
  text[2]->AddText(Form("RMS (to average) = %0.1f #mum",RMStoaverage));
  text[2]->AddText(Form("mean (to nominal) = %0.1f #mum",meantonominal));
  text[2]->AddText(Form("Max #DeltaZ (to nominal) = %0.1f #mum",deltazmax_stave_HSL*1000));
  text[2]->Draw("same");
  legStaveHSL_X->Draw("same");
  cStaveHSleft_plan->Update();

  TCanvas*& cStaveHSright_plan = fCanvas[6];
  cStaveHSright_plan->cd();
  gResToNomPlaneStaveHSright_Z->GetYaxis()->SetRangeUser(-800,800);
  gResToNomPlaneStaveHSright_Z->Draw("AP");
  fAvPlaneStaveHSright->Draw("same");
  ComputePlanarityRMSandMean(gResToNomPlaneStaveHSright_Z,planaritytonominal,RMStonominal,meantonominal);
  ComputePlanarityRMSandMean(gResToAvPlaneStaveHSright_Z,planaritytoaverage,RMStoaverage,meantoaverage);
  text[3]->AddText(Form("planarity (to average) = %0.1f #mum",planaritytoaverage));
  text[3]->AddText(Form("RMS (to average) = %0.1f #mum",RMStoaverage));
  text[3]->AddText(Form("mean (to nominal) = %0.1f #mum",meantonominal));
  text[3]->AddText(Form("Max #DeltaZ (to nominal) = %0.1f #mum",deltazmax_stave_HSR*1000));
  text[3]->Draw("same");
  legStaveHSR_X->Draw("same");
  cStaveHSright_plan->Update();

  TCanvas*& cFinalStave_HSleft_res = fCanvas[7];
  cFinalStave_HSleft_res->Divide(2,1);
  cFinalStave_HSleft_res->cd(1);
  hResToNominalStave_HSleft_X->Draw();
  cFinalStave_HSleft_res->cd(2);
  hResToNominalStave_HSleft_Y->Draw();
  cFinalStave_HSleft_res->Update();

  TCanvas*& cFinalStave_HSright_res = fCanvas[8];
  cFinalStave_HSright_res->Divide(2,1);
  cFinalStave_HSright_res->cd(1);
  hResToNominalStave_HSright_X->Draw();
  cFinalStave_HSright_res->cd(2);
  hResToNominalStave_HSright_Y->Draw();
  cFinalStave_HSright_res->Update();

  TCanvas*& cStaveMeasExtrRes = fCanvas[10];
  cStaveMeasExtrRes->Divide(3,1);
  cStaveMeasExtrRes->cd(1);
  hStaveMeasExtrapResX->Draw();
  cStaveMeasExtrRes->cd(2);
  hStaveMeasExtrapResY->Draw();
  cStaveMeasExtrRes->cd(3);
  hStaveMeasExtrapResZ->Draw();
  cStaveMeasExtrRes->Update();

  //Save output files
  //dat file with extrapolated values
  TString outfilenamedat = infilename_FinalStave;
  outfilenamedat.ReplaceAll(".dat","_MeasAndExtrap.dat");
  outfilenamedat.ReplaceAll(".txt","_MeasAndExtrap.txt");
  outfilenamedat.ReplaceAll(".csv","_MeasAndExtrap.csv");
  CreateDatFile(outfilenamedat.Data(),x_FinalStave_MeasPlusExtr,y_FinalStave_MeasPlusExtr,z_FinalStave_MeasPlusExtr);

  //root file with QA plots
  TString outfilenameroot = outfilenamedat;
  outfilenameroot.ReplaceAll(".dat",".root");
  outfilenameroot.ReplaceAll(".txt",".root");
  outfilenameroot.ReplaceAll(".csv",".root");
  TFile outfile(outfilenameroot.Data(),"RECREATE");
  cHS->Write();
  cHSleft_plan->Write();
  cHSright_plan->Write();
  cHSleft_res->Write();
  cHSright_res->Write();
  cStave->Write();
  cStaveHSleft_plan->Write();
  cStaveHSright_plan->Write();
  cFinalStave_HSleft_res->Write();
  cFinalStave_HSright_res->Write();
  cStaveMeasExtrRes->Write();
  gHSleftMeas->Write();
  gHSrightMeas->Write();
  gResToNomPlaneHSleft_Z->Write();
  gResToNomPlaneHSleft_Z_posX->Write();
  gResToNomPlaneHSleft_Z_negX->Write();
  gResToNomPlaneHSright_Z->Write();
  gResToNomPlaneHSright_Z_posX->Write();
  gResToNomPlaneHSright_Z_negX->Write();
  gResToAvPlaneHSleft_Z->Write();
  gResToAvPlaneHSleft_Z_posX->Write();
  gResToAvPlaneHSleft_Z_negX->Write();
  gResToAvPlaneHSright_Z->Write();
  gResToAvPlaneHSright_Z_posX->Write();
  gResToAvPlaneHSright_Z_negX->Write();
  hResToNominalHSleft_X->Write();
  hResToNominalHSleft_Y->Write();
  hResToNominalHSright_X->Write();
  hResToNominalHSright_Y->Write();
  gStaveMeas->Write();
  gResToNomPlaneStaveHSleft_Z->Write();
  gResToNomPlaneStaveHSleft_Z_posX->Write();
  gResToNomPlaneStaveHSleft_Z_negX->Write();
  gResToNomPlaneStaveHSright_Z->Write();
  gResToAvPlaneStaveHSleft_Z->Write();
  gResToAvPlaneStaveHSleft_Z_posX->Write();
  gResToAvPlaneStaveHSleft_Z_negX->Write();
  gResToAvPlaneStaveHSright_Z->Write();
  hResToNominalStave_HSleft_X->Write();
  hResToNominalStave_HSleft_Y->Write();
  hResToNominalStave_HSright_X->Write();
  hResToNominalStave_HSright_Y->Write();
  gStaveExtrap->Write();
  hStaveMeasExtrapResX->Write();
  hStaveMeasExtrapResY->Write();
  hStaveMeasExtrapResZ->Write();
  outfile.Close();

  //pdf file with QA plots
  TString outfilenamepdf = outfilenameroot;
  outfilenamepdf.ReplaceAll(".root",".pdf");
  cHS->Print(Form("%s[",outfilenamepdf.Data()));
  cHS->Print(outfilenamepdf.Data());
  cHSleft_plan->Print(outfilenamepdf.Data());
  cHSright_plan->Print(outfilenamepdf.Data());
  cHSleft_res->Print(outfilenamepdf.Data());
  cHSright_res->Print(outfilenamepdf.Data());
  cStave->Print(outfilenamepdf.Data());
  cStaveHSleft_plan->Print(outfilenamepdf.Data());
  cStaveHSright_plan->Print(outfilenamepdf.Data());
  cFinalStave_HSleft_res->Print(outfilenamepdf.Data());
  cFinalStave_HSright_res->Print(outfilenamepdf.Data());
  cStaveMeasExtrRes->Print(outfilenamepdf.Data());
  cStaveMeasExtrRes->Print(Form("%s]",outfilenamepdf.Data()));

  return 0;
}

//_____________________________________________________________________________________________
bool ReadDatFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double xoffset)
{
  if(!FileName.Contains("txt") && !FileName.Contains("dat") && !FileName.Contains(".csv")) {
    std::cerr << "Wrong file format. Exit." << std::endl;
    return false;
  }

  ifstream inSet(FileName.Data());
  if(!inSet) {
    std::cerr<<"Please check if "<<FileName.Data() <<" is the right path. Exit."<<std::endl;
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
    std::cerr << "Wrong file format. Exit." << std::endl;
    return false;
  }

  string valueSeparator = ";";
  ifstream inSet(FileName.Data());
  if(!inSet) {
    std::cerr<<"Please check if "<<FileName.Data() <<" is the right path. Exit."<<std::endl;
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
void CreateDatFile(TString FileName, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
{
  ofstream outSet;
  outSet.open(FileName.Data());

  for(unsigned int iEntry=0; iEntry<z.size(); iEntry++) {
    outSet << x[iEntry] << " " << y[iEntry] << " " << z[iEntry] << std::endl;
  }

  std::cout << "File " << FileName.Data() << " saved." << std::endl;
  outSet.close();
}

//______________________________________________________________________________________________
bool MatchMeasToNominal(double xmeas, double ymeas, double xnom, double ynom, double accdeltax, double accdeltay, double yshift) {

  if(TMath::Abs(xmeas-xnom)<accdeltax && TMath::Abs((ymeas-yshift)-ynom)<accdeltay) return true;

  return false;
}

//______________________________________________________________________________________________
void FillVectorsOfMatchedPositions(std::vector<double> xmeas, std::vector<double> ymeas, std::vector<double> zmeas, std::vector<double> xnom, std::vector<double> ynom, std::vector<double>& xmeasfilled, std::vector<double>& ymeasfilled, std::vector<double>& zmeasfilled, std::vector<double>& xres, std::vector<double>& yres) {

  //positions of initial markers in std::vector
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

      xres.push_back((xmeas[matchedmarkerpos]-xnom[iNomEntry])*1000);
      yres.push_back((ymeas[matchedmarkerpos]-ynom[iNomEntry])*1000);

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
int FindOppositeNominalMarkerPosition(double x, double y, std::vector<double> xvec, std::vector<double> yvec) {
  int posmarker=-1;
  for(unsigned int iEntry=0; iEntry<xvec.size(); iEntry++) {
    if(yvec[iEntry]==y && xvec[iEntry]!=x) {
      posmarker=iEntry;
      break;
    }
  }

  return posmarker;
}

//______________________________________________________________________________________________
double ComputeResidualToPlane(double x, double y, double z, double *pars) {

  double zplane = pars[0]+pars[1]*x+pars[2]*y;
  double res = z - zplane;

  return res;
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

//______________________________________________________________________________________________
bool PrintMatchedCoordinates(std::vector<double> xmeas, std::vector<double> ymeas, std::vector<double> xnom, std::vector<double> ynom, double shiftx) {

  if(xmeas.size()!=xnom.size() || ymeas.size()!=ynom.size()) {
    std::cerr << "Different number of entries between measured and nominal positions. Exit." << std::endl;
    return false;
  }

  for(unsigned int iCoord=0; iCoord<xmeas.size(); iCoord++) {
    std::cout << Form("%d:  xmeas = %0.4f  ymeas = %0.4f  xnom = %0.4f  ynom = %0.4f  resx = %0.4f resy =  %0.4f",iCoord,xmeas[iCoord]+shiftx,ymeas[iCoord],xnom[iCoord]+shiftx,ynom[iCoord],xmeas[iCoord]-xnom[iCoord],ymeas[iCoord]-ynom[iCoord])<<std::endl;
  }

  return true;
}

//______________________________________________________________________________________________
void DoInvisibleMarkerExtrapolation(int LeftOrRight, int &nextrap, std::vector<double> x_nominal, std::vector<double> y_nominal, std::vector<double> x_filled_HS, std::vector<double> y_filled_HS, std::vector<double> z_filled_HS, std::vector<double> x_filled_Stave, std::vector<double> y_filled_Stave, std::vector<double> z_filled_Stave, std::vector<double>& x_extrap, std::vector<double>& y_extrap, std::vector<double>& z_extrap, bool correctfortilt, double* parsHSplane, double* parsStaveplane) {

  //positions of initial markers in std::vector
  std::vector<unsigned int>  modposcornermarkers;
  std::vector<unsigned int>  modnegcornermarkers;

  modnegcornermarkers.push_back(14);
  modnegcornermarkers.push_back(27);
  modnegcornermarkers.push_back(42);
  modnegcornermarkers.push_back(55);
  modnegcornermarkers.push_back(70);
  modnegcornermarkers.push_back(83);
  modnegcornermarkers.push_back(98);
  modnegcornermarkers.push_back(111);
  modnegcornermarkers.push_back(126);
  modnegcornermarkers.push_back(139);
  modnegcornermarkers.push_back(154);
  modnegcornermarkers.push_back(167);
  modnegcornermarkers.push_back(182);
  modnegcornermarkers.push_back(195);

  modposcornermarkers.push_back(0);
  modposcornermarkers.push_back(13);
  modposcornermarkers.push_back(28);
  modposcornermarkers.push_back(41);
  modposcornermarkers.push_back(56);
  modposcornermarkers.push_back(69);
  modposcornermarkers.push_back(84);
  modposcornermarkers.push_back(97);
  modposcornermarkers.push_back(112);
  modposcornermarkers.push_back(125);
  modposcornermarkers.push_back(140);
  modposcornermarkers.push_back(153);
  modposcornermarkers.push_back(168);
  modposcornermarkers.push_back(181);

  std::vector<unsigned int>::iterator it;
  bool isposcorner = false;
  bool isnegcorner = false;

  nextrap=0;
  for(unsigned int iEntry=0; iEntry<x_nominal.size(); iEntry++) {
    if((LeftOrRight==MainFrame::kHSL && x_nominal[iEntry]>0) || (LeftOrRight==MainFrame::kHSR && x_nominal[iEntry]<0)) {
      int posoppmarker = FindOppositeNominalMarkerPosition(x_nominal[iEntry],y_nominal[iEntry],x_nominal,y_nominal);

      it = find(modposcornermarkers.begin(), modposcornermarkers.end(), iEntry);
      if(it!=modposcornermarkers.end()) isposcorner=true;
      if(!isposcorner) {
        it = find(modnegcornermarkers.begin(), modnegcornermarkers.end(), iEntry);
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
        if(x_filled_HS[oppmarkarray[j]]!=-10000 && x_filled_HS[iEntry]!=-10000 && x_filled_Stave[oppmarkarray[j]]!=-10000) {
          deltax[j]= x_filled_HS[iEntry]-x_filled_HS[oppmarkarray[j]];
          deltay[j] = y_filled_HS[iEntry]-y_filled_HS[oppmarkarray[j]];
     			deltaz[j] = z_filled_HS[iEntry]-z_filled_HS[oppmarkarray[j]];
          if(correctfortilt) {
            double z_HS_base = parsHSplane[1]*x_filled_HS[iEntry]+parsHSplane[2]*y_filled_HS[iEntry];
            double xtrasl_SF = 12.9;
            if(LeftOrRight==MainFrame::kHSR) xtrasl_SF=-12.9;
            double z_HS_SF = z_HS_SF=parsStaveplane[1]*x_filled_HS[iEntry]+parsStaveplane[2]*y_filled_HS[iEntry]-parsStaveplane[1]*xtrasl_SF;
            double tiltcorrection = z_HS_base-z_HS_SF;
            deltaz[j]+=tiltcorrection;
    		  }
		      rawexrap_x[j]=x_filled_Stave[oppmarkarray[j]]+deltax[j];
			    rawexrap_y[j]=y_filled_Stave[oppmarkarray[j]]+deltay[j];
			    rawexrap_z[j]=z_filled_Stave[oppmarkarray[j]]+deltaz[j];
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
      nextrap++;
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
