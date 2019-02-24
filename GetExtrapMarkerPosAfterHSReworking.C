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

#include <RQ_OBJECT.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGTab.h>
#include <TGFileDialog.h>
#include <TRootEmbeddedCanvas.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGButtonGroup.h>
#include <TG3DLine.h>

#endif

//**********************************************************************//
//                                                                      //
//        Main Function: GetExtrapMarkerPosAfterHSReworking()           //
//                                                                      //
//**********************************************************************//

//author: Fabrizio Grosa, INFN Torino
//e-mail: fabrizio.grosa@to.infn.it

using namespace std;

//_____________________________________________________________________________________________
//METHOD PROTOTYPES
void GetExtrapMarkerPosAfterHSReworking();
int ExtrapMarkerPosAfterHSReworking(int layer, int LeftOrRight, bool ismodulereplaced[7], TString file_before_rew, TString file_after_rew, TString file_after_rew_repl[7], TString initdir, TCanvas** fCanvas);
void DoInvisibleMarkerExtrapolation(int LeftOrRight, vector<double> x_nominal, vector<double> y_nominal, vector<double> x_after_rew_filled, vector<double> y_after_rew_filled, vector<double> z_after_rew_filled, vector<double> x_before_rew_filled, vector<double> y_before_rew_filled, vector<double> z_before_rew_filled, vector<double>& x_extrap, vector<double>& y_extrap, vector<double>& z_extrap);
bool MatchMeasToNominal(double xmeas, double ymeas, double xnom, double ynom, double accdeltax, double accdeltay, double yshift);
void FillVectorsOfMatchedPositions(vector<double> xmeas, vector<double> ymeas, vector<double> zmeas, vector<double> xnom, vector<double> ynom, vector<double>& xmeasfilled, vector<double>& ymeasfilled, vector<double>& zmeasfilled);
int FindOppositeNominalMarkerPosition(double x, double y, vector<double> xvec, vector<double> yvec);
bool ReadDatFile(TString FileName, vector<double>& x, vector<double>& y, vector<double>& z);
bool ReadDatFileMitutoyo(TString FileName, vector<double>& x, vector<double>& y, vector<double>& z);
void CreateDatFileMitutoyo(TString FileName, vector<double> x_after_rew_filled, vector<double> y_after_rew_filled, vector<double> z_after_rew_filled, vector<double> x_extrap, vector<double> y_extrap, vector<double> z_extrap);
void MergeVectorMeasurements(vector<double>& x_merged, vector<double>& y_merged, vector<double>& z_merged, vector<double> x_all_afterUarms, vector<double> y_all_afterUarms, vector<double> z_all_afterUarms, vector<double> x_repl_beforeUarms[7], vector<double> y_repl_beforeUarms[7], vector<double> z_repl_beforeUarms[7], bool ismodulereplaced[7], int layer);
void SetStyle();

//_____________________________________________________________________________________________
//GUI IMPLEMENTATION
const char *fileformats[] = { "All files",     "*",
                              "ROOT files",    "*.root",
                              "ROOT macros",   "*.C",
                              "Text files",    "*.[tT][xX][tT]",
                              0,               0 };

//_____________________________________________________________________________________________
class GroupBoxInFile : public TGGroupFrame {
  private:
    TGTextButton *fButton;
    TGTextEntry  *fEntry;

  public:
    GroupBoxInFile(const TGWindow *p, const char *name);
    TGTextEntry*  GetEntry()  { return fEntry; }
    TGTextButton* GetButton() { return fButton; }

    ClassDef(GroupBoxInFile, 1);
};

//______________________________________________________________________________
GroupBoxInFile::GroupBoxInFile(const TGWindow *p, const char *name) :
   TGGroupFrame(p, name)
{
   // Group frame containing text entry and text button.
   fEntry = new TGTextEntry(this);
   fEntry->SetDefaultSize(fEntry->GetWidth()*4, fEntry->GetHeight());
   AddFrame(fEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   TGHorizontalFrame *horz = new TGHorizontalFrame(this);
   AddFrame(horz, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

   fButton = new TGTextButton(horz, "&Open");
   horz->AddFrame(fButton, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 3,3,3,3));
}

//_____________________________________________________________________________________________
class GroupBoxInFileAndSel : public TGGroupFrame {
  private:
    TGCheckButton     *fButton;
    TGTextButton      *fOpenButton;
    TGTextEntry       *fEntry;
    int               fModNum;

  public:
    GroupBoxInFileAndSel(const TGWindow *p, const char *name, int modnum);
    TGCheckButton* GetButton() { return fButton; }
    TGTextButton* GetOpenButton() { return fOpenButton; }
    TGTextEntry*  GetEntry()  { return fEntry; }
    int GetModNum() const { return fModNum; }

    ClassDef(GroupBoxInFileAndSel, 1);
};

//______________________________________________________________________________
GroupBoxInFileAndSel::GroupBoxInFileAndSel(const TGWindow *p, const char *name, int modnum) :
   TGGroupFrame(p, name),
   fModNum(modnum)
{
   // Group frame containing text button, text entry, text button.
   fEntry = new TGTextEntry(this);
   fEntry->SetDefaultSize(fEntry->GetWidth(), fEntry->GetHeight());
   AddFrame(fEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   TGHorizontalFrame *horz = new TGHorizontalFrame(this);
   AddFrame(horz, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

   fOpenButton = new TGTextButton(horz, "&Open");
   horz->AddFrame(fOpenButton, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 3,3,3,3));

   fButton = new TGCheckButton(horz, "is replaced");
   horz->AddFrame(fButton, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 3,3,3,3));
}

//______________________________________________________________________________________________
class MainFrameModRepl {
  RQ_OBJECT("MainFrameModRepl")

 private:
  TGMainFrame*                fMain;
  TGTab*                      fTab;
  GroupBoxInFile*             fGboxBeforeRewBeforeUarms;
  GroupBoxInFile*             fGboxAfterRewAfterUarms;
  GroupBoxInFileAndSel*       fGboxModReplBeforeUarms[7];
  TString                     fFileNameBeforeRewBeforeUarms;
  TString                     fFileNameAfterRewAfterUarms;
  TString                     fFileNameReplModBeforeUarms[7];
  bool                        fIsModRepl[7];
  int                         fHSFlavour;
  int                         fLayerType;
  TString                     fInitDir;
  TRootEmbeddedCanvas*        fEcanvas;

 public:
  enum HStype {kHSL, kHSR};
  enum layertype {kOL, kML};

  MainFrameModRepl(const TGWindow* p, unsigned int w, unsigned int h);
  virtual ~MainFrameModRepl();

  void DoOpenBeforeRewFile();
  void DoOpenAfterRewFile();
  void DoOpenModReplFile(int);
  void DoOpenMod1ReplFile()    { DoOpenModReplFile(0); }
  void DoOpenMod2ReplFile()    { DoOpenModReplFile(1); }
  void DoOpenMod3ReplFile()    { DoOpenModReplFile(2); }
  void DoOpenMod4ReplFile()    { DoOpenModReplFile(3); }
  void DoOpenMod5ReplFile()    { DoOpenModReplFile(4); }
  void DoOpenMod6ReplFile()    { DoOpenModReplFile(5); }
  void DoOpenMod7ReplFile()    { DoOpenModReplFile(6); }

  void SetIsMod1Replaced(bool enabled) { fIsModRepl[0] = enabled; return; }
  void SetIsMod2Replaced(bool enabled) { fIsModRepl[1] = enabled; return; }
  void SetIsMod3Replaced(bool enabled) { fIsModRepl[2] = enabled; return; }
  void SetIsMod4Replaced(bool enabled) { fIsModRepl[3] = enabled; return; }
  void SetIsMod5Replaced(bool enabled) { fIsModRepl[4] = enabled; return; }
  void SetIsMod6Replaced(bool enabled) { fIsModRepl[5] = enabled; return; }
  void SetIsMod7Replaced(bool enabled) { fIsModRepl[6] = enabled; return; }

  void DoSelectLayer(int layer) { layer == 4 ? fLayerType = kOL : fLayerType = kML; return; }
  void DoSelectHS(int HS) { HS == 4 ? fHSFlavour = kHSL : fHSFlavour = kHSR; return; }

  int DoExtrapolationAndMerge();
};

//_______________________________________________________________________________
MainFrameModRepl::MainFrameModRepl(const TGWindow* p, unsigned int w, unsigned int h):
  fFileNameBeforeRewBeforeUarms(""),
  fFileNameAfterRewAfterUarms(""),
  fHSFlavour(0),
  fLayerType(0)
{
  for(int iMod=0; iMod<7; iMod++) {
    fFileNameReplModBeforeUarms[iMod] = "";
    fIsModRepl[iMod] = false;
  }

  fInitDir = gSystem->WorkingDirectory();

  fMain = new TGMainFrame(p,w,h);
  fMain->DontCallClose(); // to avoid double deletions.
  // use hierarchical cleaning
  fMain->SetCleanup(kDeepCleanup);
  //
  TGLayoutHints* fL1=new TGLayoutHints(kLHintsTop|kLHintsLeft,2, 2, 2, 2);

   // hFrameInFilesRepl    
  TGHorizontalFrame *hFrameChoice = new TGHorizontalFrame(fMain);
  fMain->AddFrame(hFrameChoice, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,6));

  TGHorizontalFrame *hFrameInFiles = new TGHorizontalFrame(fMain);
  fMain->AddFrame(hFrameInFiles, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5));

  TGHorizontalFrame *hFrameInFilesRepl = new TGHorizontalFrame(fMain);
  fMain->AddFrame(hFrameInFilesRepl, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5));

  TGVerticalFrame *hFrameDoMerge = new TGVerticalFrame(fMain);
  fMain->AddFrame(hFrameDoMerge, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,5,5,5));

  fGboxBeforeRewBeforeUarms = new GroupBoxInFile(hFrameInFiles, "Before reworking");  
  hFrameInFiles->AddFrame(fGboxBeforeRewBeforeUarms, fL1);
  fGboxBeforeRewBeforeUarms->GetButton()->Connect("Clicked()","MainFrameModRepl",this,"DoOpenBeforeRewFile()");
  fGboxAfterRewAfterUarms = new GroupBoxInFile(hFrameInFiles, "After reworking");
  hFrameInFiles->AddFrame(fGboxAfterRewAfterUarms, fL1);
  fGboxAfterRewAfterUarms->GetButton()->Connect("Clicked()","MainFrameModRepl",this,"DoOpenAfterRewFile()");  
  for(int iMod=0; iMod<7; iMod++) {
    fGboxModReplBeforeUarms[iMod] = new GroupBoxInFileAndSel(hFrameInFilesRepl, Form("Mod %d",iMod+1),iMod);
    hFrameInFilesRepl->AddFrame(fGboxModReplBeforeUarms[iMod], fL1);
    fGboxModReplBeforeUarms[iMod]->GetOpenButton()->Connect("Clicked()","MainFrameModRepl",this,Form("DoOpenMod%dReplFile()",iMod+1));
    fGboxModReplBeforeUarms[iMod]->GetButton()->Connect("Toggled(bool)","MainFrameModRepl",this,Form("SetIsMod%dReplaced(bool)",iMod+1));
  }

  TGButtonGroup *layerButton = new TGButtonGroup(hFrameChoice, "Layer type");
  layerButton->SetTitlePos(TGGroupFrame::kCenter);
  new TGRadioButton(layerButton, "OL", kTextCenterX);
  new TGRadioButton(layerButton, "ML", kTextLeft);
  layerButton->SetButton(kTextCenterX);
  layerButton->Connect("Pressed(int)", "MainFrameModRepl", this,"DoSelectLayer(int)");
  hFrameChoice->AddFrame(layerButton, new TGLayoutHints(kLHintsTop|kLHintsExpandY));

  TGButtonGroup *HSButton = new TGButtonGroup(hFrameChoice, "HS type");
  HSButton->SetTitlePos(TGGroupFrame::kCenter);
  new TGRadioButton(HSButton, "HS-left", kTextCenterX);
  new TGRadioButton(HSButton, "HS-right", kTextLeft);
  HSButton->SetButton(kTextCenterX);
  HSButton->Connect("Pressed(int)", "MainFrameModRepl", this,"DoSelectHS(int)");
  hFrameChoice->AddFrame(HSButton, new TGLayoutHints(kLHintsTop|kLHintsExpandY,5));

  TGHorizontalFrame* hFrameBot = new TGHorizontalFrame(hFrameDoMerge);
  hFrameDoMerge->AddFrame(hFrameBot, new TGLayoutHints(kLHintsTop|kLHintsRight,2, 2, 5, 1));

  TGTextButton* MetrologyButton=new TGTextButton(hFrameBot,"&Extrapolate and merge files");
  MetrologyButton->Connect("Clicked()","MainFrameModRepl",this,"DoExtrapolationAndMerge()");
  hFrameBot->AddFrame(MetrologyButton, fL1);

  TGTextButton* exit = new TGTextButton(hFrameBot,"&Exit","gApplication->Terminate(0)");
  hFrameBot->AddFrame(exit, fL1);

  //--------- create Tab widget and some composite frames
  fTab = new TGTab(fMain, 300, 300);
  const int kNtab=4;
  TRootEmbeddedCanvas* fEcanvas[kNtab];
  TString tabnames[kNtab] = {"HS before and after reworking (3D)","HS before and after reworking (with extrap)","Residuals after-before reworking","Extrapolation goodness"};
  for(int iTab=0; iTab<kNtab; iTab++) {
    TGCompositeFrame* tf = fTab->AddTab(tabnames[iTab].Data());
  }
  fMain->AddFrame(fTab, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,2, 5, 5, 5));

  // Set a name to the main frame
  fMain->SetWindowName("Extrapolation of invisible positions after module replacement");

  // Map all subwindows of main frame
  fMain->MapSubwindows();

  // Initialize the layout algorithm
  fMain->Resize(fMain->GetDefaultSize());

  // Map main frame
  fMain->MapWindow();
}

//_______________________________________________________________________________
MainFrameModRepl::~MainFrameModRepl() {

  // Clean up used widgets: frames, buttons, layout hints
  fMain->Cleanup();
  delete fMain;
}

//_______________________________________________________________________________
void MainFrameModRepl::DoOpenBeforeRewFile() {

  TGFileInfo fi;
  fi.fFileTypes = fileformats;
  fi.fIniDir = StrDup(fFileNameBeforeRewBeforeUarms);
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  fFileNameBeforeRewBeforeUarms=fi.fFilename;
  TGTextEntry* entry = fGboxBeforeRewBeforeUarms->GetEntry();
  entry->SetText(fFileNameBeforeRewBeforeUarms.Data());
  cout << "\n\n\nOpen before reworking (before u-arms) file " << fFileNameBeforeRewBeforeUarms << endl;

  return;
}

//_______________________________________________________________________________
void MainFrameModRepl::DoOpenAfterRewFile() {

  TGFileInfo fi;
  fi.fFileTypes = fileformats;
  fi.fIniDir = StrDup(fFileNameAfterRewAfterUarms);
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  fFileNameAfterRewAfterUarms=fi.fFilename;
  TGTextEntry* entry = fGboxAfterRewAfterUarms->GetEntry();
  entry->SetText(fFileNameAfterRewAfterUarms.Data());
  cout << "Open after reworking (after u-arms) file " << fFileNameAfterRewAfterUarms << endl;

  return;
}

//_______________________________________________________________________________
void MainFrameModRepl::DoOpenModReplFile(int iMod) {

  TGFileInfo fi;
  fi.fFileTypes = fileformats;
  fi.fIniDir = StrDup(fFileNameReplModBeforeUarms[iMod]);
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  fFileNameReplModBeforeUarms[iMod]=fi.fFilename;
  TGTextEntry* entry = fGboxModReplBeforeUarms[iMod]->GetEntry();
  entry->SetText(fFileNameReplModBeforeUarms[iMod].Data());
  cout << "Open replaced module " << iMod+1 << " (before u-arms) file " << fFileNameReplModBeforeUarms[iMod] << endl;

  return;
}

//_______________________________________________________________________________
int MainFrameModRepl::DoExtrapolationAndMerge()
{
  SetStyle();

  const int kNtab=4;
  TRootEmbeddedCanvas* fEcanvas[kNtab];
  TCanvas* fCanvas[kNtab];
  for(int iTab=0; iTab<kNtab; iTab++) {
    TGCompositeFrame* tb=fTab->GetTabContainer(iTab);
    tb->RemoveAll();
    fEcanvas[iTab] = new TRootEmbeddedCanvas(Form("efc%d",iTab), tb, 1200, 800);
    tb->AddFrame(fEcanvas[iTab], new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY,1, 1, 1, 1));
    fEcanvas[iTab]->GetCanvas()->SetBorderMode(0);
    fCanvas[iTab]=fEcanvas[iTab]->GetCanvas();
  }

  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();

  return ExtrapMarkerPosAfterHSReworking(fLayerType,fHSFlavour,fIsModRepl,fFileNameBeforeRewBeforeUarms,fFileNameAfterRewAfterUarms,fFileNameReplModBeforeUarms,fInitDir,fCanvas);
}

//_______________________________________________________________________________
//METHOD IMPLEMENTATIONS
void GetExtrapMarkerPosAfterHSReworking() {
  new MainFrameModRepl(gClient->GetRoot(), 200, 200);
}

//_____________________________________________________________________________________________
int ExtrapMarkerPosAfterHSReworking(int layer, int LeftOrRight, bool ismodulereplaced[7], TString file_before_rew, TString file_after_rew, TString file_after_rew_repl[7], TString initdir, TCanvas** fCanvas) {

  vector<double> x_before_rew, y_before_rew, z_before_rew;
  vector<double> x_after_rew, y_after_rew, z_after_rew;
  vector<double> x_after_rew_repl[7], y_after_rew_repl[7], z_after_rew_repl[7];
  vector<double> x_after_rew_merged, y_after_rew_merged, z_after_rew_merged;
  vector<double> x_nominal, y_nominal, z_nominal;

  bool read_before_rew = ReadDatFileMitutoyo(file_before_rew,x_before_rew,y_before_rew,z_before_rew);
  if(!read_before_rew) return 1;
  cout << "Loaded before reworking (before u-arms) modules" << endl;
  bool read_after_rew = ReadDatFileMitutoyo(file_after_rew,x_after_rew,y_after_rew,z_after_rew);
  if(!read_after_rew) return 2;
  cout << "Loaded after reworking (after u-arms) modules" << endl;
  for(int iMod=0; iMod<7; iMod++) {
    if(ismodulereplaced[iMod] && file_after_rew_repl[iMod]!="") {
      bool read_after_rew_repl = ReadDatFileMitutoyo(file_after_rew_repl[iMod],x_after_rew_repl[iMod],y_after_rew_repl[iMod],z_after_rew_repl[iMod]);
      if(!read_after_rew_repl) return 11+iMod;
      cout << "Loaded after reworking (before u-arms) replaced module" <<  iMod+1 << endl;
    }
  }

  MergeVectorMeasurements(x_after_rew_merged,y_after_rew_merged,z_after_rew_merged,x_after_rew,y_after_rew,z_after_rew,x_after_rew_repl,y_after_rew_repl,z_after_rew_repl,ismodulereplaced,layer);

  TString nominalposfilename = "";
  if(layer==MainFrameModRepl::kOL)
    nominalposfilename = Form("%s/OL_marker_nominal_positions_HS.dat",initdir.Data());
  else if(layer==MainFrameModRepl::kML)
    nominalposfilename = Form("%s/ML_marker_nominal_positions_HS.dat",initdir.Data());

  bool read_nominal_pos = ReadDatFile(nominalposfilename,x_nominal,y_nominal,z_nominal);

  vector<double> x_after_rew_filled, y_after_rew_filled, z_after_rew_filled;
  vector<double> x_before_rew_filled, y_before_rew_filled, z_before_rew_filled;
  FillVectorsOfMatchedPositions( x_after_rew_merged, y_after_rew_merged, z_after_rew_merged, x_nominal, y_nominal, x_after_rew_filled, y_after_rew_filled, z_after_rew_filled);
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

  TGraph2D* g_after_rew = new TGraph2D(x_after_rew_merged.size());
  g_after_rew->SetMarkerColor(kBlue);
  g_after_rew->SetMarkerStyle(kFullSquare);
  for(unsigned int iMarker=0; iMarker<x_after_rew_merged.size(); iMarker++) {
    g_after_rew->SetPoint(iMarker,x_after_rew_merged[iMarker],y_after_rew_merged[iMarker],z_after_rew_merged[iMarker]);
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
  for(unsigned int iMarker=0; iMarker<x_after_rew_merged.size(); iMarker++) {
    if(x_after_rew_merged[iMarker]>0) {
      g_after_rew_posX->SetPoint(iPosMarker,y_after_rew_merged[iMarker],z_after_rew_merged[iMarker]*1000);
      iPosMarker++;
    }
    else {
      g_after_rew_negX->SetPoint(iNegMarker,y_after_rew_merged[iMarker],z_after_rew_merged[iMarker]*1000);
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

  TH1F* h_beforeafter_res_X = new TH1F("h_beforeafter_res_X","after-before reworking marker position X;x_{after}-x_{before} (#mum);Entries",100,-50,50);
  TH1F* h_beforeafter_res_Y = new TH1F("h_beforeafter_res_Y","after-before reworking marker position Y;y_{after}-y_{before} (#mum);Entries",100,-50,50);
  TH1F* h_beforeafter_res_Z = new TH1F("h_beforeafter_res_Z","after-before reworking marker position Z;z_{after}-z_{before} (#mum);Entries",100,-50,50);
  h_beforeafter_res_X->SetLineWidth(2);
  h_beforeafter_res_Y->SetLineWidth(2);
  h_beforeafter_res_Z->SetLineWidth(2);
  for(unsigned int iMarker=0; iMarker<x_nominal.size(); iMarker++) {
    if(x_after_rew_filled[iMarker]!=-10000.) {
      if(x_before_rew_filled[iMarker]!=-10000.) h_beforeafter_res_X->Fill((x_after_rew_filled[iMarker]-x_before_rew_filled[iMarker])*1000);
    }
    else if(x_extrap[iMarker]!=-10000.) {
      if(x_before_rew_filled[iMarker]!=-10000.) h_beforeafter_res_X->Fill((x_extrap[iMarker]-x_before_rew_filled[iMarker])*1000);
    }
    if(y_after_rew_filled[iMarker]!=-10000.) {
      if(y_before_rew_filled[iMarker]!=-10000.) h_beforeafter_res_Y->Fill((y_after_rew_filled[iMarker]-y_before_rew_filled[iMarker])*1000);
    }
    else if(y_extrap[iMarker]!=-10000.) {
      if(y_before_rew_filled[iMarker]!=-10000.) h_beforeafter_res_Y->Fill((y_extrap[iMarker]-y_before_rew_filled[iMarker])*1000);
    }
    if(z_after_rew_filled[iMarker]!=-10000.) {
      if(z_before_rew_filled[iMarker]!=-10000.) h_beforeafter_res_Z->Fill((z_after_rew_filled[iMarker]-z_before_rew_filled[iMarker])*1000);
    }
    else if(z_extrap[iMarker]!=-10000.) {
      if(z_before_rew_filled[iMarker]!=-10000.) h_beforeafter_res_Z->Fill((z_extrap[iMarker]-z_before_rew_filled[iMarker])*1000);
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
  TH3F* hFrameHS = new TH3F("hFrame",";x (mm);y (mm);z (mm)",100,-50.,50.,100,-1000,1000,100,0.,1.);
  hFrameHS->SetStats(0);

  TLegend* leg_before_after = new TLegend(0.6,0.75,0.9,0.9);
  leg_before_after->SetTextSize(0.045);
  leg_before_after->AddEntry(g_before_rew,"Before reworking","p");
  leg_before_after->AddEntry(g_after_rew,"After reworking","p");

  TLegend* leg_pos_neg = new TLegend(0.4,0.75,0.9,0.9);
  leg_pos_neg->SetTextSize(0.05);
  leg_pos_neg->AddEntry(g_before_rew_posX,"X = 15 mm","p");
  leg_pos_neg->AddEntry(g_before_rew_negX,"X = -15 mm","p");

  TCanvas*& cHS3D = fCanvas[0];
  cHS3D->cd();
  hFrameHS->Draw();
  g_before_rew->Draw("Psame");
  g_after_rew->Draw("Psame");
  leg_before_after->Draw();
  cHS3D->Update();

  TCanvas*& cHS = fCanvas[1];
  cHS->cd();
  cHS->Divide(3,1);
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
  cHS->Update();

  TCanvas*& cResBeforeAfter = fCanvas[2];
  cResBeforeAfter->cd();
  cResBeforeAfter->Divide(3,1);
  cResBeforeAfter->cd(1);
  h_beforeafter_res_X->Draw();
  cResBeforeAfter->cd(2);
  h_beforeafter_res_Y->Draw();
  cResBeforeAfter->cd(3);
  h_beforeafter_res_Z->Draw();
  cResBeforeAfter->Update();

  TCanvas*& cExtrapQuality = fCanvas[3];
  cExtrapQuality->cd();
  cExtrapQuality->Divide(3,1);
  cExtrapQuality->cd(1);
  h_extrap_res_X->Draw();
  cExtrapQuality->cd(2);
  h_extrap_res_Y->Draw();
  cExtrapQuality->cd(3);
  h_extrap_res_Z->Draw();
  cExtrapQuality->Update();

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
  cResBeforeAfter->Print(Form("%s",outfile_pdf.Data()));
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
    if((LeftOrRight==MainFrameModRepl::kHSL && x_nominal[iMarker]>0) || (LeftOrRight==MainFrameModRepl::kHSR && x_nominal[iMarker]<0)) {
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
void MergeVectorMeasurements(vector<double>& x_merged, vector<double>& y_merged, vector<double>& z_merged, vector<double> x_all_afterUarms, vector<double> y_all_afterUarms, vector<double> z_all_afterUarms, vector<double> x_repl_beforeUarms[7], vector<double> y_repl_beforeUarms[7], vector<double> z_repl_beforeUarms[7], bool ismodulereplaced[7], int layer) {

  int nMod = 0;
  double ymins[7] = {-9999, -9999, -9999, -9999, -9999, -9999, -9999};
  double ymaxs[7] = {-9999, -9999, -9999, -9999, -9999, -9999, -9999};
  
  if(layer == MainFrameModRepl::kOL) {
    nMod = 7;

    ymins[0] = -738.575;
    ymins[1] = -527.475; 
    ymins[2] = -316.375; 
    ymins[3] = -105.275; 
    ymins[4] = 105.825; 
    ymins[5] = 316.925; 
    ymins[6] = 528.025;
    
    ymaxs[0] = -528.025; 
    ymaxs[1] = -316.925; 
    ymaxs[2] = -105.825; 
    ymaxs[3] = 105.275; 
    ymaxs[4] = 316.375; 
    ymaxs[5] = 527.475; 
    ymaxs[6] = 738.575;
  }
  else if(layer == MainFrameModRepl::kML) {
    nMod = 4;
    
    ymins[0] = -421.925;
    ymins[1] = -210.825; 
    ymins[2] = 0.275; 
    ymins[3] = 211.375; 
    
    ymaxs[0] = -211.375; 
    ymaxs[1] = -0.275; 
    ymaxs[2] = 210.825; 
    ymaxs[3] = 421.925; 
  }

  for(int iMod=0; iMod<nMod; iMod++) {
    if(ismodulereplaced[iMod]) {
      for(unsigned int iMarker=0; iMarker<y_repl_beforeUarms[iMod].size(); iMarker++) {
        if(y_repl_beforeUarms[iMod][iMarker] > ymins[iMod]-0.2 && y_repl_beforeUarms[iMod][iMarker] < ymaxs[iMod]+0.2) {
          x_merged.push_back(x_repl_beforeUarms[iMod][iMarker]);
          y_merged.push_back(y_repl_beforeUarms[iMod][iMarker]);
          z_merged.push_back(z_repl_beforeUarms[iMod][iMarker]);
        }  
      }
    }
    else {
      for(unsigned int iMarker=0; iMarker<y_all_afterUarms.size(); iMarker++) {
        if(y_all_afterUarms[iMarker] > ymins[iMod]-0.2 && y_all_afterUarms[iMarker] < ymaxs[iMod]+0.2) {
          x_merged.push_back(x_all_afterUarms[iMarker]);
          y_merged.push_back(y_all_afterUarms[iMarker]);
          z_merged.push_back(z_all_afterUarms[iMarker]);
        }  
      }
    }
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