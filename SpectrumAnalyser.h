//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 31 15:30:10 2023 by ROOT version 6.26/10
// from TTree ft10/SIDDHARTA-2 TTree
// found on file: 20220611_0024_0611_0315_xray_25kv_50ua_tube1_cal.root
//
// Modified Fri Mar 31 2023 by Aleksander Khreptak
//////////////////////////////////////////////////////////

#ifndef SPECTRUM_ANALYSER_H
#define SPECTRUM_ANALYSER_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

// Header file for the class stored in the TTree

class SpectrumAnalyser {
public :
  TTree          *fChain;     // pointer to the analyzed TTree or TChain
  Int_t           fCurrent;   // current Tree number in a TChain

  // Fixed size dimensions of arrays
  static const Int_t MAXHITS  = 40000;
  static const Int_t MAXEDGES = 280;

  // Declaration of leaf types
  Int_t           buf;
  Short_t         kt[4];
  Int_t           nhits;
  Int_t           bus[MAXHITS];     //[nhits]
  UShort_t        evnr[MAXHITS];    //[nhits]
  UShort_t        ht[MAXHITS];      //[nhits]
  UShort_t        trigg[MAXHITS];   //[nhits]
  UShort_t        sdd[MAXHITS];     //[nhits]
  Short_t         adc[MAXHITS];     //[nhits]
  Short_t         drift[MAXHITS];   //[nhits]
  Int_t           date;
  Short_t         ie;
  Short_t         ip;
  Short_t         dum;
  Int_t           edges_veto;
  UShort_t        v_edge[MAXEDGES]; //[edges_veto]
  UShort_t        v_val[MAXEDGES];  //[edges_veto]
  UShort_t        v_ch[MAXEDGES];   //[edges_veto]
  Int_t           v_found;
  UShort_t        v_ch_r[100];      //[v_found]
  Int_t           v_tdc[100];       //[v_found]
  Int_t           v_tot[100];       //[v_found]

  // List of branches
  TBranch        *b_buf;
  TBranch        *b_kt;
  TBranch        *b_nhits;
  TBranch        *b_bus;
  TBranch        *b_evnr;
  TBranch        *b_ht;
  TBranch        *b_trigg;
  TBranch        *b_sdd;
  TBranch        *b_adc;
  TBranch        *b_drift;
  TBranch        *b_date;
  TBranch        *b_ie;
  TBranch        *b_ip;
  TBranch        *b_dum;
  TBranch        *b_edges_veto;
  TBranch        *b_v_edge;
  TBranch        *b_v_val;
  TBranch        *b_v_ch;
  TBranch        *b_v_found;
  TBranch        *b_v_ch_r;
  TBranch        *b_v_tdc;
  TBranch        *b_v_tot;

  UInt_t runtime = 0;

  static const int NUM_BUSES  = 6;    // Number of buses
  static const int NUM_SDDS   = 64;   // Number of SDDs
  
  // Binning for the ADC
  int num_adc_bins = 10000;
  const int MIN_ADC = 0;
  const int MAX_ADC = 10000;

  enum class CanvasFormat { Default, A4, A5 };
  enum class CanvasOrientation { Default, Landscape, Portrait };

  SpectrumAnalyser(TTree *tree = 0);
  virtual ~SpectrumAnalyser();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  void             initHistograms(int rebin_factor = 1);
  void             writeHistograms(const std::string& filename);
  void             drawADCSpectra(const std::string& filename);
  void             drawSDDMap(const std::string& filename);
    
private:
  // Histograms
  TH1D *h_adc[NUM_BUSES][NUM_SDDS], *h_adc_raw[NUM_BUSES][NUM_SDDS];
  TH1D *h_xtalk[NUM_BUSES][NUM_SDDS];
  TH2D *h_sdd_map, *h_sdd_rate, *h_sdd_rate_sig, *h_sdd_rate_noise;
  
  std::string convertTime(time_t t);
  bool        crossTalkTiming(Short_t drift, Short_t drift_pre);
  void        sddHitMap(int sddnumber, int busnumber, int &column, int &row);
  int         SFERAnumber(int sdd);
  TStyle*     setHistogramStyle();
  TCanvas*    createCanvas(
                  CanvasFormat format = CanvasFormat::Default,
                  CanvasOrientation orientation = CanvasOrientation::Default,
                  int width = 800, int height = 600);
  void        drawSpectrum(TCanvas* canvas);
  void        draw2DHistogram(TH2D* hist, const std::string& title);
  
};

#endif  // SPECTRUM_ANALYSER_H

#ifdef SPECTRUM_ANALYSER_CXX
SpectrumAnalyser::SpectrumAnalyser(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("${XRAYDATA}/20220611/20220611_0024_0611_0315_xray_25kv_50ua_tube1_cal.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("${XRAYDATA}/20220611/20220611_0024_0611_0315_xray_25kv_50ua_tube1_cal.root");
    }
    f->GetObject("ft10", tree);

  }
  Init(tree);
}

SpectrumAnalyser::~SpectrumAnalyser()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t SpectrumAnalyser::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t SpectrumAnalyser::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void SpectrumAnalyser::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("buf", &buf, &b_buf);
  fChain->SetBranchAddress("kt", kt, &b_kt);
  fChain->SetBranchAddress("nhits", &nhits, &b_nhits);
  fChain->SetBranchAddress("bus", bus, &b_bus);
  fChain->SetBranchAddress("evnr", evnr, &b_evnr);
  fChain->SetBranchAddress("ht", ht, &b_ht);
  fChain->SetBranchAddress("trigg", trigg, &b_trigg);
  fChain->SetBranchAddress("sdd", sdd, &b_sdd);
  fChain->SetBranchAddress("adc", adc, &b_adc);
  fChain->SetBranchAddress("drift", drift, &b_drift);
  fChain->SetBranchAddress("date", &date, &b_date);
  fChain->SetBranchAddress("ie", &ie, &b_ie);
  fChain->SetBranchAddress("ip", &ip, &b_ip);
  fChain->SetBranchAddress("dum", &dum, &b_dum);
  fChain->SetBranchAddress("edges_veto", &edges_veto, &b_edges_veto);
  fChain->SetBranchAddress("v_edge", v_edge, &b_v_edge);
  fChain->SetBranchAddress("v_val", v_val, &b_v_val);
  fChain->SetBranchAddress("v_ch", v_ch, &b_v_ch);
  fChain->SetBranchAddress("v_found", &v_found, &b_v_found);
  fChain->SetBranchAddress("v_ch_r", v_ch_r, &b_v_ch_r);
  fChain->SetBranchAddress("v_tdc", v_tdc, &b_v_tdc);
  fChain->SetBranchAddress("v_tot", v_tot, &b_v_tot);
  Notify();
}

Bool_t SpectrumAnalyser::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void SpectrumAnalyser::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t SpectrumAnalyser::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  return 1;
}
#endif  // SPECTRUM_ANALYSER_CXX
