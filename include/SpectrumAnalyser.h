//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 31 15:30:10 2023 by ROOT version 6.26/10
// from TTree ft10/SIDDHARTA-2 TTree
// found on file: 20220611_0024_0611_0315_xray_25kv_50ua_tube1_cal.root
//
// Modified on Fri Mar 31 2023 by Aleksander Khreptak
// Last updated on Tue Apr 11 2023
//////////////////////////////////////////////////////////

#ifndef SPECTRUM_ANALYSER_H
#define SPECTRUM_ANALYSER_H

#include <iostream>
#include <map>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpectrum.h>

// Header file for the class stored in the TTree

class SpectrumAnalyser {
public :
  TTree          *fChain;     // pointer to the analyzed TTree or TChain
  Int_t           fCurrent;   // current Tree number in a TChain

  // Fixed size dimensions of arrays
  static const Int_t MAX_HITS  = 40000;
  static const Int_t MAX_EDGES = 280;

  // Declaration of leaf types
  Int_t           buf;
  Short_t         kt[4];
  Int_t           nhits;
  Int_t           bus[MAX_HITS];     //[nhits]
  UShort_t        evnr[MAX_HITS];    //[nhits]
  UShort_t        ht[MAX_HITS];      //[nhits]
  UShort_t        trigg[MAX_HITS];   //[nhits]
  UShort_t        sdd[MAX_HITS];     //[nhits]
  Short_t         adc[MAX_HITS];     //[nhits]
  Short_t         drift[MAX_HITS];   //[nhits]
  Int_t           date;
  Short_t         ie;
  Short_t         ip;
  Short_t         dum;
  Int_t           edges_veto;
  UShort_t        v_edge[MAX_EDGES]; //[edges_veto]
  UShort_t        v_val[MAX_EDGES];  //[edges_veto]
  UShort_t        v_ch[MAX_EDGES];   //[edges_veto]
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

  const double PI = 3.14159265358979323846;

  UInt_t runtime = 0;

  static const int num_buses  = 6;    // Number of buses
  static const int num_sdds   = 64;   // Number of SDDs per bus
  
  // Binning for the ADC
  const int num_adc_bins  = 10000;
  const int adc_min = 0;
  const int adc_max = 10000;

  const int rebin_factor = 4;

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

  void             InitHistograms();
  void             WriteHistograms(const std::string& filename);
  void             DrawADCSpectra(const std::string& filename);
  void             DrawSDDMap(const std::string& filename);
  void             FindADCPeaks(
                      const float& x_min, const float& x_max,
                      const std::string& filename);

private:
  // Histograms
  TH1D *h_adc[num_buses][num_sdds], *h_adc_raw[num_buses][num_sdds];
  TH1D *h_xtalk[num_buses][num_sdds];
  TH2D *h_sdd_map, *h_sdd_rate, *h_sdd_rate_sig, *h_sdd_rate_noise;

  static std::map<std::string, double> xray_line_energies;

  const std::vector<double> triad_peak_energies = {xray_line_energies["TiKa"],
                                                   xray_line_energies["CuKa"],
                                                   xray_line_energies["CuKb"]};
  const std::vector<std::string> triad_peak_names = {"TiKa", "CuKa", "CuKb"};

  std::vector<double> peak_energies;
  std::vector<std::string> peak_names;
  int num_peaks = 0;

  std::pair<std::vector<std::string>, std::vector<double>>
    AddLines(const std::vector<std::string>& lines);

  const int kMaxPeaks = 20;
  const int kMaxNumPeaksPF  = 15; // Maximum number of peaks that the Peak Finder can detect
  const int kNumPeaksPF     = 3;  // Desired number of peaks to find
  
  std::string ConvertTime(time_t t);
  bool        CrossTalkTiming(Short_t drift, Short_t drift_pre);
  void        SDDHitMap(int sddnumber, int busnumber, int &column, int &row);
  int         SFERAnumber(int sdd);
  TStyle*     SetHistogramStyle();
  TCanvas*    CreateCanvas(
                  CanvasFormat format = CanvasFormat::Default,
                  CanvasOrientation orientation = CanvasOrientation::Default,
                  int width = 800, int height = 600);
  void        DrawSpectrum(TCanvas* canvas);
  void        Draw2DHistogram(TH2D* hist, const std::string& title);

  void        SetPreCalibFitFunction(
                  TString& fit_func, int n_peaks, const int kNumGaussParams,
                  const int kNumBkgParams, int* gauss_amp,
                  int* gauss_mean, int* gauss_sigma,
                  int bkg_p0, int bkg_p1, int bkg_p2);  
  void        FindPeakCandidates(
                  TH1D* histo, const int sigma_pf, 
                  std::vector<std::pair<float, float>>& peak_pos,
                  TSpectrum& spectrum, const float init_threshold);
  void        GetCalibrationParameters(
                  const std::vector<std::pair<float, float>>& peak_pos,
                  const int& num_found, const float& energy_ratio,
                  float& slope, float& offset, float& highest_peak);
  void        FitSpectrum(
                  TH1D* histo, const std::vector<double> peak_energies,
                  const float slope, const float offset, const int num_peaks,
                  const float highest_peak, const int bus_idx, const int sdd_idx,
                  const std::string& filename);
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

std::map<std::string, double> SpectrumAnalyser::xray_line_energies = 
{
  {"TiKa",  4508.83},
  {"TiKb",  4931.81},
  {"CuKa",  8041.05},
  {"CuKb",  8905.29},
  {"MnKa",  5895.23},
  {"FeKa",  6399.47},
  {"FeKb",  7058.0},
  {"WKa",   8631.10},
  {"WKg",   11285.9},
  {"BrKa",  11908.26},
  {"BrKb",  13292.0},
  {"BiKa",  10828.1},
  {"BiKb",  13011.6},
  {"BiKg",  15247.7},
  {"SrKa",  14142.04},
  {"SrKb",  15836.0},
  {"PbKa",  10541.39},
  {"PbKb",  12616.14},
  {"PbKg",  14764.4}
};

std::pair<std::vector<std::string>, std::vector<double>> SpectrumAnalyser::
    AddLines(const std::vector<std::string>& lines)
{
  std::vector<double> add_peak_energies = triad_peak_energies;
  std::vector<std::string> add_peak_names = triad_peak_names;

  // Add lines to names and energies vectors
  for (const auto& line : lines) {
    const auto& iter = xray_line_energies.find(line);
    if (iter != xray_line_energies.end()) {
      // If the line exists in the energies map,
      // add it to peak names and peak energies vectors
      add_peak_names.push_back(line);
      add_peak_energies.push_back(iter->second);
    } else {
      // If the line does not exist in the energies map, print an error message
      std::cerr << "ERROR: Line \"" << line
                << "\" does not exist in energies map!" << std::endl;
      return std::make_pair(std::vector<std::string>(), std::vector<double>());
    }
  }

  // Sort the lines by energy
  std::vector<std::pair<double, std::string>> sorted_lines;
  for (std::size_t i = 0; i < add_peak_energies.size(); ++i) {
    sorted_lines.emplace_back(add_peak_energies[i], add_peak_names[i]);
  }
  std::sort(sorted_lines.begin(), sorted_lines.end());

  // Clear the previous peak data
  peak_energies.clear();
  peak_names.clear();

  // Add the sorted lines to the peak data
  for (const auto& pair : sorted_lines) {
    peak_energies.push_back(pair.first);
    peak_names.push_back(pair.second);
  }

  num_peaks = peak_energies.size();

  // Output the sorted peaks to the terminal
  std::cout << "Number of lines: " << num_peaks << std::endl;
  for (Int_t i = 0; i < num_peaks; ++i) {
    std::cout << "Line#" << i+1 << ": line " << peak_names[i]
              << "; energy = " << peak_energies[i] << std::endl;
  }

  return { peak_names, peak_energies };
}
#endif  // SPECTRUM_ANALYSER_CXX
