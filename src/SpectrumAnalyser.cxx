/*
* File:         SpectrumAnalyser.cxx
* Author:       Aleksander Khreptak <aleksander.khreptak@lnf.infn.it>
* Created:      31 Mar 2023
* Last updated: 07 Apr 2023
*
* Description:
* This file implements a SpectrumAnalyser class that is used
* to analyse the SDD data from the SIDDHARTA-2 experiment
*/

#define SPECTRUM_ANALYSER_CXX
#include "../include/SpectrumAnalyser.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <algorithm> // for std::sort
#include <vector>
#include <utility> // for std::pair
#include <cmath>

#include <TBranch.h>
#include <TLeaf.h>
#include <TStyle.h>
#include <TPDF.h>
#include <TCanvas.h>
#include <TError.h>
#include <TSpectrum.h>
#include <TF1.h>

void SpectrumAnalyser::Loop()
{
  int column, row;
  bool ISGOODtime;
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  // The total runtime of the data acquisition process
  auto* date_leaf = static_cast<TLeaf*>(fChain->GetBranch("date")
                                              ->GetLeaf("date"));
  fChain->GetEntry(0);
  time_t time_start = static_cast<time_t>(date_leaf->GetValue(0));
  fChain->GetEntry(nentries-1);
  time_t time_end = static_cast<time_t>(date_leaf->GetValue(0));
  UInt_t time_diff = time_end - time_start;
  runtime += time_diff;
  std::cout << "Start time: " << ConvertTime(time_start) << std::endl;
  std::cout << "End time: " << ConvertTime(time_end) << std::endl;
  std::cout << "Runtime: " << runtime/3600 << " h, "
            << (runtime%3600)/60 << " m, " << runtime%60 << " s" << std::endl;
  
  // Loop over each entry in the TTree, skipping entries with zero hits
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if (jentry%1000 == 0) {
      std::cout << "Reading file entry " << jentry
                << "/" << nentries << std::endl;
    }

    // Skip the loop if nhits == 0
    if (nhits < 1) {
      std::cout << "Skip entry " << jentry<< " with " 
                << nhits << " hits" << std::endl;
      continue;
    }

    // Fill histograms
    for (Int_t ihit = 0; ihit < nhits; ++ihit) {
      h_adc_raw[bus[ihit]-1][sdd[ihit]-1]->Fill(adc[ihit]);
      ISGOODtime = true;
      if (ihit > 0) {
        ISGOODtime = CrossTalkTiming(drift[ihit], drift[ihit-1]);
      }
      if (ISGOODtime) {
        h_adc[bus[ihit]-1][sdd[ihit]-1]->Fill(adc[ihit]);
      } else {
        h_xtalk[bus[ihit]-1][sdd[ihit]-1]->Fill(adc[ihit]);
      }
      SDDHitMap(sdd[ihit], bus[ihit], column, row);
      h_sdd_map->SetBinContent(column, row, sdd[ihit]);
      h_sdd_rate->Fill(column, row);
      if (adc[ihit] > 1200) { h_sdd_rate_sig->Fill(column, row); }
      if (adc[ihit] <= 1200) { h_sdd_rate_noise->Fill(column, row); }
    }
  }

  h_sdd_rate->Scale(1.0/runtime);
  h_sdd_rate_sig->Scale(1.0/runtime);
  h_sdd_rate_noise->Scale(1.0/runtime);

  for (int i = 1; i <= 49; i++) {
    for (int j = 1; j < 9; j++) {
      if (h_sdd_rate->GetBinContent(i, j) < 0.001) {
        h_sdd_rate->SetBinContent(i, j, 0);
        h_sdd_rate_sig->SetBinContent(i, j, 0);
        h_sdd_rate_noise->SetBinContent(i, j, 0);
      }
    }
  }
  //// Peak Finder ////
  std::vector<std::string> lines = {"TiKb"};
  AddLines(lines);
  
  // Output the precalibration fit function
  //std::cout << "Precalibration fit function is: " << fit_function << std::endl;
  
  FindADCPeaks(1700,3700,8);  
}

std::string SpectrumAnalyser::ConvertTime(time_t t)
{
    std::stringstream ss;
    tm* timeinfo = localtime(&t);
    ss << std::setw(4) << std::setfill('0') << timeinfo->tm_year + 1900 << '-';
    ss << std::setw(2) << std::setfill('0') << timeinfo->tm_mon + 1 << '-';
    ss << std::setw(2) << std::setfill('0') << timeinfo->tm_mday << ' ';
    ss << std::setw(2) << std::setfill('0') << timeinfo->tm_hour -1 << ':';
    ss << std::setw(2) << std::setfill('0') << timeinfo->tm_min << ':';
    ss << std::setw(2) << std::setfill('0') << timeinfo->tm_sec;
    return ss.str();
}

void SpectrumAnalyser::SetPreCalibFitFunction(
    TString& fit_func, int n_peaks, const int kNumGaussParams, 
    const int kNumBkgParams, int* gauss_amp, 
    int* gauss_mean, int* gauss_sigma,
    int bkg_p0, int bkg_p1, int bkg_p2) 
{  
  // Set precalibration fit function: Gaussians + background
  fit_func.Clear();
  for (int i = 0; i < n_peaks; ++i) {
    gauss_amp[i] = kNumGaussParams * i;
    gauss_mean [i]  = kNumGaussParams * i + 1;
    gauss_sigma [i] = kNumGaussParams * i + 2;
    if (i != 0) { fit_func += "+"; }
    fit_func += TString::Format(
        "[%d]*exp(-0.5*pow(x-[%d],2)/pow([%d],2))/sqrt(2.0*%f)/[%d]",
        gauss_amp[i], gauss_mean[i], gauss_sigma[i], PI, gauss_sigma[i]);
  }

  bkg_p0 = kNumGaussParams * n_peaks;
  bkg_p1 = kNumGaussParams * n_peaks + 1;
  bkg_p2 = kNumGaussParams * n_peaks + 2;
  fit_func += TString::Format("+exp([%d]+[%d]*x)+[%d]", 
                            bkg_p0, bkg_p1, bkg_p2);
}

void SpectrumAnalyser::FindADCPeaks(
    const float& x_min, const float& x_max, const int& factor)
{
  const float energy_ratio = (triad_peak_energies[2] - triad_peak_energies[1]) /
                             (triad_peak_energies[1] - triad_peak_energies[0]);

  std::cout << "energy_ratio = " << energy_ratio << std::endl;

  float init_threshold = 0.01; // initial threshold parameter for the Peak Finder (std in TSpectrum is 0.05)
  float init_tolerance = 0.05; // tolerance to check that the peak assumption is correct (5%)
  
  const int sigma_pf = 20 / factor;  // sigma for the Peak Finder (in units of ADC channels)

  const int kMinStats = 1000;

  TSpectrum spectrum{kMaxNumPeaksPF};

  for(int bus_idx = 0; bus_idx < num_buses; ++bus_idx) {
    for(int sdd_idx = 0; sdd_idx < num_sdds; ++sdd_idx) {
      if (!h_adc[bus_idx][sdd_idx]) continue;

      auto histo = static_cast<TH1D*>(h_adc[bus_idx][sdd_idx]->Clone("histo"));

      std::cout << std::endl << "--- STARTING PEAK FINDER ---" << std::endl;
      std::cout << "for [bus#" << bus_idx + 1 << ", sdd#" << sdd_idx + 1 
                << "] with a statistics of " << histo->GetEntries() << std::endl 
                << std::endl;
      
      if (histo->GetEntries() < kMinStats) {
        std::cout << "NOT ENOUGH STATISTICS (" << histo->GetEntries() 
                  << "). SKIPPING..." << std::endl;
        continue;
      }

      histo->GetXaxis()->SetRangeUser(x_min, x_max);

      std::vector<std::pair<float, float>> peak_pos;

      FindPeakCandidates(histo, sigma_pf, peak_pos, spectrum, init_threshold);

      if (peak_pos.empty()) {
        std::cout << "PEAK FINDER COULD NOT FIND ANY PEAKS. SKIPPING..." 
                  << std::endl;
        continue;
      }

      //int num_found = peak_pos.size();

      std::cout << peak_pos.size() << " candidate peaks have been identified:" 
                << std::endl;
      
      // Sort the peaks based on their x positions
      std::sort(peak_pos.begin(), peak_pos.end(),
          [](const auto& a, const auto& b) { return a.first < b.first; });

      // Print found peaks
      for (int i = 0; i < peak_pos.size(); ++i) {
        std::cout << "- peak#" << i + 1 << ": position = "
                  << peak_pos[i].first << " [ADC], intensity = "
                  << peak_pos[i].second << std::endl;
      }

      float slope = 0.0;
      float offset = 0.0;
      float highest_peak = 0.0;
      std::cout << "Checking peak triads for optimal calibration..." << std::endl;
      GetCalibrationParameters(
          peak_pos, peak_pos.size(), energy_ratio, slope, offset, highest_peak);

      FitSpectrum(histo, peak_energies, slope, offset, peak_energies.size(),
          highest_peak, bus_idx, sdd_idx, "output_file");

      histo->Delete();
    }
  }
}

void SpectrumAnalyser::FindPeakCandidates(
    TH1D* histo, const int sigma_pf, std::vector<std::pair<float, float>>& peak_pos,
    TSpectrum& spectrum, const float init_threshold)
{
  spectrum.Search(histo, sigma_pf, "", init_threshold);

  int num_tries = 0;
  int num_found = 0;
  float peak_threshold = init_threshold;

  while(num_found < kNumPeaksPF && num_tries < 10) {
    num_found = spectrum.Search(histo, sigma_pf, "", peak_threshold);
    //num_found = spectrum.GetNPeaks();
    peak_threshold *= 0.1;  // Initial = 0.01. It changes until it finds peaks
    num_tries++;
  }

  peak_pos.reserve(num_found);
  double* peak_x = spectrum.GetPositionX();
  double* peak_y = spectrum.GetPositionY();

  for (int i = 0; i < num_found; ++i) {
    peak_pos.emplace_back(static_cast<float>(peak_x[i]), static_cast<float>(peak_y[i]));
  }

}

void SpectrumAnalyser::GetCalibrationParameters(
    const std::vector<std::pair<float, float>>& peak_pos,
    const int& num_found, const float& energy_ratio,
    float& slope, float& offset, float& highest_peak) 
{
  // Make a triad from all the peaks found:
  int peak_idx[3] = {-1, -1, -1};  
 
  for (int i0 = 0; i0 < num_found && peak_idx[2] == -1; i0++) {
    for (int i1 = i0 + 1; i1 < num_found && peak_idx[2] == -1; i1++) {
      for (int i2 = i1 + 1; i2 < num_found && peak_idx[2] == -1; i2++) {        
        float adc_ratio = (peak_pos[i2].first - peak_pos[i1].first) / 
                          (peak_pos[i1].first - peak_pos[i0].first);        

        float peak_tolerance = 0.05;  // 5%
        bool tolerance_pass = fabs(1.0 - (energy_ratio / adc_ratio)) 
                            < peak_tolerance;
        if (!tolerance_pass) continue;

        // Get the Peak Finder calibration offset and slope
        float adc_diff = peak_pos[i0].first - peak_pos[i1].first;
        float energy_diff = triad_peak_energies[0] - triad_peak_energies[1];
        slope  = energy_diff / adc_diff;
        offset = -1.0 * peak_pos[i0].first * slope + triad_peak_energies[0];
        
        // define an acceptable slope and offset
        float min_slope = 2.9, max_slope = 3.9;
        float min_offset = -3000, max_offset = -1000;
        if ((slope > max_slope) || (slope < min_slope) || 
            (offset > max_offset) || (offset < min_offset)) {
            continue;
        }        
        peak_idx[0] = i0;
        peak_idx[1] = i1;
        peak_idx[2] = i2;
      }
    }
  }

  if (peak_idx[2] == -1) {
    std::cout << "No valid triad has been found. Unable to calculate calibration parameters." 
              << std::endl;
    return;
  }

  std::cout << "Optimal peak triad has been found: " << peak_idx[0] + 1 << " "
            << peak_idx[1] + 1 << " " << peak_idx[2] + 1 << std::endl;
  std::cout << "Calculating calibration parameters..." << std::endl;
  std::cout << "Offset = " << offset << ", slope = " << slope << std::endl;

  // Find also the highest height among the selected ones  
  if (peak_idx[0] > -1) {
    if (peak_pos[peak_idx[0]].second > highest_peak)
      highest_peak = peak_pos[peak_idx[0]].second;
    if (peak_pos[peak_idx[1]].second > highest_peak)
      highest_peak = peak_pos[peak_idx[1]].second;
    if (peak_pos[peak_idx[2]].second > highest_peak)
      highest_peak = peak_pos[peak_idx[2]].second;
  }
}

void SpectrumAnalyser::FitSpectrum(
    TH1D* histo, const std::vector<double> peak_energies,
    const float slope, const float offset, const int num_peaks,
    const float highest_peak, const int bus_idx, const int sdd_idx,
    const std::string& filename) 
{
  float sigma = 20.f;   // Start with an initial peak resolution in ADC of 20
  std::vector<float> x_min_pre(num_peaks), x_max_pre(num_peaks);

  // Calculate fit limits using a range-based loop
  for (int i = 0; i < num_peaks; ++i) {   
    const auto x_min = (peak_energies[i] - offset) / slope - 2.5f * sigma;
    const auto x_max = (peak_energies[i] - offset) / slope + 3.0f * sigma;
    x_min_pre[i] = x_min;
    x_max_pre[i] = x_max;
  }

  // Set bin errors to zero outside the fit limits using a range-based loop
  for (int bin_idx = 1; bin_idx <= histo->GetNbinsX(); ++bin_idx) {
    bool keep = false;
    for (int i = 0; i < num_peaks; ++i) {
      if (histo->GetBinCenter(bin_idx) > x_min_pre[i] 
          && histo->GetBinCenter(bin_idx) < x_max_pre[i]) {
        keep = true;
        break;
      }
    }
    if (!keep) histo->SetBinError(bin_idx, 0.);
  }

  // Define Gaussian and background functions
  const int kNumGaussParams = 3;  // Number of parameters for the Gaussian function
  const int kNumBkgParams   = 3;  // Number of parameters for the background: p0 + exp(p1 + p2*x)
  int gauss_amp[kMaxNumPeaksPF]   = {-1};   // Gaussian parameters
  int gauss_mean[kMaxNumPeaksPF]  = {-1};
  int gauss_sigma[kMaxNumPeaksPF] = {-1};
  int bkg_p0 = -1;  // Background parameters
  int bkg_p1 = -1;
  int bkg_p2 = -1;

  TString fit_function = "";
  SetPreCalibFitFunction(fit_function, num_peaks, kNumGaussParams, kNumBkgParams,
                         gauss_amp, gauss_mean, gauss_sigma, 
                         bkg_p0, bkg_p1, bkg_p2);

  auto fit_func = std::make_unique<TF1>(
      Form("fit_func_bus%d_sdd%d", bus_idx, sdd_idx),
      fit_function, x_min_pre[0], x_max_pre[num_peaks - 1]);

  for (int i = 0; i < num_peaks; ++i) {  //set initial parameters for Gaussian function
    fit_func->SetParameter(gauss_amp[i], highest_peak * sigma);
    fit_func->SetParName(gauss_amp[i], Form("gauss%i_Amp", i));
    fit_func->SetParameter(gauss_mean[i], (peak_energies[i] - offset) / slope);
    fit_func->SetParName(gauss_mean[i], Form("gauss%i_Mea", i));
    fit_func->SetParameter(gauss_sigma[i], sigma);
    fit_func->SetParName(gauss_sigma[i], Form("gauss%i_Sig", i));
  }

  //set parameters for background function
  TF1* bkg_func = new TF1(Form("bkg_func_bus%d_sdd%d", bus_idx, sdd_idx), 
                          "expo(0) + pol0(2)", 1500, 4500);
  histo->Fit(Form("bkg_func_bus%d_sdd%d", bus_idx, sdd_idx), "", "", 1500, 4500);

  fit_func->SetParameter(num_peaks*3, bkg_func->GetParameter(0));
  fit_func->SetParameter(num_peaks*3+1, bkg_func->GetParameter(1));
  fit_func->SetParameter(num_peaks*3+2, bkg_func->GetParameter(2));

  histo->Fit(fit_func.get(), "R", "", x_min_pre[0], x_max_pre[num_peaks-1]);

  //put back original errors:
  for (int bin_idx = 1; bin_idx <= histo->GetNbinsX(); ++bin_idx) {
    histo->SetBinError(
        bin_idx, h_adc[bus_idx][sdd_idx]->GetBinError(bin_idx));
  }

  std::ofstream  file;
  file.open("output/files/peaks_" + filename + ".dat", std::ios::app);
  file << Form("%d  %d  %f  %f  %f  %f", bus_idx + 1, sdd_idx + 1, 
      fit_func->GetParameter(1),
      fit_func->GetParError(1),
      fit_func->GetParameter(7),
      fit_func->GetParError(7)) << std::endl;
  file.close();
}


bool SpectrumAnalyser::CrossTalkTiming(Short_t drift, Short_t drift_pre) 
{
  bool ISGOODtime = true;
  int t1 = 0, t2 = 0, timediff = 0;
  t1 = drift + 32768;
  t2 = drift_pre + 32768.;
  if (t1 > t2) { timediff = t1 - t2; }
  if (t2 > t1) { timediff = t1 + (32768.*2) - t2; }
  if ((timediff > 0. && timediff < 625.)) { ISGOODtime = false; } // 625 -> 5 microseconds
  return ISGOODtime;
}

void SpectrumAnalyser::SDDHitMap(
    int sddnumber, int busnumber, int &column, int &row) 
{
  row = 0; column = 0;
  int SFERA = 0;

  // Back side view
  if (sddnumber <= 16) { sddnumber = sddnumber; SFERA = 0; }
  if (sddnumber >= 17 && sddnumber <= 32) { sddnumber = sddnumber - 16; SFERA = 2; }
  if (sddnumber >= 33 && sddnumber <= 48) { sddnumber = sddnumber - (16*2); SFERA = 4; }
  if (sddnumber >= 49 && sddnumber <= 64) { sddnumber = sddnumber - (16*3); SFERA = 6; }

  if (sddnumber == 1)   { row = 4; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 2)   { row = 3; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 3)   { row = 2; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 4)   { row = 1; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 5)   { row = 1; column = 2 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 6)   { row = 2; column = 2 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 7)   { row = 3; column = 2 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 8)   { row = 4; column = 2 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 9)   { row = 8; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 10)  { row = 7; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 11)  { row = 6; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 12)  { row = 5; column = 1 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 13)  { row = 5; column = 2 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 14)  { row = 6; column = 2 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 15)  { row = 7; column = 2 + SFERA + ((busnumber - 1)*8); }
  if (sddnumber == 16)  { row = 8; column = 2 + SFERA + ((busnumber - 1)*8); }
}

int SpectrumAnalyser::SFERAnumber(int sdd) 
{
  int SFERA = 0;
  if (sdd <= 16) { SFERA = 1; }
  if (sdd >= 17 && sdd <= 32) { SFERA = 2; }
  if (sdd >= 33 && sdd <= 48) { SFERA = 3; }
  if (sdd >= 49 && sdd <= 64) { SFERA = 4; }
  return SFERA;
}
