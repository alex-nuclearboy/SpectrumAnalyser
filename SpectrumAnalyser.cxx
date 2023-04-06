/*
* File:         SpectrumAnalyser.cxx
* Author:       Aleksander Khreptak <aleksander.khreptak@lnf.infn.it>
* Created:      31 Mar 2023
* Last updated: 05 Apr 2023
*
* Description:
* This file implements a SpectrumAnalyser class that is used
* to analyse the SDD data from the SIDDHARTA-2 experiment
*/

#define SPECTRUM_ANALYSER_CXX
#include "SpectrumAnalyser.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <algorithm> // for std::sort
#include <vector>
#include <utility> // for std::pair

#include <TBranch.h>
#include <TLeaf.h>
#include <TStyle.h>
#include <TPDF.h>
#include <TCanvas.h>
#include <TError.h>
#include <TSpectrum.h>

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
  TString fit_function = "";
  SetPreCalibFitFunction(fit_function,num_peaks);
  // Output the precalibration fit function
  std::cout << "Precalibration fit function is: " << fit_function << std::endl;
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

void SpectrumAnalyser::InitHistograms(int rebin_factor) 
{
  num_adc_bins /= rebin_factor;   // divide number of bins by rebin factor
  //factor = rebin_factor;
  
  // Loop over all the buses and SDDs to create the histograms
  for (int bus_idx = 0; bus_idx < kNumBuses; ++bus_idx) {
    for (int sdd_idx = 0; sdd_idx < kNumSDDs; ++sdd_idx) {

      // ADC histograms
      h_adc[bus_idx][sdd_idx] = new TH1D(
          Form("h_adc_bus%i_sdd%i", bus_idx + 1, sdd_idx + 1),
          Form("h_adc_bus%i_sdd%i", bus_idx + 1, sdd_idx + 1),
          num_adc_bins, kMinADC, kMaxADC
      );
      h_adc[bus_idx][sdd_idx]->GetXaxis()->SetTitle("ADC [channel]");
      h_adc[bus_idx][sdd_idx]->GetYaxis()->SetTitle(
          Form("counts / %i channel", rebin_factor));
      h_adc_raw[bus_idx][sdd_idx] = static_cast<TH1D*>(
          h_adc[bus_idx][sdd_idx]->Clone(
              Form("h_adc_raw_bus%i_sdd%i", bus_idx + 1, sdd_idx + 1)));

      // Crosstalk histogram
      h_xtalk[bus_idx][sdd_idx] = new TH1D(
          Form("h_xtalk_bus%i_sdd%i", bus_idx + 1, sdd_idx + 1),
          Form("h_xtalk_bus%i_sdd%i", bus_idx + 1, sdd_idx + 1),
          num_adc_bins, kMinADC, kMaxADC);
    }
  }

  // Create the SDD maps and rate histograms
  h_sdd_map = new TH2D("h_sdd_map", "h_sdd_map", 48, 1, 49, 8, 1, 9);
  h_sdd_rate = new TH2D("h_sdd_rate", "h_sdd_rate", 48, 1, 49, 8, 1, 9);
  h_sdd_rate_sig = new TH2D(
      "h_sdd_rate_sig", "h_sdd_rate_sig", 48, 1, 49, 8, 1, 9);
  h_sdd_rate_noise = new TH2D(
      "h_sdd_rate_noise", "h_sdd_rate_noise", 48, 1, 49, 8, 1, 9);
}

void SpectrumAnalyser::WriteHistograms(const std::string& filename)
{
  // Create the output file
  TFile* output_file = new TFile(filename.c_str(), "RECREATE");
  output_file->cd();
  // Write the histograms to the file
  for (int bus_idx = 0; bus_idx < kNumBuses; ++bus_idx) {
    for (int sdd_idx = 0; sdd_idx < kNumSDDs; ++sdd_idx) {
      if (h_adc[bus_idx][sdd_idx]->GetEntries() > 0) {
        h_adc_raw[bus_idx][sdd_idx]->Write();
        h_adc[bus_idx][sdd_idx]->Write();
        h_xtalk[bus_idx][sdd_idx]->Write();
      }
    }
  }
  h_sdd_map->Write();
  h_sdd_rate->Write();
  h_sdd_rate_sig->Write();
  h_sdd_rate_noise->Write();
  // Save and close the output file
  output_file->Save();
  output_file->Close();
  // Print output file name
  std::cout << "Output file: " << filename << std::endl;
}

TStyle* SpectrumAnalyser::SetHistogramStyle()
{
  // Global histogram style settings
  static TStyle* style = []{
    auto style = new TStyle("style", "");
    style->SetOptStat(kFALSE);
    style->SetPadGridX(kTRUE);
    style->SetPadGridY(kTRUE);
    style->SetPadLeftMargin(0.12);
    style->SetPadRightMargin(0.12);
    style->SetPadBottomMargin(0.12);
    style->SetPadTopMargin(0.1);
    style->SetTitleSize(0.05,"XY");
    style->SetTitleOffset(1,"X");
    style->SetTitleOffset(1,"Y");
    style->SetLabelSize(0.045,"XY");
    style->SetTitleFont(42,"XYZ");
    style->SetLabelFont(42,"XYZ");
    style->SetTextFont(42);
    gROOT->SetStyle("style");
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE);
    return style;
  }();
  return style;
}

TCanvas* SpectrumAnalyser::CreateCanvas(
    CanvasFormat format, CanvasOrientation orientation, 
    int width, int height)
{
  SetHistogramStyle();
  int canvas_width, canvas_height;
  switch (format) {
  case CanvasFormat::A4:
    canvas_width = (orientation == CanvasOrientation::Landscape) ? 2339 : 1654;
    canvas_height = (orientation == CanvasOrientation::Landscape) ? 1654 : 2339;
    break;
  case CanvasFormat::A5:
    canvas_width = (orientation == CanvasOrientation::Landscape) ? 1654 : 1169;
    canvas_height = (orientation == CanvasOrientation::Landscape) ? 1169 : 1654;
    break;
  case CanvasFormat::Default:
    canvas_width = width;
    canvas_height = height;
    break;
  default:
    std::cerr << "Invalid canvas format" << std::endl;
    return nullptr;
  }
  auto canvas = std::unique_ptr<TCanvas>(new TCanvas("canvas", "canvas", 
                                         canvas_width, canvas_height));  
  return canvas.release();
}

void SpectrumAnalyser::DrawADCSpectra(const std::string& filename) 
{   
  auto canvas = CreateCanvas(CanvasFormat::A4, CanvasOrientation::Portrait);
  auto pdf = std::make_unique<TPDF>(
      Form("output/plots/hADC_spectra_%s.pdf", filename.c_str()));
  DrawSpectrum(canvas);
  pdf->Close();
  canvas->Close();
  std::cout << "ADC spectra have been saved to output/plots/hADC_spectra_" 
            << filename.c_str() << ".pdf" << std::endl;
}

void SpectrumAnalyser::DrawSpectrum(TCanvas* canvas)
{
  const int NUM_COLUMNS = 2;
  const int NUM_ROWS    = 4;
  const int kMinX = 1500;
  const int kMaxX = 4500;
  
  for (int bus_idx = 0; bus_idx < kNumBuses; ++bus_idx) {
    for (int sdd_idx = 0; sdd_idx < kNumSDDs; ++sdd_idx) {
      // Clone histograms
      auto hist1 = dynamic_cast<TH1D*>(h_adc_raw[bus_idx][sdd_idx]->Clone());
      auto hist2 = dynamic_cast<TH1D*>(h_adc[bus_idx][sdd_idx]->Clone());
      auto hist3 = dynamic_cast<TH1D*>(h_xtalk[bus_idx][sdd_idx]->Clone());

      // Create a new page for every 8 SDDs
      if (sdd_idx % 8 == 0) {
        canvas->Clear();
        canvas->Divide(NUM_COLUMNS, NUM_ROWS);
      }
      // Add histograms to the canvas
      canvas->cd((sdd_idx % 8) + 1);
      gPad->SetLogy();
      if (hist1->GetEntries() <= 100) continue; // skip empty histo
      hist1->SetTitle(
          Form("Fluorescence x-ray spectrum (BUS: %d, SDD: %d)", 
                bus_idx+1, sdd_idx+1));
      hist1->GetXaxis()->SetRangeUser(kMinX, kMaxX);
      hist1->UseCurrentStyle();
      hist1->SetLineColor(kBlue);
      hist2->SetLineColor(kRed);
      hist3->SetLineColor(kBlack);
      hist1->Draw("HIST");
      hist2->Draw("HIST SAME");
      hist3->Draw("HIST SAME");
      canvas->Update();
    }
  }
}

void SpectrumAnalyser::DrawSDDMap(const std::string& filename)
{
  // Disable "Info in <TCanvas::Print>" message
  gErrorIgnoreLevel = kWarning;

  // Clone histograms
  auto hist1 = dynamic_cast<TH2D*>(h_sdd_map->Clone());
  auto hist2 = dynamic_cast<TH2D*>(h_sdd_rate->Clone());
  auto hist3 = dynamic_cast<TH2D*>(h_sdd_rate_sig->Clone());
  auto hist4 = dynamic_cast<TH2D*>(h_sdd_rate_noise->Clone());

  // Create a canvas and draw histograms
  static auto canvas = CreateCanvas(CanvasFormat::A4, 
                                    CanvasOrientation::Landscape);
  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  Draw2DHistogram(hist1,"SDD positions around the target");
  canvas->cd(2);
  Draw2DHistogram(hist2,"Map of SDD counts around the target");
  canvas->cd(3);
  Draw2DHistogram(hist3,"Map of SDD signal counts around the target");
  canvas->cd(4);
  Draw2DHistogram(hist4,"Map of SDD noise counts around the target");

  // Save canvas as image
  canvas->SaveAs(Form("output/plots/hSDD_map_%s.pdf", filename.c_str()));

  std::cout << "2D histograms have been saved to output/plots/hSDD_map_" 
            << filename.c_str() << ".pdf" << std::endl;
}

void SpectrumAnalyser::Draw2DHistogram(TH2D* hist, const std::string& title) 
{
  hist->SetTitle(title.c_str());
  hist->GetXaxis()->SetTitle("column");
  hist->GetYaxis()->SetTitle("row");
  hist->UseCurrentStyle();
  hist->Draw("COLZ");
}

void SpectrumAnalyser::SetPreCalibFitFunction(
    TString& fit_func, int n_peaks) 
{
  const int kNumGaussParams = 3;  // Number of parameters for the Gaussian function
  const int kNumBkgParams   = 3;  // Number of parameters for the background: p0 + exp(p1 + p2*x)

  int gauss_amp[kMaxNumPeaksPF]   = {-1};   // Gaussian parameters
  int gauss_mean[kMaxNumPeaksPF]  = {-1};
  int gauss_sigma[kMaxNumPeaksPF] = {-1};
  int bkg_p0 = -1;  // Background parameters
  int bkg_p1 = -1;
  int bkg_p2 = -1;

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

  float init_threshold = 0.01; // initial threshold parameter for the Peak Finder (std in TSpectrum is 0.05)
  float init_tolerance = 0.05; // tolerance to check that the peak assumption is correct (5%)
  
  const int sigma_pf = 20 / factor;  // sigma for the Peak Finder (in units of ADC channels)

  const int kMinStats = 1000;

  TSpectrum spectrum{kMaxNumPeaksPF};

  for(int bus_idx = 0; bus_idx < kNumBuses; ++bus_idx) {
    for(int sdd_idx = 0; sdd_idx < kNumSDDs; ++sdd_idx) {
      if (!h_adc[bus_idx][sdd_idx]) continue;

      auto histo = static_cast<TH1F*>(h_adc[bus_idx][sdd_idx]->Clone("histo"));
      if (histo->GetEntries() < kMinStats) continue;

      std::vector<std::pair<float, float>> peak_pos;

      std::cout << std::endl << "--- PEAK FINDER: BUS#" 
                << bus_idx + 1 << " SDD#" 
                << sdd_idx + 1<< ". Statistics: " 
                << histo->GetEntries() << " ---" << std::endl << std::endl;

      histo->GetXaxis()->SetRangeUser(x_min, x_max);

      int num_found = 0;
      int num_tries = 0;
      float peak_threshold = init_threshold;

      // Use TSpectrum to find the peak candidates      
      while(num_found < kNumPeaksPF && num_tries < 10) {
        num_found = spectrum.Search(histo, sigma_pf, "", peak_threshold);        
        peak_threshold *= 0.1;  // Initial = 0.01. It changes until it finds peaks
        num_tries++;
      }

      std::cout << "Found " << num_found << " candidate peaks:" << std::endl;

      if (num_found == 0) {
        std::cout << "PEAK FINDER COULD NOT FIND ANY PEAKS. CONTINUING..." 
                  << std::endl;
        continue;
      }
      if (num_tries >= 10) {
        std::cout << "PEAK FINDER DOES NOT WORK. CONTINUING..." << std::endl;
        continue;
      }

      for (int i = 0; i < num_found; ++i) {
        peak_pos.emplace_back(
            spectrum.GetPositionX()[i], spectrum.GetPositionY()[i]);
      }

      // Sort the peaks based on their x positions
      std::sort(peak_pos.begin(), peak_pos.end(),
          [](const auto& a, const auto& b) { return a.first < b.first; });

      // Print found peaks
      for (int i = 0; i < num_found; ++i) {
        std::cout << "Peak#" << i + 1 << ". Position [ADC]: "
                  << peak_pos[i].first << ". Intensity: "
                  << peak_pos[i].second << std::endl;
      }
      
      // Check if the peaks found are compatible with the assumption
      float slope  = 0.0;
      float offset = 0.0;

      // Make a triad from all the peaks found:
      int peak_idx[3];
      
      for (int i0 = 0; i0 < num_found; i0++) {
        for (int i1 = i0 + 1; i1 < num_found; i1++) {
          for(int i2 = i1 + 1; i2 < num_found; i2++) {
            std::cout << std::endl << "-> Trying the triad: " 
                      << i0 << " " << i1 << " " << i2 << std::endl;
            float adc_ratio = (peak_pos[i2].first - peak_pos[i1].first) / 
                              (peak_pos[i1].first - peak_pos[i0].first);

            std::cout << "Checking assumption: Energy relation " 
                      << energy_ratio << " vs ADC relation " 
                      << adc_ratio << std::endl;

            float peak_tolerance = init_tolerance;  // 5%
            bool tolerance_pass = true;
            if (abs(1.0 - (energy_ratio / adc_ratio)) > peak_tolerance) {
              tolerance_pass = false;
            }

            // Get the Peak Finder calibration offset and slope
            float adc_diff = peak_pos[i0].first - peak_pos[i1].first;
            float energy_diff = triad_peak_energies[0] - triad_peak_energies[1];
            slope  = energy_diff / adc_diff;
            offset = -1.0 * peak_pos[i0].first * slope + triad_peak_energies[0];

            std::cout << "Offset = " << offset 
                      << ". Slope = " << slope << std::endl;

            //define an acceptable gain and offset; check if conditions (tolerance, G and G0) are met
            float min_slope = 2.9;  float max_slope = 3.9;
            float min_offset = -3000;   float max_offset = -1000;
            if ((slope < max_slope) && (slope > min_slope) && 
                (offset < max_offset) && (offset > min_offset) && 
                tolerance_pass) {
              std::cout << " -- TEST PASSED! --" << std::endl;              
              peak_idx[0] = i0;
              peak_idx[1] = i1;
              peak_idx[2] = i2;
            }
          }
        }
      }

      // Find also the highest height among the selected ones
      float highest_peak = 0.0;
      if (peak_idx[0] > -1) {
        if (peak_pos[peak_idx[0]].second > highest_peak) 
          highest_peak = peak_pos[peak_idx[0]].second;
        if (peak_pos[peak_idx[1]].second > highest_peak) 
          highest_peak = peak_pos[peak_idx[1]].second;
        if (peak_pos[peak_idx[2]].second > highest_peak) 
          highest_peak = peak_pos[peak_idx[2]].second;
      }
      
      std::cout << "The highest peak is " << highest_peak << std::endl;
      
      histo->Delete();
    }
  }
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
