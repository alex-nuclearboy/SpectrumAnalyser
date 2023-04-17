#include "../include/SpectrumAnalyser.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include <TH1D.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPDF.h>
#include "TError.h"

void SpectrumAnalyser::SearchPeaks(const std::string& file_name)
{
  gErrorIgnoreLevel = kError;

  std::ofstream file;

  double threshold = 0.01;
  // Create a TSpectrum objects to search for peaks in the histogram
  auto spectrum_Cu = std::make_unique<TSpectrum>(4);
  auto spectrum_Ti = std::make_unique<TSpectrum>(2);

  std::map<std::string, float> CuKas;
  std::map<std::string, float> TiKas; 
  std::map<std::string, std::vector<float>> centers;

  TF1* fit_Cu; 
  TF1* fit_Ti;

  for (int bus_id = 0; bus_id < num_buses; ++bus_id) {
    for (int sdd_id = 0; sdd_id < num_sdds; ++sdd_id) {
      // Clone the histogram to avoid modifying the original
      
      h_adc_pf[bus_id][sdd_id] =
          dynamic_cast<TH1D*>(h_adc[bus_id][sdd_id]->Clone());
      if (h_adc_pf[bus_id][sdd_id]->GetEntries() <= 1000) continue;
      h_adc_pf[bus_id][sdd_id]->Rebin(rebin_factor);
      
      // CuKa
      h_adc_pf[bus_id][sdd_id]->GetXaxis()->SetRangeUser(1700, 3700);

      const int num_peaks_Cu =
          spectrum_Cu->Search(h_adc_pf[bus_id][sdd_id], 5, "", threshold);
      std::map<float, float> peaks_bin_Cu;
      for (int i = 0; i < num_peaks_Cu; i++) {
        float bin_pos = spectrum_Cu->GetPositionX()[i];
        if (bin_pos > 2500 && bin_pos < 3500) {
          float bin_pos_cont = h_adc_pf[bus_id][sdd_id]->GetBinContent(
              h_adc_pf[bus_id][sdd_id]->GetXaxis()->FindBin(bin_pos));
              peaks_bin_Cu[bin_pos] = bin_pos_cont;
        } else {
          peaks_bin_Cu[bin_pos] = 0;
        }
      }
      if (!peaks_bin_Cu.empty()) {
        // Sort the peaks by their intensity
        std::multimap<float, float, std::greater<float>> sorted_peaks;
        for (const auto& peak : peaks_bin_Cu) {
          sorted_peaks.insert({peak.second, peak.first});
        }
        // Get the highest peak
        float CuKa_center = sorted_peaks.begin()->second;
        
        //TF1 fit_Cu("m1", "gaus(0)+pol0(3)", CuKa_center - 70, CuKa_center + 70);
        fit_Cu = new TF1("m1", "gaus(0)+pol0(3)", CuKa_center - 70, CuKa_center + 70);
        fit_Cu->SetLineColor(2);
        fit_Cu->SetParameter(0, h_adc_pf[bus_id][sdd_id]->GetMaximum());
        fit_Cu->SetParLimits(0, 0, 1.5 * h_adc_pf[bus_id][sdd_id]->GetMaximum());
        fit_Cu->SetParameter(1, CuKa_center);
        fit_Cu->SetParLimits(1, CuKa_center - 30, CuKa_center + 30);
        fit_Cu->SetParLimits(2, 10, 120);
        fit_Cu->SetParameter(2, 80);        
        h_adc_pf[bus_id][sdd_id]->Fit(fit_Cu, "RQ");
        float fit_Cu_chi2 = fit_Cu->GetChisquare() / fit_Cu->GetNDF();
        std::string hist_name = "bus_" + std::to_string(bus_id+1) + "_sdd_" + std::to_string(sdd_id+1);
        CuKas[hist_name] = fit_Cu->GetParameter(1);
        std::vector<double> center = { fit_Cu->GetParameter(1) };
        centers[hist_name].push_back(fit_Cu->GetParameter(1));

        if (fit_Cu_chi2 > 2) {
          std::cout << "WARNING: Chi^2 / NDF for Cu Ka [bus" << bus_id + 1 
                    << ", sdd" << sdd_id << "] > 2 (" << fit_Cu_chi2 
                    << ")" << std::endl;
        }
      }

      // TiKa
      h_adc_pf[bus_id][sdd_id]->GetXaxis()->SetRangeUser(1700, 2400);

      const int num_peaks_Ti =
          spectrum_Ti->Search(h_adc_pf[bus_id][sdd_id], 5, "nobackground", threshold);

      std::map<float, float> peaks_bin_Ti;
      float TiKa_center = 0;
      for (int i = 0; i < num_peaks_Ti; i++) {
        float bin_pos = spectrum_Ti->GetPositionX()[i];
        if (bin_pos > 1700 && bin_pos < 2300) {
          float bin_pos_cont = h_adc_pf[bus_id][sdd_id]->GetBinContent(
              h_adc_pf[bus_id][sdd_id]->GetXaxis()->FindBin(bin_pos));
              peaks_bin_Ti[bin_pos] = bin_pos_cont;
        } else {
          peaks_bin_Ti[bin_pos] = 0;
        }
      }
      if (!peaks_bin_Ti.empty()) {
        // Sort the peaks by their intensity
        std::multimap<float, float, std::greater<float>> sorted_peaks;
        for (const auto& peak : peaks_bin_Ti) {
          sorted_peaks.insert({peak.second, peak.first});
        }
        // Get the highest peak
        TiKa_center = sorted_peaks.begin()->second;
        //float highest_peak_cont = sorted_peaks.begin()->first;
        //TF1 fit_Ti("m1", "gaus(0)+pol0(3)", TiKa_center - 70, TiKa_center + 70);
        fit_Ti = new TF1("m2", "gaus(0)+pol0(3)", TiKa_center - 70, TiKa_center + 70);
        fit_Ti->SetLineColor(2);
        fit_Ti->SetParameter(0, h_adc_pf[bus_id][sdd_id]->GetMaximum());
        fit_Ti->SetParLimits(0, 0, 1.5 * h_adc_pf[bus_id][sdd_id]->GetMaximum());
        fit_Ti->SetParameter(1, TiKa_center);
        fit_Ti->SetParLimits(1, TiKa_center - 30, TiKa_center + 30);
        fit_Ti->SetParLimits(2, 10, 120);
        fit_Ti->SetParameter(2, 80);        
        h_adc_pf[bus_id][sdd_id]->Fit(fit_Ti, "RQ+");
        float fit_Ti_chi2 = fit_Ti->GetChisquare() / fit_Ti->GetNDF();
        std::string hist_name = "bus_" + std::to_string(bus_id+1) + "_sdd_" + std::to_string(sdd_id+1);
        TiKas[hist_name] = fit_Ti->GetParameter(1);
        std::vector<double> center = { fit_Ti->GetParameter(1) };
        centers[hist_name].push_back(fit_Ti->GetParameter(1));

        if (fit_Ti_chi2 > 2) {
          std::cout << "WARNING: Chi^2 / NDF for Ti Ka [bus" << bus_id + 1 
                    << ", sdd" << sdd_id << "] > 2 (" << fit_Ti_chi2 
                    << ")" << std::endl;
        }
      }

      file.open(
          "output/files/peaks_TiKa_CuKa" + file_name + ".dat", std::ios::app);
      file << Form("%d  %d  %f  %f", bus_id + 1, sdd_id + 1, 
          fit_Ti->GetParameter(1), fit_Cu->GetParameter(1)) << std::endl;
      file.close();
    }
  }
  DrawPFSpectra(file_name);
}

void SpectrumAnalyser::DrawPFSpectra(const std::string& file_name) 
{   
  auto canvas = CreateCanvas(CanvasFormat::A4, CanvasOrientation::Portrait);
  auto pdf = std::make_unique<TPDF>(
      Form("output/plots/PeakFinder_spectra_%s.pdf", file_name.c_str()));
  canvas->SetName("canvas_PF");
  DrawPFSpectrum(canvas);
  pdf->Close();
  canvas->Close();
  std::cout << "Peak Finder spectra have been saved to output/plots/PeakFinder_spectra_" 
            << file_name.c_str() << ".pdf" << std::endl;
}

void SpectrumAnalyser::DrawPFSpectrum(TCanvas* canvas)
{
  const int num_columns = 2;
  const int num_rows    = 4;
  const int x_min = 1500;
  const int x_max = 4500;

  for (int bus_id = 0; bus_id < num_buses; ++bus_id) {
    for (int sdd_id = 0; sdd_id < num_sdds; ++sdd_id) {     
      h_adc_pf[bus_id][sdd_id]->GetXaxis()->SetRangeUser(x_min, x_max);
      if (h_adc_pf[bus_id][sdd_id]->GetEntries() <= 100) continue; // skip empty histo

      // Create a new page for every 8 SDDs
      if (sdd_id % 8 == 0) {
        canvas->Clear();
        canvas->Divide(num_columns, num_rows);
      }
      // Add histograms to the canvas
      canvas->cd((sdd_id % 8) + 1);
      gPad->SetLogy();      
      h_adc_pf[bus_id][sdd_id]->SetTitle(
          Form("Fluorescence x-ray spectrum (BUS: %d, SDD: %d)", 
                bus_id+1, sdd_id+1));      
      h_adc_pf[bus_id][sdd_id]->GetYaxis()->SetTitle(
          Form("counts / %d channels", rebin_factor));    
      h_adc_pf[bus_id][sdd_id]->UseCurrentStyle();
      //h_adc_pf[bus_id][sdd_id]->SetLineColor(kBlue);
      h_adc_pf[bus_id][sdd_id]->Draw();      
      canvas->Update();
    }
  }
}
