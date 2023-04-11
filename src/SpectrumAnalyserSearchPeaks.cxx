#include "../include/SpectrumAnalyser.h"

#include <iostream>
#include <vector>

#include <TH1D.h>
#include <TSpectrum.h>
#include <TF1.h>

constexpr int MIN_STATS = 1000;
constexpr int NUM_PEAKS = 4;
constexpr double threshold = 0.01;
constexpr double min_position = 2500;
constexpr double max_position = 3500;
constexpr double fit_range = 70;
constexpr double fit_width = 80;
constexpr double sigma_min = 10;
constexpr double sigma_max = 120;

std::vector<Double_t> SpectrumAnalyser::SearchPeakCuKa() {
  std::vector<Double_t> CuKa_centers;

  for (int bus_id = 0; bus_id < num_buses; ++bus_id) {
    for (int sdd_id = 0; sdd_id < num_sdds; ++sdd_id) {
      if (!h_adc[bus_id][sdd_id]) continue;

      TH1D* hist = static_cast<TH1D*>(h_adc[bus_id][sdd_id]->Clone("hist"));
      std::cout << std::endl << "--- STARTING PEAK FINDER ---" << std::endl
                << "for [bus#" << bus_id + 1 << ", sdd#" << sdd_id + 1
                << "] with a statistics of " << hist->GetEntries() << std::endl
                << std::endl;

      if (hist->GetEntries() < MIN_STATS) {
        std::cout << "NOT ENOUGH STATISTICS (" << hist->GetEntries()
                  << "). SKIPPING..." << std::endl;
        continue;
      }

      hist->Rebin(rebin_factor);

      // Set the range of the histogram
      const int bin_min = hist->GetXaxis()->FindBin(1500);
      const int bin_max = hist->GetXaxis()->FindBin(4500);
      hist->GetXaxis()->SetRange(bin_min, bin_max);

      // Perform peak search with TSpectrum
      TSpectrum spectrum(NUM_PEAKS);
      const Int_t num_found = spectrum.Search(hist, 50, "", threshold);

      // Select the CuKa peak
      double max_content = 0;
      double CuKa_center = -1;
      for (Int_t i = 0; i < num_found; ++i) {
        const double bin_pos = spectrum.GetPositionX()[i];
        const int bin = hist->GetXaxis()->FindBin(bin_pos);
        const double bin_pos_cont = hist->GetBinContent(bin);
        if (bin_pos > min_position && bin_pos < max_position &&
            bin_pos_cont > max_content) {
          max_content = bin_pos_cont;
          CuKa_center = bin_pos;
        }
      }

      if (CuKa_center < 0) {
        std::cout << "CANNOT FIND CuKa PEAK. SKIPPING..." << std::endl;
        continue;
      }

      std::cout << "CuKa PEAK FOUND AT " << CuKa_center << " [ADC]" << std::endl;

      // Fit the CuKa peak
      const double fit_min = CuKa_center - fit_range / 2;
      const double fit_max = CuKa_center + fit_range / 2;
      TF1* fit_func =
          new TF1("fit_func", "gaus(0)+pol0(3)", fit_min, fit_max);
      fit_func->SetParameter(0, max_content);
      fit_func->SetParameter(1, CuKa_center);
      fit_func->SetParameter(2, fit_width);
      hist->Fit(fit_func, "QR");

      const double CuKa_sigma = fit_func->GetParameter(2);
      if (CuKa_sigma < sigma_min || CuKa_sigma > sigma_max) {
        std::cout << "CANNOT FIND PROPER CuKa PEAK. SKIPPING..." << std::endl;
        continue;
      }
      
      std::cout << "CuKa PEAK FITTED SUCCESSFULLY" << std::endl
                << "FIT PARAMETERS: " << std::endl
                << "AMPLITUDE: " << fit_func->GetParameter(0) << " [counts]" 
                << std::endl
                << "MEAN: " << fit_func->GetParameter(1) << " [ADC]" 
                << std::endl
                << "SIGMA: " << CuKa_sigma << " [ADC]" << std::endl;
                
      CuKa_centers.push_back(fit_func->GetParameter(1));
      
      delete fit_func;
      delete hist;
    }
  }
  
  return CuKa_centers;
}
