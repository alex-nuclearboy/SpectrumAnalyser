#include "../include/SpectrumAnalyser.h"

#include <TSpectrum.h>

void SpectrumAnalyser::FindADCPeaks(
    const float& x_min, const float& x_max, const std::string& filename)
{
  const float energy_ratio = (triad_peak_energies[2] - triad_peak_energies[1]) /
                             (triad_peak_energies[1] - triad_peak_energies[0]);

  std::cout << "energy_ratio = " << energy_ratio << std::endl;

  float init_threshold = 0.01; // initial threshold parameter for the Peak Finder (std in TSpectrum is 0.05)
  float init_tolerance = 0.05; // tolerance to check that the peak assumption is correct (5%)
  
  const int sigma_pf = 20 / rebin_factor;  // sigma for the Peak Finder (in units of ADC channels)

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
          highest_peak, bus_idx, sdd_idx, filename);

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
