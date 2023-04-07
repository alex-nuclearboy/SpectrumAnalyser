#include "../include/SpectrumAnalyser.h"

#include <TH1D.h>
#include <TH2D.h>

void SpectrumAnalyser::InitHistograms(int rebin_factor) 
{
  num_adc_bins /= rebin_factor;   // divide number of bins by rebin factor
  //factor = rebin_factor;
  
  // Loop over all the buses and SDDs to create the histograms
  for (int bus_id = 0; bus_id < num_buses; ++bus_id) {
    for (int sdd_id = 0; sdd_id < num_sdds; ++sdd_id) {

      // ADC histograms
      h_adc[bus_id][sdd_id] = new TH1D(
          Form("h_adc_bus%i_sdd%i", bus_id + 1, sdd_id + 1),
          Form("ADC distribution for bus %i, SDD %i", bus_id + 1, sdd_id + 1),
          num_adc_bins, adc_min, adc_max
      );
      h_adc[bus_id][sdd_id]->GetXaxis()->SetTitle("ADC [channel]");
      h_adc[bus_id][sdd_id]->GetYaxis()->SetTitle(
          Form("counts / %i channel", rebin_factor));
      h_adc_raw[bus_id][sdd_id] = static_cast<TH1D*>(
          h_adc[bus_id][sdd_id]->Clone(
              Form("h_adc_raw_bus%i_sdd%i", bus_id + 1, sdd_id + 1)));

      // Crosstalk histogram
      h_xtalk[bus_id][sdd_id] = new TH1D(
          Form("h_xtalk_bus%i_sdd%i", bus_id + 1, sdd_id + 1),
          Form("Crosstalk distribution for bus %i, SDD %i", 
                bus_id + 1, sdd_id + 1),
          num_adc_bins, adc_min, adc_max);
    }
  }

  // Create the SDD maps and rate histograms
  h_sdd_map = new TH2D("h_sdd_map", "SDD occupancy map", 48, 1, 49, 8, 1, 9);
  h_sdd_rate = new TH2D("h_sdd_rate", "SDD rate", 48, 1, 49, 8, 1, 9);
  h_sdd_rate_sig = new TH2D(
      "h_sdd_rate_sig", "SDD signal rate", 48, 1, 49, 8, 1, 9);
  h_sdd_rate_noise = new TH2D(
      "h_sdd_rate_noise", "SDD noise rate", 48, 1, 49, 8, 1, 9);
}
