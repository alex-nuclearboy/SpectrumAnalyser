#include "../include/SpectrumAnalyser.h"

#include <TFile.h>

void SpectrumAnalyser::WriteHistograms(const std::string& file_name)
{
  // Create the output file
  TFile* output_file = new TFile(file_name.c_str(), "RECREATE");
  output_file->cd();
  // Write the histograms to the file
  for (int bus_idx = 0; bus_idx < num_buses; ++bus_idx) {
    for (int sdd_idx = 0; sdd_idx < num_sdds; ++sdd_idx) {
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
  std::cout << "Output file: " << file_name << std::endl;
}
