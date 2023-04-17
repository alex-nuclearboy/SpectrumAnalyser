#include "../include/SpectrumAnalyser.h"

#include <fstream>
#include <TF1.h>

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
      fit_func->GetParameter(0),
      fit_func->GetParameter(7),
      fit_func->GetParameter(6)) << std::endl;
  file.close();
}
