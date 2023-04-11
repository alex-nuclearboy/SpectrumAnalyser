#include "../include/SpectrumAnalyser.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPDF.h>
#include <TError.h>

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

void SpectrumAnalyser::DrawADCSpectra(const std::string& file_name) 
{   
  auto canvas = CreateCanvas(CanvasFormat::A4, CanvasOrientation::Portrait);
  auto pdf = std::make_unique<TPDF>(
      Form("output/plots/hADC_spectra_%s.pdf", file_name.c_str()));
  DrawSpectrum(canvas);
  pdf->Close();
  canvas->Close();
  std::cout << "ADC spectra have been saved to output/plots/hADC_spectra_" 
            << file_name.c_str() << ".pdf" << std::endl;
}

void SpectrumAnalyser::DrawSpectrum(TCanvas* canvas)
{
  const int NUM_COLUMNS = 2;
  const int NUM_ROWS    = 4;
  const int kMinX = 1500;
  const int kMaxX = 4500;
  
  for (int bus_idx = 0; bus_idx < num_buses; ++bus_idx) {
    for (int sdd_idx = 0; sdd_idx < num_sdds; ++sdd_idx) {
      // Clone histograms
      auto hist1 = dynamic_cast<TH1D*>(h_adc_raw[bus_idx][sdd_idx]->Clone());
      auto hist2 = dynamic_cast<TH1D*>(h_adc[bus_idx][sdd_idx]->Clone());
      auto hist3 = dynamic_cast<TH1D*>(h_xtalk[bus_idx][sdd_idx]->Clone());

      hist1->Rebin(rebin_factor);
      hist2->Rebin(rebin_factor);
      hist3->Rebin(rebin_factor);

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
      hist1->GetYaxis()->SetTitle(Form("counts / %d channels", rebin_factor));    
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

void SpectrumAnalyser::DrawSDDMap(const std::string& file_name)
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
  canvas->SaveAs(Form("output/plots/hSDD_map_%s.pdf", file_name.c_str()));

  std::cout << "2D histograms have been saved to output/plots/hSDD_map_" 
            << file_name.c_str() << ".pdf" << std::endl;
}

void SpectrumAnalyser::Draw2DHistogram(TH2D* hist, const std::string& title) 
{
  hist->SetTitle(title.c_str());
  hist->GetXaxis()->SetTitle("column");
  hist->GetYaxis()->SetTitle("row");
  hist->UseCurrentStyle();
  hist->Draw("COLZ");
}
