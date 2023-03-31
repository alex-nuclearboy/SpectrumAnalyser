/*
* File:         SpectrumAnalyser.cxx
* Author:       Aleksander Khreptak <aleksander.khreptak@lnf.infn.it>
* Created:      31 Mar 2023
* Last updated: 01 Apr 2023
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

#include <TBranch.h>
#include <TLeaf.h>

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
  std::cout << "Start time: " << convertTime(time_start) << std::endl;
  std::cout << "End time: " << convertTime(time_end) << std::endl;
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
      if (ihit > 0) {ISGOODtime = crossTalkTiming(drift[ihit], drift[ihit-1]);}
      if (ISGOODtime) {
        h_adc[bus[ihit]-1][sdd[ihit]-1]->Fill(adc[ihit]);
      } else {
        h_xtalk[bus[ihit]-1][sdd[ihit]-1]->Fill(adc[ihit]);
      }
      sddHitMap(sdd[ihit], bus[ihit], column, row);
      h_sdd_map->SetBinContent(column, row, sdd[ihit]);
      h_sdd_rate->Fill(column, row);
      if (adc[ihit] > 1200) {h_sdd_rate_sig->Fill(column, row);}
      if (adc[ihit] <= 1200) {h_sdd_rate_noise->Fill(column, row);}
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
}

std::string SpectrumAnalyser::convertTime(time_t t)
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

void SpectrumAnalyser::initHistograms(int rebin_factor) 
{
  num_adc_bins /= rebin_factor;   // divide number of bins by rebin factor
  
  // Loop over all the buses and SDDs to create the histograms
  for (int bus_idx = 0; bus_idx < NUM_BUSES; ++bus_idx) {
    for (int sdd_idx = 0; sdd_idx < NUM_SDDS; ++sdd_idx) {

      // ADC histograms
      h_adc[bus_idx][sdd_idx] = new TH1D(
          Form("h_adc_bus%i_sdd%i", bus_idx + 1, sdd_idx + 1),
          Form("h_adc_bus%i_sdd%i", bus_idx + 1, sdd_idx + 1),
          num_adc_bins, MIN_ADC, MAX_ADC
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
          num_adc_bins, MIN_ADC, MAX_ADC);
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

void SpectrumAnalyser::writeHistograms(const std::string& filename) {
  // Create the output file
  TFile* output_file = new TFile(filename.c_str(), "RECREATE");
  output_file->cd();
  // Write the histograms to the file
  for (int bus_idx = 0; bus_idx < NUM_BUSES; ++bus_idx) {
    for (int sdd_idx = 0; sdd_idx < NUM_SDDS; ++sdd_idx) {
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

bool SpectrumAnalyser::crossTalkTiming(Short_t drift, Short_t drift_pre) {
  bool ISGOODtime = true;
  int t1 = 0, t2 = 0, timediff = 0;
  t1 = drift + 32768;
  t2 = drift_pre + 32768.;
  if (t1 > t2) {timediff = t1 - t2;}
  if (t2 > t1) {timediff = t1 + (32768.*2) - t2;}
  if ((timediff > 0. && timediff < 625.)) {ISGOODtime = false;}   // 625 -> 5 microseconds
  return ISGOODtime;
}

void SpectrumAnalyser::sddHitMap(int sddnumber, int busnumber,
                                    int &column, int &row) {
  row = 0; column = 0;
  int SFERA = 0;

  // Back side view
  if (sddnumber <= 16) {sddnumber = sddnumber; SFERA = 0;}
  if (sddnumber >= 17 && sddnumber <= 32) {sddnumber = sddnumber - 16; SFERA = 2;}
  if (sddnumber >= 33 && sddnumber <= 48) {sddnumber = sddnumber - (16*2); SFERA = 4;}
  if (sddnumber >= 49 && sddnumber <= 64) {sddnumber = sddnumber - (16*3); SFERA = 6;}

  if (sddnumber == 1)   {row = 4; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 2)   {row = 3; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 3)   {row = 2; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 4)   {row = 1; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 5)   {row = 1; column = 2 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 6)   {row = 2; column = 2 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 7)   {row = 3; column = 2 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 8)   {row = 4; column = 2 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 9)   {row = 8; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 10)  {row = 7; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 11)  {row = 6; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 12)  {row = 5; column = 1 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 13)  {row = 5; column = 2 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 14)  {row = 6; column = 2 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 15)  {row = 7; column = 2 + SFERA + ((busnumber - 1)*8);}
  if (sddnumber == 16)  {row = 8; column = 2 + SFERA + ((busnumber - 1)*8);}
}

int SpectrumAnalyser::SFERAnumber(int sdd) {
  int SFERA = 0;
  if (sdd <= 16) {SFERA = 1;}
  if (sdd >= 17 && sdd <= 32) {SFERA = 2;}
  if (sdd >= 33 && sdd <= 48) {SFERA = 3;}
  if (sdd >= 49 && sdd <= 64) {SFERA = 4;}
  return SFERA;
}
