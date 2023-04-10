/*
* File:         SpectrumAnalyser.cxx
* Author:       Aleksander Khreptak <aleksander.khreptak@lnf.infn.it>
* Created:      31 Mar 2023
* Last updated: 10 Apr 2023
*
* Description:
* This file implements a SpectrumAnalyser class that is used
* to analyse the SDD data from the SIDDHARTA-2 experiment
*/

#define SPECTRUM_ANALYSER_CXX
#include "../include/SpectrumAnalyser.h"

#include <iostream>
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
