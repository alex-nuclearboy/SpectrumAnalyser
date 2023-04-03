/*
* File:         run_analysis.cpp
* Author:       Aleksander Khreptak <aleksander.khreptak@lnf.infn.it>
* Created:      31 Mar 2023
* Last updated: 01 Apr 2023
*
* Description:
* This program performs analysis on SDD data using the SpectrumAnalyser class
* 
* Usage:
* ./run_analysis [input_file_name]
*/

#include <iostream>
#include <fstream>
#include <memory>
#include <string>

#include <TFile.h>
#include <TTree.h>

#include "SpectrumAnalyser.h"
#include "version.h"

std::string checkExtension(const std::string& filename);

int main(int argc, char* argv[]) {
  // Check the number of command-line arguments
  if (argc < 2) {
    std::cerr << "ERROR: No input file specified." << std::endl;
    std::cerr << "Usage: " << argv[0] << " [input file name]" << std::endl;
    return EXIT_FAILURE;
  }

  printWelcomeMessage();

  // Set the input directory and input file name
  const std::string data_directory  = "${XRAYDATA}";
  const std::string input_file_name = checkExtension(argv[1]);
  const std::string file_date       = input_file_name.substr(0, 8);  
  const std::string input_file_path = data_directory + "/" + file_date 
                                  + "/" + input_file_name;

  // Open an input file in read mode
  std::unique_ptr<TFile> input_file_handle(
      TFile::Open(input_file_path.c_str(), "READ"));
  // Check if the input file was opened successfully
  if (!input_file_handle || input_file_handle->IsZombie()) {
    std::cerr << "ERROR: Could not open file " << input_file_name << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Reading events from file: " << input_file_name << std::endl;

  // Retrieve the TTree from the input file
  TTree* tree = nullptr;
  input_file_handle->GetObject("ft10", tree);
  // Check if the TTree was retrieved successfully
  if (!tree) {
    std::cerr << "ERROR: Could not retrieve TTree from file "
              << input_file_name << std::endl;    
    return EXIT_FAILURE;
  }

  // Create an instance of SpectrumAnalyser and analyse the events in the TTree
  SpectrumAnalyser analyser(tree);
  
  analyser.initHistograms(8);

  analyser.Loop();
  
  // Set the output file name
  std::string output_file_name = input_file_name;
  std::string output_file_path = "output/rootfiles/histos_" + output_file_name;
  analyser.writeHistograms(output_file_path);
  analyser.myFunction();
  // Deallocate memory used by the TFile object
  input_file_handle->Close();

  return 0;
}

std::string checkExtension(const std::string& filename) {
  if (filename.substr(filename.size() - 5) != ".root") {
    return filename + ".root";
  }
  return filename;
}
