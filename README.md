# SDD Spectrum Analyser

The SDD Spectrum Analyser is a software tool developed for the analysis and visualisation of data obtained from Silicon Drift Detectors (SDDs) used in the <a href="https://iopscience.iop.org/article/10.1088/1402-4896/ac7fc0">SIDDHARTA-2 experiment</a> at the DAΦNE collider. This experiment aims to perform precision spectroscopy measurements of kaonic deuterium x-ray transitions to the fundamental 1s level, which will provide insight into low-energy quantum chromodynamics (QCD).

To achieve this challenging goal, novel <a href="https://iopscience.iop.org/article/10.1088/1361-6501/ac777a/meta">large-area Silicon Drift Detectors</a> (SDDs) have been developed with a special geometry, field configuration, and readout electronics that ensure excellent performance in terms of linearity and stability. The x-ray detection system of the SIDDHARTA-2 apparatus consists of 48 arrays of SDDs, each featuring eight cells, resulting in a total of 384 read-out channels. The signals coming from the detectors are buffered by a CMOS low-noise charge-sensitive preamplifier (CUBE), then processed by a dedicated analog SDD front-end readout ASIC (SFERA) before being digitised by an analogue-to-digital converter (ADC, model <a href="https://www.ni.com/it-it/support/model.pci-6115.html">NI PCI-6115</a>).

The analyser includes modules for drawing ADC histograms, automatic peak searches, peak identification, energy calibration, and background subtraction. These modules are designed to facilitate the analysis of SDD data and ensure that the results obtained are reliable and accurate.

## Required software

A data analysis framework

      ROOT 

Cross-platform build-system generator

      CMake

## Necessary environment variables

Path where data from experiment are located:

      XRAYDATA

## Installation & Compilation

      git clone https://github.com/alex-nuclearboy/SpectrumAnalyser.git
      cd SpectrumAnalyser
      mkdir SDD-build
      cd SDD-build
      cmake ..
      make

## Usage

To perform an analysis of the spectra, use the following command:

      ./run_analysis <ROOT file name>

## The result of analysis

The Spectral Analyser will generate a ROOT file containing histograms for each individual detector (including all events, crosstalk, and events after crosstalk subtraction). 
It also saves all histograms in PDF format. Additionally, it will generate a text file containing position of Ti Kα and Cu Kα peaks.

The ROOT files can then be used for calibration employing the software available <a href="https://github.com/alex-nuclearboy/SDDCalibration">here</a>.
