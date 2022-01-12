# README for Code Underlying "The Recreational Value of Rare Species: Causal Evidence from the Cassia Crossbill"

This repository contains code for data cleaning, simulation, data analysis, and manuscrupt preparation for the paper "The Recreational Value of Rare Species: Causal Evidence from the Cassia Crossbill." 

Running the code in this repository requires downloading data and assigning it to the correct filepaths. Below, each code file is described along with its necessary inputs.

data_preparation.R
This file prepares raw data for analysis. Note that this code involves a number of computationally intense calculations (such as calculating a large distance matrix). It was originally run on departmental server and took overnight. It requires:
- A Census API key for use by package tidycensus. This can be input near the top of the file.
- The "Sampling Event Data" version of the eBird Basic Dataset. This dataset requires a very easy and short application to the Cornell Lab of Ornithology in order to access. Once access is granted, it can be accessed from https://ebird.org/data/download/ebd. The filepath of this data can also be modified near the top of the file. 
- Forest Service Ranger District boundaries. These can be downloaded from https://data.fs.usda.gov/geodata/edw/datasets.php?dsetCategory=boundaries. The file assumes the filepath "S_USA.RangerDistrict/S_Usa.RangerDistrict.shp"
- Canadian Census District Boundaries. These can be downloaded from https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2016-eng.cfm. The file assumes the filepath "lcd_000b16a_e/lcd_000b16a_e.shp"
- Canadian Census District Populations. These can be downloaded from https://www12.statcan.gc.ca/census-recensement/2011/dp-pd/prof/details/download-telecharger/comprehensive/comp-csv-tab-dwnld-tlchrgr.cfm?Lang=E#tabs2011. The file assumes the filepath "98-316-XWE2011001-701_CSV/98-316-XWE2011001-701.csv"

simulation_parameters.R
This short file defines the parameters that are inputs to the paper's behavioral model simulations. Changes to this file are automatically reflected in the text of the manuscript (with the exception of figure labels and notes).

simulation.R
This file takes the parameters defined in simulation_parameters.R and runs simulations as described in the text of the manuscript. Results are saved because the file can take hours to run.

manuscript.RMD
This file produces the manuscript pdf, but uses RMarkdown to natively conduct much of the analysis as well. This file should compile in a reasonable amount of time (<10 mins) if all the above steps have been completed properly. The file requires each of the above files to have already been run. In addition, it requires the following data:
- Custom Download from the eBird Basic Dataset site linked above for all sightings of the Cassia Crossbill. The filepath of this data can be modified in the first chunk of the file.
- References bibtex file. This is included in the repository, named references.bib.
