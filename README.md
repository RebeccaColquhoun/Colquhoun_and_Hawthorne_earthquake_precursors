# Basic code repository for Colquhoun and Hawthorne -- Investigating earthquake precursors: can they help to constrain earthquake nucleation processes?

Basic example of codes used in the calculation of Colquhoun, Hawthorne & Kamranzad (in prep.)

## Requirements
A number of packages are required to run the code. Distributed packages are listed in the yml file, from which a conda environment can be created as follows:
```bash
conda env create -f phasecoherencedependencies.yml
```

Other requirements are provided in this depository.

## Getting to grips with phase coherence
A jupyter notebook is provided which allows you to explore how different parameters affect the final results. 

## Use of code:
1. calculatePC_multiprocessing --> this creates a phase coherence object for each earthquake, and calculates the phase coherence. Uses PC_calculator and PhaseCoherence
2. consolidate_detections --> this takes phase coherence detections from each earthquakes objects, and concatenates them into lists, based on detection limit. 

## Data:
Expects data to be arranged in folders, by earthquake. 
