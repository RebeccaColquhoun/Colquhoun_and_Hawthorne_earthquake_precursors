# Basic code repository for Colquhoun and Hawthorne -- Investigating earthquake precursors: can they help to constrain earthquake nucleation processes?

Basic example of codes used in the calculation of Colquhoun and Hawthorne (submitted to EPSL, 2022)

## Use of code:
1. calculatePC_multiprocessing --> this creates a phase coherence object for each earthquake, and calculates the phase coherence. Uses PC_calculator and PhaseCoherence
2. consolidate_detections --> this takes phase coherence detections from each earthquakes objects, and concatenates them into lists, based on detection limit. 
