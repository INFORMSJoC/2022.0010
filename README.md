[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Data and Codes for A Column Generation Scheme for Distributionally Robust Multi-Item Newsvendor Problems

This repository includes the data and the codes for the following paper:
Wang S, Delage E. A Column Generation Scheme for Distributionally Robust Multi-Item Newsvendor Problems. INFORMS Journal on Computing, 2023.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.
```
@article{contdecbeh,
  author =        {S. Wang and E. Delage},
  publisher =     {INFORMS Journal on Computing},
  title =         {A column generation scheme for distributionally robust multi-item newsvendor problems},
  year =          {2023},
  url =           {https://github.com/INFORMSJoC/2022.0010},
} 
```

## Data files
Collection of multi-item newsvendor instances in text format. 
Each file contains 10 instances.
Each instance name is in the oNm_k format, where N is the number of scenarios, m denotes the instances number, and k is optional and represents the number of items. There are 10 items without k.
The first row is unit selling revenue.
The second row is the upper bound of demands.
The third row to the end is the scenarios of demands.

## Code files
There are four folders. 
ew_mad: it includes all the different algorithms described in "A Column Generation Scheme for Distributionally Robust Multi-Item Newsvendor Problems" to solve the multi-item newsvendor problem under the EW-MAD ambiguity set.

wasserstein: it includes all the different algorithms to solve the multi-item newsvendor problem under the Wasserstein ambiguity set.

oos_performance: it solves the DR and SAA multi-item newsvendor models in a data-driven environment described in Section 5.2.

functions: it includes all the functions used in folders "ew_mad", "wasserstein", and "oos_performance".

