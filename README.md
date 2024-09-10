# SPRFittingPaper2023 - Tools for fitting bivalent antibody SPR assays
[![Documentation](https://img.shields.io/badge/documentation-blue.svg)](https://isaacsas.github.io/SPRFittingPaper2023.jl/dev/)
[![DOI](https://zenodo.org/badge/687753865.svg)](https://zenodo.org/doi/10.5281/zenodo.13743460)

Library for estimating kinetic parameters and molecular reach for surface plasmon resonance (SPR) assays via fitting to a spatial, stochastic, particle-based model. Methodology and results are available in the accompanying manuscript [here](https://doi.org/10.1101/2023.09.06.556503 ) (citation info below).
Please see the [documentation](https://isaacsas.github.io/SPRFittingPaper2023.jl/dev/) for how to 
1. Run forward simulations.
2. Construct surrogate models.
3. Fit surrogate models to SPR data to estimate kinetic parameters and molecular reach.


## Supporting and citing SPRFittingPaper2023
The software in this ecosystem was developed as part of academic research. If you would like to help
support it, please star the repository as such metrics may help us secure funding in the future. If
you use it as part of your research, teaching, or other activities, we would be grateful if you
could cite our work:
```
@article {Huhn2023.09.06.556503,
	author = {Huhn, Anna and Nissley, Daniel and Wilson, Daniel B. and Kutuzov, Mikhail and Donat, Robert and Tan, Tiong Kit and Zhang, Ying and Barton, Michael I. and Liu, Chang and Dejnirattisai, Wanwisa and Supasa, Piyada and Mongkolsapaya, Juthathip and Townsend, Alain and James, William and Screaton, Gavin and van der Merwe, P. Anton and Deane, Charlotte M. and Isaacson, Samuel A. and Dushek, Omer},
	title = {Analysis of emergent bivalent antibody binding identifies the molecular reach as a critical determinant of SARS-CoV-2 neutralisation potency},
	elocation-id = {2023.09.06.556503},
	year = {2024},
	doi = {10.1101/2023.09.06.556503},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Key functions of antibodies, such as viral neutralisation, depend on bivalent binding but the factors that influence it remain poorly characterised. Here, we develop and employ a new bivalent model to mechanistically analyse binding between \&gt;45 patient-isolated IgG1 antibodies interacting with SARS-CoV-2 RBD surfaces. Our method reproduces the monovalent on/off-rates and enables measurements of the bivalent on-rate and the molecular reach: the maximum antigen separation that supports bivalent binding. We find large variations in these parameters across antibodies, including variations in reach (22-46 nm) that exceed the physical antibody size (\~{}15 nm) due to the antigen size. The bivalent model integrates all parameters, including reach and antigen density, to predict an emergent binding potency for each antibody that matches their neutralisation potency. Indeed, antibodies with similar monovalent affinities to the same RBD-epitope but with different reaches display differences in emergent bivalent binding that match differences in their neutralisation potency. Together, our work highlights that antibodies within an isotype class binding the same antigen can display differences in molecular reach that can substantially modulate their emergent binding and functional properties.Lay Summary Antibodies are soluble proteins that can neutralise pathogens by sticking to them. They contain two identical {\textquoteleft}arms{\textquoteright} that allow them to simultaneously bind two identical {\textquoteleft}antigen{\textquoteright} molecules on pathogen surfaces. Although we know that bivalent binding is important for neutralisation, we don{\textquoteright}t know how different antibodies achieve it. We developed a new model to analyse the mechanism of bivalent binding and used it to study over 45 antibodies from COVID-19 patients that bind the RBD antigen of SARS-CoV-2. Unexpectedly, we found that the molecular reach of an antibody, which is the maximum antigen separation that supports bivalent binding, varied widely between antibodies and exceeded their physical size. We show how antibody binding emerges from the interplay of multiple factors, including reach, and that this emergent binding predicts their neutralisation function. The ability to analyse and predict bivalent binding should improve our understanding and exploitation of antibodies.Competing Interest StatementGRS is on the GSK Vaccines Scientific Advisory Board, a founder shareholder of RQ biotechnology, and a Jenner investigator. Oxford University holds intellectual property related to the Oxford-AstraZeneca vac- cine.},
	URL = {https://www.biorxiv.org/content/early/2024/04/15/2023.09.06.556503},
	eprint = {https://www.biorxiv.org/content/early/2024/04/15/2023.09.06.556503.full.pdf},
	journal = {bioRxiv}
}

```
