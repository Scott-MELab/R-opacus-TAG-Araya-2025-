**This repository contains:**

- **An Excel file** (CoreModel_Summary.xlsx) summarizing the reactions involved in the central metabolism model for TAG accumulation in *R. opacus* DSM43205. For each reaction, the GPR rule constructed from the genome available at NCBI (RefSeq assembly GCF_001646735.1, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001646735.1/), the reaction directionality, and the enzyme commission (EC) number are provided.

- **The Python script** (CoreModel_RO_Araya2025.py) developed for the simulation of FBA (Flux Balance Analysis).
Note: By default, the synthesis of TAGs from glucose is maximized, assuming a specific glucose uptake rate of 1.94 mmol/gCDW·h (experimentally determined).

- **The neutral model** built in JSON format (Neutral MODEL DSM4305.json) is presented in this repository. This version does not include any objective function for TAG maximization or substrate uptake constraints used in this work. Additionally, the corresponding map (Neutral Core Model Map.json) is also provided. This model can be used to evaluate other substrates by loading both the JSON model and the map into Escher (https://sbrg.github.io/escher-fba/#/).
Note: To perform this analysis, both files must be uploaded into Escher, and the desired TAG synthesis function can be maximized by setting the appropriate uptake rate for the chosen substrate.
