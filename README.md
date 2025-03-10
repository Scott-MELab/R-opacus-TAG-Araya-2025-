*Rhodococcus opacus*  DSM 43205 as a bacterial chassis for triacylglycerols accumulation using carbon dioxide, hydrogen, and derived organic substrates derived.

Blanca Araya<sup>1</sup>, Paz Torres-Praderio<sup>1</sup>, Raúl Conejeros<sup>2</sup>, Alberto Vergara-Fernández<sup>1</sup>, Felipe Scott<sup>1,*</sup> 

<sup>1</sup>Green Technologies Research Group, Facultad de Ingeniería y Ciencias Aplicadas, Universidad de los Andes, Chile. Av. Mons. Álvaro del Portillo 12.455. Las Condes, Santiago, Chile.

<sup>2</sup>School of Biochemical Engineering, Pontificia Universidad Católica de Valparaíso, Av. Brasil 2085, Valparaíso, Chile

<sup>*</sup> Corresponding author: Felipe Scott. Tel.+ 562 2618 1909. E-mail: fscott@uandes.cl.


**This repository contains:**

- **An Excel file** (CoreModel_Summary.xlsx) summarizing the reactions involved in the central metabolism model for TAG accumulation in *R. opacus* DSM43205. For each reaction, the GPR rule constructed from the genome available at NCBI (RefSeq assembly GCF_001646735.1, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001646735.1/), the reaction directionality, and the enzyme commission (EC) number are provided.

- **The Python script** (CoreModel_RO_Araya2025.py) developed for the simulation of FBA (Flux Balance Analysis).
Note: By default, the synthesis of TAGs from glucose is maximized, assuming a specific glucose uptake rate of 1.94 mmol/gCDW·h (experimentally determined).

- **The neutral model** built in JSON format (Neutral MODEL DSM4305.json) is presented in this repository. This version does not include any objective function for TAG maximization or substrate uptake constraints used in this work. Additionally, the corresponding map (Neutral Core Model Map.json) is also provided. This model can be used to evaluate other substrates by loading both the JSON model and the map into Escher (https://sbrg.github.io/escher-fba/#/).
Note: To perform this analysis, both files must be uploaded into Escher, and the desired TAG synthesis function can be maximized by setting the appropriate uptake rate for the chosen substrate.
- **Maps folder** contains the metabolic flux distribution maps for TAG accumulation for each substrate evaluated in this study.
