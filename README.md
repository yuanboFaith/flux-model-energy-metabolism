<h1>Flux Modeling of Mammalian Energy Metabolism</h1>


This repository contains *all* source code and raw data to reproduce the figures and computation in the paper ***An Organism-Level Quantitative Flux Model of Mammalian Energy Metabolism***. 🌟

The source code is prefixed with numbers **1** to **12**. Run each script sequentially by the prefix. Running each file will generate an associated **.RData** file, which will be imported in the following script.

### 13CO2 data analysis

The raw data from the metabolic cage-isotope analyzer is in folder **data-13CO2**. Experimental parameters are in **data_13CO2_mouseID_paper.xlsx**. 

- **1_core_13CO2_import_data.R:** imports and clean up all 13CO2 data.
- **2_core_13CO2_bolus_injection.R:** analyzes bolus-based 13CO2 data.
- **3_core_13CO2_infusion_physiology_basics.R:** analyzes infusion-based 13CO2 data, and integrates all other respirometry data, such as total O2 consumption, CO2 production, RER, energy expenditure, etc. 

### Serum labeling and integrated analysis

The serum labeling raw data is in **data_serum labeling.xlsm**, and the experimental parameters are in **data_serum labeling mouse ID.xlsx**.

- **4_core_serum_labeling_import data.R:** imports and cleanup all 13C-serum labeling data, and performs natural abundance correction. It employs a homemade mini-package **vivoFlux** (leveraging Xiaoyang's **[Accucor](https://github.com/XiaoyangSu/AccuCor)**)

- **5_core_labeling_analysis.R:** analyzes serum labeling data, calculate arterial-venous labeling difference, and computes circulatory turnover fluxes.
- **6_core_consumption_flux.R:** calculates the consumption (sink) fluxes, including the oxidation and storing fluxes, and the *direct* fluxes. 
- **7_core_production_flux.R:** calculates the storage-releasing fluxes, and inter-converting fluxes. 
- **8_core_flux extrapolation.R:** an integrated analysis, and performs statistical test of metabolic fluxes between phenotypes.
- **9_core_ATP_estimate.R:** calculates the ATP wasted by the futile cycles. 
- **10_core_draw flux network.R:** draws the metabolic flux network for all three phenotypes. The variant **11** maps flux values (significantly different between phenotypes) to color scale; and **12** creates *interactive* networks. 

---

**All scripts have been checked to run smoothly. In case of error:**

- Check path of file import
- Check if you have generated the needed .RData file
- Restart R session to clear package conflict
- Check you have the needed package
- Did you change the raw data? This may cause disruption. 


<br>
