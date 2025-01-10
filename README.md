<h1><strong>An Organism-Level Quantitative Flux Model of Energy Metabolism in Mice</strong></h1>

This repository contains the source code and raw data to reproduce the figures and computation in our paper [***An Organism-Level Quantitative Flux Model of Energy Metabolism in Mice***](https://www.biorxiv.org/content/10.1101/2024.02.11.579776v2). 🌟

![](https://github.com/yuanboFaith/flux-model-energy-metabolism/blob/main/github%20summary%20figure.png?raw=true)


The source code is prefixed with numbers **1** to **14**. Run each script sequentially by the prefix. Execution of the so-named core scripts (**1-9**) as well as script **13** generates the associated **.RData** file, which will be imported in the following script.

# **13CO2 data analysis**

The raw data from the metabolic cage-isotope analyzer is in folder **data-13CO2**. Experimental parameters are in **data_13CO2_mouseID_paper.xlsx**. 

- **1_core_13CO2_import_data.R** imports and clean up all 13CO2 data.
- **2_core_13CO2_bolus_injection.R** analyzes bolus-based 13CO2 data.
- **3_core_13CO2_infusion_physiology_basics.R** analyzes infusion-based 13CO2 data, and integrates all other respirometry data, such as total O2 consumption, CO2 production, RER, energy expenditure, etc. 

# **Serum labeling analysis**

The serum labeling raw data is in **data_serum labeling.xlsm**, and the experimental parameters are in **data_serum labeling mouse ID.xlsx**. The following scripts cleans up the labeling data, and performs some simple calculations (e.g., the circulatory turnover flux).  

- **4_core_serum_labeling_import data.R** imports and tidies up all 13C-serum labeling data, and performs natural abundance correction. It employs a homemade mini-package **vivoFlux** (wrapper of **[Accucor](https://github.com/XiaoyangSu/AccuCor)**)

- **5_core_labeling_analysis.R** analyzes serum labeling data, calculate arterial-venous labeling difference, and computes circulatory turnover fluxes.


# **Integrated flux analysis**

The following scripts performs more advanced computations, such as the metabolites interconversion fluxes, and direct oxidation and futile storing fluxes, etc. 

- **6_core_consumption_flux.R** calculates the consumption (sink) fluxes, including the oxidation and storing fluxes, and the *direct* fluxes. 
- **7_core_production_flux.R** calculates the storage-releasing fluxes, and inter-converting fluxes. 
- **8_core_flux extrapolation.R** an integrated analysis, and performs statistical test of metabolic fluxes between phenotypes.
- **9_core_ATP_estimate.R** calculates the ATP wasted by the futile cycles. 
- **10_core_draw flux network.R** draws the metabolic flux network for all three phenotypes. The script **11** maps flux values (with significant difference from the control) to color scale; and **12** creates *interactive* networks. 

# **More analysis...**

The following data analysis is for hypothesis check, and protocol validation. They are not directly related with the model construction per se. 

- **13_validation_experiment.R** performs additional analysis, most of which requested by the peer review.
- **14_export_source_data.R** outputs the data underlying each figure panel into Excel (Data S1 - Source Data). 

----

**All scripts have been checked to run smoothly. In case of error:**

- Check path of file import;
- Check if you have generated the needed .RData file;
- Restart R session to clear package conflict;
- Check you have the needed package;
- Did you change the raw data? This may cause disruption. 


----

**Reference of R code for data wrangling and data visualization:**

- [DataBrewer.co](https://databrewer.co/) 
- [R for Data Science (2e)](https://r4ds.hadley.nz/)


