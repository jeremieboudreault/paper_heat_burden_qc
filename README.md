Estimating the heat-related mortality and morbidity burden in Quebec, Canada
================================================================================

This is a code example for the implementation presented in the paper "*Estimating the heat-related mortality and morbidity burden in Quebec, Canada*" to compute the heat burden for multiple health outcomes.

- By [Jérémie Boudreault](https://jeremieboudreault.github.io/), Éric Lavigne, Céline Campagna and Fateh Chebana

---

### Data

Datasets are located in `data/` folder :

- `health_data_synthetic.csv` : Daily count in each health region (RSS) for each health outcome (HO)*
- `weather_data.csv` : Daily lagged value of temperature, humidity (and air pollution) in each health region (RSS)
- `meta_predictors.csv` : Meta-predictors for each health region (RSS) and each health outcome (HO)
- `heat_thresholds.csv` : Extreme temperature thresholds for the heat burden quantification by health region (RSS)

*Real health data cannot be shared. So we replaced the count of each health outcomes with simulations from a Poisson distribution with λ=100. Obtained results will not make sence, but at least the user can run the code and understand the data struture.

---

### Analysis

The analyses are all performed in the `main.R` script. The script is separated in the three main steps :

- **Step 0** : Data preparation
- **Step 1.1** : DLNM fitting by region
- **Step 1.2** : Heat burden computation with regional DLNM
- **Step 2.1** : BLUP from meta-regression
- **Step 2.2** : Heat burden computation with regional BLUP
- **Step 3.1** : Pooled effect across the province
- **Step 3.2** : Heat burden computation with pooled effect

Helper functions are located in the `R/` folder.

---


***Enjoy !***
