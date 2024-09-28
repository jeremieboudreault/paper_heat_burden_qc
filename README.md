Estimating the heat-related mortality and morbidity burden in Quebec, Canada
================================================================================

This is a code example to compute the heat burden for multiple health outcomes presented in the paper "[***Estimating the heat-related mortality and morbidity burden in Quebec, Canada***](https://doi.org/10.1016/j.envres.2024.119347)" published in *Environmental Research* in 2024.

- By [Jérémie Boudreault](https://jeremieboudreault.github.io/), Éric Lavigne, Céline Campagna and Fateh Chebana

---

### Data

Datasets are located in `data/` folder :

- `health_data_synthetic.csv` : Daily count* in each health region (RSS) for each health outcome (HO)
- `weather_data.csv` : Daily lagged values of weather and air pollution in each health region (RSS)
- `meta_predictors.csv` : Meta-predictors for each health region (RSS) and each health outcome (HO)
- `heat_thresholds.csv` : Extreme temperature thresholds (Q95) by health region (RSS)

*Real health data cannot be shared. So we replaced the count of each health outcomes with simulations from a Poisson distribution with λ=100. Obtained results will not make sense, but at least the user can run the code and understand the data structure.

---

### Analysis

The analyses are all performed in the `main.R` script. The script is separated in three main steps :

- **Step 1** : DLNM fitting by region and heat burden computation
- **Step 2** : BLUP from meta-regression and heat burden computation
- **Step 3** : Pooled effect across the province  and heat burden computation

Helper functions are located in the `R/` folder.

---

***Enjoy !***
