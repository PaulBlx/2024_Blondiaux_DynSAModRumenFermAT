
This folder contains the R implementation of a dynamic sensitivity analysis of a mathematical model describing the effect of the macroalgae Asparagopsis taxiformis on rumen fermentation and methane production under in vitro continuous conditions.

The mathematical model was modified based on the in vitro batch condition model described in  Muñoz-Tamayo et al. (2021). 

The sensitivity analysis methods implemented were the Shapley effects (Owen, 2014; Song et al., 2016) using the sensitivity R package (Iooss et al., 2023), and the full and independent Sobol indices (Mara et al., 2015) using the sensobol R package (Puy et al., 2022). For both methods, the scripts were runned using the MESO@LR-Platform by splitting the sampling matrix into 20 parts due to the high computation time. The computational time for one dietary scenario was of 24h for the Shapley effects and 5h for the Sobol indices using the MESO@LR-Platform. The scripts in the folder do not consider this split due to the numerous number of scripts considered (one script for each part of the splitted sampling matrix). 

Very importantly: To use these implementations, you must respect the conditions detailed in the document 2024_Agreement_INRAE.pdf

References

Bertrand Iooss, A., Da Veiga, S., Janon, A., Pujol, G., contribu-tions from Baptiste Broto,  with, Boumhaout, K., Clouvel, L., De-lage, T., El Amri, R., Fruth, J., Gilquin, L., Guillaume, J., Herin, M., Il Idrissi, M., Le Gratiet, L., Lemaitre, P., Marrel, A., Mey-naoui, A., Nelson, B.L., Monari, F., Oomen, R., Rakovec, O., Ramos, B., Rous-tant, O., Sarazin, G., Song, E., Staum, J., Sueur, R., Touati, T., Verges, V., Weber, F., 2023. Package “sensitivity” Title Global Sensitivity Analysis of Model Outputs.

Mara, T.A., Tarantola, S., Annoni, P., 2015. Non-parametric methods for global sensitivity analysis of model output with dependent inputs. Environ. Model. Softw. 72. https://doi.org/10.1016/j.envsoft.2015.07.010

Muñoz-Tamayo, R., Chagas, J.C., Ramin, M., Krizsan, S.J., 2021b. Modelling the impact of the macroalgae Asparagopsis taxiformis on rumen microbial fermentation and methane production. Peer Community J. 1. https://doi.org/10.24072/pcjournal.11

Owen, A.B., 2014. Sobol’ indices and shapley value. SIAM-ASA J. Uncertain. Quantif. 2. https://doi.org/10.1137/130936233

Puy, A., Piano, S. Lo, Saltelli, A., Levin, S.A., 2022. sensobol : An R Package to Compute Variance-Based Sensitivity Indices. J. Stat. Softw. 102. https://doi.org/10.18637/jss.v102.i05

Song, E., Nelson, B.L., Staum, J., 2016. Shapley effects for global sensitivity analysis: Theory and computation. SIAM-ASA J. Uncertain. Quantif. 4. https://doi.org/10.1137/15M1048070




