# Analysis of bone strength growth curves in the Fels longitudinal study (FLS)
## Code for reproducing analyses

This repository contains all code for reproducing the analyses in "name of paper" by [list authors]. To request the FLS data, please contact [insert name and contact info].

The code is organized as follows:

* Main folder
    + prep_data.R: preparing data
    + FLS_stats.R: obtaining basic sample statistics
    + predict*.R: functions used to obtain predicted trends, derivatives, and percent change
* bodysize: code for determining the regression model for centering bodysize. Parts of this code are incorporated into prep_data.R in the main folder. The files in the bodysize folder are not used in subsequent analysis.
* joR
    + code (listed in order that files needs to be run)
        + exploratory_joR.R: assessing how well the gamma distribution describes BSI
        + joR_modelFit*.R: fitting models with centered and scaled bodysize
        + joR_modelEval*.R: evaluating models with centered and scaled bodysize. These files produce the figures shown in the paper.
        + impute_phvage_single.R: fitting and evaluating single imputation model (produced plots)
        + impute_phvage_multiple_post.R: evaluate the results from the multiple imputations (produces plots)
    + imputations: batch code for running imputation models. This code must be run before running joR/code/impute_phvage_multiple_post.R
* cti
    + code (listed in order that files needs to be run)
        + cti_modelFits_cent.R: fitting models
        + cti_modelEval_cent.R: evaluating models. These files produce the figures shown in the paper.
        + impute_phvage_multiple_post.R: evaluate the results from the multiple imputations (produces plots)
    + imputations: batch code for running imputation models. This code must be run before running cti/code/impute_phvage_multiple_post.R
* ttar
    + code (listed in order that files needs to be run)
        + ttar_modelFits_cent.R: fitting models
        + ttar_modelEval_cent.R: evaluating models. These files produce the figures shown in the paper.
        + impute_phvage_multiple_post.R: evaluate the results from the multiple imputations (produces plots)
    + imputations: batch code for running imputation models. This code must be run before running ttar/code/impute_phvage_multiple_post.R

