
# Usage
On this page, we go over the most important user-facing functions that the Redshift-Claasifier (X-ray) web app has and what each of these functions does. Read this page carefully to understand the various arguments/options that can be set while using these functions.


As discussed in [Getting Started](./Getting_Started.md), copy the files and run the app.py file in an IDE (preferably, VSCode) and write the following in the terminal: 
```streamlit run app.py```

Below is the illustration of our web app interface showing how it works.

![Introductory Image](./../../images/Classifier_home.png)


```{Attention}
To start using the Redshift-Classifier and for first-hand testing, the user only needs to use the ``NEW-XRAY_DATA_RAW_w_errorbar_WITHOUT-M-estimator.csv`` data file from the [data](https://github.com/Milind018/Redshift-Classifier/tree/main/data) folder. However, we have already provided all the datasets in the [data](https://github.com/Milind018/Redshift-Classifier/tree/main/data) folder that we have used in this Redshift-Classifier, which user can also use to test each code separately in their local machine.
```

```{Note}
Our Web-App allows users to customize the classification process based on their specific requirements. Users can select whether they want to remove outliers or not by using the M-estimator based on their desired threshold. Users can also opt to impute the missing variables using MICE imputation and data balancing using SMOTE. Additionally, user can select their preferred redshift cutoff for a flexible GRB redshift classification.
```

## References
* LASSO feature selection: [here](https://javatpoint.com/feature-selection-techniques-in-machine-learning#:~:text=Feature%20selection%20is%20a%20way,%2C%20irrelevant%2C%20or%20noisy%20features), [here](https://www.analyticsvidhya.com/blog/2020/10/feature-selection-techniques-in-machine-learning/), [here](https://medium.com/@agrawalsam1997/feature-selection-using-lasso-regression-10f49c973f08), [here](https://towardsdatascience.com/feature-selection-in-machine-learning-using-lasso-regression-7809c7c2771a), and [here](https://tahera-firdose.medium.com/lasso-regression-a-comprehensive-guide-to-feature-selection-and-regularization-2c6a20b61e23).
* M-estimator: [here](https://www.statisticshowto.com/m-estimator/) and [here](https://ui.adsabs.harvard.edu/abs/2022FrASS...936215G/abstract).
* MICE: [here](https://medium.com/@brijesh_soni/topic-9-mice-or-multivariate-imputation-with-chain-equation-f8fd435ca91#:~:text=MICE%20stands%20for%20Multivariate%20Imputation,produce%20a%20final%20imputed%20dataset), [here](https://cran.r-project.org/web/packages/midastouch/midastouch.pdf), and [here](https://ui.adsabs.harvard.edu/abs/2022FrASS...936215G/abstract).
* SMOTE: [here](https://arxiv.org/abs/1106.1813) and [here](https://github.com/dalpozz/unbalanced/tree/master/R).
* SuperLearner: [here](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html#:~:text=SuperLearner%20is%20an%20algorithm%20that,using%20the%20test%20data%20performance.) and [here](https://cran.r-project.org/web/packages/SuperLearner/index.html).
* Cross-validation: [here](https://towardsdatascience.com/cross-validation-in-machine-learning-72924a69872f) and [here](https://www.geeksforgeeks.org/cross-validation-machine-learning/).
* Confusion Matrix: [here](https://changjunlee.com/blogs/posts/4_confusion_mat_and_roc), [here](https://www.digitalocean.com/community/tutorials/confusion-matrix-in-r), and [here](https://www.analyticsvidhya.com/blog/2020/09/precision-recall-machine-learning/).
* AUC: [here](https://developers.google.com/machine-learning/crash-course/classification/roc-and-auc#:~:text=AUC%20stands%20for%20%22Area%20under,Area%20under%20the%20ROC%20Curve), [here](https://towardsdatascience.com/understanding-auc-roc-curve-68b2303cc9c5), and [here](https://lexjansen.com/nesug/nesug10/hl/hl07.pdf).
* ROC: [here](https://changjunlee.com/blogs/posts/4_confusion_mat_and_roc), [here](https://www.digitalocean.com/community/tutorials/plot-roc-curve-r-programming), [here](https://developers.google.com/machine-learning/crash-course/classification/roc-and-auc#:~:text=AUC%20stands%20for%20%22Area%20under,Area%20under%20the%20ROC%20Curve), [here](https://search.r-project.org/CRAN/refmans/AUC/html/plot.AUC.html), and [here](https://lexjansen.com/nesug/nesug10/hl/hl07.pdf).


## Working of the Web-App

1. First, the user has to select if they want to use the M-estimator or not to remove the outliers.

   a) If selected ``No``, outliers will not be removed.

   b) If selected ``Yes``, outliers will be removed. Then, the user has to select the weight cutoff to remove the outliers.

2. Second, the user has to select if they want to impute the missing variables or not using MICE.

   a) If selected ``No``, missing data will not be imputed.

   b) If selected ``Yes``, missing data will be imputed.

3. Then, the user has to select if they want to use SMOTE to balance the dataset or not.
   
   a) If selected ``No``, the data will not be balanced.

   b) If selected ``Yes``, the will be balanced.
   
4. Then, the user needs to choose the redshift cutoff according to their desired interest to define the high-redshift and low-redshift GRBs.

5. Finally, the user needs to upload the raw data file ``NEW-XRAY_DATA_RAW_w_errorbar_WITHOUT-M-estimator.csv`` data file from the [data](https://github.com/Milind018/Redshift-Classifier/tree/main/data) folder to perform the GRB redshift classification.
