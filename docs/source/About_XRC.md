# About

## Why was Redshift-Classifier (X-Ray) developed?
Gamma-ray bursts (GRBs) are essential astrophysical objects that can be used to understand the evolution of the early Universe since they can be detected up to a redshift ~ 20. However, the afterglows of GRBs dim extremely rapidly; hence, performing high−redhsift measurements poses challenges because even observations of optically bright GRBs are hindered by reduced telescope time and the limited number of follow-up programs for GRBs.

Hence, the requirement is to develop a system that can quickly identify whether a newly detected GRB event is observed at high redshift to allow follow-up observations in other wavelengths.

Previously, supervised machine learning (ML) techniques were also used to identify high−z GRBs, utilizing only a single algorithm to train the model. However, they did not include the plateau emission phase while training their model. Hence, this framework includes X-ray plateaus and utilizes an ensemble method, combining multiple models within the SuperLearner framework to improve the prediction accuracy and power of the ML framework.



In order to address these above challenges, we developed the Redshift-Classifier (X-Ray).

:::{admonition} Redshift-Classifier (X-Ray) Feature Summary:-
:class: note

1. Redshift-Classifier (X-Ray) is a redshift-based GRB classifier to classify GRBs as high-redshift or low-redshift based on the user-defined redshift.

    * It generates the scatter matrix plot for the input dataset.
    
    * It takes into account outlier removal (bad data points) using the M-estimator technique. It generates the histogram plot showing the weights assigned by the M-estimator. It also generates the scatter matrix plot for the filtered data, showing the data points that are outliers.

    * It also takes into account feature selection to select the best features using the feature selection method, Least Absolute Shrinkage and Selection Operator (LASSO).  It generates the plot showing the weights assigned by LASSO to select the best features based on the user-defined cutoff.
      
    * It also incorporates missing data imputation using the Multiple Imputation by Chained Equations (MICE) technique. It generates the distribution plot showing the missing data points (what number of data points are missing and for what features). It also generates the scatter matrix plot for the imputed data, showing the data points that are imputed using MICE. It also creates a histogram plot showing the imputed and the original data.

    * It further incorporates balancing the dataset using the Synthetic Minority Over-sampling Technique (SMOTE) technique. It saves the new synthetically generated dataset and creates a redshift distribution plot.

2. It performs the ML training using the SuperLearner framework to combine multiple models at the same time. The SuperLearner allows us to see which algorithm works the best and lets us choose the algorithm based on our desired cutoff. It works on a 100-nested loop 10-fold cross-validation technique.

    * It generates the plots showing the weights assigned by SuperLearner to each algorithm to select the best algorithms based on user-defined cutoff.
  
    * It creates the confusion matrix, a classification matrix used to evaluate the results and the classifier's performance.

    * It generates the Receiver Operating Characteristic (ROC) curves showing the Area Under the Curve (AUC) for each algorithm and the combined SuperLearner model.
      
    * It also generates Precision/Recall and Accuracy plots.
:::


 

Don't hesitate to contact us if you want our help/advice in using Redshift-Classifier (X-Ray)! 



## Publications
The Redshift-Classifier(X-Ray) was initially introduced in [Dainotti et al. 2025](https://arxiv.org/abs/2408.08763). Please cite this publication if you make use of the Redshift-Classifier (X-ray) Web app or some code herein.


## Attribution Info.

Please cite the below mentioned publication if you make use of Redshift-Classifier (X-Ray) or some code herein.

``` tex
@misc{dainotti2025grbredshiftclassifierfollowup,
      title={GRB Redshift Classifier to Follow-up High-Redshift GRBs Using Supervised Machine Learning}, 
      author={Maria Giovanna Dainotti and Shubham Bhardwaj and Christopher Cook and Joshua Ange and Nishan Lamichhane and Malgorzata Bogdan and Monnie McGee and Pavel Nadolsky and Milind Sarkar and Agnieszka Pollo and Shigehiro Nagataki},
      year={2025},
      eprint={2408.08763},
      archivePrefix={arXiv},
      primaryClass={astro-ph.HE},
      url={https://arxiv.org/abs/2408.08763}, 
}
```


## Getting Help/Contributing

If you have a question, please send me an e-mail at this
`mariagiovannadainotti@xxxxx.it` yahoo address.

If you have spotted a bug in the code/documentation or you want to
propose a new feature, please feel free to open an issue/a pull request
on GitHub.
