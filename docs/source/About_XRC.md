# About

## Why was Redshift-Estimator (X-Ray) developed?
Gamma-ray bursts (GRBs) are essential astrophysical objects that can be used to understand the evolution of the early Universe since they can be detected up to a redshift ~ 20. However, the afterglows of GRBs dim extremely rapidly; hence, performing high−redhsift measurements poses challenges because even observations of optically bright GRBs are hindered by reduced telescope time and the limited number of follow-up programs for GRBs.

Hence, the requirement is to develop a system that can provide reliable redshift estimates of GRBs and enable follow-up observations in other wavelengths.

Previously, supervised machine learning (ML) techniques were also used to predict GRB redshifts, utilizing only a single algorithm to train the model. However, they did not include the plateau emission phase while training their model. Hence, this framework includes X-ray plateaus and utilizes an ensemble method, combining multiple models within the SuperLearner framework to improve the prediction accuracy and power of the ML framework.



In order to address these above challenges, we developed the Redshift-Estimator (X-Ray).

:::{admonition} Redshift-Estimator (X-Ray) Feature Summary:-
:class: note

1. Redshift-Estimator (X-Ray) is a ensemble machine learning model designed to predict redshift of GRBs based using photometric GRB properties, obtained from the prompt, plateau and afterglow emissions.

    * It generates the scatter matrix plot for the generalisation data provided by the user.
    
    * It allows the user to enter GRB data manually.

    * It enables MICE based data imputation.
  
    * It also has the option to use a custom trained model of the user.
  
    *  It presents the data in a tabular format.
  
    *  It also presents the box plot with all the redshift predictions.

    

2.  Data Visualisation : This enables the user to visualise their data. 

    * It displays the scatter matrix plot.
  
    * It displays the missing data plot.
  
    * It provides the user the option to upload training or generalisation data.


  
3.  SuperLearner : This section enables the user to generate their own trained models.

    * It enables the user to either apply MICE imputation or not on their training set.
  
    * It enables the user to apply upsampling or not on their training data.
  
    * It enables the user to remove outliers using M-Estimator or not from their training data.
  
    * It also enables the user to use custom models within superlearner framework.
  
    * It can present results with or without catastrophic outliers.
  
4.  GAM/GLM formula generator : It provides the user the ability to generate their own formulas.

    * It enables the user to set a custom correlation cutoff.
  
    * It enables the user to set a custom RMSE cutoff.
  
    * It also enables the user to create formulas either using GAM or GLM.
  
5.  Model Comparison : It enables the user to upload their own data and compare the performance of different machine learning models using the SuperLearner Framework.

    * It enables the user to apply MICE on the data.
  
    * It enables the user to set the number of loops to be executed.
  
    * It enables the user to select their custom models.
  
6.  Relative Importance : This section shows the Relative importance of the features in the user's training data. By default, it shoes the feature importance of our training data.

    
:::


 

Don't hesitate to contact us if you want our help/advice in using Redshift-Estimator (X-Ray)! 



## Publications
The Redshift-Estimator(X-Ray) was initially introduced in [Narendra et al. 2025](https://arxiv.org/abs/2410.13985). Please cite this publication if you make use of the Redshift-Estimator (X-ray) Web app or some code herein.


## Attribution Info.

Please cite the below mentioned publication if you make use of Redshift-Estimator (X-Ray) or some code herein.

``` tex
@misc{narendra2024grbredshiftestimationusing,
      title={GRB Redshift Estimation using Machine Learning and the Associated Web-App}, 
      author={Aditya Narendra and Maria Dainotti and Milind Sarkar and Aleksander Lenart and Malgorzata Bogdan and Agnieszka Pollo and Bing Zhang and Aleksandra Rabeda and Vahe Petrosian and Iwasaki Kazunari},
      year={2024},
      eprint={2410.13985},
      archivePrefix={arXiv},
      primaryClass={astro-ph.HE},
      url={https://arxiv.org/abs/2410.13985}, 
}
```


## Getting Help/Contributing

If you have a question, please send me an e-mail at this
`mariagiovannadainotti@xxxxx.it` yahoo address.

If you have spotted a bug in the code/documentation or you want to
propose a new feature, please feel free to open an issue/a pull request
on GitHub.
