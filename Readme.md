# Redshift-Classifier (X-ray)

<!--[![Status of Build and Tests Workflow](https://github.com/aritraghsh09/GaMPEN/actions/workflows/main.yml/badge.svg)](https://github.com/aritraghsh09/GaMPEN/actions/workflows/main.yml)-->
[![Documentation Status](https://readthedocs.org/projects/gampen/badge/?version=latest)](https://redshift-classifier.readthedocs.io/en/latest/)
[![R Version 4.3](https://img.shields.io/badge/R-4.3-blue)](https://cran.r-project.org/)
[![Python Version 3.11](https://img.shields.io/badge/Python-3.8-blue)](https://www.python.org/downloads/)
[![GitHub license](https://img.shields.io/github/license/Milind018/Redshift-Classifier)](https://github.com/Milind018/Redshift-Classifier/blob/main/LICENSE)
[![image](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Code DOI](https://zenodo.org/badge/841602738.svg)](https://doi.org/10.5281/zenodo.14704737)
<!--[![Publication DOI](https://img.shields.io/badge/publication%20doi-10.3847%2F1538--4357%2Fac7f9e-blue)](https://doi.org/10.3847/1538-4357/ac7f9e)-->
[![arXiv](https://img.shields.io/badge/arXiv-2408.08763-blue)](https://arxiv.org/abs/2408.08763)


Redshift-Classifier (X-ray) is a reliable ensemble machine learning (ML) framework developed to classify X-ray gamma-ray bursts (GRBs) as high-redshift or low-redshift based on the desired redshift cutoff (user-defined). We have used all the GRB light-curve (LC) phases in this framework, prompt, X-ray plateau, and X-ray afterglow, as features to train and test the models. However, our model is developed to use any set of features based on the availability and user's interest.

This ML framework is developed to quickly identify if a GRB is a high-redshift or a low-redshift as soon as a new GRB is detected, which will help other facilities to have quick and more precise follow-up observations. In this ML framework, we have also incorporated ML techniques like outlier removal (removing bad data points), feature selection (to select the best features for redshift classification and to reduce the dimensionality), imputing missing data (to increase the sample size), and also balancing the data (to not have a biased ML trained model), to improve the performance of the ML model.

Our ML model is more reliable because we use an ensemble model (combining many ML algorithms in a single framework using SuperLearner) rather than traditional methods that use a single ML algorithm. This further increases the prediction accuracy.

# Documentation

The Redshift-Classifier's (X-Ray) documentation is available in this repository and also hosted 
on [readthedocs.io](https://redshift-classifier.readthedocs.io/en/latest/). Although the documentation
is fairly complete, if you are trying to use the  Redshift-Classifier (X-ray) Web app and run into issues, 
please get in touch with us!

# Publications
The Redshift-Classifier(X-Ray) was initially introduced in the following publication. Please cite this publication if you make use of the Redshift-Classifier (X-ray) Web app or some code herein.

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
# License

Copyright 2025 Maria Giovanna Dainotti, Malgorzata Bogdan, Aditya Narendra, Agnieszka Pollo, Shubham Bhardwaj, Milind Sarkar, Joshua Ange, Christopher Cook & contributors

Made available under a [GNU GPL
v3.0](https://github.com/Milind018/Redshift-Classifier/blob/main/LICENSE)
license.

# Contributors

The Redshift-Classifier(X-ray) Ensemble Learning Framework and Web App was initially developed by Maria Giovanna Dainotti, Malgorzata Bogdan, Aditya Narendra, Agnieszka Pollo, Shubham Bhardwaj, Milind Sarkar, Joshua Ange, Christopher Cook & contributors.

The initial documentation was developed by [Shubham Bhardwaj](https://github.com/QGravityGRGW), [Milind Sarkar](https://milind018.github.io/), and Joshua Ange.

# Getting Help/Contributing

If you have a question, please send me an e-mail at this
`mariagiovannadainotti@xxxxx.it` yahoo address.

If you have spotted a bug in the code/documentation or you want to
propose a new feature, please feel free to open an issue/a pull request
on GitHub.
