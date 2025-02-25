# Redshift-Classifier (X-ray)
<!--[![Status of Build and Tests Workflow](https://github.com/aritraghsh09/GaMPEN/actions/workflows/main.yml/badge.svg)](https://github.com/aritraghsh09/GaMPEN/actions/workflows/main.yml)-->
[![Documentation Status]()
[![R Version 4.3](https://img.shields.io/badge/R-4.3-blue)](https://cran.r-project.org/)
[![Python Version 3.11](https://img.shields.io/badge/Python-3.8-blue)](https://www.python.org/downloads/)
[![GitHub license](https://img.shields.io/github/license/Milind018/Redshift-Classifier)](https://github.com/Milind018/Redshift-Classifier/blob/main/LICENSE)
[![image](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
<!--[![Code DOI](https://zenodo.org/badge/841602738.svg)](https://doi.org/10.5281/zenodo.14704737)-->
<!--[![Publication DOI](https://img.shields.io/badge/publication%20doi-10.3847%2F1538--4357%2Fac7f9e-blue)](https://doi.org/10.3847/1538-4357/ac7f9e)-->
[![arXiv](https://img.shields.io/badge/arXiv-2410.13985-blue)](https://arxiv.org/pdf/2410.13985)

Redshift-estimator (X-ray) is a reliable ensemble machine learning (ML) framework developed to estimate the redshift of X-ray gamma-ray bursts (GRBs) based on the GRB features. We have used all the GRB light-curve (LC) phases in this framework, prompt, X-ray plateau, and X-ray afterglow, as features to train and test the models.
This ML framework is developed to predict the redshift of a GRB as soon as a new GRB is detected, which will help other facilities to have quick and more precise follow-up observations. In this ML framework, we have also incorporated ML techniques like outlier removal (removing bad data points), feature selection (to select the best features for redshift regression and to reduce the dimensionality), imputing missing data (to increase the sample size) improve the performance of the ML model.
Our ML model is more reliable because we use an ensemble model (combining many ML algorithms in a single framework using SuperLearner) rather than traditional methods that use a single ML algorithm. This further increases the prediction accuracy.
# Documentation
The Redshift-estimator's (X-Ray) documentation is available in this repository and also hosted
on [readthedocs.io](https://grb-web-app.readthedocs.io/en/latest/).
Although the documentation is fairly complete, if you are trying to use the Redshift-estimator (X-ray) Web app and run into issues,
please get in touch with us!
# Publications
The Redshift-estimator (X-Ray) was initially introduced in the following publication. Please cite this publication if you make use of the Redshift-estimator (X-ray) Web app or some code herein.
``` tex
@ARTICLE{2024arXiv241013985N,
       author = {{Narendra}, Aditya and {Dainotti}, Maria and {Sarkar}, Milind and {Lenart}, Aleksander and {Bogdan}, Malgorzata and {Pollo}, Agnieszka and {Zhang}, Bing and {Rabeda}, Aleksandra and {Petrosian}, Vahe and {Kazunari}, Iwasaki},
        title = "{GRB Redshift Estimation using Machine Learning and the Associated Web-App}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2024,
        month = oct,
          eid = {arXiv:2410.13985},
        pages = {arXiv:2410.13985},
          doi = {10.48550/arXiv.2410.13985},
archivePrefix = {arXiv},
       eprint = {2410.13985},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv241013985N},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
# License
Copyright 2025 Maria Giovanna Dainotti, Malgorzata Bogdan, Aditya Narendra, Agnieszka Pollo, Milind Sarkar, Aleksander Lenart, Aleksandra Rabeda & contributors
Made available under a [GNU GPL
v3.0](https://github.com/Milind018/Redshift-Classifier/blob/main/LICENSE)
license.
# Contributors
The Redshift-estimator (X-ray) Ensemble Learning Framework and Web App were initially developed by Maria Giovanna Dainotti, Malgorzata Bogdan, Aditya Narendra, Agnieszka Pollo, Milind Sarkar, Aleksander Lenart, Aleksandra Rabeda & contributors.
The initial documentation was developed by [Milind Sarkar](https://milind018.github.io/), and [Aditya Narendra](https://github.com/AdityaNarendra).
# Getting Help/Contributing
If you have a question, please send me an e-mail at this
`mariagiovannadainotti@xxxxx.it` yahoo address.
If you have spotted a bug in the code/documentation or you want to
propose a new feature, please feel free to open an issue/a pull request
on GitHub.











Message Aditya Narendra:spiral_calendar_pad:In a meeting • Outlook Calendar




