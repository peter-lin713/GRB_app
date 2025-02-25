```{toctree}
:maxdepth: 2
:hidden:

self
About_XRC
Getting_Started
Using_XRC
```

# Home

```{attention}
Note that although Redshift-Classifier (X-Ray) current documentation is fairly substantive, we are still working on some parts of the documentation. If you run into issues while trying to use Redshift-Classifier (X-Ray), please contact us! We will be more than happy to help you!
```


***

Redshift-Classifier (X-ray) is a reliable ensemble machine learning (ML) framework developed to classify X-ray gamma-ray bursts (GRBs) as high-redshift or low-redshift based on the desired redshift cutoff (user-defined). We have used all the GRB light-curve (LC) phases in this framework, prompt, X-ray plateau, and X-ray afterglow, as features to train and test the models. However, our model is developed to use any set of features based on the availability and user's interest.

This ML framework is developed to quickly identify if a GRB is a high-redshift or a low-redshift as soon as a new GRB is detected, which will help other facilities to have quick and more precise follow-up observations. In this ML framework, we have also incorporated ML techniques like outlier removal (removing bad data points), feature selection (to select the best features for redshift classification and to reduce the dimensionality), imputing missing data (to increase the sample size), and also balancing the data (to not have a biased ML trained model), to improve the performance of the ML model.

Our ML model is more reliable because we use an ensemble model (combining many ML algorithms in a single framework using SuperLearner) rather than traditional methods that use a single ML algorithm. This further increases the prediction accuracy.


## First Steps with Redshift-Classifier (X-Ray)
0. For a quick blog-esque introduction to the most important features of Redshift-Classifier(X-Ray), please check out [About](./About_XRC.md).
:::{tip}
For a deep-dive, Please refer to [Dainotti et al. 2025](https://arxiv.org/abs/2408.08763).
:::
1. Follow the installation instructions and quick-start guide in [Getting Started](./Getting_Started.md).
2. Review the [Usage](./Using_XRC.md) page to dive into the details about the various user-facing functions that our classifier provides and how to use them in a GUI environment.
