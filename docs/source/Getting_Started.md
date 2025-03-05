```{attention}
Note that all the .Rmd and .R files in the primary and ancillary subfolders of the scripts folder must be copied into the working directory where the app.py file is present!!
```

# Getting Started

Redshift-Estimator (X-Ray) is written in R and relies on the [Streamlit]([https://pytorch.org/](https://streamlit.io/)) and Python for its GUI Web app.

## Installation

To set up the project and ensure all dependencies are installed, ensure that you have access to
* **R** (version >= 4.3) with RStudio
* **Python** (version >= 3.11)

1. Locally clone the project repository from GitHub or directly download the .zip file:
```
git clone https://github.com/gammarayapp/GRB-Web-App.git 
```

2. Install the necessary R packages, open and run `Package_install.Rmd` in RStudio. Alternatively, use the `install.packages()` command in R with `packages.txt` and `https://cran.r-project.org` as the base URL of repositories.

This project also requires the Python libraries **Streamlit**, **Pandas**, and **Pillow**, which can be installed with `pip install -r requirements.txt`
