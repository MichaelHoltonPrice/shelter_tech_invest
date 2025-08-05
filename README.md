# Homeward Bound: Exploring Decisions Involving Shelter and Transport

This repository contains open-source code to generate the main results and plots for the article "Homeward Bound: Exploring Decisions Involving Shelter and Transport", by Martin H. Welker and Michael H. Price.

## Creating the results and plots

### (1) Clone the repository

```bash
git clone https://github.com/MichaelHoltonPrice/shelter_tech_invest
```

### (2) Install dependencies

Open R. Then, if necessary, install the "MASS" and "triangle" packages:

```r
install.packages(c("MASS", "triangle"))
```

### (3) Run the script to make publication results

Run the following script to generate all of the plots and other results:

```r
source('make_publication_results.R')
```