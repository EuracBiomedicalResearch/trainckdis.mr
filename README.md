Practical MR session using R
================

This repository contains material of the Mendelian Randomization (MR)
session within the [Online Summer School 2021: Genetic Epidemiology of
Kidney Function and Chronic Kidney
Disease](https://geneticepisummerschool.eurac.edu/).

The session will cover practical aspects regarding the implementation of
Mendelian Randomization (MR) analyses. The course is focused on
two-sample MR applications using summary-level data from genome-wide
association studies (GWAS). References on methodology and software for
one-sample MR methods will be also provided.

To download all the session material, click on `Code -> Download ZIP`.
Then save the *.zip* file locally and extract all the material of the
course.

The main folders are organized as follows:

-   **data** folder contains the raw GWAS summary-level data that are
    used in the case study. The files are in *.txt* format.

-   **scripts** folder contains the R script that is used for the case
    study during the session

-   **slides** folder contains all the files that are used to create the
    slides for the presentation

-   **vignettes** folder contains all the files that are used to create
    the tutorial document.

The tutorial document contains a summary of the case study used to show
MR analysis workflow. It can be found [here](vignettes/tutorial.md).

A PDF version of the presentation can be found
[here](slides/session_slides.pdf). Please use Google Chrome browser for
a correct visualization of the PDF file.

All the other folders and files should not be moved or deleted for
reproducing the material of the session.

**Note**: the material will be continuously updated until the day of the
session, i.e. July 14th, 2021. After the session, updates will be
performed based on requests or suggestions from the users. If you need
to make a request, please open an issue
[here](https://github.com/EuracBiomedicalResearch/trainckdis.mr/issues).

## Program

1.  **Data preparation**:

-   Collecting summary-level genome-wide association studies (GWAS) data
    from online repositories
-   Data harmonization with `TwosampleMR` R package

2.  **Running MR analysis**:

-   Estimating causal effects using `MendelianRandomization` and
    `TwoSampleMR` R packages
-   Displaying and interpreting results

3.  **Sensitivity analyses**:

-   Evaluating if MR assumptions are supported by the data
-   Estimating causal effects with robust MR methods

## Prerequisites

The session will be taught in `R` (see
[here](https://www.r-project.org/) for more information). Previous
knowledge of `R` programming can be helpful although not mandatory.

## Software setup

-   Please dowload and install the latest `R` version from
    [here](https://www.r-project.org/)

-   Although it is not mandatory, it is recommended to download and
    install the latest `RStudio` version from
    [here](https://www.rstudio.com/products/rstudio/download/) (free
    version) since `RStudio` will be used during the session

-   Please download and install the latest version of the following `R`
    packages from CRAN:

``` r
install.packages(
  c(
    "remotes",
    "devtools",
    "MendelianRandomization",
    "ggplot2"
  ),
  dependencies = TRUE
)
```

-   Please download and install the latest version of `TwoSampleMR` from
    Github repository:

``` r
remotes::install_github("MRCIEU/TwoSampleMR")
```

## Software introduction

Introductory material on `R`, `Rstudio` and `R` packages to implement
one-sample and two-sample MR is provided below. We suggest the students
to check the material before the beginning of the session to start
practicing with `R` programming and the packages that will be used
throughout the session

### Setting-up R

-   See [Chapter 1 of
    ModernDive](https://moderndive.netlify.app/1-getting-started.html)
    for beginner-friendly guide to install `R` and `Rstudio`

-   See [BasicBasics 1 and
    2](https://rladiessydney.org/courses/ryouwithme/01-basicbasics-0/)
    for a gentle introduction to basic commands in `R`

### Two-sample MR packages

-   [Introduction to the `TwoSampleMR`
    package](https://mrcieu.github.io/TwoSampleMR/articles/index.html)

-   [Introduction to the `MendelianRandomization`
    package](https://mrcieu.github.io/TwoSampleMR/articles/index.html)

### One-sample MR packages

-   `ivreg` function in [`AER` R
    package](https://cran.r-project.org/web/packages/AER/index.html) or
    `tsls` function in [`sem` R
    package](https://cran.r-project.org/web/packages/sem/index.html) to
    implement two-stage least squares (2SLS) method

-   [`ivtools` R
    package](https://cran.r-project.org/web/packages/ivtools/index.html)
    to implement several one-sample MR methods

## Contributions

| Name                                                           | Affiliation                                                                                                        |
|:---------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------|
| [Daniele Bottigliengo](https://github.com/danielebottigliengo) | [EURAC Research, Institute for Biomedicine](https://www.eurac.edu/en/institutes-centers/institute-for-biomedicine) |

## Reading list

Before starting the session, the reading of the following papers is
suggested:

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bowden_mendelian_2015" class="csl-entry">

Bowden, Jack, George Davey Smith, and Stephen Burgess. 2015. “Mendelian
Randomization with Invalid Instruments: Effect Estimation and Bias
Detection Through Egger Regression.” *International Journal of
Epidemiology* 44 (2): 512–25. <https://doi.org/10.1093/ije/dyv080>.

</div>

<div id="ref-burgess_guidelines_2019" class="csl-entry">

Burgess, Stephen, George Davey Smith, Neil M. Davies, Frank Dudbridge,
Dipender Gill, M. Maria Glymour, Fernando P. Hartwig, et al. 2019.
“Guidelines for Performing Mendelian Randomization Investigations.”
*Wellcome Open Research* 4: 186.
<https://doi.org/10.12688/wellcomeopenres.15555.2>.

</div>

<div id="ref-burgess_review_2017" class="csl-entry">

Burgess, Stephen, Dylan S. Small, and Simon G. Thompson. 2017. “A Review
of Instrumental Variable Estimators for Mendelian Randomization.”
*Statistical Methods in Medical Research* 26 (5): 2333–55.
<https://doi.org/10.1177/0962280215597579>.

</div>

<div id="ref-davey_smith_mendelian_2003" class="csl-entry">

Davey Smith, George, and Shah Ebrahim. 2003. “‘Mendelian Randomization’:
Can Genetic Epidemiology Contribute to Understanding Environmental
Determinants of Disease?\*.” *International Journal of Epidemiology* 32
(1): 1–22. <https://doi.org/10.1093/ije/dyg070>.

</div>

<div id="ref-davey_smith_mendelian_2014" class="csl-entry">

Davey Smith, George, and Gibran Hemani. 2014. “Mendelian Randomization:
Genetic Anchors for Causal Inference in Epidemiological Studies.” *Human
Molecular Genetics* 23 (R1): R89–98.
<https://doi.org/10.1093/hmg/ddu328>.

</div>

<div id="ref-davies_reading_2018" class="csl-entry">

Davies, Neil M., Michael V. Holmes, and George Davey Smith. 2018.
“Reading Mendelian Randomisation Studies: A Guide, Glossary, and
Checklist for Clinicians.” *BMJ (Clinical Research Ed.)* 362 (July):
k601. <https://doi.org/10.1136/bmj.k601>.

</div>

<div id="ref-hartwig_two-sample_2016" class="csl-entry">

Hartwig, Fernando Pires, Neil Martin Davies, Gibran Hemani, and George
Davey Smith. 2016. “Two-Sample Mendelian Randomization: Avoiding the
Downsides of a Powerful, Widely Applicable but Potentially Fallible
Technique.” *International Journal of Epidemiology* 45 (6): 1717–26.
<https://doi.org/10.1093/ije/dyx028>.

</div>

<div id="ref-noauthor_mendelian_nodate" class="csl-entry">

“Mendelian Randomization Dictionary.” n.d.
<https://mr-dictionary.mrcieu.ac.uk/>.

</div>

<div id="ref-slob_comparison_2020" class="csl-entry">

Slob, Eric A. W., and Stephen Burgess. 2020. “A Comparison of Robust
Mendelian Randomization Methods Using Summary Data.” *Genetic
Epidemiology* 44 (4): 313–29.
https://doi.org/<https://doi.org/10.1002/gepi.22295>.

</div>

</div>
