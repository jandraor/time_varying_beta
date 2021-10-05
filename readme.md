
# Andrade & Duggan

This repository contains code for the paper:

[Jair Andrade](https://www.linkedin.com/in/jandraor/) and [Jim
Duggan](https://ie.linkedin.com/in/jduggan). *Inferring the effective
reproductive number from deterministic and semi-deterministic
compartmental models using incidence and mobility data*

The analysis in this study can be reproduced by executing the files:

-   **S1.rmd**
-   **S2.rmd**
-   **S3.rmd**
-   **S4.rmd**

## Abstract

The effective reproduction number (ℜ<sub>*t*</sub>) is a theoretical
indicator of the course of an infectious disease that allows
policymakers to evaluate whether current or previous control efforts
have been successful or whether additional interventions are necessary.
This metric, though, cannot be directly observed and must be estimated
from available data. One of the approaches to obtain such estimates is
fitting compartmental models to incidence data. We can envision these
dynamic models as the ensemble of structures that describe the disease’s
natural history and individuals’ behavioural patterns. In the context of
the response to the COVID-19 pandemic, the assumption of a constant
transmission rate is rendered unrealistic, and it is critical to
identify a mathematical formulation that accounts for the changes in
contact patterns. In this work, we leverage existing approaches to
propose three complementary formulations that yield similar estimates
for ℜ<sub>*t*</sub> during Ireland’s first COVID-19 wave. We describe
these Data Generating Processes (DGP) in terms of State-Space models.
Two (DGP1 and DGP2) correspond to stochastic process models whose
transmission rate is modelled as Brownian motion processes (Geometric
and Cox-Ingersoll-Ross). These DGPs share a measurement model that
accounts for incidence and transmission rates, where mobility data is
assumed as a proxy of the transmission rate. We perform inference on
these structures using iterated filtering and the particle filter. The
final DGP (DGP3) is built from a pool of deterministic models that
describe the transmission rate as information delays. We calibrate the
pool of models to incidence reports using Hamiltonian Monte Carlo. By
following this complementary approach, we assess the tradeoffs involved
in each formulation and reflect on the interpretation of the estimated
uncertainty. We anticipate this work will help to evaluate the
implications of choosing a particular formulation for the time-varying
transmission rate
