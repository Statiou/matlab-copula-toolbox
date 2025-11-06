---
title: "Copula Toolbox — All-in-One: clustering & copula modeling for MATLAB"
tags:
  - copulas
  - dependence modeling
  - clustering
  - financial econometrics
  - MATLAB
authors:
  - name: D. A. Statiou
    orcid: 0009-0005-9581-2064
    affiliation: 1
  - name: P. Hatzopoulos
    
    affiliation: 2
affiliations:
  - name: Department of Economics, University of ______
    index: 1
  - name: Department of Statistics, Institute of ______
    index: 2
date: 2025-01-01
bibliography: references.bib
---

# Summary

**Copula Toolbox — All-in-One** is a MATLAB package and GUI for end-to-end dependence modeling with copulas. It provides: (i) fitting of Gaussian, Student-t, and 2-D Archimedean families (Clayton, Frank, Gumbel); (ii) automatic marginal modeling and Probability Integral Transform (PIT) diagnostics; (iii) clustering on raw or copula (U) space with *k*-means, Gaussian mixtures, hierarchical clustering, and DBSCAN; (iv) per-cluster copula selection via information criteria; (v) diagnostic plots (rank heatmaps, PIT, empirical copula, goodness-of-fit via Rosenblatt transform); (vi) parametric bootstrap and *K*-fold cross-validation; (vii) dynamic Kendall’s τ forecasting using a Generalized Autoregressive Score (GAS) model; and (viii) systemic-risk measures CoVaR/CoES.

The toolbox targets researchers and practitioners who need a reproducible, scriptable, and GUI-assisted workflow to explore high-dimensional dependence structures, learn regime-specific copulas, and produce publishable diagnostics in finance, risk management, and applied statistics.

# Statement of need

Copula modeling is widely used to separate marginal behavior from dependence and to capture non-linear, tail, and asymmetric relationships beyond linear correlation [@Sklar1959; @Nelsen2006; @Joe2014]. Despite mature theory, researchers face practical obstacles: preparing marginals and PITs, switching among copula families, validating fits, and—critically—handling heterogeneity where a single global copula is inadequate. Existing MATLAB examples are fragmented, and many open-source tools emphasize either estimation **or** clustering, not both in one workflow.

**Copula Toolbox — All-in-One** fills this gap by coupling standard copula estimation with clustering that can operate in raw space or on copula ranks (*U*-space), followed by per-cluster copula selection. It adds routine but often tedious tasks—diagnostics, bootstrap/CV, and GAS dynamics—in a single, reproducible interface that works out-of-the-box for MATLAB users.

# State of the field

The package builds on standard copula theory [@Sklar1959; @Nelsen2006; @Joe2014] and practice in financial econometrics [e.g., @Patton2012]. Goodness-of-fit leverages the Rosenblatt transform [@Rosenblatt1952]. Clustering options include *k*-means, Gaussian mixture modeling with EM [@Dempster1977], hierarchical clustering, and DBSCAN for density-based structure [@Ester1996]. Validity diagnostics report silhouette [@Rousseeuw1987] and gap statistics [@Tibshirani2001]. Dynamic τ modeling follows the GAS framework [@Creal2013]. For independence testing we expose HSIC [@Gretton2005]. Model selection uses AIC/BIC [@Akaike1974; @Schwarz1978].

While Python and R ecosystems offer rich copula packages, MATLAB users often lack an integrated solution that combines clustering, per-cluster copulas, PIT diagnostics, and forecasting within a unified GUI and scriptable API. This toolbox addresses that niche.

# Functionality

**Core features**

- **Marginals & PIT:** Normal, Lognormal, Exponential, Gamma, *t*-location-scale, or empirical ranks; PIT histograms/diagnostics.
- **Copulas:** Gaussian, Student-*t* (with ν), and 2-D Archimedeans (Clayton, Frank, Gumbel). Empirical copula contours; PDF/CDF on *U*.
- **Model selection:** AIC/BIC across families; *K*-fold CV of held-out log-likelihood.
- **Diagnostics:** Rank-based heatmaps (Kendall, Spearman), Kendall plots, Rosenblatt GOF (2-D), tail-dependence estimates.
- **Clustering:** *k*-means, GMM, hierarchical, DBSCAN; operates on raw *X* or *U*-space. Automated *K* heuristics (elbow/silhouette/gap). Per-cluster copula fitting with BIC ranking and visual diagnostics.
- **Dynamics:** Rolling Kendall’s τ and one-step-ahead τ forecasts via a GAS(1,1) score update with optional exogenous features.
- **Risk:** CoVaR / CoES (2-D) under the fitted copula and chosen marginals.
- **Uncertainty:** Parametric bootstrap for R (or θ) and ν; simulation from fitted copulas.
- **Reproducibility:** Export/Save sessions; quick text/PNG reports.

**User interface & API**

A single `fitcopulatoolbox_allinone()` function opens the GUI. All internal routines are also accessible as functions (e.g., `cluster_copula_module_inline`, `forecastTauGAS_core`) for scripted pipelines and CI tests. Examples and a small CSV dataset are provided in `examples/` for immediate replication.

# Quality control

We include smoke tests in `tests/` that (i) load the example dataset, (ii) transform marginals, (iii) fit at least one copula family, (iv) run a minimal clustering + per-cluster fitting pass, and (v) call the GAS τ forecaster on a short rolling window, asserting finite outputs. Continuous integration via GitHub Actions executes `runtests` on every commit. Randomness is controlled via fixed seeds where feasible. Plots are generated programmatically from the returned structures to ensure consistency between GUI and API.

# Example usage

```matlab
% Load data (matrix or table), pick columns
T = readtable('examples/stocks_small.csv');
X = T{:, {'AAPL','MSFT','GOOGL'}};

% Cluster in U-space and fit per-cluster copulas
CC = cluster_copula_module_inline(X, ...
      'Marginals','Normal','Algorithm','kmeans', ...
      'UseUspace',true,'Kmax',8,'Verbose',true);

% Inspect best family per cluster
for k = 1:CC.K
    disp(CC.clusters(k).fit);           % family, R/theta, nu
    disp(CC.clusters(k).criteria.BIC);  % information criterion
end

% Forecast dynamic Kendall's tau for a pair
out = forecastTauGAS_core(X(:,1:2),'Normal',60,'t',0.10,0.20);
plot(out.tau_real); hold on; plot(out.tau_pred,'--');
