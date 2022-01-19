
# ELCIC: Empirical Likelihood-based Consistent Information Criterion
[![codecov](https://codecov.io/gh/chencxxy28/ELCIC/branch/master/graph/badge.svg?token=MP1P4H0OHP)](https://nam11.safelinks.protection.outlook.com/?url=https%3A%2F%2Fapp.codecov.io%2Fgh%2Fchencxxy28%2FELCIC&amp;data=04%7C01%7CChixiang.Chen%40som.umaryland.edu%7Caac18e7607274a30deb208d9db27dc2a%7C717009a620de461a88940312a395cac9%7C0%7C0%7C637781786524712369%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000&amp;sdata=5Qr4dvXJsoD9lxRjMSzHpoTDiZHGet61d%2BuYIJXvKVI%3D&amp;reserved=0)
[![R-CMD-check](https://github.com/chencxxy28/ELCIC/workflows/R-CMD-check/badge.svg)](https://github.com/chencxxy28/ELCIC/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/ELCIC)](https://cran.r-project.org/package=ELCIC)

Conventional likelihood-based information criteria for model selection
rely on the distribution assumption of data. However, for complex data
that are increasingly available in many scientific fields, the
specification of their underlying distribution turns out to be
challenging, and the existing criteria may be limited and are not
general enough to handle a variety of model selection problems. We
proposed a robust and consistent model selection criterion, named as
ELCIC, based upon the empirical likelihood function which is
data-driven. In particular, this framework adopts plug-in estimators
that can be achieved by solving external estimating equations, not
limited to the empirical likelihood, which avoids potential
computational convergence issues and allows versatile applications, such
as generalized linear models, generalized estimating equations,
penalized regressions, and so on. The formulation of our proposed
criterion is initially derived from the asymptotic expansion of the
marginal likelihood under the variable selection framework, but more
importantly, the consistent model selection property is established
under a general context.

ELCIC offers a robust model assessment and can be applied to address
more complicated situations where existing methods fail to work.

# How to cite ELCIC

Please cite the following publication: Chixiang Chen, Ming Wang,
Rongling Wu, and Runze, Li, A Robust Consistent Information Criterion
for Model Selection based on Empirical Likelihood
<https://arxiv.org/pdf/2006.13281.pdf>

# Installation

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("chencxxy28/ELCIC")
```

# Vignettes

Please visit [Tutorial](https://chencxxy28.github.io/ELCIC/articles/ELCIC.html)
