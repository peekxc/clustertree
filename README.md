# clustertree
`clustertree` is a standalone, scalable, and extensible [R](https://www.r-project.org/package) package for estimating the empirical cluster tree, a hierarchical representation of *high-density clusters*, defined (recursively) as the connected subsets of: 

$$\{ x : f(x) \geq \lambda \} \text{ for some } \lambda > 0$$

## Installation 
The package currently only exists on github. 

1. It can be download and built manually, or, 
2. Alternatively installed with help from the [devtools](https://github.com/hadley/devtools) package, i.e. 

```R
install.packages("devtools")
devtools::install_github("peekxc/clustertree")
```
### Development note 
The package as it stands is in very early stage of development, and should be regarded as such. A release candidate for [CRAN](https://cran.r-project.org/) is planned for approximately sometime around September 5$^{\mathrm{th}}$ 2017. 

## Usage 
Usage section to be continued.... 

### Estimators implemented
There are multiple algorithms for approximating the empirical cluster tree. Below is a (growing) list of estimators provided by this package:

1. The **Robust Single Linkage (RSL)** algorithm, which can be found in: 

> Chaudhuri, Kamalika, and Sanjoy Dasgupta. "Rates of convergence for the cluster tree." Advances in Neural Information Processing Systems. 2010.


## Additional References
The cluster tree theory itself has a long history. Specifically, see section  **11.13** of: 
> Hartigan, John A., and J. A. Hartigan. Clustering algorithms. Vol. 209. New York: Wiley, 1975.
for an overview of what Hartigan refers to as the *density-contour tree,* and brief overview of the background and associated concepts of formal high-density clustering. 

Established notions of the consistency of estimators of the cluster tree can be found in: 
> Hartigan, John A. "Consistency of single linkage for high-density clusters." Journal of the American Statistical Association 76.374 (1981): 388-394.

## License 
The package and its associated software 

## Acknowledgements 
This package is being developed as part of the [Google Summer of Code 2017](https://summerofcode.withgoogle.com/dashboard/project/5111030546956288/overview/) under the [R Project](https://www.r-project.org/). 