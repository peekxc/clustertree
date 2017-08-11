# clustertree
[![CRAN](http://www.r-pkg.org/badges/version/clustertree)](#)
[![Travis-CI Build Status](https://travis-ci.org/peekxc/clustertree.svg?branch=master)](https://travis-ci.org/peekxc/clustertree)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/peekxc/clustertree?branch=master&svg=true)](https://ci.appveyor.com/project/peekxc/clustertree)

`clustertree` is a fast and extensible [R package](https://www.r-project.org/package) for estimating the empirical cluster tree, a hierarchical representation of *high-density clusters*, defined (recursively) as the connected subsets of: 
<div style = "text-align:center" align="center"> <img src="clustertree.png" width = "278"/> </div>
From a high-level perspective, the cluster tree provides a highly interpretable, multi-resolution, and statistically sound summary of the underlying density of a [finite] sample. The package includes both tools and complete implementations of the following (growing) list of estimators: 

--
1. The **Robust Single Linkage (RSL)** algorithm from: 
	
	> Chaudhuri, Kamalika, and Sanjoy Dasgupta. "Rates of convergence for the cluster tree." Advances in Neural Information Processing Systems. 2010.
	
--

For more information regarding the utility of this package and of the cluster tree itself, see the **Usage** and **Additional References** sections, respectively.

<!--The applications are many—density-based clustering is one such application. The benefits of density-based clustering are numerous, including the ability to capture clusters of arbitrary or non-convex shapes, they do not require *a priori* knowledge concerning number of clusters to find, and they are more often than not robust to varying amounts noise. Akin to some density-based clustering approaches, the cluster tree shares another benefit relatively absent in other clustering approaches: the definition of what constitutes a cluster and its overall object of inference, the hierarchical tree of high-density clusters, is clearly and formally stated. -->

## Installation 
The package currently only exists on github. The installation options are as follows: 

1. Installed directly from the repo with help from the [devtools](https://github.com/hadley/devtools) package, i.e. 

	```R
	install.packages("devtools")
	devtools::install_github("peekxc/clustertree")
	```
2. Download the most recent successful build from [AppVeyor](https://ci.appveyor.com/api/projects/peekxc/clustertree/artifacts/bin/debug.zip?branch=master
) 
#### Development note 

The package is in very early stage of development, and should be regarded as such. A release candidate for [CRAN](https://cran.r-project.org/) is planned for approximately sometime around 09-5-2017. Also note that much of the preliminary code is being developed within **separate branches**,
which are very unstable and highly experimental. 

## Usage 
Usage section to be continued.... 

## Additional References
The cluster tree theory itself has a long history. For a brief overview of the definition, see section  **11.13** of: 
> Hartigan, John A., and J. A. Hartigan. Clustering algorithms. Vol. 209. New York: Wiley, 1975.
for an overview of what Hartigan refers to as the *density-contour tree,* and brief overview of the background and associated concepts of formal high-density clustering. 

Established notions of consistency of the above class of estimators can be found in: 
> Hartigan, John A. "Consistency of single linkage for high-density clusters." Journal of the American Statistical Association 76.374 (1981): 388-394.


## Acknowledgements 
This package is part of the [Google Summer of Code 2017](https://summerofcode.withgoogle.com/dashboard/project/5111030546956288/overview/) initiative, under the [R Project](https://www.r-project.org/). 