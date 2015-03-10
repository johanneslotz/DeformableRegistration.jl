# Image Registration

[![wercker status](https://app.wercker.com/status/6e79cb4e56aefa2c386d4f14cd2e3a0f/m "wercker status")](https://app.wercker.com/project/bykey/6e79cb4e56aefa2c386d4f14cd2e3a0f)

A first step towards an image registration package written in Julia.

Some of the code is inspired by the [FAIR](http://www.mic.uni-luebeck.de/de/people/jan-modersitzki/software/fair.html)  software package by Jan Modersitzki.

## Getting started

```
Pkg.clone("https://github.com/lruthotto/KrylovMethods.jl")
Pkg.clone("https://github.com/johanneslotz/ImageRegistration.jl")
```

To get started, compute a first registration based on the example code in

```
tests/{Parametric,Nonparametric}Registration.jl
```

## Features

- Disance Measures: SSD, NGF (only nonparametric so far)
- Regularizer: Diffusive, Elastic
- Optimization: Gauss-Newton/Armijo

## Module organization

Besides the mentioned packages the **ImageRegistration** module itself consists of a few submodules:

* **.Distance**: Sum of Squared Differences (SSD), masked SSD, Normalized Gradient Fields (NGF)

* **.Regularizer**: Functions to regularize the deformation field. So far an elastic and a diffusive regularizer have been implemented (both matrix-based).

* **.Optimization**: Gauss-Newton optimization with Armijo line search is provided by this submodule. A wrapper for the NLOpt package is in development.

* **.ImageProcessing**: Different functions to create and load images and handle their properties.

* **.Transformation**: Convenience functions to create cell-centered and staggered grids, to convert between them and to transform them.

* **.Visualization**

There are still some open issues in the bug tracker. Code contributions of additional features, bugfixes and optimizations are very welcome.
