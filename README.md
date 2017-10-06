# DeformableRegistration.jl

[![Build Status](https://travis-ci.com/johanneslotz/DeformableRegistration.jl.svg?token=rpbV4sPrj6BdxqGJ84cq&branch=master)](https://travis-ci.com/johanneslotz/DeformableRegistration.jl)


This is a first step towards an image registration package written in Julia. Much of this code is still in a very early stage, lacks a lot of features and even has some bugs. Basic registration is possible though and pull-requests adding features or fixing bugs are explicitly welcome. :)

Some of the code is inspired by the [FAIR](http://www.mic.uni-luebeck.de/de/people/jan-modersitzki/software/fair.html)  software package by Jan Modersitzki.

## Getting started

```
Pkg.clone("https://github.com/lruthotto/KrylovMethods.jl")
Pkg.clone("https://github.com/johanneslotz/DeformableRegistration.jl")
```

To get started, compute a first registration based on the example code in

```
test/{Parametric,Nonparametric}Registration.jl
```

## Features

- Image Registration in 2D
- Disance Measures: SSD, NGF
- Regularizer: Diffusive, Elastic
- Optimization: Gauss-Newton/Armijo

## Module organization

Besides the mentioned packages the **DeformableRegistration** module itself consists of a few submodules:

* **.Distance**: Sum of Squared Differences (SSD), masked SSD, Normalized Gradient Fields (NGF)

* **.Regularizer**: Functions to regularize the deformation field. So far an elastic and a diffusive regularizer have been implemented (both matrix-based).

* **.Optimization**: Gauss-Newton optimization with Armijo line search is provided by this submodule. A wrapper for the NLOpt package is in development.

* **.ImageProcessing**: Different functions to create and load images and handle their properties.

* **.Interpolation**: Linear interpolation, a specialized implementation and a generic based on the Grid package.

* **.Transformation**: Convenience functions to create cell-centered and staggered grids, to convert between them and to transform them.

* **.Visualization**
