# DeformableRegistration.jl

[![Build Status](https://travis-ci.org/johanneslotz/DeformableRegistration.jl.svg?branch=master)](https://travis-ci.org/johanneslotz/DeformableRegistration.jl)


This is a first step towards an image registration package written in Julia. The current goal is to build a platform to experiment with registration and explore algorithms. Basic registration is possible though and pull-requests adding features or fixing bugs are explicitly welcome. :)

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
- Regularizer: Diffusive, Curvature
- Optimization: Gauss-Newton/Armijo
- some unit and integration tests

## Module organization

**DeformableRegistration** is organized in modules:

* **.Distance**: Sum of Squared Differences (SSD), masked SSD, Normalized Gradient Fields (NGF)

* **.Regularizer**: Functions to regularize the deformation field. So far, curvature and a diffusive regularization have been implemented (both matrix-based).

* **.Optimization**: Gauss-Newton optimization with Armijo line search.

* **.ImageProcessing**: Different functions to create and load images and handle their properties.

* **.Interpolation**: Linear interpolation, a specialized implementation and a generic based on the Interpolations package.

* **.Transformation**: Convenience functions to create and transform cell-centered grids.

* **.Visualization**
