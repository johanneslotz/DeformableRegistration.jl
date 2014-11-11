# Image Registration

[![Build Status](https://magnum.travis-ci.com/johanneslotz/ImageRegistration.jl.svg?token=rpbV4sPrj6BdxqGJ84cq&branch=master)](https://magnum.travis-ci.com/johanneslotz/ImageRegistration.jl)

A first step towards an image registration package written in Julia. 

Some of the code is inspired by the [FAIR](http://www.mic.uni-luebeck.de/de/people/jan-modersitzki/software/fair.html)  software package by Jan Modersitzki. ImageRegistration.jl relies on the Julia packages
- [Images](https://github.com/timholy/Images.jl) for image loading and management of image properties
- [Grid](https://github.com/timholy/Grid.jl) for the interpolation
- [KrylovMethods](https://github.com/lruthotto/KrylovMethods.jl) to solve the linear system in the Gauss-Newton optimization
- [PyPlot](https://github.com/stevengj/PyPlot.jl) 


Besides the mentioned packages the **ImageRegistration** module itself consists of a few submodules:

* **.Distance**: Sum of Squared Differences (SSD), masked SSD, Normalized Gradient Fields (NGF)

* **.Regularizer**: Functions to regularize the deformation field. So far an elastic and a diffusive regularizer have been implemented (both matrix-based).

* **.Optimization**: Gauss-Newton optimization with Armijo line search is provided by this submodule. A wrapper for the NLOpt package is in development.

* **.ImageProcessing**: Different functions to create and load images and handle their properties.

* **.Transformation**: Convenience functions to create cell-centered and staggered grids, to convert between them and to transform them. 

* **.Visualization**

There are still some open issues in the bug tracker. Code contributions of additional features, bugfixes and optimizations are very welcome.