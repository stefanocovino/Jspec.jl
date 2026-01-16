# JSpecAstro

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stefanocovino.github.io/JSpecAstro.jl/stable/)
[![Build Status](https://github.com/stefanocovino/JSpecAstro.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stefanocovino/JSpecAstro.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package to allow to read and analyse spectra obtained from multi-channel instruments (e.g., Swift-XRT) with data from any other source (e.g. optical/NIR observations). Altough several features are in common, no attempt to mimic the full funtionalities offered by [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) was tried. In addition, at present, only fits with a Gaussian likelihood are implemented. Fits in a full Poissonian regime will possibly be included in a future version (or never...). 


## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/stefanocovino/JSpecAstro.jl.git")
```

or just

```julia
using Pkg
Pkg.add("JSpecAstro")
```

will install this package, with the latter enabled once (if ever) the package is registered.


### Instruments

JSpecAstro currently handles data from [Swift-BAT](https://science.nasa.gov/mission/swift/), [Swift-XRT](https://science.nasa.gov/mission/swift/), [SVOM-MXT](https://www.svom.eu/en/the-svom-mission/), [NuSTAR-FPM](https://heasarc.gsfc.nasa.gov/docs/nustar/), [XMM-EMOS](https://www.cosmos.esa.int/web/xmm-newton) and [XMM-EPN](https://www.cosmos.esa.int/web/xmm-newton).


### Documentation and examples

Documentation for the package can be found [here](https://stefanocovino.github.io/JSpecAstro.jl/stable/).


### Similar tools

If you are interested in similar (or better) capabilities you may also check the "[SpectralFitting.jl](https://github.com/fjebaker/SpectralFitting.jl?tab=readme-ov-file)" and "[LibXSPEC_jll.jl](https://github.com/astro-group-bristol/LibXSPEC_jll.jl)".


## Getting Started

The purpose of the package is to provide tools to mode data from multi-channel instruments togeter, if needed, with data from any other surce. The package compute the needed response matrices that can then used for creating models, carry out fits, etc.

No attempt has been tried, on purpose, to mimic the simplified XSPEC syntax to create models, etc. Therefore models, etc. will be coded according to plain Julia syntax.

The present version of the package uses Gaussian statistics and data has to be, if needed, adequately rebinned. A future version "might" offer analyses based on Poissonian statistics.

The package has been used for real scientific analyses as in, e.g., [Brivio et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025A%26A...695A.239B/abstract).



### A toy session


The instruments currently covered can be obtained with:

```julia
GetKnownInstruments()
```

Ad the first step is to create a new dataset. For instance, assuming we want to model 'Swift-XRT' data and data from an optical telescope, we might write:

```julia
XRTdt = CreateDataSet("XRTdata","Swift-XRT")
Optdt = CreateDataSet("Optdata","Other")
```

'XRTdt' and 'Optdt' are dictionaris that are going to include all the needed information.

Assuming we have the following 'ex.rmf', 'ex.arf', 'exsrc.pi' and 'exbck.pi' XRT files we can import them as follows:

```julia
ImportData(XRTdt, rmffile="ex.rmf", arffile="ex.arf", srcfile="exsrc.pi", bckfile="exbck.pi")
```

Optical data, but also data from any other source where a non-diagonal response matrix is not needed, should be converted to energy, in keV, and photon flux density in photons ``cm^{-2}~s{^-1}~keV{^-1}``. Alternatively, data can be represented by flux but in such a case the bandwidth, again in ``keV``, must be provided too.

```julia
ImportOtherData(Optdt, energy=[1.,2.,3.,4], phflux=[0.1,0.2,0.3,0.4], ephflux=[0.01,0.02,0.03,0.04])
```

It is possible to visualize the imported data with, e.g.:
```julia
PlotRaw(XRTdt)
```

or:
```julia
PlotRaw(Optdt,ylbl=L"Photons s$^{-1}$ cm$^{-2}$ keV$^{-1}$")
```

Often, for multi-channel instruments, channels can (or need to be) ignored. This can be achieved with, e.g. (please be aware that the first channels is the 0-channel, at variance with the julia convention):

```julia
IgnoreChannels(XRTdt,[0:30,1000:2047])
```


And, data must often be rebinned to make the analysis based on a Gaussian likelihood meaningful. This can be achieved easily choosing the minimum S/N per bin with, e.g.:

```julia
RebinData(XRTdt,minSN=7)
```

In case no rebinning is needed the step should be executed anyway with 'minSN=0'.

Once a rebinning schema has been defined, the ancillary data should also be properly rebinned:

```julia
RebinAncillaryData(XRTdt)
```

Now, it is also possible to visualize the rebinned data with:

```julia
PlotRebinned(XRTdt)
```

And, finally, a response matrix properly rebinned following the rebin schema identified above can be generated:

```julia
GenResponseMatrix(XRTdt)
```

At this point, we need to define a model for our data. This can be expressed by regular `Julia` syntax, although the function 
arguments should be set up to allow JSpecAstro to convolve the function results with the response matrices and be able to compare the
 model prediction with the observation.
 
For instance, with XRT and optical data, it might be something as a simple power-law with local optical and X-ray extinction and absorption:

```julia
OptExt(E,EBV) = Extinction(12.4 ./ E,EBV,gal="SMC",Rv=FFGals["SMC"],z=0.)

function MyModel(pars,E)
    N,β,NH,EBV = pars
    return PL(E,N,β) .* ifelse.(E .< 0.01, OptExt(E,EBV), XAbs(E,NH=NH,z=0.))
end
```

where we have used functions defined in the [`FittingFunction.jl`](https://github.com/stefanocovino/FittingFunction.jl.git) package. It is important to check that the defined function can be used for all the energy ranges covered by the imported datasets.

Vectors with the observatins, uncertainties and input energies for all the imported datased can be obtained with:


```julia
obs,eobs,engy = GenFullObsData([Optdt,XRTdt])
```

Finally, having defined a theoretical model, we can convolve it with the various, possiby rebinned, responce matrices of the imported datasets to obtain model predictions to be compared with the observed data. This can be obtained by means of the `JSpecAstroFunc` function as:

```julia
modpreds = JSpecAstroFunc(pars,[Optdt,XRTdt],MyModel)
```

`JSPECFunc` with the right parameters can be used for any optimization problem or for any Baysian analysis involving sampling, etc.


For instance, a [`Turing`](https://turinglang.org/) model might look like:

```julia
@model function SEDModel(f,ef)
        lN ~ Uniform(-10., 10.)
        β ~ Normal(-2.,2)
        lNH ~ Uniform(log(1e17), log(1e22))
        EBV ~ Uniform(0,2)
        #
        ris = JSpecAstroFunc([exp(lN),β,exp(lNH),EBV],[Optdt,XRTdt],MyModel)
        #
        f ~ MvNormal(ris1,ef)
    end

model = SEDModel(obs,eobs);
```

And, having defined a model you can get the maximum likelihod, carry out a sampling, etc.
