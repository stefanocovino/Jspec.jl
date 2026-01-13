module Jspec

using CairoMakie
using DataFrames
using FITSIO
using LaTeXStrings
using LinearAlgebra



export Angstrom2KeV
export CreateDataSet
export GenFullObsData
export GenResponseMatrix
export GetKnownInstruments
export IgnoreChannels
export ImportData
export ImportOtherData
export JspecFunc
export Jy2PhFlux
export KeV2Angstrom
export KeV2Channel
export PlotRaw
export PlotRebinned
export RebinAncillaryData
export RebinData





"""
    Angstrom2KeV(wave)

Convert wavelengths in Angstrom (``10^{-8}~m``) to ``keV``.

#Arguments

- `wave` the input wavelengths.


# Examples
```jldoctest

Angstrom2KeV(5000:5005)

# output

6-element Vector{Float64}:
 0.00248
 0.002479504099180164
 0.0024790083966413435
 0.0024785128922646415
 0.0024780175859312552
 0.0024775224775224775
```
"""
function Angstrom2KeV(wave)
    return 12.4 ./ wave
end




KnownInstruments = ["Swift-XRT", "Swift-BAT", "SVOM-MXT", "Other", "NuSTAR-FPM", "XMM-EMOS", "XMM-EPN"]




"""
    CreateDataSet(Name::String, Instrument::String; verbose=true)::Dict

Create a JspecDataSent entry.

#Arguments

- `Name` is the arbitrary name of the dataset.
- `Instrument` is one of supperted instrument by the package.
- `verbose` enables warning messages.


# Examples
```jldoctest

newdataset = CreateDataSet("XRTTest","Swift-XRT")

# output

Dict{Any, Any} with 3 entries:
  "Name"       => "XRTTest"
  "Instrument" => "Swift-XRT"
  "Created"    => true
```
"""
function CreateDataSet(Name::String, Instrument::String; verbose=true)::Dict
    nd = Dict()
    if uppercase(Instrument) in [uppercase(i) for i in Jspec.KnownInstruments]
        nd["Name"]=Name
        nd["Instrument"]=Instrument
        nd["Created"] = true
    else
        if verbose
            println("Warning! Unknown instrument: "*Instrument)
        end
        nd["Name"]=Name
        nd["Instrument"]=""
        nd["Created"] = false
    end
    return nd
end







"""
    FindRebinSchema(x,ey;minSN=5)::AbstractVector{Real}

Compute the rebin schema to guarantee that the S/N is at least `minSN` in each bin (or channel).

# Arguments

- `x` input array.
- `ex` uncertainties.

# Examples
```jldoctest

x = [1.,2.,3.,4.,]
ex = [0.1,0.5,0.6,0.05]

Jspec.FindRebinSchema(x,ex)

# output

3-element Vector{Real}:
 1
 3
 4

```
"""
function FindRebinSchema(x::AbstractVector{Float64},ex::AbstractVector{Float64};minSN=5)::AbstractVector{Real}
    sbin = []
    i = 1
    while i <= length(x)
        for l in i:length(x)
            c = sum(x[i:l])
            b = sqrt(sum(ex[i:l].^2))
            if abs(c)/b >= minSN || l == length(x)
                push!(sbin,l)
                i = l+1
                break
            end
        end
    end
    return sbin
end




"""
    GenFullObsData(datasets;verbose=true)

Gneerate input data basing on the available datasets.

# Arguments

- `datasets` array of Jspec dictionaries.
- `verbose` enable warning message.

It reports three arrays: inputdata, error on inputdata, input energies.


# Examples
```julia
GenFullObsData([dataset1, dataset2])
```
"""
function GenFullObsData(datasets;verbose=true)
    inp = Vector{Float64}()
    einp = Vector{Float64}()
    eninp = Vector{Float64}()
    for d in datasets
        if d["Instrument"] == "Other"
            if "ImportedData" in keys(d) && d["ImportedData"]
                append!(inp,d["PhFlux"])
                append!(einp,d["PhFluxErr"])
                append!(eninp,d["Energy"][!,"E"])
            else
                if verbose
                    println("Warning! Dataset ",d["Name"]," not properly processed yet.")
                end
            end
        else
            if "RebinnedData" in keys(d) && d["RebinnedData"] && "RebinnedAncillaryData" in keys(d) && d["RebinnedAncillaryData"]
                append!(inp,d["RebinnedMaskedInputData"])
                append!(einp,d["RebinnedMaskedInputDataErr"])
                append!(eninp,d["RebinnedMaskedEnergy"])
            else
                if verbose
                    println("Warning! Dataset ",d["Name"]," not properly processed yet.")
                end
            end
        end
    end
    return inp,einp,eninp
end








"""
    GenRebin(x,rebs)::AbstractVector{Real}

Rebin input data following a given rebin schema.

# Arguments

- `x` input array.
- `rebs` array with rebin schema.


# Examples
```jldoctest

x = [1.,2.,3.,4.,]
rbs = [1,3,4]

Jspec.GenRebin(x,rbs)

# output

3-element Vector{Real}:
 1.0
 2.5
 4.0
```
"""
function GenRebin(x,rebs)::AbstractVector{Real}
    newa = zeros(Real,length(rebs))
    old = 1
    for l in enumerate(rebs)
        newa[l[1]] = sum(x[old:l[2]])/Float64(length(x[old:l[2]]))
        old = l[2]+1
    end
    return newa
end



"""
    GenResponseMatrix(ds::Dict; verbose=true)

Generate a rebinned response matrix following the rebin schema identified for input data.

# Arguments

- `ds` Jspec data set dictionary.


# Examples
```julia

GenResponseMatrix(newdataset)

```
"""
function GenResponseMatrix(ds::Dict; verbose=true)
    if "RebinnedData" in keys(ds) && ds["RebinnedData"] && "RebinnedAncillaryData" in keys(ds) && ds["RebinnedAncillaryData"]
        if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
            matx = ds["MaskedRMF"][!,"MATRIX"] .* ds["ARF"][!,"SPECRESP"]
        elseif uppercase(ds["Instrument"]) == uppercase("Swift-BAT")
            matx = ds["MaskedRMF"][!,"MATRIX"]
        elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
            matx = ds["MaskedRMF"][!,"MATRIX"] .* ds["ARF"][!,"SPECRESP"]
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
            matx = ds["MaskedRMF"][!,"MATRIX"] .* ds["ARF"][!,"SPECRESP"]
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
            matx = ds["MaskedRMF"][!,"MATRIX"] .* ds["ARF"][!,"SPECRESP"]
        end
        #
        vt = zeros(length(ds["RebinSchema"]),length(matx))
        for i in 1:size(vt)[2]
            vt[:,i] = Jspec.GenRebin(matx[i],ds["RebinSchema"])
        end
        ds["RebinnedMaskedRMF"] = vt
        ds["RebinnedResponseMatrix"] = true
    else
        if verbose
            println("Warning! Data not fully rebinned yet!")
        end
    end
end



"""
    GetKnownInstruments()

Returns the instruments currently supported by the Jspec package.


# Examples
```julia
@show GetKnownInstruments()
```
"""
function GetKnownInstruments()
    return KnownInstruments
end





"""
    IgnoreChannels(ds:Dict, chns; verbose=true)

Ignore channels in the input data.

# Arguments

- `ds` Jspec data set dictionary.
- `chns` vector of channels to be ignored,
    e.g. [0,1,2,3] or even [0:4, 1000:1023]. Pay attention that channel numbering starts at 0.
- `verbose` enables warning messages.


# Examples
```julia

IgnoreChannels(newdataset,[0,1,2,3])
```
"""
function IgnoreChannels(ds::Dict,chns; verbose=true)
    if !ds["ImportedData"]
        if verbose
            println("Warning! Data not imported yet.")
        end
    elseif uppercase(ds["Instrument"]) == uppercase("Other")
        if verbose
            println("Warning! Data are not from a multi-channel instrument.")
        end
        ds["IgnoredChannels"] = false
    else
        mask = ones(Bool, length(ds["InputData"]))
        for e in chns
            mask[intersect(1:end, e .+ 1)] .= false
        end
        ds["MaskedInputData"] = ds["InputData"][mask]
        ds["MaskedInputDataErr"] = ds["InputDataErr"][mask]
        if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
            ds["MaskedInputSrcData"] = ds["InputSrcData"][mask]
            ds["MaskedInputBckData"] = ds["InputBckDataCorr"][mask]
        elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
            ds["MaskedInputSrcData"] = ds["InputSrcData"][mask]
            ds["MaskedInputBckData"] = ds["InputBckDataCorr"][mask]
        elseif uppercase(ds["Instrument"]) == uppercase("SVOM-MXT")
            ds["MaskedInputSrcData"] = ds["InputSrcData"][mask]
            ds["MaskedInputBckData"] = ds["InputBckDataCorr"][mask]
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
            ds["MaskedInputSrcData"] = ds["InputSrcData"][mask]
            ds["MaskedInputBckData"] = ds["InputBckDataCorr"][mask]
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
            ds["MaskedInputSrcData"] = ds["InputSrcData"][mask]
            ds["MaskedInputBckData"] = ds["InputBckDataCorr"][mask]
        end
        #ds["MaskedChanNumber"] = ds["ChanNumber"][mask]
        #
        ds["MaskedChannels"] = copy(ds["Channels"])
        deleteat!(ds["MaskedChannels"],findall(iszero,mask))
        ds["MaskedChannels"][!,"CHANNEL"] .= 0:nrow(ds["MaskedChannels"])-1
        #
        ds["MaskedRMF"] = copy(ds["RMF"])
        for i in 1:nrow(ds["RMF"])
            ds["MaskedRMF"][i,"MATRIX"] = ds["RMF"][i,"MATRIX"][mask]
        end
        ds["IgnoredChannels"] = true
    end
end




"""
    ImportData(ds::Dict; rmffile::String="", arffile::String="", srcfile::String="", bckfile::String="", verbose=true)

Import data from "multi-channel" instruments (e.g., Swift-XRT).

# Arguments

- `ds`` Jspec data set dictionary.
- `rmfile`` RMF response matrix.
- `arffile` effective area matrix.
- `srcfile` source counts (or rate).
- `bckfile` background counts (or rate).
- `verbose` enables warning messages.


# Examples

```julia
fnrmf = joinpath("data","wt.rmf")
fnarf = joinpath("data","wt.arf")
fnpisrc = joinpath("data","wtsource.pi")
fnpibck = joinpath("data","tback.pi");

ImportData(newdataset, rmffile=fnrmf,arffile=fnarf,srcfile=fnpisrc,bckfile=fnpibck)
```
"""
function ImportData(ds::Dict; rmffile::String="", arffile::String="", srcfile::String="", bckfile::String="", verbose=true)
    if haskey(ds,"Created") && !ds["Created"]
        if verbose
            println("Warning! Dataset not created yet.")
        end
    elseif ds["Instrument"] == "Other"
        if verbose
            println("Warning! Use ImportOtherData instead.")
        end
    elseif uppercase(ds["Instrument"]) ∉ [uppercase(i) for i in Jspec.KnownInstruments]
        if verbose
            println("Warning! Unknown intrument.")
        end
    else
        if rmffile != ""
            rmf = FITS(rmffile)
            if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
                ds["RMF"] = DataFrame(rmf[3])
                ds["Channels"] = DataFrame(rmf[2])
                ds["ChanNumber"] = ds["Channels"][!,"CHANNEL"]
                rmfid = 3
            elseif uppercase(ds["Instrument"]) == uppercase("Swift-BAT")
                ds["RMF"] = DataFrame(rmf[2])
                ds["Channels"] = DataFrame(rmf[3])
                ds["ChanNumber"] = ds["Channels"][!,"CHANNEL"]
                rmfid = 2
            elseif uppercase(ds["Instrument"]) == uppercase("SVOM-MXT")
                ds["RMF"] = DataFrame(rmf[2])
                ds["Channels"] = DataFrame(rmf[3])
                ds["ChanNumber"] = ds["Channels"][!,"CHANNEL"]
                rmfid = 2
            elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
                ds["RMF"] = DataFrame(rmf[3])
                ds["Channels"] = DataFrame(rmf[2])
                ds["ChanNumber"] = ds["Channels"][!,"CHANNEL"]
                rmfid = 3
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
                ds["RMF"] = DataFrame(rmf[2])
                ds["Channels"] = DataFrame(rmf[3])
                ds["ChanNumber"] = ds["Channels"][!,"CHANNEL"]
                rmfid = 2
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
                ds["RMF"] = DataFrame(rmf[2])
                ds["Channels"] = DataFrame(rmf[3])
                ds["ChanNumber"] = ds["Channels"][!,"CHANNEL"]
                rmfid = 2
            end
            #
            ds["Channels"][!,"E"] = (ds["Channels"][!,"E_MIN"] + ds["Channels"][!,"E_MAX"])/2.
            #
            en = (ds["RMF"][!,"ENERG_LO"] .+ ds["RMF"][!,"ENERG_HI"]) ./ 2
            de = (ds["RMF"][!,"ENERG_HI"] .- ds["RMF"][!,"ENERG_LO"])
            ds["Energy"] = DataFrame(E=en,ΔE=de, MinE=ds["RMF"][!,"ENERG_LO"], MaxE=ds["RMF"][!,"ENERG_HI"])
            #
            if !("MATRIX" in names(ds["RMF"]))
                if verbose
                    println("Warning! Data appear to be grouped or rebinned already.")
                end
                mtr = read(rmf[rmfid],"MATRIX")
                if !("N_CHAN" in names(ds["RMF"]))
                    nch = read(rmf[rmfid],"N_CHAN")
                else
                    nch = ds["RMF"][!,:N_CHAN]
                end
                if !("F_CHAN" in names(ds["RMF"]))
                    fch = read(rmf[rmfid],"F_CHAN")
                else
                    fch = ds["RMF"][!,:F_CHAN]
                end
                #
                headr = read_header(rmf[rmfid])
                #
                cols = []
                for i in 1:nrow(ds["RMF"])
                    vtrx = zeros(headr["TLMAX4"]+1)
                    idx = 1
                    for g in 1:ds["RMF"][i,:N_GRP]
                        chstart = fch[i][g]
                        for l in idx:idx+nch[i][g]-1
                            vtrx[chstart+l-idx+1] = mtr[i][l]
                        end
                        idx = idx+nch[i][g]
                    end
                    push!(cols,vtrx)
                end
                #
                ds["RMF"][!,"MATRIX"] = cols
                ds["RMF"][!,"F_CHAN"] .= 0
                ds["RMF"][!,"N_CHAN"] .= 1
                ds["RMF"][!,"N_GRB"] .= 1
                ds["DataOriginallyGrouped"] = true
                if verbose
                    println("Warning! Data reprocessed.")
                end
            end
        end
        end
        if arffile != ""
            arf = FITS(arffile)
            if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
                ds["ARF"] = DataFrame(arf[2])
            elseif uppercase(ds["Instrument"]) == uppercase("SVOM-MXT")
                ds["ARF"] = DataFrame(arf[2])
            elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
                ds["ARF"] = DataFrame(arf[2])
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
                ds["ARF"] = DataFrame(arf[2])
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
                ds["ARF"] = DataFrame(arf[2])
            end
        end
        if srcfile != ""
            pisrc = FITS(srcfile)
            if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
                ds["SrcCnt"] = DataFrame(pisrc[2])
                ds["InputSrcData"] = ds["SrcCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("Swift-BAT")
                ds["SrcRate"] = DataFrame(pisrc[2])
                ds["InputSrcRate"] = ds["SrcRate"][!,"RATE"]
                ds["InputSrcRateErr"] = ds["SrcRate"][!,"STAT_ERR"]
                ds["InputSrcRateSysErr"] = ds["SrcRate"][!,"SYS_ERR"]
            elseif uppercase(ds["Instrument"]) == uppercase("SVOM-MXT")
                ds["SrcCnt"] = DataFrame(pisrc[2])
                ds["InputSrcData"] = ds["SrcCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
                ds["SrcCnt"] = DataFrame(pisrc[2])
                ds["InputSrcData"] = ds["SrcCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
                ds["SrcCnt"] = DataFrame(pisrc[2])
                ds["InputSrcData"] = ds["SrcCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
                ds["SrcCnt"] = DataFrame(pisrc[2])
                ds["InputSrcData"] = ds["SrcCnt"][!,"COUNTS"]
            end
            heasrc = read_header(pisrc[2])
            ds["SrcExpTime"] = heasrc["EXPOSURE"]
            ds["SrcBackScal"] = heasrc["BACKSCAL"]
        end
        if bckfile != ""
            pibck = FITS(bckfile)
            if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
                ds["BckCnt"] = DataFrame(pibck[2])
                ds["InputBckData"] = ds["BckCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("SVOM-MXT")
                ds["BckCnt"] = DataFrame(pibck[2])
                ds["InputBckData"] = ds["BckCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
                ds["BckCnt"] = DataFrame(pibck[2])
                ds["InputBckData"] = ds["BckCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
                ds["BckCnt"] = DataFrame(pibck[2])
                ds["InputBckData"] = ds["BckCnt"][!,"COUNTS"]
            elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
                ds["BckCnt"] = DataFrame(pibck[2])
                ds["InputBckData"] = ds["BckCnt"][!,"COUNTS"]
            end
            #
            heabck = read_header(pibck[2])
            ds["BckExpTime"] = heabck["EXPOSURE"]
            ds["BckBackScal"] = heabck["BACKSCAL"]
        else
            ds["BckExpTime"] = 1.
            ds["BckBackScal"] = 1.
        end
        ds["BackScaleRatio"] = ds["SrcBackScal"]/ds["BckBackScal"]
        ds["ExposureRatio"] = ds["SrcExpTime"]/ds["BckExpTime"]
        #
        if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
            ds["InputBckDataCorr"] = ds["InputBckData"]*ds["BackScaleRatio"]*ds["ExposureRatio"]
            ds["InputData"] = ds["InputSrcData"] .- ds["InputBckDataCorr"]
            ds["InputDataErr"] = sqrt.(ds["InputSrcData"] .+ ds["InputBckDataCorr"])
        elseif uppercase(ds["Instrument"]) == uppercase("Swift-BAT")
            ds["InputData"] = ds["InputSrcRate"]
            ds["InputDataErr"] = ds["InputSrcRateErr"]
        elseif uppercase(ds["Instrument"]) == uppercase("SVOM-MXT")
            ds["InputBckDataCorr"] = ds["InputBckData"]*ds["BackScaleRatio"]*ds["ExposureRatio"]
            ds["InputData"] = ds["InputSrcData"] .- ds["InputBckDataCorr"]
            ds["InputDataErr"] = sqrt.(ds["InputSrcData"] .+ ds["InputBckDataCorr"])
        elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
            ds["InputBckDataCorr"] = ds["InputBckData"]*ds["BackScaleRatio"]*ds["ExposureRatio"]
            ds["InputData"] = ds["InputSrcData"] .- ds["InputBckDataCorr"]
            ds["InputDataErr"] = sqrt.(ds["InputSrcData"] .+ ds["InputBckDataCorr"])
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
            ds["InputBckDataCorr"] = ds["InputBckData"]*ds["BackScaleRatio"]*ds["ExposureRatio"]
            ds["InputData"] = ds["InputSrcData"] .- ds["InputBckDataCorr"]
            ds["InputDataErr"] = sqrt.(ds["InputSrcData"] .+ ds["InputBckDataCorr"])
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
            ds["InputBckDataCorr"] = ds["InputBckData"]*ds["BackScaleRatio"]*ds["ExposureRatio"]
            ds["InputData"] = ds["InputSrcData"] .- ds["InputBckDataCorr"]
            ds["InputDataErr"] = sqrt.(ds["InputSrcData"] .+ ds["InputBckDataCorr"])
        end
        #
    ds["ImportedData"] = true
end



"""
    ImportOtherData(ds::Dict, energy, phflux, ephflux; bandwidth=1., verbose=true)

Import data already in physical units.

# Arguments

- `ds` Jspec dictionary.
- `energy` input energy (KeV).
- `phflux` photon flux density (``photons~cm{^-2}~s{^-1}~KeV{^-1})``.
- `ephflux` photon flux density uncertainty.
- `bandwidth` band width (KeV).
- `verbose` enable warning message.

Bandwidth is needed only in case photon flux (``photons~cm{^-2}~s{^-1})``, rather then photon flux
density, is provided.


# Examples
```julia
ImportOtherData(newdataset, [1.,2.,3.,4], [0.1,0.2,0.3,0.4], [0.01,0.02,0.03,0.04])
```
"""
function ImportOtherData(ds::Dict, energy, phflux, ephflux; bandwidth=1., verbose=true)
    if haskey(ds,"Created") && !ds["Created"]
        if verbose
            println("Warning! Dataset not created yet.")
        end
    elseif uppercase(ds["Instrument"]) == uppercase("Other")
        ds["Energy"] = DataFrame(E=energy,ΔE=bandwidth)
        #ds["Energy"] = energy
        ds["PhFlux"] = phflux
        ds["PhFluxErr"] = ephflux
        #ds["BandWidth"] = bandwidth
        ds["RMF"] = I
        ds["ImportedData"] = true
        # Steps added to performace reasons
        ds["RebinnedMaskedRMF"] = I
        ds["RebinnedResponseMatrix"] = true
    else
        if verbose
            println("This function can be used only for data declared as 'Other'.")
        end
        ds["ImportedData"] = false
    end
end



"""
JspecFunc(pars,dts,inpfnc)

Convolve the output of the `inpfnc` function with the responce matrices of the imported data.


# Arguments

- `pars` parameters for the `inpfnc` function.
- `dts` array of Jspec dictionaries.
- `inpfnc` function to model the imported data.

`inpfnc` can be any legal `Julia` function. THe function should be declared as in the following example.


# Examples
```julia

function Myfunc(pars,Energy)
    N, λ = pars
    return anyfunc(E,N,λ)
end

JspecFunc([N,λ], [Optdt,XRTdt], Myfunc)
```
"""
function JspecFunc(pars,dts,inpfnc)
    res = map(d -> d["RebinnedMaskedRMF"] * (inpfnc(pars,d["Energy"][!,"E"]) .* d["Energy"][!,"ΔE"]), dts)
    return collect(Iterators.flatten(res))
end



"""
    Jy2PhFlux(energy,jyspectrum)

Convert an input sectrum in ``Jy`` to ``ph~s^{-1}~cm^{-2}~KeV^{-1}``.

# Arguments

- `energy` input energy of the spectrum in ``KeV``.
- `jyspectrum` input spectrum in ``Jy``.



# Examples
```jldoctest

e = [1.,2.,3.,4.,]
sp = [1e-3,3e-3,4e-3,5e-3]

Jy2PhFlux(e,sp)

# output

4-element Vector{Float64}:
 1.51
 2.265
 2.013333333333333
 1.8875
"""
function Jy2PhFlux(energy,jyspectrum)
    return jyspectrum .* 1.51e3 ./ energy
end




"""
    KeV2Angstrom(energy)

Convert photon energy (``KeV``) to wavelengths in Angstrom (``10^{-8}~m``).

#Arguments

- `energy` the input photon energy.


# Examples
```jldoctest

KeV2Angstrom(1:3)

# output

0.08064516129032258:0.08064516129032258:0.24193548387096775
```
"""
function KeV2Angstrom(energy)
    return energy ./ 12.4
end




"""
    KeV2Channel(ds::Dict, energy)

Convert photon energy (``KeV``) to original detector channel.

#Arguments

- `ds` Jspec dictionary.
- `energy` the input photon energy.

It returns `-1` if the energy is not in the covered range.


# Examples
```julia

KeV2Channel(ds,1.2)

# output

11
```
"""
function KeV2Channel(ds::Dict,energy)
    minE = minimum(ds["Channels"][!,:E_MIN])
    maxE = maximum(ds["Channels"][!,:E_MAX])
    #
    if !(minE <= energy <= maxE)
      return -1
    else
      fdt = filter(row -> row.E_MIN >= energy, ds["Channels"])
      return fdt[!,:CHANNEL][1]
    end
end







"""
    PlotRaw(ds:Dict; xlbl="Channels", ylbl="Counts", tlbl=ds.Name, verbose=true)::Figure

Draw a plot of the raw input data.

# Arguments

- `ds` Jspec data set dictionary.
- `xlbl` x-axis label.
- `ylbl` y-axis label.
- `tlbl` plot title.
- `verbose` enables warning messages.

# Examples
```julia
figraw = PlotRaw(newdataset)
```
"""
function PlotRaw(ds::Dict; xlbl="Channels", ylbl=L"Counts ch$^{-1}$", tlbl=ds["Name"], verbose=true)::Figure
    fig = Figure(fontsize=30)
    #
    ax = Axis(fig[1, 1],
        spinewidth=3,
        xlabel = xlbl,
        ylabel = ylbl,
        title = tlbl,
        #yscale=log10,
        #xscale=log10,
        )
    if !ds["ImportedData"]
        if verbose
            println("Warning! Data not imported yet.")
        end
    elseif uppercase(ds["Instrument"]) == uppercase("Other")
        scatter!(ds["Energy"][!,"E"],ds["PhFlux"],color=(:orange,0.2))
        errorbars!(ds["Energy"][!,"E"],ds["PhFlux"],ds["PhFluxErr"],color=(:orange,0.2))
    else
        scatter!(ds["ChanNumber"],ds["InputData"],color=(:orange,0.2))
        errorbars!(ds["ChanNumber"],ds["InputData"],ds["InputDataErr"],color=(:orange,0.2))
    end
    #
    return fig
end




"""
    PlotRebinned(ds:Dict; xlbl="Channels", ylbl="Counts", tlbl=ds.Name, verbose=true)::Figure

Draw a plot of the rebinned input data.

# Arguments

- `ds` Jspec data set dictionary.
- `xlbl` x-axis label.
- `ylbl` y-axis label.
- `tlbl` plot title.
- `verbose` enables warning messages.


# Examples
```julia
figreb = PlotRebinned(newdataset)
```
"""
function PlotRebinned(ds::Dict; xlbl="Channels", ylbl=L"Counts ch$^{-1}$ s$^{-1}$", tlbl=ds["Name"],verbose=true)::Figure
    fig = Figure(fontsize=30)
    #
    ax = Axis(fig[1, 1],
        spinewidth=3,
        xlabel = xlbl,
        ylabel = ylbl,
        title = tlbl,
        #yscale=log10,
        #xscale=log10,
        )
    if uppercase(ds["Instrument"]) == uppercase("Other")
        scatter!(ds["Energy"],ds["PhFlux"],color=(:orange,0.2))
        errorbars!(ds["Energy"],ds["PhFlux"],ds["PhFluxErr"],color=(:orange,0.2))
    elseif "RebinnedData" in keys(ds) && ds["RebinnedData"] && "RebinnedAncillaryData" in keys(ds) && ds["RebinnedAncillaryData"]
      scatter!(ds["RebinnedMaskedChannel"],ds["RebinnedMaskedInputData"],color=(:orange,0.2))
      errorbars!(ds["RebinnedMaskedChannel"],ds["RebinnedMaskedInputData"],ds["RebinnedMaskedInputDataErr"],color=(:orange,0.2))
    else
      if verbose
        println("Warning! Data not fully rebinned yet.")
      end
    end
    #
    return fig
end



"""
    RebinAncillaryData(ds::Dict; verbose=true)

Rebin ancillary data (channels, channel energy, etc.) with the rebin schema identified for input data.

# Arguments

- `ds` Jspec data set disctionary.
- `verbose` enables warning messages.


# Examples
```julia

RebinAncilaryData(newdataset)

```
"""
function RebinAncillaryData(ds::Dict; verbose=true)
    if "RebinnedData" in keys(ds) && !ds["RebinnedData"]
        if verbose
            println("Warning! Channels not rebinned yet.")
        end
        ds["RebinnedAncillaryData"] = false
    else
        ds["RebinnedMaskedEnergy"] = Jspec.GenRebin(ds["MaskedChannels"][!,"E"],ds["RebinSchema"])
        ds["RebinnedMaskedChannel"] = Jspec.GenRebin(ds["MaskedChannels"][!,"CHANNEL"],ds["RebinSchema"])
        ds["RebinnedAncillaryData"] = true
    end
end



"""
    RebinData(ds::Dict;minSN=5,verbose=true)

Rebin input data with a mininum S/N per bin.

# Arguments

- `ds` Jspec data set disctionary.
- `minSN` minimum S/N per bin.
- `verbose` enables warning messages.


# Examples
```julia

RebinData(newdataset)

```
"""
function RebinData(ds::Dict;minSN=5,verbose=true)
    if !ds["IgnoredChannels"]
        if verbose
            println("Warning! Channels not ignored yet.")
        end
        ds["RebinnedData"] = false
    else
        ds["RebinSchema"] = Jspec.FindRebinSchema(ds["MaskedInputData"],ds["MaskedInputDataErr"],minSN=minSN)
        ncts = zeros(Real,length(ds["RebinSchema"]))
        encts = zeros(Real,length(ds["RebinSchema"]))
        old = 1
        for l in enumerate(ds["RebinSchema"])
            c = sum(ds["MaskedInputData"][old:l[2]])
            b = sqrt(sum(ds["MaskedInputDataErr"][old:l[2]].^2))
            ncts[l[1]] = c/length(ds["MaskedInputData"][old:l[2]])
            encts[l[1]] = b/length(ds["MaskedInputDataErr"][old:l[2]])
            old = l[2]+1
        end
        if uppercase(ds["Instrument"]) == uppercase("Swift-XRT")
            ds["RebinnedMaskedInputData"] = ncts/ds["SrcExpTime"]
            ds["RebinnedMaskedInputDataErr"] = encts/ds["SrcExpTime"]
        elseif uppercase(ds["Instrument"]) == uppercase("Swift-BAT")
            ds["RebinnedMaskedInputData"] = ncts
            ds["RebinnedMaskedInputDataErr"] = encts
        elseif uppercase(ds["Instrument"]) == uppercase("NuSTAR-FPM")
            ds["RebinnedMaskedInputData"] = ncts/ds["SrcExpTime"]
            ds["RebinnedMaskedInputDataErr"] = encts/ds["SrcExpTime"]
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EMOS")
            ds["RebinnedMaskedInputData"] = ncts/ds["SrcExpTime"]
            ds["RebinnedMaskedInputDataErr"] = encts/ds["SrcExpTime"]
        elseif uppercase(ds["Instrument"]) == uppercase("XMM-EPN")
            ds["RebinnedMaskedInputData"] = ncts/ds["SrcExpTime"]
            ds["RebinnedMaskedInputDataErr"] = encts/ds["SrcExpTime"]
        end
        #
        ds["RebinnedData"] = true
    end
end






end
