using CairoMakie
using Jspec
using Test

@testset "Jspec.jl" begin
    # Write your tests here.
    # CreateDataSet
    newdt = CreateDataSet("XRTTest","Swift-XRT")
    @test newdt["Instrument"] == "Swift-XRT"
    #
    # ImportData
    frmf = joinpath("testdata","xrtwt.rmf")
    farf = joinpath("testdata","xrtwt.arf")
    fsrc = joinpath("testdata","xrtwtsource.pi")
    fbck = joinpath("testdata","xrtwtback.pi")
    ImportData(newdt, rmffile=frmf, arffile=farf, srcfile=fsrc, bckfile=fbck)
    @test newdt["ImportedData"] == true
    #
    # ImportOtherData
    newodt = CreateDataSet("OptData","Other")
    ImportOtherData(newodt, [1.,2.,3.,4], [0.1,0.2,0.3,0.4], [0.01,0.02,0.03,0.04])
    @test newodt["ImportedData"] == true
    #
    # PlotRaw
    @test typeof(PlotRaw(newdt)) == Figure
    #
    # IgnoreChannels
    IgnoreChannels(newdt,[0,1,2,3])
    @test newdt["IgnoredChannels"] == true
    #
    # FindRebinSchema
    @test Jspec.FindRebinSchema([1.,2.,3.,4.,],[0.1,0.5,0.6,0.05]) == [1, 3, 4]
    #
    # GenRebin
    @test Jspec.GenRebin([1.,2.,3.,4.],[1,3,4]) == [1.0,2.5,4.0]
    #
    # RebinData
    RebinData(newdt)
    @test newdt["RebinnedData"] == true
    #
    # RebinAncillaryData
    RebinAncillaryData(newdt)
    @test newdt["RebinnedAncillaryData"] == true
    #
    # PlotRebinned
    @test typeof(PlotRebinned(newdt)) == Figure
    #
    # GenResponseMatrix
    GenResponseMatrix(newdt)
    @test newdt["RebinnedResponseMatrix"] == true
    #
    # GenFullObsData
    @test typeof(GenFullObsData([newodt,newdt])) == Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    #
    # JspecFunc
    function Myfunc(pars,E)
        A,B = pars
        return (A.+B).*E
    end
    @test JspecFunc([1.,2.],[newodt,],Myfunc) == [3.0,6.0,9.0,12.0]
    #
    # Angstrom2KeV
    @test Angstrom2KeV(5000:5004) == [0.00248,0.002479504099180164,0.0024790083966413435,0.0024785128922646415,0.0024780175859312552]
    #
    # Jy2PhFlux
    e = [1.,2.,3.,4.,]
    sp = [1e-3,3e-3,4e-3,5e-3]
    @test Jy2PhFlux(e,sp) == [1.51,2.265,2.013333333333333,1.8875]
    #
    # KeV2Angstrom
    @test KeV2Angstrom(1:3) == 0.08064516129032258:0.08064516129032258:0.24193548387096775
    #
    # KeV2Channel
    @test KeV2Channel(newdt,1.2) == 120
    #
end
