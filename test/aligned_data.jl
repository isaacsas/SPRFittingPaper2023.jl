using SPRFittingPaper2023, CSV, Test

datadir = joinpath(@__DIR__,"..", "data")

# test csv loading to aligned data
fname = joinpath(datadir, "tester-13.8_aligned.csv")
times = [[5.0,6.0,7.0,8.0,9.0,10.0,11.0],[5.0,6.0,7.0,8.0]]
refdata = [[.821,1.12,1.22,1.52,1.42,1.82,1.92],[3.80,4.30,4.80,5.50]]
antigenconcen = 13.8
antibodyconcens = [25.0,100.0]
ad = get_aligned_data(fname)
@test times == ad.times
@test refdata == ad.refdata
@test antigenconcen == ad.antigenconcen
@test antibodyconcens == ad.antibodyconcens