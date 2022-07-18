using SafeTestsets

@time begin
    @time @safetestset "Aligned Data" begin include("aligned_data.jl") end
    @time @safetestset "Callbacks" begin include("callbacks.jl") end
    @time @safetestset "Monovalent Fitting" begin include("monovalent_fit.jl") end
end
