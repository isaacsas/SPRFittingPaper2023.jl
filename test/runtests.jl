using SafeTestsets

@time begin
    @time @safetestset "Aligned Data" begin include("aligned_data.jl") end
    @time @safetestset "Callbacks" begin include("callbacks.jl") end
end
