using Test
using Logging
using ThreeBodyDecays

@testset "Recipes expand without warnings" begin
    # Expand the actual definitions even when the package was loaded from cache.
    source = read(joinpath(pkgdir(ThreeBodyDecays), "src", "dalitz.jl"), String)
    expressions = Meta.parseall(source)
    recipes = filter(expressions.args) do expression
        expression isa Expr &&
            expression.head == :macrocall &&
            expression.args[1] == Symbol("@recipe")
    end
    @test !isempty(recipes)
    for recipe in recipes
        @test_logs min_level = Logging.Warn macroexpand(ThreeBodyDecays, recipe)
    end
end

# Supply the backend hook without requiring a plotting backend.
ThreeBodyDecays.RecipesBase.is_key_supported(::Symbol) = true

@testset "Dalitz projection keywords" begin
    ms = ThreeBodyMasses(0.1, 0.2, 0.3; m0 = 2.0)
    intensity = σs -> 1.0
    integrator = (f, a, b) -> ((b-a) * f((a+b)/2), 0.0)
    recipe = DalitzProjection((ms, intensity, integrator))
    apply_recipe = ThreeBodyDecays.RecipesBase.apply_recipe
    @test_throws ArgumentError apply_recipe(Dict{Symbol,Any}(), recipe)
    for k = 1:3, bins in (8, 100)
        attributes = Dict{Symbol,Any}(:k => k)
        bins == 8 && (attributes[:bins] = bins)
        result = only(apply_recipe(attributes, recipe))
        xs, ys = result.args
        @test length(xs) == length(ys) == bins
        @test all(isfinite, ys)
        @test result.plotattributes[:seriestype] == :path
    end
end
