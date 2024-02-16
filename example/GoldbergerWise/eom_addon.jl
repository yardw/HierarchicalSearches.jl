using ForwardDiff
function φtest(F′, F, y, l², k, γ², m²)
    return -3ϕP^2/(2u*l²*ϕ0(y)) * (F′ - 2A′(y, l², k, γ²)*F)  #(3.12)
end
φtest(F′, F, l², k, γ², m²) = y->φtest(F′, F, y, l², k, γ², m²)
φtest(F′, F, params) = φtest(F′, F, params...)
φ′test(F′, F, l², k, γ², m²) = y->ForwardDiff.derivative(φtest(F′, F, l², k, γ², m²), y)
φ′test(F′, F, params) = φ′test(F′, F, params...)