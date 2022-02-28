```@meta
CurrentModule = MoleculeSpectrum
```

```@docs
State(N, mₙ, I₁, mᵢ₁, I₂, mᵢ₂)
```

```@docs
    plot_transition_strengths(spectrum::Spectrum, g::State, frequency_range; polarization::SphericalUnitVector=Unpolarized())
```

```@docs
    transition_strengths(spectrum::Spectrum, g::State, frequency_range; polarization::SphericalUnitVector=Unpolarized())
```