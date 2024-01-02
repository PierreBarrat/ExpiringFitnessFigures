# ExpiringFitnessFigures

## Install

- install [julia v1.9.4](https://julialang.org/downloads/). *Note*: [juliaup](https://github.com/JuliaLang/juliaup) makes managing julia versions easier.
- `git clone` this somewhere. Navigate to the repository. 
- start julia and enter the following

```
using Pkg; 
Pkg.activate(".") # activates the local environment
Pkg.instantiate() # download and precompile all packages
```

  Note that precompilation could take some time, especially `PartialSweepSIR` which depends on `OrdinaryDiffEq`.

- start `Pluto` by typing `using Pluto; Pluto.run()`. Once the Pluto server is ready, use it to navigate to the `notebooks/` folder and open any notebook. Running each notebook will generate a `.png` file that will be saved to `figures/`. Notebooks in the `intermediate_notebooks/` subfolder perform simulations necessary to generate the data used in the main notebooks. 
- Alternatively, you can run notebooks as julia scripts: `julia notebook_name.jl`. 

*Note*: Some hints about how to use Pluto can be found on the [github page](https://github.com/fonsp/Pluto.jl), in the [wiki/faq](https://github.com/fonsp/Pluto.jl/wiki) and in the featured notebooks shown on the welcome screen. 




