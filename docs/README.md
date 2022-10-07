Build docs:
    * Install Documenter in Julia: `] add Documenter.jl`
    * Compile source files: `julia make.jl`(this will generate html files in `build/`)

Launching docs page (`build/`) locally:
    * Install LiveServer in Julia: `] add LiveServer`
    * Launch live server: `julia -e 'using LiveServer; serve(dir="docs/build")'`
    * Open the given URL in your browser
