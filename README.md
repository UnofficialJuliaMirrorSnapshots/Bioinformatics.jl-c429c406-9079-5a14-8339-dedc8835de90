# Bioinformatics.jl

![GitHub release (latest by date)](https://img.shields.io/github/v/release/mrtkp9993/Bioinformatics.jl)
[![Build Status](https://travis-ci.org/mrtkp9993/Bioinformatics.jl.svg?branch=master)](https://travis-ci.org/mrtkp9993/Bioinformatics.jl)
[![codecov](https://codecov.io/gh/mrtkp9993/Bioinformatics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mrtkp9993/Bioinformatics.jl)
[![GitHub license](https://img.shields.io/github/license/mrtkp9993/Bioinformatics.jl.svg)](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/LICENSE)

## Usage

Just include the `Bioinformatics.jl` file with

```julia
include("Bioinformatics.jl")
```

or install the library with

```julia
(v1.1) pkg> add "https://github.com/mrtkp9993/Bioinformatics.jl"
```

and you're ready to use!

## Features

* Global alignment (Needleman-Wunsch algorithm with linear gap penalty)
* Dotplots
* Transcription, Reverse Complement, Translation, Open reading frames (ORFs)
* GC-content, minimum skew
* Protein mass, Extinction coefficient, Instability index, Isoelectric point, Grand average of hydropathicity (GRAVY) statistics of protein sequences

## Documentation

Please refer [here](https://mrtkp9993.github.io/Bioinformatics.jl/latest/).

## Contribute

If you find any bugs/issues, please create an issue or fix and create a pull request.

## License

GNU General Public License v3.0 © [Murat Koptur](https://github.com/mrtkp9993)
