# ChemicalIdentifiers.jl

[![Build Status](https://github.com/longemen3000/ChemicalIdentifiers.jl/workflows/CI/badge.svg)](https://github.com/longemen3000/ChemicalIdentifiers.jl/actions)
[![codecov](https://codecov.io/gh/longemen3000/ChemicalIdentifiers.jl/graph/badge.svg?token=7OF7XX2TVJ)](https://codecov.io/gh/longemen3000/ChemicalIdentifiers.jl)

A chemical identifiers search package, using the databases present at CalebBell/chemicals.

## Instalation:
```
using Pkg
Pkg.add("ChemicalIdentifiers.jl")
```
The databases are downloaded, parsed, processed and stored as Apache Arrow files at the first package usage, so the first usage may take some time.

## Usage
This package exports `search_chemical`, that, given a search string, performs a search on a database of over 70000 compounds, returning a Named Tuple with the identifiers of the substance in question.

```
julia>using ChemicalIdentifiers
julia> res = search_chemical("water")
(pubchemid = 962, CAS = (7732, 18, 5), formula = "H2O", MW = 18.01528, smiles = "O", InChI = "H2O/h1H2", InChI_key = "XLYOFNOQVPJJNP-UHFFFAOYSA-N", iupac_name = "oxidane", common_name = "water")

#worst case scenario, not found on the present databases
julia> @btime search_chemical("dimethylpyruvic acid22",nothing)
  273.700 μs (264 allocations: 15.05 KiB)
missing

#common compound found in the short database
julia> @btime search_chemical("methane",nothing)
  7.075 μs (57 allocations: 2.97 KiB)
(pubchemid = 297, CAS = (74, 82, 8), formula = "CH4", MW = 16.04246, smiles = "C", InChI = "CH4/h1H4", InChI_key = "VNWKTOKETHGBQD-UHFFFAOYSA-N", iupac_name = "methane", common_name = "methane")
```
A query is usually a string and its type is detected automatically when possible. the supported query types are:

- **PubChemID** : by using any ``<:Integer` (or a string containing an Integer)
```
julia> search_chemical(8003)
(pubchemid = 8003, CAS = (109, 66, 0), formula = "C5H12", MW = 72.14878, smiles =
"CCCCC", InChI = "C5H12/c1-3-5-4-2/h3-5H2,1-2H3", InChI_key = "OFBQJSOFQDEBGM-UHFFFAOYSA-N", iupac_name = "pentane", common_name = "pentane")
```

- **CAS registry number** : by using a Tuple of integers or a string with the digits separated by `-` :

```
julia> search_chemical((67,56,1))
(pubchemid = 887, CAS = (67, 56, 1), formula = "CH4O", MW = 32.04186, smiles = "CO", InChI = "CH4O/c1-2/h2H,1H3", InChI_key = "OKKJLVBELUTLKV-UHFFFAOYSA-N", iupac_name = "methanol", common_name = "methanol")

 search_chemical((67,56,1),nothing) == search_chemical("67-56-1",nothing) #true
```

- **SMILES** : by using a string starting with `SMILES=`, or by passing the keyword argument (`by = :smiles`):
```
julia> search_chemical("SMILES=N")
(pubchemid = 222, CAS = (7664, 41, 7), formula = "H3N", MW = 17.03052, smiles = "N", InChI = "H3N/h1H3", InChI_key = "QGZKDVFQNNGYKY-UHFFFAOYSA-N", iupac_name = "azane", common_name = "ammonia")

julia> search_chemical("O",by = :smiles)
(pubchemid = 962, CAS = (7732, 18, 5), formula = "H2O", MW = 18.01528, smiles = "O", InChI = "H2O/h1H2", InChI_key = "XLYOFNOQVPJJNP-UHFFFAOYSA-N", iupac_name = "oxidane", common_name = "water")
```

- **InChI** : by using a string starting with `InChI=1/` or `InChI=1S/` :
```
julia> search_chemical("InChI=1/C2H4/c1-2/h1-2H2")
(pubchemid = 6325, CAS = (74, 85, 1), formula = "C2H4", MW = 28.05316, smiles = "C=C", InChI = "C2H4/c1-2/h1-2H2", InChI_key = "VGGSQFUCUMXWEO-UHFFFAOYSA-N", iupac_name = "ethene", common_name = "ethene")
```

- **InChI key** : by using a string with the pattern `XXXXXXXXXXXXXX-YYYYYYYYFV-P`:
```
julia> search_chemical("IMROMDMJAWUWLK-UHFFFAOYSA-N")
(pubchemid = 11199, CAS = (9002, 89, 5),
formula = "C2H4O", MW = 44.05256, smiles
= "C=CO", InChI = "C2H4O/c1-2-3/h2-3H,1H2", InChI_key = "IMROMDMJAWUWLK-UHFFFAOYSA-N", iupac_name = "ethenol", common_name
= "ethenol")
```
Searches by CAS and PubChemID are a little bit faster thanks to being encoded as native numeric types, other properties are stored as strings.

The package stores each query in `ChemicalIdentifiers.SEARCH_CACHE` as a `Dict{String,Any}`, so subsequent queries on the same (or similar) strings, dont pay the cost of searching in the database.

If you don't want to store the query, you could use `search_chemical(query,nothing)`, or, if you want your own cache to be used, pass your own cache via `search_chemical(query,mycache)`.

## Custom Databases
If you want to add your own databases, you could use the (unexported) data utilities to do so. lets say we also want to add the inorganic database located at https://github.com/CalebBell/chemicals/blob/master/chemicals/Identifiers/Inorganic%20db.tsv. we could do:
```julia
using ChemicalIdentifiers
inorganic_url = "https://github.com/CalebBell/chemicals/blob/master/chemicals/Identifiers/Inorganic%20db.tsv"
ChemicalIdentifiers.load_data!(:inorganic,url = inorganic_url)
ChemicalIdentifiers.load_db!(:inorganic)
```
or if you already have a local database:

```julia
using ChemicalIdentifiers
filepath = "path/to/my/db.tsv"
ChemicalIdentifiers.load_data!(:custom,file = filepath)
ChemicalIdentifiers.load_db!(:custom)
```
`ChemicalIdentifiers.load_data!` will generate a named tuple of file paths (stored in `ChemicalIdentifiers.DATA_INFO`), and `ChemicalIdentifiers.load_db!` will use that data to generate the corresponding Apache Arrow files and store those in a [scratchspace](https://github.com/JuliaPackaging/Scratch.jl) (stored in `ChemicalIdentifiers.download_cache`). This download cache can be cleaned (in case a download goes wrong) with `ChemicalIdentifiers.clear_download_cache!()`

The raw databases are then stored in `ChemicalIdentifiers.DATA_DB`. if the data was already processed, then the arrow files are read directly, saving significant loading time.

In case of adding user databases, those are searched first, so there is a possibility of collision.
