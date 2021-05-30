# ChemicalIdentifiers

A chemical identifiers search package, using the databases present at CalebBell/chemicals.
## Warning
Work in progress, many things could change.

## Instalation:
```
using Pkg
Pkg.add(https://github.com/longemen3000/ChemicalIdentifiers.jl)
```
The databases are downloaded, parsed,processed and stored as Apache Arrow files at the first package usage, so the first usage may take some time.

## Usage
This package exports `search_chemical`, that, given a search string, performs a search on a database of over 70000 compounds, returning a Named Tuple with the identifiers of the substance in question. 

```
julia> using ChemicalIdentifiers
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

The package stores each query in `ChemicalIdentifiers.SEARCH_CACHE` as a `Dict{String,Any}`, so subsequent queries on the same (or similar) strings, dont pay the cost of searching in the database.

If you don't want to store the query, you could use `search_chemical(query,nothing)`, or, if you want your own cache to be used, pass your own cache via `search_chemical(query,mycache)`. 

#Custom Databases
If you want to add your own databases, you could use the (unexported) data utilities to do so. lets say we also want to add the inorganic database located at https://github.com/CalebBell/chemicals/blob/master/chemicals/Identifiers/Inorganic%20db.tsv. we could do:
```
using ChemicalIdentifiers
inorganic_url = "https://github.com/CalebBell/chemicals/blob/master/chemicals/Identifiers/Inorganic%20db.tsv"
ChemicalIdentifiers.load_data!(:inorganic,url = inorganic_url)
ChemicalIdentifiers.load_db!(:inorganic)
```
or if you already have a local database:

```
using ChemicalIdentifiers
filepath = "path/to/my/db.tsv"
ChemicalIdentifiers.load_data!(:custom,file = filepath)
ChemicalIdentifiers.load_db!(:custom)
```
`ChemicalIdentifiers.load_data!` will generate a named tuple of file paths (stored in `ChemicalIdentifiers.DATA_INFO`), and `ChemicalIdentifiers.load_db!` will use that data to generate the corresponding Apache Arrow files and store those in a [scratch](https://github.com/JuliaPackaging/Scratch.jl) space (`ChemicalIdentifiers.download_cache`). 

The raw databases are then stored in `ChemicalIdentifiers.DATA_DB`. if the data was already processed, then the arrow files are read directly, saving significant loading time.

In case of adding user databases, those are searched first, so there is a possibility of collision. 

