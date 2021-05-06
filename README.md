# ChemicalIdentifiers

A chemical identifiers search package, using the databases present at CalebBell/chemicals.
# Warning
Work in progress, many things could change.

## Instalation:
```
using Pkg
Pkg.add(https://github.com/longemen3000/ChemicalIdentifiers.jl)
```
The databases are downloaded, parsed,processed and stored as Apache Arrow files at the first package usage, so the first usage may take some time.

## usage
This package exports `search_chemical`, that, given a search string, performs a search on a database of over 70000 compounds, returning a Named Tuple with the identifiers of the substance in question. 

```
julia> using ChemicalIdentifiers
julia> res = search_chemical("water")
(pubchemid = 962, CAS = (7732, 18, 5), formula = "H2O", MW = 18.01528, smiles = "O", InChI = "H2O/h1H2", InChI_key = "XLYOFNOQVPJJNP-UHFFFAOYSA-N", iupac_name = "oxidane", common_name = "water")
```

The package stores each query in `ChemicalIdentifiers.SEARCH_CACHE` as a `Dict{String,Any}`, so subsequent queries on the same (or similar) strings, dont pay the cost of searching in the database.

If you don't want to store the query, you could use `search_chemical(query,nothing)`, or, if you want your own cache to be used, pass your own cache via `search_chemical(query,mycache)`

This package is slower than the original python package, because we don't build a dictionary of all existent synonyms (about 700K synonyms), but perform a sorted search. this is in exchange for a quicker loading time on rarer queries and a lower memory footprint.
