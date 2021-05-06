# ChemicalIdentifiers

a chemical identifiers search package, using the databases present at CalebBell/chemicals.
# Warning
Work in progress, many things could change

## Instalation:
```
using Pkg
Pkg.add(https://github.com/longemen3000/ChemicalIdentifiers.jl)
```

## usage
this package exports `search_chemical`, that, given a search string, performs a search on a database of over 70000 compounds, returning a Named Tuple with the identifiers of the substance in question. the databases are downloaded, parsed,processed and stored as arrow files at the first package usage, so the first instalation may take some time.

```
julia> using ChemicalIdentifiers
julia> res = search_chemical("water")
(pubchemid = 962, CAS = (7732, 18, 5), formula = "H2O", MW = 18.01528, smiles = "O", InChI = "H2O/h1H2", InChI_key = "XLYOFNOQVPJJNP-UHFFFAOYSA-N", iupac_name = "oxidane", common_name = "water")
```

the package stores the query in `ChemicalIdentifiers.SEARCH_CACHE` as a `Dict{String,Any}`, so subsequent queries on the same (or similar) strings, dont pay the cost of searching in the database.

This package is slower than the original python package, because we don't build a dictionary of all existent synonyms (about 700K synonyms), but perform a sorted search. but it should have better memory footprint for the same reason.



