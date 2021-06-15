module ChemicalIdentifiers

    const DATA_DB = Dict{Symbol,Any}()   
    const DATA_INFO = Dict{Symbol,Any}()
    
    import Unicode,Downloads,Arrow
    import Scratch, Preferences
    export search_chemical
    
    download_cache = ""    



"""
    function search_chemical(query,cache=cache=ChemicalIdentifiers.SEARCH_CACHE)


Given a query, performs a search on a database of over 70000 common compounds from CalebBell/Chemicals, returning a Named Tuple with the identifiers of the substance in question. 

## Examples
```
julia>using ChemicalIdentifiers
julia>res = search_chemical("water")
(pubchemid = 962, CAS = (7732, 18, 5), formula = "H2O", MW = 18.01528, smiles = "O", InChI = "H2O/h1H2", InChI_key = "XLYOFNOQVPJJNP-UHFFFAOYSA-N", iupac_nam e = "oxidane", common_name = "water")

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

### PubChemID: 
By using any ``<:Integer` (or a string containing an Integer)
```julia
julia> search_chemical(8003)
(pubchemid = 8003, CAS = (109, 66, 0), formula = "C5H12", MW = 72.14878, smiles = 
"CCCCC", InChI = "C5H12/c1-3-5-4-2/h3-5H2,1-2H3", InChI_key = "OFBQJSOFQDEBGM-UHFFFAOYSA-N", iupac_name = "pentane", common_name = "pentane")
```
### CAS registry number:
By using a Tuple of integers or a string with the digits separated by `-` :

```
julia> search_chemical((67,56,1))        
(pubchemid = 887, CAS = (67, 56, 1), formula = "CH4O", MW = 32.04186, smiles = "CO", InChI = "CH4O/c1-2/h2H,1H3", InChI_key = "OKKJLVBELUTLKV-UHFFFAOYSA-N", iupac_name = "methanol", common_name = "methanol")

 search_chemical((67,56,1),nothing) == search_chemical("67-56-1",nothing) #true  
```
### SMILES: 
By using a string starting with `SMILES=` :
```
julia> search_chemical("SMILES=N")       
(pubchemid = 222, CAS = (7664, 41, 7), formula = "H3N", MW = 17.03052, smiles = "N", InChI = "H3N/h1H3", InChI_key = "QGZKDVFQNNGYKY-UHFFFAOYSA-N", iupac_name = "azane", common_name = "ammonia")
```
### InChI:
By using a string starting with `InChI=1/` or `InChI=1S/` :

```
julia> search_chemical("InChI=1/C2H4/c1-2/h1-2H2")     
(pubchemid = 6325, CAS = (74, 85, 1), formula = "C2H4", MW = 28.05316, smiles = "C=C", InChI = "C2H4/c1-2/h1-2H2", InChI_key = "VGGSQFUCUMXWEO-UHFFFAOYSA-N", iupac_name = "ethene", common_name = "ethene") 
```

### InChI key: 
By using a string with the pattern `XXXXXXXXXXXXXX-YYYYYYYYFV-P`:
```
julia> search_chemical("IMROMDMJAWUWLK-UHFFFAOYSA-N")
(pubchemid = 11199, CAS = (9002, 89, 5), 
formula = "C2H4O", MW = 44.05256, smiles 
= "C=CO", InChI = "C2H4O/c1-2-3/h2-3H,1H2", InChI_key = "IMROMDMJAWUWLK-UHFFFAOYSA-N", iupac_name = "ethenol", common_name 
= "ethenol")
```
Searches by CAS and PubChemID are a little bit faster thanks to being encoded as native numeric types, other properties are stores as strings.

The package stores each query in `ChemicalIdentifiers.SEARCH_CACHE` as a `Dict{String,Any}`, so subsequent queries on the same (or similar) strings, dont pay the cost of searching in the database.

If you don't want to store the query, you could use `search_chemical(query,nothing)`, or, if you want your own cache to be used, pass your own cache via `search_chemical(query,mycache)`. 
    
"""
    function search_chemical end
    include("data_script.jl")
    include("search_types.jl")
    include("search.jl")

  function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(search_id_impl,(Int,Symbol))
    Base.precompile(search_id_impl,(String,Symbol))
    Base.precompile(search_id_impl,(Tuple{Int32,Int16,Int16},Symbol))
    Base.precompile(search_chemical_id,(AnyQuery,))
  end

  function __init__()
      global download_cache = Scratch.@get_scratch!("databases")
      

      url_short = "https://github.com/CalebBell/chemicals/raw/master/chemicals/Identifiers/chemical%20identifiers%20pubchem%20small.tsv"
      url_long = "https://github.com/CalebBell/chemicals/raw/master/chemicals/Identifiers/chemical%20identifiers%20pubchem%20large.tsv"
      load_data!(:short,url= url_short)
      load_data!(:long,url = url_long)

      load_db!(:short)
      load_db!(:long)
  end
end

