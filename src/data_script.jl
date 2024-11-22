
"""
    cas_parse(str)

given a CAS number string, returns a tuple of 3 Int32 containing the CAS numbers, it performs no validation on the data

"""
function cas_parse(str)
    a1,a2,a3 = split(str,'-')
    n1 = parse(Int32,a1)
    n2 = parse(Int16,a2)
    n3 = parse(Int16,a3)
    return (n1,n2,n3)
end

function unique_idxs_sorted(vec)
    n = length(vec)
    res = zeros(Int,n)
    k = 1
    res[1] = k
    @inbounds for i = 2:n
        if !isequal(vec[i],vec[i-1])
            k+=1
            res[k] = i
        end
    end
    resize!(res,k)
    return res
end

"""
    load_db!(dbtype::Symbol)

Downloads, processes and stores a database corresponding to the one with the same key stored in DATA_INFO
"""
function load_db!(dbtype::Symbol)
    data = DATA_INFO[dbtype]
    if !isfile(data.textdb)
        @info ":" * string(dbtype) * " database file not found, downloading from " * data.url
        url  = data.url
        fname = data.textdb
        path = Downloads.download(url,fname)
        @info ":" * string(dbtype) *  " database file downloaded."
    end

    path = data.textdb
    if !isfile(data.db)
        @info ":" * string(dbtype)  * " arrow file not generated, processing..."
        arrow_db,arrow_synonym_db,arrow_sort_db =parse_and_write_db!(dbtype)
    else
        arrow_db = Arrow.Table(data.db)
        arrow_synonym_db = Arrow.Table(data.symsdb)
        arrow_sort_db = Arrow.Table(data.sorteddb)

    end
    DATA_DB[dbtype] = (arrow_db,arrow_synonym_db,arrow_sort_db)
    return arrow_db,arrow_synonym_db,arrow_sort_db
end

function parse_and_write_db!(dbtype::Symbol)
    data = DATA_INFO[dbtype]
    path = data.textdb
    i = 0
    for _ in eachline(path)
        i +=1
    end

    pubchemid = zeros(Int64,i)
    CAS = Vector{Tuple{Int32,Int16,Int16}}(undef,i)
    formula = Vector{String}(undef,i)
    MW = Vector{Float64}(undef,i)
    smiles = Vector{String}(undef,i)
    InChI = Vector{String}(undef,i)
    InChI_key = Vector{String}(undef,i)
    iupac_name = Vector{String}(undef,i)
    common_name = Vector{String}(undef,i)
    _synonyms  = Vector{Vector{String}}(undef,i)

    i = 0
    for line in eachline(path)
        i += 1
        strs = line |> z->rstrip(z) |> z->split(z,'\t',limit=10)
        pubchemid[i] = parse(Int64,strs[1])
        CAS[i] = cas_parse(strs[2])
        formula[i] = strs[3]
        MW[i] = parse(Float64,strs[4])
        smiles[i] = strs[5]
        InChI[i] = strs[6]
        InChI_key[i] = strs[7]
        iupac_name[i] = strs[8]
        common_name[i] = strs[9]
        if length(strs) >= 10 #not any synonyms
            sym_i = split(strs[10],('\t',';'))
            push!(sym_i,strs[8])
            _synonyms[i]  = sym_i
        else
            _synonyms[i] = String[strs[8]]
        end
    end

    syms_i = mapreduce(length,+,_synonyms)
    synonyms_list = Vector{String}(undef,syms_i)
    synonyms_index = Vector{Int}(undef,syms_i)

    #for some reason,some empty strings are generated as synonyms.
    #those are eliminated here.

    k = 0
    for (ii,sym_vec) in pairs(_synonyms)

        for (jj,sym) in pairs(sym_vec)
            if !isempty(sym)
                k+=1
                synonyms_list[k] = sym
                synonyms_index[k] = ii
            end
        end
    end

    resize!(synonyms_list,k)
    resize!(synonyms_index,k)

    pubchemid_sort =sortperm(pubchemid)
    CAS_sort = sortperm(CAS)
    formula_sort =sortperm(formula)
    MW_sort =sortperm(MW)
    smiles_sort =sortperm(smiles)
    InChI_sort =sortperm(InChI)
    InChI_key_sort =sortperm(InChI_key)
    iupac_name_sort =sortperm(iupac_name)
    common_name_sort =sortperm(common_name)
    synonyms_sort = sortperm(synonyms_list)

    list = synonyms_list[synonyms_sort]
    index = synonyms_index[synonyms_sort]
    #there is the posibility of repeated elements.
    #this is not present in CalebBell/chemicals because it uses a dict, so equal keys
    #store the same value
    list_unique_idx = unique_idxs_sorted(list)
    list = list[list_unique_idx]
    index = index[list_unique_idx]

    db = (;pubchemid, CAS, formula, MW, smiles, InChI, InChI_key, iupac_name, common_name)
    synonym_db = (;list,index)
    sort_db = (;pubchemid_sort, CAS_sort, formula_sort, MW_sort, smiles_sort, InChI_sort, InChI_key_sort, iupac_name_sort, common_name_sort)

    Arrow.write(data.db,db)
    Arrow.write(data.symsdb,synonym_db)
    Arrow.write(data.sorteddb,sort_db)

    arrow_db = Arrow.Table(data.db)
    arrow_synonym_db = Arrow.Table(data.symsdb)
    arrow_sort_db = Arrow.Table(data.sorteddb)
    return arrow_db,arrow_synonym_db,arrow_sort_db
end

"""
    load_data!(key::Symbol;url=nothing,file=nothing)

Generates and adds to the global DATA_INFO dict a new database. download and process this database with ´load_db(key)´
"""
function load_data!(key::Symbol;url=nothing,file=nothing)
    if url == file == nothing
        throw(ArgumentError("a file or a url must be provided."))
    elseif (url !== nothing) & (file !== nothing)
        throw(ArgumentError("a file or a url must be provided."))
    elseif (url !== nothing)
        fname = joinpath(download_cache, "pubchem_" * string(key))
        farrow = fname * ".arrow"
        fsyms = fname * "_synonyms.arrow"
        sorted = fname * "_sorted.arrow"
        data = (url=url,textdb=fname,db=farrow,symsdb=fsyms,sorteddb=sorted)

    else #file
        fname = file
        url = fname
        farrow = fname * ".arrow"
        fsyms = fname * "_synonyms.arrow"
        sorted = fname * "_sorted.arrow"
        data = (url=url,textdb=fname,db=farrow,symsdb=fsyms,sorteddb=sorted)
    end
    DATA_INFO[key] = data
end
