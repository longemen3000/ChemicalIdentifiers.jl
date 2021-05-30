
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


"""
    load_db!(dbtype::Symbol)
    downloads, processes and stores a database corresponding to the one with the same key stored in DATA_INFO

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

        #@show i
        i = 0
        for line in eachline(path)
            i += 1
            strs = line |> z->rstrip(z,'\n') |> z->split(z,'\t')
            pubchemid[i] = parse(Int64,strs[1])
            CAS[i] = cas_parse(strs[2])
            formula[i] = strs[3]
            MW[i] = parse(Float64,strs[4])
            smiles[i] = strs[5]
            InChI[i] = strs[6]
            InChI_key[i] = strs[7]
            iupac_name[i] = strs[8]
            common_name[i] = strs[9]
            _synonyms[i]  = strs[10:end]
        end

        syms_i = mapreduce(length,+,_synonyms)
        synonyms_list = Vector{String}(undef,syms_i)
        synonyms_index = Vector{Int}(undef,syms_i)
        
        k = 0
        for (ii,sym_vec) in pairs(_synonyms)
            
            for (jj,sym) in pairs(sym_vec)
                k+=1
                synonyms_list[k] = sym
                synonyms_index[k] = ii
            end
        end


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
        db = (;pubchemid, CAS, formula, MW, smiles, InChI, InChI_key, iupac_name, common_name)
        synonym_db = (;list,index)
        sort_db = (;pubchemid_sort, CAS_sort, formula_sort, MW_sort, smiles_sort, InChI_sort, InChI_key_sort, iupac_name_sort, common_name_sort)
        
        Arrow.write(data.db,db)
        Arrow.write(data.symsdb,synonym_db)
        Arrow.write(data.sorteddb,sort_db)

        arrow_db = Arrow.Table(data.db)
        arrow_synonym_db = Arrow.Table(data.symsdb)
        arrow_sort_db = Arrow.Table(data.sorteddb)
    else 
        arrow_db = Arrow.Table(data.db)
        arrow_synonym_db = Arrow.Table(data.symsdb)
        arrow_sort_db = Arrow.Table(data.sorteddb)

    end
    DATA_DB[dbtype] = (arrow_db,arrow_synonym_db,arrow_sort_db) 
    return arrow_db,arrow_synonym_db,arrow_sort_db
end

"""
    load_data!(key::Symbol;url=nothing,file=nothing)

generates and adds to the global DATA_INFO dict a new database. download and process this database with ´load_db(key)´
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





#3 * x^0.7 - 2 * x + 1
#a,b,c,a
#3,2,1,0.7
#=
ax = bx-c
c = (b-a)x
=#