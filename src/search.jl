const SEARCH_CACHE = Dict{String,Any}()

function build_result(idx,key)
    if idx == -1
        return missing
    end
    db,sdb = DATA_DB[key]
    pubchemid= db.pubchemid[idx]
    CAS= db.CAS[idx]
    formula= db.formula[idx]
    MW= db.MW[idx]
    smiles= db.smiles[idx]
    InChI= db.InChI[idx]
    InChI_key= db.InChI_key[idx]
    iupac_name= db.iupac_name[idx]
    common_name= db.common_name[idx]
    return (;pubchemid, CAS, formula, MW, smiles, InChI, InChI_key, iupac_name, common_name)
end

function detect_query(id::String)
    if is_element(id)
        return ElementQuery(id)
    elseif is_cas(id)
        return CASQuery(id)
    elseif is_pubchemid(id)
        return PubChemIDQuery(id)
    elseif is_inchi(id)
        return InChIQuery(id)
    elseif is_inchikey(id)
        return InChIKeyQuery(id)
    elseif is_smiles(id)
        return SMILESQuery(id)
    else
        return AnyQuery(id)
    end
end

function detect_query(id::Int)
        return PubChemIDQuery(id)
end

function detect_query(id::Tuple{Integer,Integer,Integer})
    return CASQuery(id)
end

function db_iteration_order(DB)
    #generating keys iteration order: user, short database, long database
    dbnames = [:short,:long]
    for k in keys(DATA_DB)
        if !(k in (:short,:long))
            push!(dbnames,k)
        end
    end
    return dbnames
end

function search_chemical(query,cache=SEARCH_CACHE)
    if cache !== nothing
        #the base is that two chemicals with different casing should be the same.
        if query isa NTuple{3,Integer}
            ID = tuple_to_casstr(query)
        else
            ID = string(query)
        end
        
        if haskey(cache,ID)
            return cache[ID]
        end
        
        normalized_id = Unicode.normalize(ID,casefold = true,stripmark=true)
        if haskey(cache,normalized_id)
            return cache[normalized_id]
        end
        
        compound_id,key = search_chemical_id(detect_query(query))
        res = build_result(compound_id,key)      #return db.common_name[compound_id]
        cache[ID] = res
        cache[normalized_id] = res
    else
        compound_id,key = search_chemical_id(detect_query(query))
        return build_result(compound_id,key)      #return db.common_name[compound_id]
    end
end

function search_chemical_id(ID::AnyQuery;skip_common_name = false,try_strategies = true)::Tuple{Int,Symbol}
    #skip_common_name skips search on the common_name col
    key = :not_found
    compound_id = -1
    id = value(ID)
    search_done = false
    _keys = db_iteration_order(DATA_DB)
    for key in _keys
    db,sdb,sortdb = DATA_DB[key]
        if !skip_common_name
            searchvec = view(db.common_name,sortdb.common_name_sort)
            idx_sort = searchsorted(searchvec,id)
            #@show idx_sort
            if length(idx_sort) == 1 #found and element
                idx = only(sortdb.common_name_sort[idx_sort])
                compound_id =idx
                search_done = true
            end
        end
        if !search_done
            idx = searchsorted(sdb.list,id)
            if length(idx) == 1 #found and element
                compound_id = only(sdb.index[idx])
                search_done = true
            end
        end
        #found in db,returning
        if search_done
            return compound_id,key 
        end
    end

    if !search_done
        if !try_strategies #bail out here if requested
            return -1,:not_found
        end
    end
    #result not found, trying same strategies as present in CalebBell/Chemicals
    #strategy 1: trying without spaces and dashs
    _ids = Vector{String}(undef,5)
    _ids[1] = Unicode.normalize(id,casefold = true,stripmark=true)
    _ids[2] = replace(_ids[1]," "=>"")
    _ids[3] = replace(_ids[2],"-"=>"")
    _ids[4] = replace(id," "=>"")
    _ids[5] = replace(_ids[4],"-"=>"")

    _ids = unique!(_ids)
    _ids = setdiff!(_ids,[id])

    for _id in _ids
        compound_id,key = search_chemical_id(AnyQuery(_id))
        if compound_id !== -1
            search_done = true
            break
        end
    end
    #strategy 2: trying to match in the form 'water (H2O)'
    if !search_done 
        re = r"\w+\s+\([\s\w]+\)"     
        if occursin(re,id)
            _id = id |> z->replace(z,")"=>"") |> z->split(z,"(") .|> strip
            id1,id2 = first(id1),last(id2)
            compound_id1,key1 = search_chemical_id(AnyQuery(id1))
            compound_id2,key2 = search_chemical_id(AnyQuery(id2))
            if (compound_id1 == compound_id2) & (key1==key2)
                search_done = true
                compound_id = compound_id2
                key = key2
            end
        end
    end

    #if something worked, return here, else, return not found
    if search_done
       return compound_id,key
    else
        return -1,:not_found
    end
end

function search_chemical_id(ID::CASQuery)::Tuple{Int,Symbol}
    id = cas(ID)  
    compound_id,key = search_id_impl(id,:CAS)
    if compound_id != -1 
        return compound_id,key
    end
    return search_chemical_id(AnyQuery(value(ID)),skip_common_name=true,try_strategies=false)
end

function search_chemical_id(ID::InChIKeyQuery)::Tuple{Int,Symbol}
    id = value(ID) 
    return search_id_impl(id,:InChI_key)
end

function search_chemical_id(ID::PubChemIDQuery)::Tuple{Int,Symbol}
    id = id_num(ID)
    return search_id_impl(id,:pubchemid)
end

function search_chemical_id(ID::InChIQuery)::Tuple{Int,Symbol}
    id = value(ID)
    search_id_impl(id,:InChI)
end

function search_chemical_id(ID::SMILESQuery)::Tuple{Int,Symbol}
    id = value(ID) 
    search_id_impl(id,:smiles)
end

arrowtype(::Type{Int}) = Arrow.Primitive{Int,Vector{Int}}
arrowtype(::Type{String}) = Arrow.List{String, Int32, Vector{UInt8}}
arrowtype(::Type{Tuple{Int32,Int16,Int16}}) = Arrow.Struct{Tuple{Int32, Int16, Int16}, Tuple{Arrow.Primitive{Int32, Vector{Int32}}, Arrow.Primitive{Int16, Vector{Int16}}, Arrow.Primitive{Int16, Vector{Int16}}}}

function search_id_impl(id::T,k::Symbol)::Tuple{Int,Symbol} where {T}
    arrowT = arrowtype(T)
    return search_id_impl(id,k,arrowT)
end

function search_id_impl(id::T,sym::Symbol,::Type{A})::Tuple{Int,Symbol} where {T,A}
    compound_id::Int = -1
    search_done = false
    dbnames = db_iteration_order(DATA_DB)::Vector{Symbol}
    for dbname in dbnames
        db,_,sortdb = DATA_DB[dbname]
        dbcol = getproperty(db,sym)::A
        sort_sym = Symbol(sym,:_sort)
        dbidx = getproperty(sortdb,sort_sym)::Arrow.Primitive{Int,Vector{Int}}
        searchvec = view(dbcol,dbidx)
        idxs = searchsorted(searchvec,id)
        #@show idx_sort
        if length(idxs) == 1 #found and element
            compound_id = only(dbidx[idxs])::Int
            search_done = true
        else
            throw("Search is not unique, multiple matches found for $id in database $dbname, on the $sym column")
        end        
        if search_done
            return compound_id,dbname 
        end
    end
    return compound_id,sym #to identify where it fails
end

