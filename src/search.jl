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
    pkg_dbs = [:short,:long]
    allkeys = collect(keys(DATA_DB))
    user_dbs = setdiff(allkeys,pkg_dbs)
    return append!(user_dbs,pkg_dbs)
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

function search_chemical_id(ID::AnyQuery;skip_common_name = false,try_strategies = true)
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

function search_chemical_id(ID::CASQuery)
    id = cas(ID)
    key = :cas_not_found
    compound_id = -1
    search_done = false
    
    _keys = db_iteration_order(DATA_DB)

    for key in _keys
        db,sdb,sortdb = DATA_DB[key]
        searchvec = view(db.CAS,sortdb.CAS_sort)
        idx_sort = searchsorted(searchvec,id)
        #@show idx_sort
        if length(idx_sort) == 1 #found and element
            idx = only(sortdb.CAS_sort[idx_sort])
            compound_id =idx
            search_done = true
        end
    #found in db,returning
        if search_done
            return compound_id,key 
        end
    end
    #try looking for cas on the synonims DB
    return search_chemical_id(AnyQuery(value(ID)),skip_common_name=true,try_strategies=false)
end

function search_chemical_id(ID::InChIKeyQuery)
    #symbol: InChI_key
    id = value(ID) 
    key = :inchikey_not_found
    compound_id = -1
    search_done = false
   
    _keys = db_iteration_order(DATA_DB)
    
    for key in _keys
        db,sdb,sortdb = DATA_DB[key]
        searchvec = view(db.InChI_key,sortdb.InChI_key_sort)
        idx_sort = searchsorted(searchvec,id)
        #@show idx_sort
        if length(idx_sort) == 1 #found and element
            idx = only(sortdb.InChI_key_sort[idx_sort])
            compound_id =idx
            search_done = true
        end
    #found in db,returning
        if search_done
            return compound_id,key 
        end
    end
    return compound_id,key
end

function search_chemical_id(ID::PubChemIDQuery)
    #pubchemid
    id = id_num(ID)
    key = :pubchem_not_found
    compound_id = -1
    search_done = false
    
    _keys = db_iteration_order(DATA_DB)

    for key in _keys
        db,sdb,sortdb = DATA_DB[key]
        searchvec = view(db.pubchemid,sortdb.pubchemid_sort)
        idx_sort = searchsorted(searchvec,id)
        #@show idx_sort
        if length(idx_sort) == 1 #found and element
            idx = only(sortdb.pubchemid_sort[idx_sort])
            compound_id =idx
            search_done = true
        end
    #found in db,returning
        if search_done
            return compound_id,key 
        end
    end
    return compound_id,key
end

function search_chemical_id(ID::InChIQuery)
    id = value(ID)
    key = :inchi_not_found
    compound_id = -1
    search_done = false

    _keys = db_iteration_order(DATA_DB)

    for key in _keys
        db,sdb,sortdb = DATA_DB[key]
        searchvec = view(db.InChI,sortdb.InChI_sort)
        idx_sort = searchsorted(searchvec,id)
        #@show idx_sort
        if length(idx_sort) == 1 #found and element
            idx = only(sortdb.InChI_sort[idx_sort])
            compound_id =idx
            search_done = true
        end
    #found in db,returning
        if search_done
            return compound_id,key 
        end
    end
    return compound_id,key
end

function search_chemical_id(ID::SMILESQuery)
    #symbol: InChI_key
    id = value(ID) 
    key = :inchikey_not_found
    compound_id = -1
    search_done = false

    _keys = db_iteration_order(DATA_DB)

    for key in _keys
        db,sdb,sortdb = DATA_DB[key]
        searchvec = view(db.smiles,sortdb.smiles_sort)
        idx_sort = searchsorted(searchvec,id)
        #@show idx_sort
        if length(idx_sort) == 1 #found and element
            idx = only(sortdb.smiles_sort[idx_sort])
            compound_id =idx
            search_done = true
        end
    #found in db,returning
        if search_done
            return compound_id,key 
        end
    end
    return compound_id,key
end
