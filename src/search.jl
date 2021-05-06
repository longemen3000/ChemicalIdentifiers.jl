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

function build_result(idx,key,symbol)
    db,sdb = DATA_DB[key]
    col= db[symbol]
    return col[idx]
end

function detect_query(id::String)
    if is_element(id)
        return ElementQuery(id)
    elseif is_cas(id)
        return CASQuery(id)
    elseif is_inchi(id)
        return InChIQuery(id)
    elseif is_inchikey(id)
        return InChIKeyQuery(id)
    else
        return AnyQuery(id)
    end
end

function search_chemical(ID::String,cache=SEARCH_CACHE)
    if cache !== nothing
        #the base is that two chemicals with different casing should be the same.
        normalized_id = Unicode.normalize(ID,casefold = true,stripmark=true)
        if haskey(cache,normalized_id)
            return cache[normalized_id]
        end
        get!(cache,ID) do 
            compound_id,key = search_chemical_id(detect_query(ID))
            build_result(compound_id,key)      #return db.common_name[compound_id]
        end
        cache[normalized_id] = cache[ID]
    else
        compound_id,key = search_chemical_id(detect_query(ID))
        build_result(compound_id,key)      #return db.common_name[compound_id]
    end
end

function search_chemical_id(ID::AnyQuery;skip_common_name = false,try_strategies = true)
    #skip_common_name skips search on the common_name col
    key = :not_found
    compound_id = -1
    id = value(ID)
    search_done = false
    #generating keys iteration order: user, short database, long database
    pkg_dbs = [:short,:long]
    allkeys = collect(keys(DATA_DB))
    user_dbs = setdiff(allkeys,pkg_dbs)
    _keys = append!(user_dbs,pkg_dbs)
    for key in _keys
    db,sdb = DATA_DB[key]
    
        if !skip_common_name
            idx_from_name = findfirst(isequal(id),db.common_name)
            if idx_from_name !== nothing
                compound_id = only(idx_from_name)
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
    id = cas_parse(value(ID))  
    key = :cas_not_found
    compound_id = -1
    #generating keys iteration order: user, short database, long database
    search_done = false
    pkg_dbs = [:short,:long]
    allkeys = collect(keys(DATA_DB))
    user_dbs = setdiff(allkeys,pkg_dbs)
    _keys = append!(user_dbs,pkg_dbs)
    for key in _keys
    db,sdb = DATA_DB[key]
        idx_from_prop = findall(isequal(id),db.CAS)
        if length(idx_from_prop) == 1 #should be unique
        compound_id = only(idx_from_prop)
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
    id = value(ID) 
    key = :inchikey_not_found
    compound_id = -1
    search_done = false
    #generating keys iteration order: user, short database, long database
    pkg_dbs = [:short,:long]
    allkeys = collect(keys(DATA_DB))
    user_dbs = setdiff(allkeys,pkg_dbs)
    _keys = append!(user_dbs,pkg_dbs)
    for key in _keys
        db,sdb = DATA_DB[key]
        idx_from_prop = findall(isequal(id),db.InChI_key)
        if length(idx_from_prop) == 1 #should be unique
        compound_id = only(idx_from_prop)
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
    id = parse(int64,value(ID)) 
    key = :pubchem_not_found
    compound_id = -1
    search_done = false
    #generating keys iteration order: user, short database, long database
    pkg_dbs = [:short,:long]
    allkeys = collect(keys(DATA_DB))
    user_dbs = setdiff(allkeys,pkg_dbs)
    _keys = append!(user_dbs,pkg_dbs)
    for key in _keys
        db,sdb = DATA_DB[key]
        idx_from_prop = findall(isequal(id),db.pubchemid)
        if length(idx_from_prop) == 1 #should be unique
        compound_id = only(idx_from_prop)
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
    id_lower = lowercase(value(ID))
    re1 = r"^inchi=1s/"
    re2 =  r"^inchi=1/"
    t1 = occursin(re1,id_lower)
    t2 = occursin(re2,id_lower)
    if t1
        id = chop(value(ID),head=9,tail=0)
    elseif t2
        id = chop(value(ID),head=8,tail=0)
    else
        throw("incorrect InChI passed")
    end

    key = :inchi_not_found
    compound_id = -1
    search_done = false
    #generating keys iteration order: user, short database, long database
    pkg_dbs = [:short,:long]
    allkeys = collect(keys(DATA_DB))
    user_dbs = setdiff(allkeys,pkg_dbs)
    _keys = append!(user_dbs,pkg_dbs)
        for key in _keys
            db,sdb = DATA_DB[key]
        idx_from_prop = findall(isequal(id),db.pubchemid)
        if length(idx_from_prop) == 1 #should be unique
        compound_id = only(idx_from_prop)
        search_done = true
        end 
    #found in db,returning
        if search_done
            return compound_id,key 
        end
    end
    return compound_id,key
end