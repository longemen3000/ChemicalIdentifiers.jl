
const RAW_DATA_DIR = ["",""]
const ARROW_DATA_DIR = ["","","",""]
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
pub_data_download()
download a fresh copy of the chemical identifiers database located at
https://github.com/CalebBell/chemicals. gives a named tuple of all columns in the database
"""


function pub_data_download(dbtype=:short)
    key = 0
    if dbtype == :long
        key = 2
        
        if !isfile(RAW_DATA_DIR[2])
            @info "long database file not found, downloading..."
            url = "https://github.com/CalebBell/chemicals/raw/master/chemicals/Identifiers/chemical%20identifiers%20pubchem%20large.tsv"
            fname = joinpath(download_cache, "pubchem_long")
            path = Downloads.download(url,fname)
            @info "long database file downloaded."
            RAW_DATA_DIR[key] = path
        end
    elseif dbtype ==:short
        key = 1
        if !isfile(RAW_DATA_DIR[2])
            @info "short database file not found, downloading..."
            url  = "https://github.com/CalebBell/chemicals/raw/master/chemicals/Identifiers/chemical%20identifiers%20pubchem%20small.tsv"
            fname = joinpath(download_cache, "pubchem_short")    
            path = Downloads.download(url,fname)
            @info "short database file downloaded."
            RAW_DATA_DIR[key] = path
        end
    else
        throw(error("incorrect symbol, only :short and :large are accepted."))
    end
    path = RAW_DATA_DIR[key]
    
    if !isfile(ARROW_DATA_DIR[key])
        @info String(dbtype) * " arrow file not generated, processing..."
        i = 0
        for line in eachline(path)
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
        synonyms  = Vector{Vector{Int}}(undef,i)

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
            synonyms[i] = zeros(Int,length(_synonyms[i]))
        end

        syms_i = mapreduce(length,+,synonyms)
        synonyms_list = Vector{String}(undef,syms_i)
        synonyms_index = Vector{Int}(undef,syms_i)
        
        k = 0
        #@show _synonyms[1]
        for (ii,sym_vec) in pairs(_synonyms)
            
            for (jj,sym) in pairs(sym_vec)
                k+=1
                synonyms_list[k] = sym
                synonyms_index[k] = ii
                synonyms[ii][jj] = k
            end
        end
        sort_order = sortperm(synonyms_list)
        
        list = synonyms_list[sort_order]
        index = synonyms_index[sort_order]
        #=
        for (ii,sym_vec) in pairs(synonyms)
            for (jj,sym) in pairs(sym_vec)
                synonyms[ii][jj] = sort_order[synonyms[ii][jj]]
            end
        end

        list = Vector{String}(undef,length(synonyms_index))
        index = Vector{Int}(undef,length(synonyms_index))
        for kk in 1:length(index)
            list[kk] = synonyms_list[sort_order[kk]]
            index[kk] = synonyms_index[sort_order[kk]]
        end
        @show index[1]
        @show synonyms_index[1]
        =#
        db = (;pubchemid, CAS, formula, MW, smiles, InChI, InChI_key, iupac_name, common_name)
        synonym_db = (list=list,index=index)
        io1 = IOBuffer()
        io2 = IOBuffer()
        Arrow.write(ARROW_DATA_DIR[key],db)
        Arrow.write(ARROW_DATA_DIR[key+2],synonym_db)
        db = Arrow.Table(ARROW_DATA_DIR[key])
        sdb = Arrow.Table(ARROW_DATA_DIR[key+2])

    else
        db = Arrow.Table(ARROW_DATA_DIR[key])
        sdb = Arrow.Table(ARROW_DATA_DIR[key+2])
    end
    DATA_DB[dbtype] = (db,sdb) 
    return db,sdb
end




#3 * x^0.7 - 2 * x + 1
#a,b,c,a
#3,2,1,0.7
#=
ax = bx-c
c = (b-a)x
=#