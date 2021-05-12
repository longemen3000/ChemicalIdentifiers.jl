module ChemicalIdentifiers

    const DATA_DB = Dict{Symbol,Any}()   
    const DATA_INFO = Dict{Symbol,Any}()
    
    import Unicode,Downloads,Arrow,ChemEquations
    import Scratch, Preferences
    export search_chemical
    
    download_cache = ""    

    include("data_script.jl")
    include("search_types.jl")
    include("search.jl")


    #const FAST_INDEX = Preferences.@load_preference("fast_index","false")
    const FAST_INDEX = false
 #=
    function fast_index!(val::Bool) 
        Preferences.@set_preferences!("fast_index" => string(val))
        val = Preferences.@load_preference("fast_index",false)
        if val
            @info("Fast index set; restart your Julia session for this change to take effect.")
        else
            @info("fast index unset; restart your Julia session for this change to take effect.")
        end
    end
=#
#=
    function fast_index()
    return Preferences.@load_preference("fast_index",false)
    end
    =#

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

