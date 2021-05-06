module ChemicalIdentifiers

    const DATA_DB = Dict{Symbol,Any}()   
    import Unicode,Downloads,Arrow,Scratch 
    #using ,Tables
    export search_chemical

    
    download_cache = ""    

    include("data_script.jl")
    include("search_types.jl")
    include("search.jl")


    function __init__()
        global download_cache = Scratch.@get_scratch!("databases")
        url_long = "https://github.com/CalebBell/chemicals/raw/master/chemicals/Identifiers/chemical%20identifiers%20pubchem%20large.tsv"
        url_short = "https://github.com/CalebBell/chemicals/raw/master/chemicals/Identifiers/chemical%20identifiers%20pubchem%20small.tsv"
        fname_long = joinpath(download_cache, "pubchem_long")
        fname_short = joinpath(download_cache, "pubchem_short")
        RAW_DATA_DIR[1] = fname_short
        RAW_DATA_DIR[2] = fname_long
        ARROW_DATA_DIR[1] = fname_short * ".arrow"
        ARROW_DATA_DIR[2] = fname_long * ".arrow"
        ARROW_DATA_DIR[3] = fname_short * "_synonyms.arrow"
        ARROW_DATA_DIR[4] = fname_long * "_synonyms.arrow"
        pub_data_download(:short)
        pub_data_download(:long)
    end
end

