abstract type AbstractSearchQuery end

struct CASQuery <: AbstractSearchQuery
    ID::String
    cas::Tuple{Int32,Int16,Int16}
end

function CASQuery(str::String)
    cas = cas_parse(str)
    return CASQuery(str,cas)
end

function CASQuery(val::NTuple{3,Integer})
    n1 = val[1]
    n2 = val[2]
    n3 = val[3]
    @assert n1 >=0
    @assert n2 >=0
    @assert n3 >=0
    v1 = convert(Int32,n1)
    v2 = convert(Int16,n2)
    v3 = convert(Int16,n3)
    cas = (v1,v2,v3)
    return CASQuery("",cas)
end

function tuple_to_casstr(n1,n2,n3)

    return string(n1) * "-" * string(n2) * "-" * string(n3)
end

tuple_to_casstr(val::NTuple{3,Integer}) = tuple_to_casstr(val...)
cas(id::CASQuery) = id.cas

struct InChIQuery <: AbstractSearchQuery
    ID::String
end

struct InChIKeyQuery <: AbstractSearchQuery
    ID::String
end

struct PubChemIDQuery <: AbstractSearchQuery
    ID::String
    val::Int
end

function PubChemIDQuery(str::String)
    val = parse(Int,str)
    @assert val >0
    return PubChemIDQuery(str,val)
end

function PubChemIDQuery(val::Integer)
    @assert val >0
    return PubChemIDQuery("",val)
end

id_num(ID::PubChemIDQuery) = ID.val

struct SMILESQuery <: AbstractSearchQuery
    ID::String
end

struct FormulaQuery <: AbstractSearchQuery
    ID::String
end

struct ElementQuery <: AbstractSearchQuery
    ID::String
end

struct AnyQuery <: AbstractSearchQuery
    ID::String
end

struct MissingQuery <: AbstractSearchQuery end

value(id::AbstractSearchQuery)::String = strip(id.ID)


function value(id::SMILESQuery)::String
    id = strip(id.ID)
    return replace(id,"SMILES=" =>"")
end

function value(id::InChIQuery)::String
    id_raw = strip(id.ID)
    id_lower = lowercase(id_raw)
    re1 = r"^inchi=1s/"
    re2 =  r"^inchi=1/"
    t1 = occursin(re1,id_lower)
    t2 = occursin(re2,id_lower)
    if t1
        return chop(id_raw,head=9,tail=0)
    elseif t2
        return chop(id_raw,head=8,tail=0)
    else
        throw("incorrect InChI passed")
    end
end

"""
    is_cas(str)::Bool
check if a given string is a cas number and returns true or false accordingly.

"""
function is_cas(str)
    str = strip(str)
    #regex from https://gist.github.com/KhepryQuixote/00946f2f7dd5f89324d8#file-pypubchemxtractor-py-L22
    cas_regex = r"^[1-9][0-9]{1,6}\\-[0-9]{2}\\-[0-9]$"
    return occursin(cas_regex,str)
end


"""
    is_inchi(str)::Bool
check if a given string is an InChI name and returns true or false accordingly.

"""
function is_inchi(str)
    str = strip(str)
    #regex from https://chemistry.stackexchange.com/a/86892
    inchi_regex = r"^InChI\=1S?\/[^\s]+(\s|$)"
    return occursin(inchi_regex,str)
end

"""
    is_inchikey(str)::Bool
check if a given string is an InChI key and returns true or false accordingly.

"""
function is_inchikey(str)
    str = strip(str)
    #14(A-Z)-10(A-Z)-1(A-Z)
    inchikey_regex  =r"[A-Z]{14}-[A-Z]{10}-[A-Z]{1}"
    return occursin(inchikey_regex,str)
end

"""
    is_element(str)::Bool
check if a given string is an element symbol and returns true or false accordingly.

"""
function is_element(str)
    return false
end

"""
    is_pubchemid(str)::Bool
check if a given string is a PubChem ID and returns true or false accordingly.
A PubChem ID is a positive integer.
"""
function is_pubchemid(str)
    str = strip(str)
    res = tryparse(Int,str)
   if res === nothing
    return false
   elseif res <= 0
    return false
   else
    return true
   end
end

"""
    is_smiles(str)::Bool
check if a given string is an element symbol and returns true or false accordingly.

"""
function is_smiles(str)
    str = strip(str)
    str = uppercase(str)
    smiles_regex  =r"^SMILES="
    return occursin(smiles_regex,str)
end
