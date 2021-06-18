include("StackyFan.jl")

function remove!(a, item)
    return deleteat!(a, findall(x->x==item, a))
end
    
function getIndex(ray::Array{Int64,1},rayMatrix::Array{Int64,2})
    slice=slicematrix(rayMatrix)
    index=findall(x->x==ray,slice)
    return index[1]
end
    
"""
    isIndependent(::Int64,::Array{Int64,1},::Array{Int64,2})

    Takes a ray matrix, a list of indices representing a cone, and an index represeting a ray of that cone. Determines whether the given ray is independent in the cone (i.e. does not contribute to the multiplicity of the cone).

# Examples
```jldoctest AlgC
julia> isIndependent(3,[1,2,3],[1 0 0; 0 1 0; 1 2 3])
false

julia> isIndependent(3,[1,2,3],[1 0 0; 0 1 0; 1 1 1])
true
```
"""
function isIndependent(rayIndex::Int64,cone::Array{Int64,1},rayMatrix::Array{Int64,2})
    scone=copy(cone)
    subcone=remove!(scone,rayIndex)
    mult=getMultiplicity(cone,rayMatrix)
    submult=getMultiplicity(subcone,rayMatrix)
    if mult==submult
        return true
    else
        return false
    end
end
    
"""
    independencyIndex(::Array{Int64,1},::Array{Int64,2})

    Returns the number of non-independent rays in a cone. Input in indices-ray matrix format.

# Examples
```jldoctest
julia> independencyIndex([1,2,3],[1 0 0 ; 1 2 0; 2 0 3; 0 0 5])
2
```
"""
function independencyIndex(cone::Array{Int64,1},rayMatrix::Array{Int64,2})
    index=0
    for elt in cone
        if isIndependent(elt,cone,rayMatrix)==false
            index+=1
        end
    end
    return index
end
    
"""
    isRelevant(::Array{Int64,1},::Array{Int64,1},::StackyFan)

    Determines if the given ray, relative to the given cone, either has a stacky value greater than 1 or is not independent.

# Examples
```jldoctest AlgC
julia> F=makeStackyFan([1 0 0; 1 2 0; 0 0 1],[[0,1,2]],[1,1,2]);

julia> isRelevant([1,2,0],[1,2,3],F)
true

julia> F=makeStackyFan([1 0 0; 0 1 0; 0 0 1],[[0,1,2]],[1,1,2]);

julia> isRelevant([0,1,0],[1,2,3],F)
false
"""
function isRelevant(ray::Array{Int64,1},cone::Array{Int64,1},F::StackyFan)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS)))
    rayStack=F.stacks[encode(ray)]
    rayIndex=getIndex(ray,rayMatrix)
    rayIndependent=isIndependent(rayIndex,cone,rayMatrix)
    if rayStack != 1 || rayIndependent == false
        return true
    else
        return false
    end
end

function toroidalIndex(cone::Array{Int64,1},F::StackyFan,div::Array{Int64,1})
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(F.fan.RAYS)))
    s=count(x->div[x]==1,cone)
    slice=slicematrix(rayMatrix)
    flipt=0
    for i in cone
        if div[i]==0
            if isRelevant(slice[i],cone,F)==false
                flipt+=1
            end
        end
    end
    t=size(cone,1)-flipt
    return t-s
