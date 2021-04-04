using Oscar
using Polymake

struct stackyFan
    fan::Polymake.BigObjectAllocated
    stacks::Dict{String, Int64}
end

function makeStackyFan(
    rays::Array{Int64,2},
    cones::Array{Array{Int64,1},1},
    stacks::Array{Int64,1})

    fan = fulton.NormalToricVariety(INPUT_RAYS=rays, INPUT_CONES=cones)
    stack_rays = mapslices(encode, fan.RAYS, dims=2)
    pairs = map((x,y) -> (x,y), stack_rays, stacks)
    stacks = Dict(pairs)

    return(stackyFan(fan, stacks))
end

function getMultiplicities(sf::stackyFan)
    cones = getCones(sf)
    return(map(coneMultiplicity, cones))
end

function getBoxPoints(sf::stackyFan)

end

function getBoxPoints(cone::Polymake.BigObjectAllocated)
    gcds = mapslices((x) -> foldl(gcd, x), cone.RAYS, dims=2)
    generators = mapslices(
        (x,y) -> div(x, repeat([y], length(x))),
        cone.RAYS,
        gcds,
        dim=2)

end

function encode(objects::Polymake.VectorAllocated{Polymake.Rational})
    return(foldl((x,y) -> string(x, ',', y), objects))
end

function ugly_gcd(first::Polymake.RationalAllocated, second::Polymake.RationalAllocated)
    
end

function getCones(sf::stackyFan)
    formatted = convertIncidenceMatrix(sf.fan.MAXIMAL_CONES)
    cones = map((x) -> Polymake.polytope.Cone(
        INPUT_RAYS=sf.fan.RAYS[x,:]), formatted)
    return(cones)
end

function convertIncidenceMatrix(A)
    A=Array(A)
    dim1=size(A,1)
    dim2=size(A,2)
    out=[]
    for i in 1:dim1
        members=[]
        for j in 1:dim2
            if A[i,j]==true
                append!(members,j)
            end
        end
        append!(out,[members])
    end
    return out
end

function coneMultiplicity(C)
    A=Polymake.common.primitive(C.RAYS)
    M=matrix(ZZ,[fmpz.(y) for y in A])
    SNF=Nemo.snf(M)
    mult=1
    for i in 1:size(SNF,1)
        mult*=SNF[i,i]
    end
    return mult
end
