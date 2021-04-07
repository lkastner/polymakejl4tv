using Oscar
using Polymake
using InvertedIndices

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
    return(map((x,y) -> (x,y), cones, map(coneMultiplicity, cones))
end

function rootConstruction(
    sf::stackyFan,
    rays::Array{Int64,2},
    scalars::Array{Int64,1})

    encoded_rays = mapslices(encode, rays, dims=2)
    for i in length(encoded_rays):
        ray = encoded_rays[i]
        scalar = scalars[i]
        current_stack = sf.stacks[ray]
        sf.stacks[ray] = current_stack * scalar
    end
    return(sf)
    # old_rays = keys(sf.stacks)
    # old_stacks = values(sf.stacks)
    # to_replace = mapslices(encode, rays, dims=2)
    # indices_to_replace = findall(in(to_replace), old_rays)
    #
    # rays_to_keep = old_rays[Not(indices_to_update),:]
    # stacks_to_keep = old_stacks[Not(indices_to_update),:]
    # new_rays = scale(scalars, rays)
    # updated_stacks = map(
    #     (x,y) -> div(x,y),
    #     old_stacks[indices_to_update,:],
    #     scalars)
    #
    # all_rays = decode([rays_to_keep; mapslices(encode, new_rays, dims=2)])
    # all_stacks = [stacks_to_keep; updated_stacks]
    #
    # old_cones = convertIncidenceMatrix(sf.fan.INPUT_CONES)
    # for cone in old_cones:
    #     temp = mapslices(encode, sf.fan.INPUT_RAYS[cone,:], dims=2)
    #     for i in length(cone):
    #         if cone[i] in(indices_to_replace):
    #             cone[i] = findall(x->x==temp[i], to_replace)[1]
    #                 +length(rays_to_keep)
    #
    # return(makeStackyFan(all_rays, old_cones, all_stacks))
end

function stackyBlowup(sf::stackyFan, cone::Array{Int64,1}, ray::Array{Int64,1})
    blowup = toric_blowup(cone, sf.fan, ray)
    sf.stacks[encode(ray)] = 1

    return(stackyFan(blowup, sf.stacks))
end




function getCones(sf::stackyFan)
    formatted = convertIncidenceMatrix(sf.fan.CONES)
    cones = map((x) -> Polymake.polytope.Cone(
        INPUT_RAYS=sf.fan.RAYS[x,:]), formatted)
    return(cones)
end

function encode(objects::Polymake.VectorAllocated{Polymake.Rational})
    return(foldl((x,y) -> string(x, ',', y), objects))
end

function decode(object::Array{String,2})
    return(map((x) -> parse(Int64, x), object))
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
