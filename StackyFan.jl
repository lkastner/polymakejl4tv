using Oscar
using Polymake
using InvertedIndices
using Combinatorics

struct StackyFan
    fan::Polymake.BigObjectAllocated
    scalars::Array{Int64, 1}
    stacks::Dict{String, Int64}
    # Constructors for the StackyFan object
    StackyFan(
        fan::Polymake.BigObjectAllocated,
        scalars::Array{Int64, 1},
        stacks::Dict{String, Int64}) = new(fan, scalars, stacks)
    StackyFan(
        rays::Array{Int64, 2},
        cones::Array{Array{Int64, 1}, 1},
        scalars::Array{Int64, 1}) = makeStackyFan(rays, cones, scalars)
end

## Helper functions

function makeStackyFan(
    rays::Array{Int64,2},
    cones::Array{Array{Int64,1},1},
    scalars::Array{Int64,1})

    fan = fulton.NormalToricVariety(INPUT_RAYS=rays, INPUT_CONES=cones)
    stack_rays = mapslices(encode, fan.RAYS, dims=2)
    pairs = map((x,y) -> (x,y), stack_rays, scalars)
    stacks = Dict(pairs)

    return(StackyFan(fan, scalars, stacks))
end

function encode(objects::Polymake.VectorAllocated{Polymake.Rational})
    return(foldl((x,y) -> string(x, ',', y), objects))
end

function decode(object::Array{String,2})
    return(map((x) -> parse(Int64, x), object))
end

## API functions

"""
    getRayStack(::StackyFan, ::Array{Int64, 1})

    Get the scalar associated with a ray in the given stacky fan structure

# Examples
"""
function getRayStack(sf::StackyFan, ray::Array{Int64, 1})
    return sf.stacks[encode(ray)]
end

function getMultiplicities(sf::StackyFan)
    cones = getCones(sf)
    return(map((x,y) -> (x,y), cones, map(coneMultiplicity, cones)))
end

function rootConstruction(
    sf::StackyFan,
    rays::Array{Int64,2},
    scalars::Array{Int64,1})

    encoded_rays = mapslices(encode, rays, dims=2)
    for i in length(encoded_rays)
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

function stackyBlowup(sf::StackyFan, cone::Array{Int64,1}, ray::Array{Int64,1})
    blowup = toric_blowup(cone, sf.fan, ray)
    sf.stacks[encode(ray)] = 1

    return(StackyFan(blowup, sf.stacks))
end



function getCones(sf::StackyFan)
    formatted = convertIncidenceMatrix(sf.fan.CONES)
    cones = map((x) -> Polymake.polytope.Cone(
        INPUT_RAYS=sf.fan.RAYS[x,:]), formatted)
    return(cones)
end


"""
    slicematrix(::AbstractMatrix{<:Number})

    Take a two-dimensional matrix and output a list of its row vectors.

# Examples
```jldoctest makeSmoothWithDependencies
julia> A=[1 2; 3 4]

julia> slicematrix(A)
[[ 1 ,  2 ], [ 3 ,  4 ]]
"""
function slicematrix(A::AbstractMatrix{<:Number})
    return [A[i, :] for i in 1:size(A,1)]
end

"""
    rowMinors(::AbstractMatrix{<:Number},::Union{AbstractSet,AbstractVector})

    Identical to slicematrix, except only returns row vectors indexed by a set S.

# Examples
```jldoctest makeSmoothWithDependencies
julia> A=[1 2 3;4 5 6; 7 8 9]

julia> S=Set([1,3])

julia> rowMinors(A,S)
2×3 LinearAlgebra.Transpose{Int64,Array{Int64,2}}:
 1  2  3
 7  8  9
"""
function rowMinors(A::AbstractMatrix{<:Number},S::Union{AbstractSet,AbstractVector})
    outList=[]
    slices=slicematrix(A)
    for i in 1:size(slices,1)
        if i in S
            append!(outList,[slices[i]])
        end
    end
    return transpose(hcat(outList...))
end

"""
    convertIncidenceMatrix(::Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric})

    Takes a Polymake incidence matrix (e.g., the output of X.MAXIMAL_CONES for a toric variety X) and outputs a list of vectors,
    with each vector recording the indices marked on a given row of the incidence matrix.

# Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]])

julia> M=X.MAXIMAL_CONES

julia> convertIncidenceMatrix(M)
[[ 1 ,  2 ], [ 2 ,  3 ]]
"""
function convertIncidenceMatrix(A::Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric})
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
    return convert.(Array{Int64, 1}, out)
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

function coneConvert(cone,rayMatrix)
    coneRays=rowMinors(rayMatrix,cone)
    C=Polymake.polytope.Cone(RAYS=coneRays)
    return C
end
    
    

function getCones(X)
    lattice=X.HASSE_DIAGRAM
    faces=@Polymake.convert_to Array{Set{Int}} lattice.FACES
    out=[]
    for i in 2:(size(faces,1)-1)
        push!(out,Array(@Polymake.convert_to Array{Int} faces[i]))
    end
    return out
end

function distinguishedAndMultiplicity(cone,rayMatrix,dist)
    l=size(rayMatrix,1)
    if dot(convertToIncidence(cone,l),dist) > 0
        C=coneConvert(cone,rayMatrix)
        mult=coneMultiplicity(C)
        if mult > 1
            return true
        else
            return false
        end
    else
        return false
    end
end

function convertToIncidence(v,l)
    out=[]
    for j in 1:l
        if j in v
            append!(out,1)
        else
            append!(out,0)
        end
    end
    return out
end

function compareCones(cone1, cone2, rayMatrix, distinguished)
    l=size(rayMatrix,1)
    c1=convertToIncidence(cone1,l)
    c2=convertToIncidence(cone2,l)
    # Calculate the number of non-distinguished rays
    nondist1 = size(cone1,1) - dot(c1, distinguished)
    nondist2 = size(cone2,1) - dot(c2, distinguished)
    if (nondist1 - nondist2 != 0)
        return nondist1 - nondist2
    else
        # Need to use the method for calculating multiplicity of cone
        mult1 = coneMultiplicity(coneConvert(cone1,rayMatrix))
        mult2 = coneMultiplicity(coneConvert(cone2,rayMatrix))
        return mult1 - mult2
    end
end

function extremalCones(S, rayMatrix, distinguished)
    # Finds the extremal cones according to # distinguished rays and multiplicity
    # distinguished is a boolean vector whose size is equal to the number of rays
    # The i-th index is 1 if the i-th ray (in rayMatrix) is distinguished
    maxCones = [S[1]]
    for i in 2:size(S,1)
        cone = S[i]
        # Compare the cone with the first element of the maximal cone list
        comp = compareCones(cone, maxCones[1], rayMatrix, distinguished)
        if (comp > 0)
            maxCones = [cone]
        elseif (comp == 0)
            push!(maxCones, cone)
        end
    end
    return maxCones
end

function interiorPoints(C)
    rayMatrix=Array(Polymake.common.primitive(C.RAYS))
    l=size(rayMatrix,1)
    dim=size(rayMatrix,2)
    if rank(rayMatrix)<l
        error("Input cone is not simplicial.")
    end
    subsets=collect(powerset([1:l;]))
    vertices=[]
    for elt in subsets
        vert=zeros(Polymake.Rational,1,dim)
        for i in 1:l
            if i in elt
                vert+=rayMatrix[[i],:]
            end
        end
        append!(vertices,[vert])
    end
    V=vcat(vertices...)
    VH=hcat(ones(Polymake.Rational,size(V,1)),V)
    P=Polymake.polytope.Polytope(POINTS=VH)
    #print(P.POINTS)
    if size(P.INTERIOR_LATTICE_POINTS,1)==0
        return nothing
    end
    intPoints=Array(P.INTERIOR_LATTICE_POINTS)[:,2:(dim+1)]
    validPoints=[]
    #return intPoints
    for i in 1:size(intPoints,1)
        point=intPoints[i,:]
        if gcd(point)==1
            append!(validPoints,[point])
        end
    end
    return Array(vcat(validPoints...))
end
