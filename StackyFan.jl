using Oscar
using Polymake
using InvertedIndices
using Combinatorics
using LinearAlgebra

"""

    Structure to store information of a stacky fan - this is a fan together with a dictionary assigning stacky values to each ray.

# Properties:
-> `fan` - the underlying fan, as a polymake object
-> `scalars` - the array of stacky values
-> `stacks` - a dictionary assigning stacky values to each ray.

"""
struct StackyFan
    fan::Polymake.BigObjectAllocated
    stacks::Dict{String, Int64}
    # Constructors for the StackyFan object
    StackyFan(
        fan::Polymake.BigObjectAllocated,
        stacks::Dict{String, Int64}) = new(fan, stacks)
    StackyFan(
        rays::Array{Int64, 2},
        cones::Array{Array{Int64, 1}, 1},
        scalars::Array{Int64, 1}) = makeStackyFan(rays, cones, scalars)
    StackyFan(
        fan::Polymake.BigObjectAllocated,
        scalars::Array{Int64, 1}) = addStackStructure(fan, scalars)
end

## Helper functions

"""
    makeStackyFan(::Array{Int64,2},::Array{Array{Int64,1},1},::Array{Int64,1}))

Function to generate a stacky fan from a matrix representing rays as row vectors, a vector of vectors representing the rays contained in each cone, and a vector of stacky values to be assigned the rays.
"""
function makeStackyFan(
    rays::Array{<:Number,2},
    cones::Array{Array{Int64,1},1},
    scalars::Array{Int64,1})

    # Construct a normal fan from the given rays and cones
    fan = fulton.NormalToricVariety(INPUT_RAYS=rays, INPUT_CONES=cones)
    
    # Construct the dictionary
    stack_rays = mapslices(encode, fan.RAYS, dims=2)
    pairs = map((x,y) -> (x,y), stack_rays, scalars)
    stacks = Dict(pairs)

    return(StackyFan(fan, stacks))
end

"""
    addStackStructure(::Polymake.BigObjectAllocated, ::Array{Int64, 1})

    Function to generate a stacky fan from a given fan and a set of scalars.
"""
function addStackStructure(
    fan::Polymake.BigObjectAllocated,
    scalars::Array{Int64, 1})
    
    # Construct the dictionary
    stack_rays = mapslices(encode, Polymake.common.primitive(fan.RAYS), dims=2)
    pairs = map((x,y) -> (x,y), stack_rays, scalars)
    stacks = Dict(pairs)

    return(StackyFan(fan, stacks))
end

"""
    encode(::Polymake.VectorAllocated{Polymake.Rational})

    Internal function that converts a Polymake vector, representing a ray in the fan,
to a string in order to allow for hashing for the dictionary.
"""
function encode(objects::Polymake.VectorAllocated{Polymake.Integer})
    return(foldl((x,y) -> string(x, ',', y), objects))
end

function encode(objects::Vector{Int64})
    return(foldl((x,y) -> string(x, ',', y), objects))
end

function encode(objects::Polymake.VectorAllocated{Polymake.Rational})
    return(foldl((x,y) -> string(x, ',', y), objects))
end

"""
    stackyWeights(::StackyFan)

Returns a list of the stacky weights of the rays of the given stacky fan with the same order as the rays of the fan.
"""
function stackyWeights(sf::StackyFan)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(sf.fan.RAYS)))
    #println("Stacky weights ray matrix: $rayMatrix")
    rayList=slicematrix(rayMatrix)
    out=Int64[]
    for ray in rayList
        try
            stack=sf.stacks[encode(ray)]
            push!(out,stack)
        catch e
            println("ERROR!!!")
            println(e)
            println(rayMatrix)
            println(ray)
            error("Ray matrix value reassignment issue encountered.")
        end
    end
    return out
end


"""
    decode(::Array{String,2})

    Unused
"""
function decode(object::Array{String,2})
    return(map((x) -> parse(Int64, x), object))
end

## API functions

"""
    getRayStack(::StackyFan, ::Array{Int64, 1})

    Get the scalar associated with a ray in the given stacky fan structure.

# Examples
"""
function getRayStack(sf::StackyFan, ray::Array{Int64, 1})
    return sf.stacks[encode(ray)]
end

"""

    getMultiplicities(::StackyFan)

    Get the multiplicities of the cones in a stacky fan.

# Examples
"""
function getMultiplicities(sf::StackyFan)
    cones = getCones(sf)
    return(map((x,y) -> (x,y), cones, map(coneMultiplicity, cones)))
end

"""
    rootConstruction(::StackyFan, ::Array{Int64, 1})

Given a fan and a set of scalars corresponding to the rays of the fan,
performs a root construction on the fan by multiplying the stack scalars
by the given values. 

rootConstruction returns a new StackyFan object, and does not modify the input.

# Examples
```jldoctest StackyFan

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);

julia> SX = StackyFan(X, [2,3,5,7]);

julia> stackyWeights(rootConstruction(SX, [1, 4, 2, 1]))
[ 2 ,  12 ,  10 ,  7 ]
```
"""
function rootConstruction(
    sf::StackyFan,
    scalars::Array{Int64, 1})
    
    # Multiply the scalars of the fan by the given values
    return StackyFan(sf.fan, stackyWeights(sf) .* scalars)
end

"""
    rootConstructionDistinguishedIndices(::StackyFan, ::Array{Int64, 1}, ::Array{Int64, 1})

    Given a fan, the indices of the distinguished rays in the fan rays (as an incidence matrix), and
a set of scalars corresponding to the rays of the fan, performs a root 
construction on the fan by multiplying the stack scalars by the given values. 

    rootConstructionDistinguishedIndices returns a new StackyFan object,
and does not modify the input.

# Examples
```jldoctest StackyFan
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);
julia> SX = StackyFan(X, [2,3,5,7]);
julia> stackyWeights(rootConstructionDistinguishedIndices(SX, [0,1,1,0], [4, 2, 1, 3]))
[ 2 ,  6 ,  5 ,  7 ]
```
"""
function rootConstructionDistinguishedIndices(
    sf::StackyFan,
    distIndices::Array{Int64, 1},
    scalars::Array{Int64, 1})
    
    numRays = size(sf.fan.RAYS, 1)
    fullScalars = fill(1, numRays)
    for i in 1:numRays
        if distIndices[i]==1 && scalars[i] != 0
            fullScalars[i] = scalars[i]
        end
    end
    # Multiply the scalars of the fan by the given values
    return rootConstruction(sf, fullScalars)
end

"""
    rootConstructionDistinguished(
        ::StackyFan, 
        ::Polymake.Matrix{Polymake.Rational},
        ::Array{Int64, 1})

    Given a fan, a set of distinguished rays, and a set of scalars of equal size,
performs a root construction on the fan on the distinguished rays by multiplying 
the stack scalars by the given values.

    rootConstructionDistinguished returns a new StackyFan object, 
and does not modify the input.

# Examples
```jldoctest StackyFan
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0;1 1 0;1 0 1;1 1 1],INPUT_CONES=[[0,1,2],[1,2,3]]);

julia> SX = StackyFan(X, [2,3,5,7]);

julia> distinguished = X.RAYS[[2,3],:];

julia> stackyWeights(rootConstructionDistinguished(SX, distinguished, [4, 2]))
[ 2 ,  12 ,  10 ,  7 ]
```
"""
function rootConstructionDistinguished(
    sf::StackyFan,
    rays::Polymake.Matrix{Polymake.Rational},
    scalars::Array{Int64, 1})
    
    # Check that the rays and scalars are the same size
    #if (size(rays, 1) != length(scalars))
    #    error("Inputs are not of equal size")
    #end
    
    encoded_rays = mapslices(encode, rays, dims=2)
    # Make a copy of the dictionary
    newStacks = copy(sf.stacks)
    for i in 1:length(encoded_rays)
        ray = encoded_rays[i]
        # Multiply the scalar of the corresponding ray
        newStacks[ray] *= scalars[i]
    end
    
    # Convert the dictionary to an array of scalars matching the indices
    #newScalars = mapslices(ray -> newStacks[encode(ray)], sf.fan.RAYS, dims=2)
    newScalars = Array{Int64, 1}()
    for i in 1:size(SX.fan.RAYS, 1)
        push!(newScalars, newStacks[encode(sf.fan.RAYS[i,:])])
    end
    
    return StackyFan(sf.fan, newScalars)
end

"""

    findBarycenter(::Union{AbstractSet,AbstractVector},::Polymake.BigObjectAllocated)

    Takes a normal toric variety X and a set s corresponding to a subset of rays of X, and outputs a polymake vector
    corresponding to the barycenter of those rays.

# Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]]);

julia> s=[1,2];

julia> findBarycenter(s,X)
pm::Matrix<pm::Integer>
2 1
```
"""

function findBarycenter(s::Union{AbstractSet,AbstractVector},X::Polymake.BigObjectAllocated)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.RAYS)))
    rays = rowMinors(rayMatrix, s)
    dim=size(rays,2)
    bary=zeros(Int64,dim,1)
    for i in 1:size(rays,1)
        bary+=rays[i,:]
    end
    return vec(bary)
end

"""
    findStackyBarycenter(::Union{AbstractSet,AbstractVector},::StackyFan)

    Takes a stacky fan SX and a set s corresponding to a subset of rays of SX, calculates the 'stacky rays' corresponding to those rays (the rays times their stacky values), and find the barycenter of the stacky rays.

# Examples
```jldoctest StackyFan
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[2,3]);

julia> findStackyBarycenter([1,2],F)
[ 5 ,  6 ]
```
"""

function findStackyBarycenter(s::Union{AbstractSet,AbstractVector},SX::StackyFan)
    rayMatrix=convert(Array{Int64,2}, Array(Polymake.common.primitive(SX.fan.RAYS)))
    # Multiply the rays by their stacky values
    stackMatrix = diagm(stackyWeights(SX)) * rayMatrix
    rays = rowMinors(stackMatrix, s)
    dim=size(rays,2)
    bary=zeros(Int64,dim,1)
    for i in 1:size(rays,1)
        bary+=rays[i,:]
    end
    return vec(bary)
end

"""

    toric_blowup(::Union{AbstractSet,AbstractVector},::Polymake.BigObjectAllocated,::AbstractVector)

    Takes a normal toric variety X, a set s corresponding to a subset of rays of X, and a (optional) polymake vector,
    v, blow up X at v. If v is not provided, blow up X at the barycenter of s.

"""
function toric_blowup(s, X, v)
    if size(v,2)==1
         v=transpose(v)
    end
    s = [i + 1 for i in s]
    if v==nothing
        v=findBarycenter(s,X)
    end
    coneList = convertIncidenceMatrix(X.MAXIMAL_CONES)
    # Extracting the indices of all the cones in X that contains the set of rays s
    starIndex = findall((t) -> all(((i) -> i in t).(s)), coneList)
    star = [coneList[i] for i in starIndex]
    rayMatrix = X.RAYS
    
    lattice = X.HASSE_DIAGRAM
    faces = @Polymake.convert_to Array{Set{Int}} lattice.FACES
    
    # Get all the subcones of X that is contained in one of the cones in star and has the same rank
    clStar = []
    # Iterate over star
    for t in star
        # Get the rank of t
        c = rank(Array(rowMinors(rayMatrix, t))) - 1
        # Get all the subcones of X with rank c
        rank_c_subcone_indices = @Polymake.convert_to Array{Int} Polymake.graph.nodes_of_rank(lattice,c)
        rank_c_subcones = [faces[i + 1] for i in rank_c_subcone_indices]
        # Iterate over rank_c_subcones, and put the cones that is contained in t into clStar
        for cone in rank_c_subcones
            new_cone = [i+1 for i in cone]
            if all((i -> i in t).(new_cone))
                push!(clStar, new_cone)
            end
        end
    end
    clStar = unique(clStar)
    
    n = size(rayMatrix, 1) + 1
    # Filter out the cones in star from conelist
    coneList = filter(x -> !(x in star), coneList)
    
    if length(s) == 1
        # If s consists of a single ray, find all the cones in clStar that does not contain s
        newCones = []
        for t in clStar
            if !(s[1] in t)
                push!(newCones, sort(push!(t, s[1])))
            end
        end
        # return newCones plus coneList
        finalCones = [[i - 1 for i in cone] for cone in append!(coneList, newCones)]
        return Polymake.fulton.NormalToricVariety(INPUT_RAYS = Array(X.RAYS), INPUT_CONES = finalCones)
    end
    newCones = []
    for t in clStar
        # Find all the cones in clStar that does not contain at least one ray in s
        # QUESTION: Why seperate this from the one element case? Any won't work with one element list?
        if any(((i) -> !(i in t)).(s))
            push!(newCones, push!(t, n))
        end
    end
    # return newCones plus coneList
    finalRays = vcat((X.RAYS),v)
    finalCones = [[i - 1 for i in cone] for cone in append!(coneList, newCones)]
    return Polymake.fulton.NormalToricVariety(INPUT_RAYS = finalRays, INPUT_CONES = finalCones)
end

"""

    stackyBlowup(::StackyFan,::Array{Int64,1},::Array{Int64,1})

    Takes a stacky fan sf, a ray excep, and a cone, and subdivides the stacky fan at the given ray. Crucially, the given cone should be the minimal cone containing the exceptional ray. The cone input should be zero-indexed.

#examples
```jldoctest StackyFan
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 1 2],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[2,3]);

julia> stackyWeights(stackyBlowup(F,[0,1],[1,1]))
[ 2 ,  1 ,  3 ]
```
"""

function stackyBlowup(sf::StackyFan, cone::Array{Int64,1}, excep::Array{Int64,1})
    # Express the exceptional ray as a scalar multiple of a primitive ray
    # Use this scalar as the stacky weight in the resulting stacky fan
    G=gcd(excep)
    excep=Polymake.common.primitive(excep)
    
    # Perform toric blowup at the given ray
    blowup = toric_blowup(cone, sf.fan, excep)
    sf.stacks[encode(excep)] = G

    return(StackyFan(blowup, sf.stacks))
end

"""
    getConesPolymake(sf::StackyFan)

    Returns a list of cones of a stacky fan, with the cones represented as polymake objects.

"""

function getConesPolymake(sf::StackyFan)
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
julia> A=[1 2; 3 4];

julia> slicematrix(A)
[[ 1 ,  2 ], [ 3 ,  4 ]]
```
"""
function slicematrix(A::AbstractMatrix{<:Number})
    return [A[i, :] for i in 1:size(A,1)]
end

"""
    rowMinors(::AbstractMatrix{<:Number},::Union{AbstractSet,AbstractVector})

    Identical to slicematrix, except only returns row vectors indexed by a set S.

# Examples
```jldoctest makeSmoothWithDependencies
julia> A=[1 2 3;4 5 6; 7 8 9];

julia> S=Set([1,3]);

julia> rowMinors(A,S)
2×3 LinearAlgebra.Transpose{Int64,Array{Int64,2}}:
 1  2  3
 7  8  9
```
"""
function rowMinors(A::AbstractMatrix{<:Number},S::Union{AbstractSet,AbstractVector})
    outList=[]
    slices=slicematrix(A)
    for i in 1:size(slices,1)
        if i in S
            append!(outList,[slices[i]])
        end
    end
    return Array(transpose(hcat(outList...)))
end

"""
    convertIncidenceMatrix(::Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric})

    Takes a Polymake incidence matrix (e.g., the output of X.MAXIMAL_CONES for a toric variety X) and outputs a list of vectors,
    with each vector recording the indices marked on a given row of the incidence matrix.

# Examples
```jldoctest makeSmoothWithDependencies
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0;1 1; 0 1],INPUT_CONES=[[0,1],[1,2]]);

julia> M=X.MAXIMAL_CONES;

julia> convertIncidenceMatrix(M)
[[ 1 ,  2 ], [ 2 ,  3 ]]
```
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
    return convert(Array{Array{Int64,1},1}, out)
end

"""

    coneMultiplicity(C::Polymake.BigObjectAllocated)

    Returns the multiplicity of a polyhedral cone (inputted as a Polymake object): here, the multiplicity is defined as the index of the sublattice generated by the edges of the cone, inside the full integer lattice contained in the linear subspace generated by the edges of the cone.

# Examples
```jldoctest StackyFan

julia> C=Polymake.polytope.Cone(INPUT_RAYS=[1 0; 1 2]);

julia> coneMultiplicity(C)
2
```
"""

function coneMultiplicity(C::Polymake.BigObjectAllocated)
    A=Polymake.common.primitive(C.RAYS)
    M=matrix(ZZ,[fmpz.(y) for y in A])
    SNF=Nemo.snf(M)
    mult=1
    for i in 1:size(SNF,1)
        mult*=SNF[i,i]
    end
    return mult
end

"""

    coneConvert(::abstractVector{Int64},::abstractMatrix{Int64})

    Takes a matrix where the columns represent rays, and a list of indices, and forms a Polymake cone object generated by the rays corresponding to those indices.

# Examples
```jldoctest StackyFan

julia> typeof(coneConvert([1, 2, 4],[1 0 0; 0 1 0; 0 0 1; 1 1 1]))
Polymake.BigObjectAllocated
```
"""

function coneConvert(cone::Array{Int64,1},rayMatrix::Array{Int64,2})
    coneRays=rowMinors(rayMatrix,cone)
    C=Polymake.polytope.Cone(RAYS=coneRays)
    return C
end
    
"""

    getCones(X::Polymake.BigObjectAllocated)
    
    Returns all the cones of a fan X as a list of lists, with each interior list containing the indices of the rays generating a given cone.

# Examples
```jldoctest StackyFan
julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 0; 1 1 0; 1 0 1; 1 1 ],INPUT_CONES=[[0,1,2,3]]);

julia> getCones(X)
[[ 0 ,  1 ,  2 ,  3 ], [ 0 ,  1 ], [ 0 ,  2 ], [ 2 ,  3 ], [ 1 ,  3 ], [ 0 ], [ 1 ], [ 2 ], [ 3 ]]
```
"""

function getCones(X::Polymake.BigObjectAllocated)
    lattice=X.HASSE_DIAGRAM
    faces=@Polymake.convert_to Array{Set{Int}} lattice.FACES
    out=Array{Int64,1}[]
    for i in 2:(size(faces,1)-1)
        newface=Array(@Polymake.convert_to Array{Int} faces[i])
        push!(out,[i+1 for i in newface])
    end
    return out
end

"""
        distinguishedAndMultiplicity(::Array{Int64,1},::Array{Int64,2},::Array{Int64,1})

    Calculates if the cone formed by a subset of rays in rayMatrix indexed by the entries of cone, and with a distinguished structure given by the incidence vector dist, both contains at least one distinguished ray and has multiplicity greater than 1.

# Examples
```jldoctest StackyFan
julia> distinguishedAndMultiplicity([1,2,4],[1 0 0; 1 2 0;2 1 3; 1 0 3],[1,0,0,0])
true
```
"""

function distinguishedAndIntPoint(cone::Array{Int64,1},rayMatrix::Array{Int64,2},dist::Array{Int64,1})
    l=size(rayMatrix,1)
    if dot(convertToIncidence(cone,l),dist) > 0 #check distinguished
        C=coneConvert(cone,rayMatrix)
        if interiorPoints(C)!=nothing #check interior point
            return true
        else
            return false
        end
    else
        return false
    end
end

"""

    convertToIncidence(v::Array{Int64,1},l::Int64)

Returns a vector of length l, with entries of 1 indexed by v and entries of 0 everywhere else.

# Examples
```jldoctest 
julia> convertToIncidence([2,3,5],6)
[ 0 , 1 , 1 , 0 , 1 , 0 ]
```
"""

function convertToIncidence(v::Array{Int64,1},l::Int64)
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

"""

    compareCones(::Array{Int64,1},::Array{Int64,1},::Array{Int64,2},::Array{Int64,1})

    Takes in two cones (in index vector notation), a ray matrix, and a incidence vector of distinguished rays. If the cones do not have an equal number of distinguished rays, returns the difference between the two values. Otherwise, returns the difference in the cone multiplicities.

# Examples
```jldoctest
julia> compareCones([1,2],[2,3],[1 0 0; 0 1 0; 0 0 1],[1,1,0])
1

julia> compareCones([1,2],[1,3],[1 0;1 2;1 -1],[1,1,1])
1
```
"""

function compareCones(cone1::Array{Int64,1}, cone2::Array{Int64,1}, rayMatrix::Array{Int64,2}, distinguished::Array{Int64,1})
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

"""

    extremalCones(::Array{Array{Int64,1},1},::Array{Int64,2},::Array{Int64,1})

    Takes a list of vectors representing cones in a fan, a ray matrix, and a vector representing the distinguished rays as 0 or 1 values, and calculates the cones that are maximal with respect to (first) the number of non-distinguished rays and (second) the multiplicity of the cone. In Bergh's algorithm A (where this ordering is used), the input S will consist only of those cones containing at least one distinguished ray and at least one interior point.

#Examples
```jldoctest StackyFan
julia> extremalCones([[1,2],[2,3],[3,4]],[1 0;1 2; 1 5; 1 8],[0,1,1,0])
[[ 3 ,  4 ]]
```
"""

function extremalCones(S::Array{Array{Int64,1},1}, rayMatrix::Array{Int64,2}, distinguished::Array{Int64,1})
    # Finds the extremal cones according to # distinguished rays and multiplicity
    # distinguished is a boolean vector whose size is equal to the number of rays
    # The i-th index is 1 if the i-th ray (in rayMatrix) is distinguished
    maxCones = [S[1]]
    for i in 2:size(S,1)
        cone = S[i]
        # Compare the cone with the first element of the maximal cone list
        comp = compareCones(cone, maxCones[1], rayMatrix, distinguished)
        if comp > 0
            maxCones = [cone]
        elseif comp == 0
            push!(maxCones, cone)
        end
    end
    return maxCones
end

"""

    interiorPoints(::Polymake.BigObjectAllocated)

    Finds all interior lattice points contained in the fundamental region of a given cone. When multiple interior lattice points lie along the same ray, only the point closest to the origin is returned. Notably, 

# Examples
```jldoctest StackyFan
julia> C=Polymake.polytope.Cone(INPUT_RAYS=[1 2; 2 1]);

julia> interiorPoints(C)
[[ 1 ,  1 ]]
```
"""

function interiorPoints(C::Polymake.BigObjectAllocated)
    rayMatrix=Array(Polymake.common.primitive(C.RAYS))
    l=size(rayMatrix,1)
    dim=size(rayMatrix,2)
    if rank(rayMatrix)<l
        print(rayMatrix)
        error("Input cone is not simplicial.")
    end
    subsets=collect(powerset([1:l;]))
    vertices=[]
    for elt in subsets #vertices of the fundamental region are in correspondence with subsets of the generators of the cone, by summing the generators in a subset to obtain a vertex
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
    P=Polymake.polytope.Polytope(POINTS=VH) #make a Polymake polytope object from the vertices of the fundamental region found in the last step
    if size(P.INTERIOR_LATTICE_POINTS,1)==0
        return nothing
    end
    intPoints=Array(P.INTERIOR_LATTICE_POINTS)[:,2:(dim+1)] #find all the interior lattice points
    validPoints=[]
    #return intPoints
    for i in 1:size(intPoints,1) #throw out all points that are integer multiples of other points
        point=intPoints[i,:]
        if gcd(point)==1
            append!(validPoints,[point])
        end
    end
    return validPoints
end

"""

    minimalByLex(::Array{Array{Int64,1},1})

    Given a list of vectors of equal length, returns the minimal vector with respect to lexicographic ordering.

# Examples
```jldoctest StackyFan
julia> A=[[1,1,1],[2,1,3],[0,5,4]];

julia> minimalByLex(A)
[ 0 ,  5 ,  4 ]
```
"""

function minimalByLex(A::Array{Array{Int64,1},1})
    l=size(A,1)
    minimal=A[1]
    d=size(minimal,1)
    for i in 2:l
        test=A[i]
        for j in 1:d
            if minimal[j]<test[j]
                break
            elseif minimal[j]>test[j]
                minimal=test
                break
            end
        end
    end
    return minimal
end

"""
    
    minimalByDist(::Array{Array{Int64,1},1},::Array{Int64,1})

    Given a list of vectors (representing rays as weighted sums of other rays) and a vector of 0's and 1's representing non-distinguished and distinguished indices, returns a vector from the list such that the sum of the entries corresponding to distinguished indices is minimized.

#Examples
```jldoctest StackyFan
julia> minimalByDist([[0,1,5,7],[3,3,2,2],[8,5,3,6],[2,1,1,10]],[0,1,1,0])
[ 3 , 3 , 2 , 2 ]
```
"""

function minimalByDist(A::Array{Array{Int64,1},1},D::Array{Int64,1})
    l=size(A,1)
    minimal=A[1]
    d=size(minimal,1)
    for i in 2:l
        test=A[i]
        if dot(test,D)<dot(minimal,D)
            minimal=test
        end
    end
    return minimal
end

"""
    coneRayDecomposition(::Array{Int64,1},::Array{Int64,2},::Array{Int64,1},::Array{Int64,1})

    This function takes in a cone (a vector of indices of cone generators in rayMatrix), a ray, and a stacky structure for rayMatrix. It first multiplies all generators of the cone by their stacky values, and then finds an expression for the ray as a sum of these stacky generators. The output is a vector of coefficients of the above representation in terms of the rays in rayMatrix, with zeros as coefficients for all rays not in the given cone.

# Examples
```jldoctest StackyFan
julia> coneRayDecomposition([1,2,3],[3 5 7; 8 16 9;2 1 3;1 1 1],[2,2,3],[1,1,1,1])
[ 6 ,  5 ,  52 ,  0 ]
```
"""

function coneRayDecomposition(cone,rayMatrix,ray,stack)
    stackMatrix=diagm(stack)*rayMatrix # multiply all rays by stack values
    coneRays=rowMinors(stackMatrix,cone) # find the (stackified) rays in the given cone
    if rank(coneRays)<size(coneRays,1)
        error("The given cone is not simplicial.")
    end
    B=Polymake.common.null_space(hcat(transpose(coneRays),-ray)) # Express the input ray in terms of the stackified cone generators
    N=convert(Array{Int64,1},vec(B))
    if size(N,1)==0
        error("The given ray is not in the span of the cone generators.")
    end
    if N[end]<0 #since the nullspace has arbitrary sign, fix it so the coefficients are all positive
        N*=-1
    end
    pop!(N)
    out=zeros(Int64,size(rayMatrix,1)) 
    for i in 1:size(N,1)#rewrite the coefficients vector in terms of all the rays in rayMatrix, by padding with zeros when appropriate.
        out[cone[i]]=N[i] 
    end
    return out
end

"""

    BerghA(F::StackyFan,D::Array{Int64,1})

Given a stacky fan F and a vector of booleans D representing the distinguished structure, returns a smooth stacky fan where the distinguished rays are independent.

The algorithm is adapted from Daniel Bergh's [paper on destackification](https://arxiv.org/abs/1409.5713). In brief, it identifies non-smooth cones containing at least one distinguished ray, finds interior points in those cones, and subdivides at those points through a series of stacky barycentric subdivisions.

# Examples
```jldoctest StackyFan.jl

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[1,1]);

julia> stackyWeights(BerghA(F,[1,1]))
[ 5 ,  2 ,  5 ,  10 ]

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[1,1]);

julia> stackyWeights(BerghA(F,[1,0]));
[ 5 ,  5 ,  1 ,  2 ,  5 ,  10 ]

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 2 5],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[1,1]);

julia> stackyWeights(BerghA(F,[0,1]))
[ 1 ,  5 ,  5 ,  2 ,  5 ,  1 ,  5 ,  10 ]

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[4 1; 7 9],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[1,1])l;

julia> stackyWeights(BerghA(F,[1,0]))
[ 609 ,  29 ,  1 ,  174 ,  29 ,  1740 ,  1218 ,  58 ,  145 ,  1044 ,  1044 ,  290 ,  290 ,  145 ,  406 ,  348 ,  261 ,  14616 ,  14616 ,  609 ,  3480 ,  870 ,  609 ,  174 ,  609 ,  9744 ,  14616 ,  1218 ,  58 ,  145 ,  435 ,  725 ,  1305 ,  1740 ,  6960 ,  3480 ,  870 ,  3480 ,  58464 ,  6090 ,  3480 ,  1392 ,  696 ,  1044 ,  2088 ,  261 ,  174 ,  261 ,  609 ,  406 ,  609 ,  609 ,  406 ,  609 ,  1218 ,  812 ,  1218 ,  2088 ,  261 ,  1044 ,  1160 ,  1740 ,  1740 ,  870 ,  1305 ,  1305 ,  14616 ,  10440 ,  145 ,  435 ,  609 ,  116 ,  580 ,  290 ,  580 ,  1740 ,  3480 ,  3480 ,  261 ,  522 ,  261 ,  522 ,  522 ,  1218 ,  1218 ,  1218 ,  1218 ,  2436 ,  2436 ,  4176 ,  2088 ,  3480 ,  2610 ,  29232 ,  6090 ,  1218 ,  290 ,  145 ,  290 ,  12180 ,  261 ,  522 ]

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[323 179; 44 135],INPUT_CONES=[[0,1]]);

julia> F=addStackStructure(X,[1,1]);

julia> stackyWeights(BerghA(F,[1,1]))
[ 491602456800 ,  49160245680 ,  294961474080 ,  468192816 ,  12173013216 ,  73038079296 ,  12173013216 ,  97384105728 ,  73038079296 ,  196640982720 ,  131093988480 ,  156064272 ,  936385632 ,  468192816 ,  196640982720 ,  49160245680 ,  9832049136 ,  292152317184 ,  5899229481600 ,  3932819654400 ,  589922948160 ,  393281965440 ,  12173013216 ,  8115342144 ,  12173013216 ,  5899229481600 ,  589922948160 ,  737403685200 ,  51126655507200 ,  292152317184 ,  51126655507200 ,  5899229481600 ,  5899229481600 ,  589922948160 ,  589922948160 ,  196640982720 ,  2949614740800 ,  5899229481600 ,  294961474080 ,  589922948160 ,  196640982720 ,  589922948160 ,  393281965440 ,  936385632 ,  24346026432 ,  24346026432 ,  11798458963200 ,  1179845896320 ,  737403685200 ,  1474807370400 ,  2949614740800 ,  5899229481600 ,  589922948160 ,  11798458963200 ,  1179845896320 ,  2949614740800 ,  1474807370400 ,  5899229481600 ,  102253311014400 ]

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 3; 4 5 6; 2 3 1],INPUT_CONES=[[0,1,2]]);

julia> F=addStackStructure(X,[1,1,1]);

julia> stackyWeights(BerghA(F,[1,1,1]))
[ 28 ,  21 ,  84 ,  28 ,  84 ,  84 ,  42 ,  84 ,  168 ]

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0 2; 2 1 1; 5 3 9],INPUT_CONES=[[0,1,2]]);

julia> F=addStackStructure(X,[1,1,1]);

julia> stackyWeights(BerghA(F,[1,1,1]))
[ 4 ,  4 ,  8 ,  4 ,  4 ,  8 ]

julia> X=Polymake.fulton.NormalToricVariety(INPUT_RAYS=[1 0; 4 3; 1 5],INPUT_CONES=[[0,1],[1,2]]);

julia> F=addStackStructure(X,[1,1,1]);

julia> stackyWeights(BerghA(F,[1,1,1]))
[ 4 ,  34 ,  6 ,  68 ,  34 ,  68 ,  34 ,  6 ,  12 ]
```
"""
function BerghA(F::StackyFan,D::Array{Int64,1};verbose::Bool=false)
    if verbose==true
        println("==algorithm is running in verbose mode==")
        println(" ")
        println("=======")
    end
    X=deepcopy(F)
    rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
    coneList=getCones(X.fan)
    dim=size(rayMatrix,2)
    numRays=size(rayMatrix,1)
    
    #check if the vector D has length equal to the number of rays in F
    if numRays != size(D,1)
        error("length of vector representing distinguished structure does not agree with number of rays in stacky fan.")
    end
    
    #A0: initialization
    i=0
    while(true)
        rayMatrix=convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS)))
        numRays=size(rayMatrix,1)
        coneList=getCones(X.fan)
        
       #Debugging
        #if verbose==true
            #println("Ray matrix: $rayMatrix")
            #println("Cone list: $coneList")
            #sW=stackyWeights(X)
            #println("Stacky weights: $sW")
            #println("Distinguished rays: $D")
        #end
        #coneMultiplicities=Int64[]
        #for cone in coneList
        #    C=coneConvert(cone,rayMatrix)
        #    push!(coneMultiplicities,coneMultiplicity(C))
        #end
        #if verbose==true
            #println("Cone multiplicities: $coneMultiplicities")
        #end
       #End debugging
        
        #A1: Find S the set of cones that contain a distinguised ray and an interior lattice point 
        #Note: cones in S are 1-indexed.
        S=filter(cone->distinguishedAndIntPoint(cone,rayMatrix,D),coneList)
        # If S is empty, the program terminates.
        if S==[]
            break
        end
        
        #A2 - find extremal cones
        Smax=extremalCones(S,rayMatrix,D)
        
        #Print information on the number of extremal cones, their number of non-distinguished rays, and their multiplicity
        #The algorithm is structured to first reduce the number of non-distinguished rays in extremal cones, and then reduce the multiplicity of said cones,
            #so this information can be used to track the algorithm's progress
        if verbose==true
            Smaxcount=size(Smax,1)
            println("Number of extremal cones: $Smaxcount")
            testCone=Smax[1]
            c1=convertToIncidence(testCone,numRays)
            nonDist=size(testCone,1)-dot(c1,D)
            mult=coneMultiplicity(coneConvert(testCone,rayMatrix))
            println("Maximal non-distinguished rays and multiplicity: $nonDist, $mult")
        end
        
        #A2 - find interior points in Smax
        intPoints=[]
        for cone in Smax
            #C=rowMinors(rayMatrix,cone)
            C=coneConvert(cone,rayMatrix)
            coneIntPoints=interiorPoints(C)
            for point in coneIntPoints
               push!(intPoints,(point,cone)) #the point is stored as a tuple along with the cone containing it
            end
        end
        
        #A2 - find stacky points (in terms of coefficients) derived from interior points
        P=Array{Int64,1}[]
        for (point,cone) in intPoints
            stackyPoint=coneRayDecomposition(cone,rayMatrix,point,stackyWeights(X)) #each interior point is rewritten as a string of coefficients
                #corresponding to its representation as a sum of stacky rays
            push!(P,stackyPoint)
        end
        
        #A2 - find element of P such that the sum of the coefficients corresponding to distinguished rays is minimal.
            #This invariant does not produce a unique ray, so there is a degree of arbitrary selection.
        psi=minimalByDist(P,D)
        #if verbose==true
            #println("Psi: $psi")
        #end
        
        
        #A3 - perform root construction
        X=rootConstructionDistinguishedIndices(X,D,psi)
        
        #A3 - modify psi with respect to root construction
        for i in 1:length(psi)
            if D[i]==1 && psi[i]>0
                psi[i]=1
            end
        end
        
        #A5 - perform repeated stacky barycentric star subdivision with respect to psi.
        while(count(x->x>0,psi)>1)
            #A4 - perform stacky star subdivision
            # Get the indices of the non-zero coefficients in psi - this is used to determine the cone 
                #containing the support of psi, which will be subdivided
            
            supportCone=findall(x->x!=0,psi)
            #find the stacky barycenter of that cone, which becomes the exceptional (blowup) ray
            exceptional=findStackyBarycenter(supportCone,X)
            
            #Performing a blowup may cause the rays defined in the StackyFan struct to be reordered - this property is inherited from
                #the Oscar/Polymake fan object. Since D (distinguished rays) and psi are lists, they must be reordered to match the
                #new order of the rays of X. This reordering is accomplished through defining dictionaries before the blowup is performed.
            
            code_rays = mapslices(encode, Polymake.common.primitive(X.fan.RAYS), dims=2)
            # Track the indices of distinguished rays
            D_pairs = map((x,y) -> (x,y), code_rays, D)
            D_Dict = Dict(D_pairs)
            # Track psi as a linear combination of the generators
            psiPairs = map((x,y) -> (x,y), code_rays,psi)
            psiDict = Dict(psiPairs)

            #perform the blowup
            X=stackyBlowup(X,[x-1 for x in supportCone],exceptional)
            
            G=gcd(exceptional) #since the blowup ray may not be primitive, it is made primitive and then assigned a stacky value so its stacky form is unchanged.
            primExcep=Polymake.common.primitive(exceptional)
            
            # Update the dictionaries storing fan information
            D_Dict[encode(primExcep)]=1
            psiDict[encode(primExcep)]=1
            
            #create new lists
            newRays=slicematrix(convert(Array{Int64,2},Array(Polymake.common.primitive(X.fan.RAYS))))
            newD=Int64[]
            newpsi=Int64[]
            for ray in newRays
                E=encode(ray)
                excepCode=encode(primExcep)
                push!(newD,D_Dict[E])
                #A4 - modify psi
                if E==excepCode
                    push!(newpsi,1)
                elseif psiDict[E]>1
                    push!(newpsi,psiDict[E]-1)
                else
                    push!(newpsi,0)
                end
            end
            psi=newpsi
            D=newD
        end
        if verbose==true
            println("=======")
        end
        i+=1
    end
    if verbose==true
        println("Number of steps: $i")
    end
    return X
end