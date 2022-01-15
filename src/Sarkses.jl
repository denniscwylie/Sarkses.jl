module Sarkses


## =============================================================================
module SkewSuffixArrays

## Implementation of skew algorithm for construction of suffix array
## based on code provided in:
## "Kärkkäinen, J., & Sanders, P. (2003, June). Simple linear work suffix
##  array construction. In International Colloquium on Automata, Languages,
##  and Programming (pp. 943-955). Springer, Berlin, Heidelberg."

function radixPass!(a::Array{Int64}, b::Array{Int64}, r::Array{Int64},
                    rShift::Int64, n::Int64, K::Int64)::Nothing
    Kp1 = K+1
    ## count occurrences
    c = zeros(Int64, Kp1)
    @inbounds for i = 1:n
        c[r[rShift + a[i] + 1] + 1] += 1
    end
    ## exclusive prefix sums
    sum = 0
    for i = 1:Kp1
        t = c[i]
        c[i] = sum
        sum += t
    end
    ## sort
    @inbounds for i = 1:n
        b[c[r[rShift + a[i] + 1] + 1] + 1] = a[i]
        c[r[rShift + a[i] + 1] + 1] += 1        
    end
    return nothing
end

function suffixArray!(s::Array{Int64}, SA::Array{Int64}, n::Int64, K::Int64)::Nothing
    n0=(n+2)÷3; n1=(n+1)÷3; n2=n÷3; n02=n0+n2
    s12 = zeros(Int64, n02+3)
    SA12 = zeros(Int64, n02+3)
    s0 = Array{Int64}(undef, n0)
    SA0 = Array{Int64}(undef, n0)
    ##
    ## generate positions of mod 1 and mod 2 suffixes
    ## the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
    let j = 1
        @inbounds for i = 0:(n+(n0-n1)-1)
            if i%3 != 0
                s12[j] = i
                j += 1            
            end
        end
    end
    ##
    ## lsb radix sort the mod 1 and mod 2 triples
    radixPass!(s12, SA12, s, 2, n02, K)
    radixPass!(SA12, s12, s, 1, n02, K)
    radixPass!(s12, SA12, s, 0, n02, K)
    ##
    ## find lexicographic names of triples
    let name = 0, c0 = -1, c1 = -1, c2 = -1
        @inbounds for i = 1:n02
            if s[SA12[i]+1] != c0 || s[SA12[i]+2] != c1 || s[SA12[i]+3] != c2
                name += 1
                c0 = s[SA12[i]+1]
                c1 = s[SA12[i]+2]
                c2 = s[SA12[i]+3]
            end
            if (SA12[i] % 3) == 1
                ## left half
                s12[SA12[i]÷3 + 1] = name
            else
                ## right half
                s12[SA12[i]÷3 + n0+1] = name
            end
        end
        ##
        ## recurse if names are not yet unique
        if name < n02
            suffixArray!(s12, SA12, n02, name)
            ## store unique names in s12 using the suffix array
            @inbounds for i = 1:n02
                s12[SA12[i]+1] = i
            end
        else
            ## generate the suffix array of s12 directly
            @inbounds for i = 1:n02
                SA12[s12[i]] = i-1
            end
        end
    end
    ##
    ## stably sort the mod 0 suffixes from SA12 by their first character
    let j = 1
        @inbounds for i = 1:n02
            if SA12[i] < n0
                s0[j] = 3 * SA12[i]
                j += 1
            end
        end
    end
    radixPass!(s0, SA0, s, 0, n0, K)
    ##
    ## merge sorted SA0 suffixes and sorted SA12 suffixes
    let p=0; t=n0-n1; k=0
        while k < n
            ## pos of current offset 12 suffix
            i = (SA12[t+1] < n0 ? SA12[t+1] * 3 + 1 : (SA12[t+1] - n0) * 3 + 2)
            ## pos of current offset 0 suffix
            j = SA0[p+1]
            ##
            if (SA12[t+1] < n0 ?
                (s[i+1], s12[SA12[t+1]+n0+1]) < (s[j+1], s12[j÷3+1]) :
                (s[i+1], s[i+2], s12[SA12[t+1]-n0+2]) < (s[j+1], s[j+2], s12[j÷3+n0+1]))
                ## suffix from SA12 is smaller
                SA[k+1] = i
                t += 1
                if t == n02
                    ## done --- only SA0 suffixes left
                    k += 1
                    while p < n0
                        SA[k+1] = SA0[p+1]
                        p += 1
                        k += 1
                    end
                end
            else
                SA[k+1] = j
                p += 1                
                if p == n0
                    ## done --- only SA12 suffixes left
                    k += 1
                    while t < n02
                        SA[k+1] = (SA12[t+1] < n0 ?
                                   SA12[t+1] * 3 + 1 :
                                   (SA12[t+1] - n0) * 3 + 2)
                        t += 1
                        k += 1
                    end
                end
            end
            k += 1
        end
    end
    return nothing
end

end # module SkewSuffixArrays


## =============================================================================
module FMIndices

## Implementation of Burrows-Wheeler transform for purpose of exact
## k-mer matching based on backward search algorithm described in:
## "Ferragina, P., & Manzini, G. (2000, November). Opportunistic
##  data structures with applications. In Proceedings 41st Annual
##  Symposium on Foundations of Computer Science (pp. 390-398). IEEE."

## Implementation of inexact k-mer matching done following procedure
## building on above reference described in:
## "Heng Li, Richard Durbin (2009, July). Fast and accurate short read
##  alignment with Burrows–Wheeler transform. Bioinformatics, Volume 25,
##  Issue 14 (pp. 1754–1760). Oxford University Press."

function fillFMIndexCO!(bwt::Array{UInt8},
                        c::Dict{UInt8,Int64},
                        o::Dict{UInt8,Array{Int64}},
                        seq::Array{UInt8},
                        sa::Array{Int64},
                        blockSize::Int64)::Nothing
    oSize = ((length(bwt)-1) ÷ blockSize) + 1
    @inbounds for i=1:length(bwt)
        if sa[i] > 1
            ichar = seq[ sa[i]-1 ]
        else
            ichar = b"!"[1]
        end
        bwt[i] = ichar
        if !(ichar in keys(c))
            c[ichar] = 0
        end
        c[ichar] += 1
        if !(ichar in keys(o))
            o[ichar] = zeros(Int64, oSize)
        end
        o[ichar][1 + ((i-1) ÷ blockSize)] += 1
    end
    cumul = 0
    cKeys = sort(collect(keys(c)))
    @inbounds for ichar in cKeys
        cumul += c[ichar]
        c[ichar] = cumul - c[ichar]
        ochar = o[ichar]
        for i=2:oSize
            ochar[i] += ochar[i-1]
        end
    end
    return nothing
end

mutable struct FMIndex
    seq::Array{UInt8}
    sa::Array{Int64}
    bwt::Array{UInt8}
    c::Dict{UInt8,Int64}
    o::Dict{UInt8,Array{Int64}}

    blockSize::Int64

    function FMIndex(seq::Array{UInt8}, sa::Array{Int64})
        this = new()
        this.seq = seq
        this.sa = sa
        this.bwt = Array{UInt8}(undef, length(this.seq))
        this.c = Dict{UInt8,Int64}()
        this.o = Dict{UInt8,Int64}()
        this.blockSize = 100
        fillFMIndexCO!(this.bwt, this.c, this.o,
                       this.seq, this.sa, this.blockSize)
        return this
    end
end

function handleBlocking(bwt::Array{UInt8},
                        c::Dict{UInt8,Int64},
                        o::Dict{UInt8,Array{Int64}},
                        i::Int64,
                        ichar::UInt8,
                        blockSize::Int64)::Int64
    iRounded = ((i-1) ÷ blockSize)
    iNew = c[ichar] + (iRounded > 0 ? o[ichar][iRounded] : 0) + 1
    for j=(1+(blockSize*iRounded)):(i-1)
        if bwt[j] == ichar
            iNew += 1
        end
    end
    return iNew
end

function findKmer(this::FMIndex, kmer::Array{UInt8})::Tuple{Int64,Int64}
    k = length(kmer)
    start = 1
    end_ = length(this.bwt) + 1
    for i=k:-1:1
        ichar = kmer[i]
        start = handleBlocking(
            this.bwt, this.c, this.o, start, ichar, this.blockSize
        )
        end_ = handleBlocking(
            this.bwt, this.c, this.o, end_, ichar, this.blockSize
        )
    end
    return (start, end_)
end

function calculateD(rev::FMIndex, kmer::Array{UInt8})::Array{Int64}
    k = length(kmer)
    d = zeros(k)
    start = 1
    end_ = length(rev.bwt) + 1
    z = 0
    for i=1:k
        ichar = kmer[i]
        start = handleBlocking(rev.bwt, rev.c, rev.o, start, ichar,
                               rev.blockSize)
        end_ = handleBlocking(rev.bwt, rev.c, rev.o, end_, ichar,
                              rev.blockSize)
        if end_ <= start
            start = 1
            end_ = length(rev.bwt) + 1
            z += 1
        end
        d[i] = z
    end
    return d
end

function findKmer(this::FMIndex, kmer::Array{UInt8},
                  nMismatch::Int64, d::Array{Int64},
                  i::Int64, start::Int64, end_::Int64;
                  indels::Bool=false)::Array{Tuple{Int64,Int64}}
    if i < 1
        if nMismatch >= 0
            return [(start, end_)]
        else
            return Array{Tuple{Int64,Int64}}(undef, 0)
        end
    end
    if nMismatch < d[i]
        return Array{Tuple{Int64,Int64}}(undef, 0)
    end
    I = Array{Tuple{Int64,Int64}}(undef, 0)
    if indels
        append!(I, findKmer(
            this, kmer, nMismatch-1, d, i-1, start, end_, indels=true
        ))
    end
    ichar = kmer[i]
    for bchar in sort(collect(keys(this.c)))
        bstart = handleBlocking(this.bwt, this.c, this.o, start, bchar,
                                this.blockSize)
        bend = handleBlocking(this.bwt, this.c, this.o, end_, bchar,
                              this.blockSize)
        if bend > bstart
            if indels
                append!(I, findKmer(
                    this, kmer, nMismatch-1, d, i, bstart, bend, indels=true
                ))
            end
            if bchar == ichar
                append!(I, findKmer(
                    this, kmer, nMismatch, d, i-1, bstart, bend, indels=indels
                ))
            else
                append!(I, findKmer(
                    this, kmer, nMismatch-1, d, i-1, bstart, bend, indels=indels
                ))
            end
        end
    end
    return I
end

function findKmer(this::FMIndex, rev::FMIndex, kmer::Array{UInt8},
                  nMismatch::Int64; indels::Bool=false)::Array{Tuple{Int64,Int64}}
    d = calculateD(rev, kmer)
    return findKmer(this, kmer, nMismatch, d,
                    length(kmer), 1, length(this.bwt)+1, indels=indels)
end

end # module FMIndices


using CodecZlib
using CSV
using DataFrames
using DataStructures
using LinearAlgebra
using Random
using Statistics
using TranscodingStreams

function multifind(needles::Array{<:AbstractString},
                   haystack::AbstractArray{<:AbstractString}
                   )::Array{Int64}
    out = -1 * ones(Int64, length(needles))
    ln = length(needles); lh = length(haystack)
    nsp = sortperm(needles); hsp = sortperm(haystack)
    i = 1
    j = 1
    while (i <= ln) && (j <= lh)
        if haystack[hsp[j]] > needles[nsp[i]]
            i += 1
        elseif haystack[hsp[j]] == needles[nsp[i]]
            out[nsp[i]] = hsp[j]
            j += 1
        else
            j += 1
        end
    end
    return out
end

function c1idx(df::DataFrame, idx::Array{<:AbstractString})::DataFrame
    ndx = multifind(idx, df[:, 1])
    return df[ndx[ndx .> 0], :]
end

function c1idx(od::OrderedDict, idx::Array)::OrderedDict
    return OrderedDict([k => od[k] for k in idx])
end

function padAndUnstack(stacked::DataFrame,
                       rowkey::Union{Integer,Symbol}=:id,
                       colkey::Union{Integer,Symbol}=:bin,
                       value::Union{Integer,Symbol}=:score;
                       side::Symbol=:left,
                       fill::Union{Nothing,Float64}=nothing,
                       n::Union{Nothing,Int64}=nothing
                       )::DataFrame
    if n == nothing
        n = maximum(stacked[:, colkey])
    end
    stacked = deepcopy(stacked)
    stacked = stacked[:, [rowkey, colkey, value]]
    uniqueIds = unique(stacked[:, rowkey])
    toAdds = Array{DataFrame}(undef, 0)
    idToInds = Dict{String,Array{Int64}}()
    for id in uniqueIds
        idToInds[id] = Array{Int64}(undef, 0)
    end
    for i in 1:size(stacked, 1)
        push!(idToInds[stacked[i, 1]], i)
    end
    for id in uniqueIds
        rmatch = idToInds[id]
        maxBin = maximum(stacked[rmatch, colkey])
        if maxBin < n
            if side == :left
                stacked[rmatch, colkey] .+= (n - maxBin)
            end
            toAdd = DataFrame(r=[id for i in 1:(n-maxBin)], c=1:(n-maxBin))
            if side != :left
                toAdd.c .+= maxBin
            end
            vVal = if fill == nothing
                mean(stacked[rmatch, value])
            else
                fill
            end
            toAdd[:, :v] = [vVal for i in 1:size(toAdd, 1)]
            rename!(toAdd, [names(toAdd)[i] =>
                            Symbol(names(stacked)[i]) for i in 1:3])
            push!(toAdds, toAdd)
        end
    end
    if length(toAdds) > 0
        stacked = vcat(stacked, vcat(toAdds...))
    end
    stacked = stacked[sortperm(stacked[:, colkey]), :]
    return unstack(stacked, rowkey, colkey, value)
end

function readFasta(filename::String)::OrderedDict{String,String}
    out = OrderedDict{String,String}()
    curName = ""
    sb = IOBuffer()
    fgz = (endswith(lowercase(filename), ".gz") ||
           endswith(lowercase(filename), ".gzip"))
    fh = (fgz ?
          TranscodingStream(GzipDecompressor(), open(filename, "r")) :
          open(filename, "r"))
    while !eof(fh)
        curLine = rstrip(readline(fh))
        if startswith(curLine, ">")
            if curName != ""
                out[curName] = String(take!(sb))
            end
            curName = curLine[2:length(curLine)]
        else
            write(sb, curLine)
        end
    end
    if curName != "" && !(curName in keys(out))
        out[curName] = String(take!(sb))
    end
    close(fh)
    close(sb)
    return out
end



## =============================================================================
abstract type AbSarks end


function concatenateSeqs!(this::AbSarks,
                          seqs::AbstractDict{<:AbstractString,<:AbstractString}
                          )::Nothing
    sb = IOBuffer()
    this.bounds = Array{Int64}(undef, length(this.blocks)+1)
    this.bounds[1] = 1
    let i = 2
        for block in this.blocks
            blockSeq = String(seqs[block])
            write(sb, blockSeq)
            write(sb, "\$")
            this.bounds[i] = this.bounds[i-1] + length(blockSeq) + 1
            i += 1
        end
    end
    this.catSeq = take!(sb)
    close(sb)
    return nothing
end

function calcSuffixArray!(this::AbSarks)::Nothing
    chars = convert.(Int64, [this.catSeq[i] for i in 1:length(this.catSeq)])
    append!(chars, [0, 0, 0])
    SA = zeros(Int64, length(chars))
    SkewSuffixArrays.suffixArray!(chars, SA, length(this.catSeq), maximum(chars))
    this.sa = SA[1:length(this.catSeq)] .+ 1
    this.saInv = Array{Int64}(undef, length(this.sa))
    for i = 1:length(this.sa)
        this.saInv[ this.sa[i] ] = i
    end
    return nothing
end


## -----------------------------------------------------------------------------
function reset!(this::AbSarks, halfWindow::Int64,
                spatialLength::Int64=0, smooth::Bool=true)::Nothing
    this.halfWindow = halfWindow
    this.spatialLength = spatialLength
    windowGini!(this)
    if smooth
        window!(this)
    else
        this.windowed = Array{Float32}(undef, 0)
    end
    if this.spatialLength > 1
        spatialGini!(this)
        if smooth
            spatialWindow!(this)
        else
            this.spatialWindowed = Array{Float32}(undef, 0)
        end
    else
        this.spatGini = Array{Float32}(undef, 0)
        this.spatialWindowed = Array{Float32}(undef, 0)
    end
    return nothing
end

function resetSpatial!(this::AbSarks, spatialLength::Int64,
                       smooth::Bool=true)::Nothing
    this.spatialLength = spatialLength
    if this.spatialLength > 1
        spatialGini!(this)
        if smooth
            spatialWindow!(this)
        else
            this.spatialWindowed = Array{Float32}(undef, 0)
        end
    else
        this.spatGini = Array{Float32}(undef, 0)
        this.spatialWindowed = Array{Float32}(undef, 0)
    end
    return nothing
end


## -----------------------------------------------------------------------------
function sourceBlock(this::AbSarks, s::Int64,
                     mindex::Int64=1, maxdex::Int64=-1)::Int64
    if maxdex == -1
        maxdex = length(this.bounds) + 1
    elseif maxdex == (mindex + 1)
        return mindex
    end
    middex = (mindex + maxdex) ÷ 2
    if s < this.bounds[middex]
        return sourceBlock(this, s, mindex, middex)
    else
        return sourceBlock(this, s, middex, maxdex)
    end
end

function sourceBlock(this::AbSarks, s0::Array{Int64})::Array{Int64}
    maxdex = length(this.bounds) + 1
    return [sourceBlock(this, s, 1, maxdex) for s in s0]
end

function fillSourceBlock!(b::Array{Int64},
                          blocks::Array{String},
                          bounds::Array{Int64},
                          saInv::Array{Int64})::Nothing
    @inbounds for block in 1:length(blocks)
        nextBlock = block + 1
        for s in bounds[block]:(bounds[nextBlock]-1)
            b[ saInv[s] ] = block
        end
    end
end

function sourceBlock(this::AbSarks)::Array{Int64}
    b = Array{Int64}(undef, length(this.sa))
    fillSourceBlock!(b, this.blocks, this.bounds, this.saInv)
    return b
end


## -----------------------------------------------------------------------------
function giniImpurity(counts::Array{Int64})::Float64
    total = sum(counts)
    out = sum(convert(Array{Float64}, counts) .* (total .- counts))
    return(out / (total * total))
end

function fillSpatGini!(spatGini::Array{Float32},
                       spatialLength::Int64, saLen::Int64,
                       windGini::Array{Float32}, saInv::Array{Int64},
                       runningSum::Float64, windSize::Float64)::Float64
    @inbounds for s=2:(saLen+1-spatialLength)
        runningSum -= windGini[ saInv[s-1] ]
        runningSum += windGini[ saInv[s+spatialLength-1] ]
        spatGini[ saInv[s] ] = convert(Float32, runningSum / windSize)
    end
    return runningSum
end

function spatialGini!(this::AbSarks)::Nothing
    windSize = convert(Float64, this.spatialLength)
    this.spatGini = Array{Float32}(undef, length(this.sa))
    runningSum = 0.0
    for s = 1:this.spatialLength
        runningSum += this.windGini[ this.saInv[s] ]
    end
    this.spatGini[ this.saInv[1] ] =
            convert(Float32, runningSum / windSize)
    runningSum = fillSpatGini!(this.spatGini, this.spatialLength,
                               length(this.sa), this.windGini,
                               this.saInv, runningSum, windSize)
    for s = (length(this.sa)+2-this.spatialLength):length(this.sa)
        this.spatGini[ this.saInv[s] ] =
                convert(Float32, runningSum / windSize)
    end
    return nothing
end


## -----------------------------------------------------------------------------
function fillSpatialWindowed!(spatialWindowed::Array{Float32},
                              spatialLength::Int64, saLen::Int64,
                              windowed::Array{Float32}, saInv::Array{Int64},
                              runningSum::Float64, windSize::Float64)::Float64
    @inbounds for s = 2:(saLen+1-spatialLength)
        runningSum -= windowed[ saInv[s-1] ]
        runningSum += windowed[ saInv[s+spatialLength-1] ]
        spatialWindowed[ saInv[s] ] = convert(Float32, runningSum / windSize)
    end
    return runningSum
end

function spatialWindow!(this::AbSarks)::Nothing
    windSize = convert(Float64, this.spatialLength)
    this.spatialWindowed = Array{Float32}(undef, length(this.sa))
    runningSum = 0.0
    for s = 1:this.spatialLength
        runningSum += this.windowed[ this.saInv[s] ]
    end
    this.spatialWindowed[ this.saInv[1] ] =
            convert(Float32, runningSum / windSize)
    runningSum = fillSpatialWindowed!(this.spatialWindowed, this.spatialLength,
                                      length(this.sa), this.windowed,
                                      this.saInv, runningSum, windSize)
    for s = (length(this.sa)+2-this.spatialLength):length(this.sa)
        this.spatialWindowed[ this.saInv[s] ] =
                convert(Float32, runningSum / windSize)
    end
    return nothing
end


## -----------------------------------------------------------------------------
function findKmer(this::AbSarks, kmer::Array{UInt8})::Tuple{Int64,Int64}
    return FMIndices.findKmer(this.fm, kmer)
end

function kmer(this::AbSarks, s::Int64, k::Int64, k0::Int64=0)::Array{UInt8}
    return this.catSeq[max(1, s+k0):min(length(this.catSeq), s+k-1)]
end

function kmers(this::AbSarks, s::Array{Int64}, k::Int64,
               k0::Int64=0)::Array{Array{UInt8}}
    return [kmer(this, sel, k, k0) for sel in s]
end

function prefixAgreeSum(this::AbSarks, i::Int64, kmax::Int64)::Float64
    kmer_ = kmer(this, this.sa[i], kmax)
    kmax = min(kmax, length(kmer_))
    agreeSum = 0
    wstart = i - this.halfWindow
    wend = i + this.halfWindow + 1
    let kstart=i, kend=i+1
        for k = kmax:-1:1
            kwin = Sarkses.findKmer(this, kmer_[1:k])
            kwin1 = max(wstart, kwin[1])
            kwin2 = min(wend, kwin[2])
            agreeSum += k * ((kstart - kwin1) +
                             (kwin2 - kend))
            kstart = kwin1
            kend = kwin2
        end
    end
    return convert(Float64, agreeSum) / (2.0 * this.halfWindow)
end


## -----------------------------------------------------------------------------
function filterGuts(sa::Array{Int64}, saInv::Array{Int64},
                    windowed::Array{Float32}, windGini::Array{Float32},
                    spatialWindowed::Array{Float32}, spatGini::Array{Float32},
                    theta::Float64=-Inf, minGini::Float64=-Inf,
                    spatialTheta::Float64=-Inf, minSpatialGini::Float64=-Inf;
                    peakify::Bool=true)::Array{Int64}
    if minGini >= 1.0
        minGini = 1.0 - (1.0 - median(windGini)) * minGini
    end
    if minSpatialGini >= 1.0
        minSpatialGini = 1.0 - (1.0 - median(windGini)) * minSpatialGini
    end
    pos = Array{Int64}(undef, 0)
    keep = false
    iLeft = -1
    iRight = -1
    @inbounds for i = 1:length(sa)
        if windowed[i] >= theta && windGini[i] >= minGini
            keep = true
            s = sa[i]
            if spatialTheta > -Inf
                if (spatialWindowed[i] < spatialTheta ||
                    spatGini[i] < minSpatialGini)
                    keep = false
                end
            end
            if keep && peakify
                iLeft = saInv[1]
                if s > 1
                    iLeft = saInv[s-1]
                end
                iRight = saInv[length(saInv)]
                if s < length(saInv)
                    iRight = saInv[s+1]
                end
                if ((windowed[iLeft] > windowed[i]) ||
                    (windowed[iRight] > windowed[i]))
                    keep = false
                end
            end
            if keep
                push!(pos, i)
            end
        end
    end
    return pos
end

function filter(this::AbSarks, theta::Float64=-Inf, minGini::Float64=-Inf,
                spatialTheta::Float64=-Inf, minSpatialGini::Float64=-Inf;
                peakify::Bool=true)::Array{Int64}
    return filterGuts(this.sa, this.saInv, this.windowed, this.windGini,
                      this.spatialWindowed, this.spatGini, theta, minGini,
                      spatialTheta, minSpatialGini, peakify=peakify)
end

function filter(this::AbSarks,
                filters::Array{Dict{String,Float64},1},
                thresholds::Array{Float64,2}=Array{Float64,2}(undef, 0, 0);
                peakify::Bool=true)::Array{Array{Int64,1},1}
    hw0 = this.halfWindow
    sl0 = this.spatialLength
    out = Array{Array{Int64,1},1}(undef, 0)
    for f=1:length(filters)
        filt = filters[f]
        hw = convert(Int64, filt["halfWindow"])
        ming = ("minGini" in keys(filt) ? filt["minGini"] : -Inf)
        sl = convert(Int64, filt["spatialLength"])
        minsg = ("minSpatialGini" in keys(filt) ? filt["minSpatialGini"] : -Inf)
        theta = -Inf
        spatialTheta = -Inf
        if size(thresholds)[1] > 0
            theta = thresholds[f, 1]
            spatialTheta = thresholds[f, 2]
        elseif "theta" in keys(filt)
            theta = filt["theta"]
        elseif "spatialTheta" in keys(filt)
            spatialTheta = filt["spatialTheta"]
        end
        if hw != this.halfWindow
            reset!(this, hw, sl, true)
        elseif sl != this.spatialLength
            resetSpatial!(this, sl, true)
        end
        push!(out, filter(
            this, theta, ming,
            spatialTheta, minsg, peakify=peakify
        ))
    end
    if this.halfWindow != hw0
        reset!(this, hw0, sl0, true)
    elseif this.spatialLength != sl0
        resetSpatial!(this, sl0, true)
    end
    return out
end


## -----------------------------------------------------------------------------
function kmerPeaks(this::AbSarks,
                   filters::Array{Dict{String,Float64},1},
                   thresholds::Array{Float64,2},
                   iFilt::Array{Array{Int64,1},1};
                   kmax::Int64=12)::DataFrame
    hw0 = this.halfWindow
    sl0 = this.spatialLength
    out = DataFrame(
        i=Array{Int64}(undef,0), s=Array{Int64}(undef,0), kmer=Array{String}(undef,0),
        khat=Array{Float64}(undef,0), block=Array{String}(undef,0),
        wi=Array{Int64}(undef,0), gini=Array{Float32}(undef,0),
        spatialGini=Array{Float32}(undef,0), score=Array{Float64}(undef,0),
        windowed=Array{Float32}(undef,0), spatialWindowed=Array{Float32}(undef,0),
        kmax=Array{Int64}(undef,0),
        halfWindow=Array{Int64}(undef,0), minGini=Array{Float64}(undef,0),
        theta=Array{Float64}(undef,0), spatialLength=Array{Int64}(undef,0),
        minSpatialGini=Array{Float64}(undef,0), spatialTheta=Array{Float64}(undef,0)
    )
    for f = 1:length(iFilt)
        filt = filters[f]
        hw = convert(Int64, filt["halfWindow"])
        ming = ("minGini" in keys(filt) ? filt["minGini"] : -Inf)
        sl = convert(Int64, filt["spatialLength"])
        minsg = ("minSpatialGini" in keys(filt) ? filt["minSpatialGini"] : -Inf)
        if hw != this.halfWindow
            reset!(this, hw, sl, true)
        elseif sl != this.spatialLength
            resetSpatial!(this, sl, true)
        end
        for i in iFilt[f]
            s = this.sa[i]
            khat = prefixAgreeSum(this, i, kmax)
            km = kmer(this, s, kmax)[1:convert(Int64, round(khat))]
            blockIndex = sourceBlock(this, s)
            block = this.blocks[blockIndex]
            wi = s - this.bounds[blockIndex] + 1
            gini = this.windGini[i]
            sg = (length(this.spatGini) > 0 ? this.spatGini[i] : -Inf32)
            score = this.catScores[i]
            sw = (length(this.spatialWindowed) > 0 ?
                  this.spatialWindowed[i] : -Inf32)
            append!(out, Dict(
                "i"=>i, "s"=>s, "kmer"=>String(km), "khat"=>khat,
                "block"=>block, "wi"=>wi, "gini"=>gini, "spatialGini"=>sg,
                "score"=>score, "windowed"=>this.windowed[i],
                "spatialWindowed"=>sw, "kmax"=>kmax,
                "halfWindow"=>hw, "minGini"=>ming,
                "theta"=>thresholds[f, 1], "spatialLength"=>sl,
                "minSpatialGini"=>minsg, "spatialTheta"=>thresholds[f, 2]
            ))
        end
    end
    reset!(this, hw0, sl0, true)
    return out
end

function spatialSubPeaks(this::AbSarks,
                         iFilt::Array{Int64,1},
                         theta::Float64)::Array{Int64}
    if this.spatialLength <= 1
        return [this.sa[i] for i in iFilt]
    end
    subPeaks = Set{Int64}()
    for i in iFilt
        sRangeMax = min(this.sa[i] + this.spatialLength - 1, length(this.saInv))
        union!(subPeaks, collect(this.sa[i]:sRangeMax))
    end
    filteredOut = Set{Int64}()
    for s in subPeaks
        if this.windowed[this.saInv[s]] < theta
            push!(filteredOut, s)
        elseif (!((s-1) in subPeaks)) &&
               (this.windowed[this.saInv[s-1]] >= theta)
            push!(filteredOut, s)
        end
    end
    setdiff!(subPeaks, filteredOut)
    return sort(collect(subPeaks))
end

function mergeKmerIntervals(this::AbSarks,
                            subPeaks::Array{Int64};
                            kmax::Int64=12)::Array{Tuple{Int64,Int64},1}
    subPeaks = sort(subPeaks)
    left = 0; right = 0; sLast = -1
    out = Array{Tuple{Int64,Int64},1}(undef, 0)
    for s in subPeaks
        khat = convert(Int64, round(prefixAgreeSum(this, this.saInv[s], kmax)))
        if (s-1) == sLast
            right = min(max(right, s+khat-1), length(this.catSeq))
        else
            if left >= 1
                push!(out, (left, right))
            end
            left = s
            right = min(s+khat-1, length(this.catSeq))
        end
        sLast = s
    end
    if left >= 1
        push!(out, (left, right))
    end
    return out
end

function multiMergeKmerIntervals(this::AbSarks,
                                 iFilt::Array{Array{Int64,1},1},
                                 thresholds::Array{Float64,2},
                                 filters::Array{Dict{String,Float64},1};
                                 kmax::Int64)::Array{Array{Tuple{Int64,Int64},1},1}
    out = Array{Array{Tuple{Int64,Int64},1},1}(undef, 0)
    hw0 = this.halfWindow
    sl0 = this.spatialLength
    for f = 1:length(filters)
        filt = filters[f]
        hw = convert(Int64, filt["halfWindow"])
        sl = convert(Int64, filt["spatialLength"])
        if hw != this.halfWindow
            reset!(this, hw, sl, true)
        elseif sl != this.spatialLength
            resetSpatial!(this, sl, true)
        end
        subPeaks = spatialSubPeaks(this, iFilt[f], thresholds[f, 2])
        push!(out, mergeKmerIntervals(this, subPeaks, kmax=kmax))
    end
    reset!(this, hw0, sl0, true)
    return out
end

function mergedKmerSubPeaks(
        this::AbSarks,
        filters::Array{Dict{String,Float64},1},
        thresholds::Array{Float64,2};
        peakify::Bool=true,
        kmax::Int64=12
    )::DataFrame
    eyes = filter(this, filters, thresholds, peakify=peakify)
    mergedKmerIntervals = multiMergeKmerIntervals(
        this, eyes, thresholds, filters, kmax=kmax
    )
    hw0 = this.halfWindow
    sl0 = this.spatialLength
    out = DataFrame(
        i=Array{Int64}(undef,0), s=Array{Int64}(undef,0), kmer=Array{String}(undef,0),
        khat=Array{Float64}(undef,0), block=Array{String}(undef,0),
        wi=Array{Int64}(undef,0), gini=Array{Float32}(undef,0),
        spatialGini=Array{Float32}(undef,0), score=Array{Float64}(undef,0),
        windowed=Array{Float32}(undef,0), spatialWindowed=Array{Float32}(undef,0),
        kmax=Array{Int64}(undef,0),
        halfWindow=Array{Int64}(undef,0), minGini=Array{Float64}(undef,0),
        theta=Array{Float64}(undef,0), spatialLength=Array{Int64}(undef,0),
        minSpatialGini=Array{Float64}(undef,0), spatialTheta=Array{Float64}(undef,0)
    )
    for f = 1:length(filters)
        filt = filters[f]
        hw = convert(Int64, filt["halfWindow"])
        ming = ("minGini" in keys(filt) ? filt["minGini"] : -Inf)
        sl = convert(Int64, filt["spatialLength"])
        minsg = ("minSpatialGini" in keys(filt) ? filt["minSpatialGini"] : -Inf)
        theta = -Inf
        spatialTheta = -Inf
        if size(thresholds)[1] > 0
            theta = thresholds[f, 1]
            spatialTheta = thresholds[f, 2]
        elseif "theta" in keys(filt)
            theta = filt["theta"]
        elseif "spatialTheta" in keys(filt)
            spatialTheta = filt["spatialTheta"]
        end
        if hw != this.halfWindow
            reset!(this, hw, sl, true)
        elseif sl != this.spatialLength
            resetSpatial!(this, sl, true)
        end
        for sInterval in mergedKmerIntervals[f]
            i = this.saInv[sInterval[1]]
            km = this.catSeq[sInterval[1]:sInterval[2]]
            khat = prefixAgreeSum(this, i, kmax)
            blockIndex = sourceBlock(this, sInterval[1])
            block = this.blocks[blockIndex]
            wi = sInterval[1] - this.bounds[blockIndex] + 1
            gini = this.windGini[i]
            sg = (length(this.spatGini) > 0 ? this.spatGini[i] : -Inf32)
            score = this.scores[blockIndex]
            sw = (length(this.spatialWindowed) > 0 ?
                  this.spatialWindowed[i] : -Inf32)
            append!(out, Dict(
                "i"=>i, "s"=>sInterval[1], "kmer"=>String(km), "khat"=>khat,
                "block"=>block, "wi"=>wi, "gini"=>gini, "spatialGini"=>sg,
                "score"=>score, "windowed"=>this.windowed[i],
                "spatialWindowed"=>sw, "kmax"=>kmax,
                "halfWindow"=>hw, "minGini"=>ming,
                "theta"=>thresholds[f, 1], "spatialLength"=>sl,
                "minSpatialGini"=>minsg, "spatialTheta"=>thresholds[f, 2]
            ))
        end
    end
    reset!(this, hw0, sl0, true)
    return out
end

function mergedKmers(this::AbSarks,
                     mergedKmerIntervals::Array{Tuple{Int64,Int64},1}
                     )::Array{String}
    merged = Set{String}()
    for interval in mergedKmerIntervals
        push!(merged, this.catSeq[interval[1], interval[2]])
    end
    return sort(collect(merged))
end

function mergedKmers(this::AbSarks,
                     mergedKmerIntervals::Array{Array{Tuple{Int64,Int64},1}}
                     )::Array{Array{String}}
    return [mergedKmers(this, mergedKmerInts)
            for mergedKmerInts in mergedKmerIntervals]
end



## =============================================================================
mutable struct Sarks <: AbSarks
    halfWindow::Int64
    spatialLength::Int64

    scores::Array{Float64}
    blocks::Array{String}
    blockPosition::Dict{String,Int64}
    catSeq::Array{UInt8}
    bounds::Array{Int64}
    sa::Array{Int64}
    saInv::Array{Int64}
    windGini::Array{Float32}
    spatGini::Array{Float32}
    windowed::Array{Float32}
    spatialWindowed::Array{Float32}

    fm::FMIndices.FMIndex

    origScores::Array{Float64}
    charFreqs::Matrix{Float64}
    compositionPreBeta::Matrix{Float64}

    function Sarks(fasta::String, scoreFile::String,
                   halfWindow::Int64, spatialLength::Int64=0;
                   adjustForCharFreqs::Bool=false)
        this = new()
        this.halfWindow = halfWindow
        this.spatialLength = spatialLength
        seqs = readFasta(fasta)
        scoreMap = readScores(scoreFile)
        this.blocks = sort(intersect(
            collect(keys(scoreMap)),
            collect(keys(seqs))
        ))
        this.blockPosition = Dict{String,Int64}()
        for i in 1:length(this.blocks)
            this.blockPosition[this.blocks[i]] = i
        end
        this.origScores = [scoreMap[this.blocks[i]] for i in 1:length(this.blocks)]
        concatenateSeqs!(this, seqs)
        cf = characterFrequencies(this.bounds, this.catSeq)
        for i = 1:size(cf, 2)
            cf[:, i] .-= mean(cf[:, i])
        end
        if adjustForCharFreqs
            this.charFreqs = cf
            this.compositionPreBeta = inv(transpose(cf) * cf) * transpose(cf)
            deltaScores = (this.origScores .- mean(this.origScores))
            hatScores = cf * (this.compositionPreBeta * deltaScores)
            this.scores = deltaScores - hatScores
        else
            this.charFreqs = zeros(0, 0)
            this.compositionPreBeta = zeros(0, 0)
            this.scores = this.origScores
        end
        calcSuffixArray!(this)
        this.fm = FMIndices.FMIndex(this.catSeq, this.sa)
        windowGini!(this)
        window!(this)
        if this.spatialLength > 1
            spatialGini!(this)
            spatialWindow!(this)
        else
            this.spatGini = Array{Float32}(undef, 0)
            this.spatialWindowed = Array{Float32}(undef, 0)
        end
        return this
    end

    function Sarks(seqs::AbstractDict{<:AbstractString,<:AbstractString},
                   scoreMap::AbstractDict{<:AbstractString,Float64},
                   halfWindow::Int64, spatialLength::Int64=0;
                   adjustForCharFreqs::Bool=false)
        this = new()
        this.halfWindow = halfWindow
        this.spatialLength = spatialLength
        this.blocks = sort(intersect(
            collect(String.(keys(scoreMap))),
            collect(String.(keys(seqs)))
        ))
        this.blockPosition = Dict{String,Int64}()
        for i in 1:length(this.blocks)
            this.blockPosition[this.blocks[i]] = i
        end
        this.origScores = [scoreMap[this.blocks[i]] for i in 1:length(this.blocks)]
        concatenateSeqs!(this, seqs)
        cf = characterFrequencies(this.bounds, this.catSeq)
        for i = 1:size(cf, 2)
            cf[:, i] .-= mean(cf[:, i])
        end
        if adjustForCharFreqs
            this.charFreqs = cf
            this.compositionPreBeta = inv(transpose(cf) * cf) * transpose(cf)
            deltaScores = (this.origScores .- mean(this.origScores))
            hatScores = cf * (this.compositionPreBeta * deltaScores)
            this.scores = deltaScores - hatScores
        else
            this.charFreqs = zeros(0, 0)
            this.compositionPreBeta = zeros(0, 0)
            this.scores = this.origScores
        end
        calcSuffixArray!(this)
        this.fm = FMIndices.FMIndex(this.catSeq, this.sa)
        windowGini!(this)
        window!(this)
        if this.spatialLength > 1
            spatialGini!(this)
            spatialWindow!(this)
        else
            this.spatGini = Array{Float32}(undef, 0)
            this.spatialWindowed = Array{Float32}(undef, 0)
        end
        return this
    end

    function Sarks(sarks::Sarks)
        this = new()
        this.halfWindow = sarks.halfWindow
        this.spatialLength = sarks.spatialLength
        this.blocks = sarks.blocks
        this.blockPosition = sarks.blockPosition
        this.scores = sarks.scores
        this.catSeq = sarks.catSeq
        this.bounds = sarks.bounds
        this.origScores = sarks.origScores
        this.charFreqs = sarks.charFreqs
        this.compositionPreBeta = sarks.compositionPreBeta
        this.sa = sarks.sa
        this.saInv = sarks.saInv
        this.fm = sarks.fm
        this.windGini = sarks.windGini
        this.spatGini = sarks.spatGini
        this.windowed = sarks.windowed
        this.spatialWindowed = sarks.spatialWindowed
        return this
    end

    function Sarks(sarks::Sarks, permutation::Array{Int64})
        this = new()
        this.halfWindow = sarks.halfWindow
        this.spatialLength = sarks.spatialLength
        this.blocks = sarks.blocks
        this.blockPosition = sarks.blockPosition
        this.scores = [sarks.scores[ permutation[b] ]
                       for b in 1:length(sarks.scores)]
        this.catSeq = sarks.catSeq
        this.bounds = sarks.bounds
        this.origScores = [sarks.origScores[ permutation[b] ]
                           for b = 1:length(sarks.origScores)]
        this.charFreqs = sarks.charFreqs
        this.compositionPreBeta = sarks.compositionPreBeta
        if size(this.charFreqs, 1) > 0
            deltaScores = (this.origScores .- mean(this.origScores))
            hatScores = this.charFreqs * (this.compositionPreBeta * deltaScores)
            this.scores = deltaScores - hatScores
        else
            this.scores = this.origScores
        end
        this.sa = sarks.sa
        this.saInv = sarks.saInv
        this.fm = sarks.fm
        this.windGini = sarks.windGini
        this.spatGini = sarks.spatGini
        window!(this)
        if this.spatialLength > 1
            spatialWindow!(this)
        else
            this.spatialWindowed = Array{Float32}(undef, 0)
        end
        return this
    end
end

function readScores(filename::String)::DataStructures.OrderedDict{String,Float64}
    ordKeys = Array{String}(undef, 0)
    unord = Dict{String,Float64}()
    fgz = (endswith(lowercase(filename), ".gz") ||
           endswith(lowercase(filename), ".gzip"))
    fh = (fgz ?
          TranscodingStream(GzipDecompressor(), open(filename, "r")) :
          open(filename, "r"))
    while !eof(fh)
        tokens = split(rstrip(readline(fh)), "\t")
        score = tryparse(Float64, tokens[2])
        if score !== nothing
            unord[tokens[1]] = score
            push!(ordKeys, tokens[1])
        end
    end
    out = DataStructures.OrderedDict{String,Float64}()
    for key in ordKeys
        out[key] = unord[key]
    end
    return out
end

function characterFrequencies(bounds::Array{Int64},
                              catSeq::Array{UInt8})::Matrix{Float64}
    alphabet = setdiff(unique(catSeq), UInt8('$'))
    charFreqs = zeros(length(bounds)-1, length(alphabet)-1)
    @inbounds for i = 1:size(charFreqs, 1)
        @inbounds for j = 1:size(charFreqs, 2)
            charFreqs[i, j] =
                    sum(catSeq[bounds[i]:(bounds[i+1]-1)] .== alphabet[j]) /
                    (bounds[i+1] - bounds[i] - 1)
            ## the -1 above accounts for the end-of-sequence marker ($)
        end
    end
    return charFreqs
end


## -----------------------------------------------------------------------------
function sourceScore(this::Sarks, s::Int64)::Float64
    return this.scores[sourceBlock(this, s)]
end

function sourceScore(this::Sarks, s::Array{Int64})::Array{Float64}
    return [this.scores[b] for b in sourceBlock(this, s)]
end

function fillSourceScore!(ss::Array{Float64},
                          scores::Array{Float64},
                          srcBlk::Array{Int64})::Nothing
    @inbounds for i = 1:length(srcBlk)
        ss[i] = scores[ srcBlk[i] ]
    end
    return nothing
end

function sourceScore(this::Sarks)::Array{Float64}
    srcBlk = sourceBlock(this)
    out = Array{Float64}(undef, length(srcBlk))
    fillSourceScore!(out, this.scores, srcBlk)
    return out
end


## -----------------------------------------------------------------------------
function blockCounts(this::Sarks, i::Int64, block::Array{Int64})::Array{Int64}
    out = zeros(Int64, length(this.blocks))
    for j = -this.halfWindow:this.halfWindow
        out[ block[i+j] ] += 1
    end
    return out
end

function fillWindGini!(windGini::Array{Float32},
                       halfWindow::Int64, saLen::Int64,
                       block::Array{Int64}, bCounts::Array{Int64},
                       runningGini::Float64)::Nothing
    total = (2 * halfWindow) + 1
    total2 = convert(Float64, total * total)
    dGini = 0
    oldb = -1
    newb = -1
    @inbounds for i = (halfWindow+2):(saLen-halfWindow)
        oldb = block[i - 1 - halfWindow]
        newb = block[i + halfWindow]
        dGini = -bCounts[oldb] * (total - bCounts[oldb])
        if bCounts[newb] > 0
            dGini -= bCounts[newb] * (total - bCounts[newb])
        end
        bCounts[oldb] -= 1
        bCounts[newb] += 1
        if bCounts[oldb] > 0
            dGini += bCounts[oldb] * (total - bCounts[oldb])
        end
        dGini += bCounts[newb] * (total - bCounts[newb])
        runningGini += (convert(Float64, dGini) / total2)
        windGini[i] = convert(Float32, runningGini)
    end
end

function windowGini!(this::Sarks)::Nothing
    block = sourceBlock(this)
    this.windGini = Array{Float32}(undef, length(this.sa))
    bCounts = blockCounts(this, this.halfWindow+1, block)
    runningGini = giniImpurity(bCounts)
    this.windGini[this.halfWindow+1] = convert(Float32, runningGini)
    fillWindGini!(this.windGini, this.halfWindow, length(this.sa),
                  block, bCounts, runningGini)
    for i = 1:this.halfWindow
        this.windGini[i] = this.windGini[this.halfWindow+1]
        this.windGini[length(this.sa) + 1 - i] =
                this.windGini[length(this.sa) - this.halfWindow]
    end
    return nothing
end


## -----------------------------------------------------------------------------
function fillWindowed!(windowed::Array{Float32},
                       halfWindow::Int64, saLen::Int64,
                       scores::Array{Float64}, srcBlk::Array{Int64},
                       runningSum::Float64, windSize=Float64)::Float64
    @inbounds for i = (halfWindow+2):(saLen-halfWindow)
        runningSum -= scores[ srcBlk[i-1-halfWindow] ]
        runningSum += scores[ srcBlk[i+halfWindow] ]
        windowed[i] = convert(Float32, runningSum / windSize)
    end
    return runningSum
end

function window!(this::Sarks)::Nothing
    windSize = convert(Float64, (2 * this.halfWindow) + 1)
    srcBlk = sourceBlock(this)
    this.windowed = Array{Float32}(undef, length(this.sa))
    runningSum = 0.0
    for i = 1:((2*this.halfWindow)+1)
        runningSum += this.scores[ srcBlk[i] ]
    end
    this.windowed[this.halfWindow + 1] = convert(Float32, runningSum / windSize)
    runningSum = fillWindowed!(this.windowed, this.halfWindow,
                               length(this.sa), this.scores, srcBlk,
                               runningSum, windSize)
    for i = 1:this.halfWindow
        this.windowed[i] = this.windowed[this.halfWindow+1]
        this.windowed[length(this.sa) + 1 - i] =
                convert(Float32, runningSum / windSize)
    end
    return nothing
end


## -----------------------------------------------------------------------------
function giniMask(this::Sarks,
                  minGini::Float64,
                  minSpatialGini::Float64)::AbstractArray{Bool}
    if minGini >= 1
        minGini = 1 - (1-median(this.windGini))*minGini
    end
    if minSpatialGini >= 1
        minSpatialGini = 1 - (1-median(this.windGini))*minSpatialGini
    end
    mask = trues(length(this.windGini))
    if minGini > 0
        for i = 1:length(mask)
            if this.windGini[i] < minGini
                mask[i] = false
            end
        end
    end
    if minSpatialGini > 0 && this.spatialLength > 1
        for i = 1:length(mask)
            if this.spatGini[i] < minSpatialGini
                mask[i] = false
            end
        end
    end
    return mask
end

function permutationDistribution(
        this::Sarks,
        reps::Int64,
        filters::Array{Dict{String,Float64},1};
        seed::Int64=0,
        permutations::Array{Int64,2}=Array{Int64,2}(undef, 0, 0)
    )::Array{Float32,3}
    rng = (seed != 0 ? MersenneTwister(seed) : MersenneTwister())
    if size(permutations)[1] == 0
        permutations = Array{Int64,2}(undef, reps, length(this.scores))
        for r = 1:reps
            permutations[r, :] = randperm(rng, length(this.scores))
        end
    end
    hw0 = this.halfWindow
    sl0 = this.spatialLength
    permResults = Array{Float32,3}(undef, reps, length(filters), 2)
    for f = 1:length(filters)
        filt = filters[f]
        hw = convert(Int64, filt["halfWindow"])
        ming = ("minGini" in keys(filt) ? filt["minGini"] : -Inf)
        sl = convert(Int64, filt["spatialLength"])
        minsg = ("minSpatialGini" in keys(filt) ? filt["minSpatialGini"] : -Inf)
        if hw != this.halfWindow
            reset!(this, hw, sl, false)
        elseif sl != this.spatialLength
            resetSpatial!(this, sl, false)
        end
        fmask = giniMask(this, ming, minsg)
        Threads.@threads for r = 1:reps
            rsarks = Sarks(this, permutations[r, :])
            if length(rsarks.spatialWindowed) == 0
                permResults[r, f, 1] = maximum(rsarks.windowed[fmask])
                permResults[r, f, 2] = -Inf32
            else
                permResults[r, f, 1] = -Inf32
                permResults[r, f, 2] = maximum(rsarks.spatialWindowed[fmask])
            end
        end
    end
    reset!(this, hw0, sl0, true)
    return permResults
end

function thresholdsFromPermutations(
        permDists::Array{Float32,3}, nSigma::Float64)::Array{Float64,2}
    thresholds = Array{Float64,2}(undef, size(permDists)[2], 2)
    for f in 1:size(permDists)[2]
        if any(permDists[:, f, 1] .> -Inf32)
            thresholds[f, 1] = mean(permDists[:, f, 1]) +
                               nSigma * std(permDists[:, f, 1])
            thresholds[f, 2] = -Inf
        else
            thresholds[f, 1] = -Inf
            thresholds[f, 2] = mean(permDists[:, f, 2]) +
                               nSigma * std(permDists[:, f, 2])
        end
    end
    return thresholds
end


## =============================================================================
mutable struct MicroSarks <: AbSarks
    halfWindow::Int64
    spatialLength::Int64
    subblockLength::Int64

    catScores::Array{Float32}
    blocks::Array{String}
    blockPosition::Dict{String,Int64}
    catSeq::Array{UInt8}
    bounds::Array{Int64}
    sa::Array{Int64}
    saInv::Array{Int64}
    windGini::Array{Float32}
    spatGini::Array{Float32}
    windowed::Array{Float32}
    spatialWindowed::Array{Float32}

    fm::FMIndices.FMIndex

    function MicroSarks(fasta::String, scoreFile::String,
                        halfWindow::Int64, spatialLength::Int64=0;
                        subblockLength::Int64=100,
                        shift::Int64=0)
        this = new()
        this.halfWindow = halfWindow
        this.spatialLength = spatialLength
        this.subblockLength = subblockLength
        seqs = readFasta(fasta)
        scores = DataFrame()
        if endswith(lowercase(scoreFile), r".gz(ip)?")
            scores = CSV.read(TranscodingStream(GzipDecompressor(),
                                                open(scoreFile, "r")),
                              DataFrame, delim="\t", comment="#",
                              header=false, missingstring="NA", copycols=true)
        else
            scores = CSV.read(scoreFile, DataFrame, delim="\t", comment="#",
                              header=false, missingstring="NA", copycols=true)
        end
        scores[!, 2] = Int64.(scores[:, 2]) .+ shift
        scores[!, 3] = Float64.(scores[:, 3])
        this.blocks = sort(collect(intersect(Set(scores[:, 1]), keys(seqs))))
        this.blockPosition = Dict{String,Int64}()
        for i in 1:length(this.blocks)
            this.blockPosition[this.blocks[i]] = i
        end
        toRemove = setdiff(keys(seqs), this.blocks)
        for k in toRemove
            delete!(seqs, k)
        end
        if shift != 0
            scores = scores[scores[:, 2] .> 0, :]
            for k in keys(seqs)
                kst = max(1, shift+1)
                ken = length(seqs[k]) + min(0, shift)
                seqs[k] = seqs[k][kst:ken]
            end
        end
        concatenateSeqs!(this, seqs)
        calcSuffixArray!(this)
        this.fm = FMIndices.FMIndex(this.catSeq, this.sa)
        concatenateScores!(this, scores[:, 1], scores[:, 2], scores[:, 3])
        windowGini!(this)
        window!(this)
        if this.spatialLength > 1
            spatialGini!(this)
            spatialWindow!(this)
        else
            this.spatGini = Array{Float32}(undef, 0)
            this.spatialWindowed = Array{Float32}(undef, 0)
        end
        return this
    end

    function MicroSarks(seqs::AbstractDict{<:AbstractString,<:AbstractString},
                        scores::DataFrame,
                        halfWindow::Int64, spatialLength::Int64=0;
                        subblockLength::Int64=100,
                        shift::Int64=0)
        this = new()
        this.halfWindow = halfWindow
        this.spatialLength = spatialLength
        this.subblockLength = subblockLength
        this.blocks = String.(sort(collect(intersect(Set(scores[:, 1]), keys(seqs)))))
        this.blockPosition = Dict{String,Int64}()
        for i in 1:length(this.blocks)
            this.blockPosition[this.blocks[i]] = i
        end
        toRemove = setdiff(keys(seqs), this.blocks)
        for k in toRemove
            delete!(seqs, k)
        end
        if shift != 0
            scores = deepcopy(scores)
            scores[:, 2] .+= shift
            scores = scores[scores[:, 2] .> 0, :]
            for k in keys(seqs)
                kst = max(1, shift+1)
                ken = length(seqs[k]) + max(0, shift)
                seqs[k] = seqs[k][kst:ken]
            end
        end
        concatenateSeqs!(this, seqs)
        calcSuffixArray!(this)
        this.fm = FMIndices.FMIndex(this.catSeq, this.sa)
        concatenateScores!(this, scores[:, 1], scores[:, 2], scores[:, 3])
        windowGini!(this)
        window!(this)
        if this.spatialLength > 1
            spatialGini!(this)
            spatialWindow!(this)
        else
            this.spatGini = Array{Float32}(undef, 0)
            this.spatialWindowed = Array{Float32}(undef, 0)
        end
        return this
    end

    function MicroSarks(sarks::MicroSarks)
        this = new()
        this.halfWindow = sarks.halfWindow
        this.spatialLength = sarks.spatialLength
        this.subblockLength = sarks.subblockLength
        this.blocks = sarks.blocks
        this.blockPosition = sarks.blockPosition
        this.bounds = sarks.bounds
        this.catSeq = sarks.catSeq
        this.catScores = sarks.catScores
        this.sa = sarks.sa
        this.saInv = sarks.saInv
        this.fm = sarks.fm
        this.windGini = sarks.windGini
        this.spatGini = sarks.spatGini
        this.windowed = sarks.windowed
        this.spatialWindowed = sarks.spatialWindowed
        return this
    end
end


## -----------------------------------------------------------------------------
function concatenateScores!(this::MicroSarks,
                            seq::AbstractArray{<:AbstractString},
                            pos::Array{Int64},
                            score::Array{Float64}
                            )::Nothing
    this.catScores = zeros(Float32, length(this.catSeq))
    for r in 1:length(seq)
        rb = this.blockPosition[String(seq[r])]
        rs = this.bounds[rb] + pos[r] - 1
        if rs < this.bounds[rb+1]
            this.catScores[ this.saInv[rs] ] = score[r]
        end
    end
    return nothing
end


## -----------------------------------------------------------------------------
function fillWindGini!(windGini::Array{Float32},
                       halfWindow::Int64, sa::Array{Int64},
                       bounds::Array{Int64},
                       block::Array{Int64}, bCounts::Dict{Int64,Int64},
                       subblockLength::Int64, runningGini::Float64)::Nothing
    total = (2 * halfWindow) + 1
    total2 = convert(Float64, total * total)
    dGini = 0
    oldb = -1
    newb = -1
    saLen = length(sa)
    nBlocks = length(bounds) - 1
    @inbounds for i = (halfWindow+2):(saLen-halfWindow)
        oldb = block[i - 1 - halfWindow]
        oldsb = 1 + ((sa[i - 1 - halfWindow] - bounds[oldb]) ÷ subblockLength)
        oldkey = (oldsb * nBlocks) + oldb
        newb = block[i + halfWindow]
        newsb = 1 + ((sa[i + halfWindow] - bounds[newb]) ÷ subblockLength)
        newkey = (newsb * nBlocks) + newb
        dGini = -bCounts[oldkey] * (total - bCounts[oldkey])
        if !(newkey in keys(bCounts))
            bCounts[newkey] = 0
        end
        if bCounts[newkey] > 0
            dGini -= bCounts[newkey] * (total - bCounts[newkey])
        end
        bCounts[oldkey] -= 1
        bCounts[newkey] += 1
        if bCounts[oldkey] > 0
            dGini += bCounts[oldkey] * (total - bCounts[oldkey])
        end
        dGini += bCounts[newkey] * (total - bCounts[newkey])
        runningGini += (convert(Float64, dGini) / total2)
        windGini[i] = convert(Float32, runningGini)
    end
end

function blockCounts(this::MicroSarks, i::Int64, block::Array{Int64})::Dict{Int64,Int64}
    out = Dict{Int64,Int64}()
    nBlocks = length(this.bounds) - 1
    for j = -this.halfWindow:this.halfWindow
        jb = block[i+j]
        jsb = 1 + ((this.sa[i+j] - this.bounds[jb]) ÷ this.subblockLength)
        jkey = (jsb * nBlocks) + jb
        if !(jkey in keys(out))
            out[jkey] = 0
        end
        out[jkey] += 1
    end
    return out
end

function windowGini!(this::MicroSarks)::Nothing
    block = sourceBlock(this)
    this.windGini = Array{Float32}(undef, length(this.sa))
    bCounts = blockCounts(this, this.halfWindow+1, block)
    runningGini = giniImpurity(collect(values(bCounts)))
    this.windGini[this.halfWindow+1] = convert(Float32, runningGini)
    fillWindGini!(this.windGini, this.halfWindow, this.sa, this.bounds,
                  block, bCounts, this.subblockLength, runningGini)
    for i = 1:this.halfWindow
        this.windGini[i] = this.windGini[this.halfWindow+1]
        this.windGini[length(this.sa) + 1 - i] =
                this.windGini[length(this.sa) - this.halfWindow]
    end
    return nothing
end


## -----------------------------------------------------------------------------
function fillWindowed!(windowed::Array{Float32},
                       halfWindow::Int64, saLen::Int64,
                       scores::Array{Float32},
                       runningSum::Float64, windSize=Float64)::Float64
    @inbounds for i = (halfWindow+2):(saLen-halfWindow)
        runningSum -= scores[i-(halfWindow+1)]
        runningSum += scores[i+halfWindow]
        windowed[i] = convert(Float32, runningSum / windSize)
    end
    return runningSum
end

function window!(this::MicroSarks)::Nothing
    windSize = convert(Float64, (2 * this.halfWindow) + 1)
    srcBlk = sourceBlock(this)
    this.windowed = Array{Float32}(undef, length(this.sa))
    runningSum = 0.0
    for i = 1:((2*this.halfWindow)+1)
        runningSum += this.catScores[i]
    end
    this.windowed[this.halfWindow + 1] = convert(Float32, runningSum / windSize)
    runningSum = fillWindowed!(this.windowed, this.halfWindow,
                               length(this.sa), this.catScores,
                               runningSum, windSize)
    for i = 1:this.halfWindow
        this.windowed[i] = this.windowed[this.halfWindow+1]
        this.windowed[length(this.sa) + 1 - i] =
                convert(Float32, runningSum / windSize)
    end
    return nothing
end



## =============================================================================
function rescore(sarks::Sarks,
                 newscores::AbstractDict{<:AbstractString,Float64}
                 )::Sarks
    this = Sarks(sarks)
    this.scores = [newscores[b] for b in this.blocks]
    this.origScores = [newscores[b] for b in this.blocks]
    if size(this.charFreqs, 1) > 0
        deltaScores = (this.origScores .- mean(this.origScores))
        hatScores = this.charFreqs * (this.compositionPreBeta * deltaScores)
        this.scores = deltaScores - hatScores
    else
        this.scores = this.origScores
    end
    Sarkses.window!(this)
    if this.spatialLength > 1
        Sarkses.spatialWindow!(this)
    else
        this.spatialWindowed = Array{Float32}(undef, 0)
    end
    return this
end

function multisarks(this::Sarks,
                    scores::Array{OrderedDict{<:AbstractString,Float64}}
                    )::Array{Sarks}
    out = Array{Sarks}(undef, length(scores))
    Threads.@threads for i in 1:length(scores)
        out[i] = rescore(this, scores[i])
    end
    return out
end

## 210530 DW: not clear centerDicts still needed...
function centerDicts(x::Array{OrderedDict{<:AbstractString,Float64}}
                     )::Array{OrderedDict{String,Float64}}
    out = Array{OrderedDict{String,Float64}}(undef, length(x))
    for i in 1:length(x)
        outi = OrderedDict{String,Float64}()
        xi = x[i]
        ximu = mean(values(xi))
        for k in keys(xi)
            outi[String(k)] = xi[k] - ximu
        end
        out[i] = outi
    end
    return out
end

function multiWindCovar(sarkses::Array{Sarks},
                        minGini::Float64=-Inf
                        )::Matrix{Float64}
    if minGini > 0.0
        gm = giniMask(sarkses[1], minGini, -Inf)
    end
    out = zeros(length(sarkses), length(sarkses))
    for i in 1:length(sarkses)
        for j in i:length(sarkses)
            if minGini > 0.0
                out[i, j] = sum(
                    (sarkses[i].windowed[gm] .- mean(sarkses[i].windowed[gm])) .*
                    (sarkses[j].windowed[gm] .- mean(sarkses[j].windowed[gm]))
                )
            else
                out[i, j] = sum(
                    (sarkses[i].windowed .- mean(sarkses[i].windowed)) .*
                    (sarkses[j].windowed .- mean(sarkses[j].windowed))
                )
            end
            if j > i
                out[j, i] = out[i, j]
            end
        end
    end
    return out
end

function eigenScores(sarkses::Array{Sarks},
                     minGini::Float64=-Inf;
                     retainFirst::Bool=false
                     )::Array{OrderedDict{String,Float64}}
    scovar = multiWindCovar(sarkses, minGini)
    ev = zeros(size(scovar, 1), size(scovar, 2))
    if retainFirst
        ev[1, end] = 1.0
        ev[2:end, 1:(end-1)] = eigen(scovar[2:end, 2:end]).vectors
    else
        ev = eigen(scovar).vectors
    end
    out = Array{OrderedDict{String,Float64}}(undef, length(sarkses))
    for i in 1:length(out)
        evi = ev[:, length(sarkses)-i+1]
        out[i] = OrderedDict([
            k => sum([evi[j] * sarkses[j].scores[b] for j in 1:length(evi)])
            for (b, k) in enumerate(sarkses[1].blocks)
        ]);
    end
    return out
end



## =============================================================================
function spatialMeansInner(bounds::Array{Int64},
                           saInv::Array{Int64},
                           windowed::Array{Float32},
                           filt::AbstractArray{Bool})::Array{Float64}
    blockMeans = zeros(length(bounds)-1)
    activeBlock = 1
    blockLength = 0.0
    @inbounds for s in 1:length(saInv)
        if s >= bounds[activeBlock+1]
            blockMeans[activeBlock] /= blockLength
            activeBlock += 1
            blockLength = 0.0
        end
        if filt[saInv[s]]
            blockMeans[activeBlock] += windowed[saInv[s]]
            blockLength += 1.0
        end
    end
    blockMeans[activeBlock] /= blockLength
    return blockMeans
end

function spatialMeans(this::Sarks,
                      minGini::Float64=-Inf)::OrderedDict{String,Float64}
    if minGini >= 1.0
        minGini = 1.0 - (1.0 - median(this.windGini)) * minGini
    end
    blockMeans = spatialMeansInner(this.bounds, this.saInv, this.windowed,
                                   this.windGini .>= minGini)
    out = OrderedDict{String,Float64}()
    for (b, block) in enumerate(this.blocks)
        out[block] = blockMeans[b]
    end
    return out
end

mutable struct SarksScorer
    sarks::Sarks
    suffixLen::Int64
    spatial::Bool
    function SarksScorer(sarks::Sarks,
                         suffixLen::Int64;
                         spatial::Bool=false)
        this = new()
        this.sarks = sarks
        this.suffixLen = suffixLen
        this.spatial = spatial
        return this
    end
end

function score(this::SarksScorer, seq::Array{UInt8},
               minGini::Float64=-Inf; naOmit=true)::Array{Float64}
    if minGini >= 1.0
        minGini = 1.0 - (1.0 - median(this.sarks.windGini)) * minGini
    end
    scores = zeros(min(length(seq), length(seq)-this.suffixLen+1))
    Threads.@threads for s = 1:length(scores)
        suffixEnd = min(s+this.suffixLen-1, length(seq))
        suffix = seq[s:suffixEnd]
        limits = Sarkses.findKmer(this.sarks, suffix)
        lower = min(limits[1], limits[2]-1)
        upper = min(max(limits[1], limits[2]-1), length(this.sarks.sa))
        num = 0.0
        den = 0.0
        for i = lower:upper
            if this.sarks.windGini[i] >= minGini
                if this.spatial
                    num += this.sarks.spatialWindowed[i]
                else
                    num += this.sarks.windowed[i]
                end
                den += 1.0
            end
        end
        scores[s] = num / den
    end
    if naOmit
        return scores[.!isnan.(scores)]
    else
        return scores
    end
end

function scoreWithGini(this::SarksScorer,
                       seq::Array{UInt8})::Tuple{Array{Float64},Array{Float64}}
    scores = zeros(min(length(seq), length(seq)-this.suffixLen+1))
    ginis = zeros(length(scores))
    Threads.@threads for s = 1:length(scores)
        suffixEnd = min(s+this.suffixLen-1, length(seq))
        suffix = seq[s:suffixEnd]
        limits = Sarkses.findKmer(this.sarks, suffix)
        lower = min(limits[1], limits[2]-1)
        upper = min(max(limits[1], limits[2]-1), length(this.sarks.sa))
        sScore = 0.0
        sGini = 0.0
        for i = lower:upper
            if this.spatial
                sScore += this.sarks.spatialWindowed[i]
                sGini += this.sarks.spatGini[i]
            else
                sScore += this.sarks.windowed[i]
                sGini += this.sarks.windGini[i]
            end
        end
        scores[s] = sScore / (1+upper-lower)
        ginis[s] = sGini / (1+upper-lower)
    end
    return (scores, ginis)
end

function sarksPositionalScore(this::Sarks,
                              seqs::OrderedDict{<:AbstractString,Array{UInt8,1}};
                              k::Int64=25)::DataFrame
    spat = max(this.spatialLength, 1)
    scorer = SarksScorer(this, k, spatial=(spat > 1))
    nSeqBins = OrderedDict{String,Int64}()
    for sk in keys(seqs)
        nSeqBins[String(sk)] = Int64(floor((length(seqs[sk])+1-k) / spat))
    end
    nTotalBins = sum(values(nSeqBins))
    id = Array{String}(undef, nTotalBins)
    bin = Array{Int64}(undef, nTotalBins)
    score = Array{Float64}(undef, nTotalBins)
    gini = Array{Float64}(undef, nTotalBins)
    offset = 0
    offsets = OrderedDict([sk => 0 for sk in keys(seqs)])
    for sk in keys(seqs)
        offsets[sk] = offset
        offset += nSeqBins[String(sk)]
    end
    seqKeys = collect(values(keys(seqs)))
    Threads.@threads for sk in seqKeys
        skOffset = offsets[sk]
        skScores = scoreWithGini(scorer, seqs[sk])
        for i in 1:nSeqBins[String(sk)]
            iOffset = i + skOffset
            id[iOffset] = String(sk)
            bin[iOffset] = i
            score[iOffset] = skScores[1][1 + (i-1)*spat]
            gini[iOffset] = skScores[2][1 + (i-1)*spat]
        end
    end
    out = DataFrame(id=id, bin=bin, score=score, gini=gini)
    return out
end

function sarksPositionalScore(this::Sarks,
                              seqs::OrderedDict{<:AbstractString,<:AbstractString};
                              k::Int64=25)::DataFrame
    return sarksPositionalScore(
        this,
        OrderedDict([k => Array{UInt8}(v*"\$") for (k, v) in seqs]),
        k = k
    )
end

function sarksPositionalScoreMat(this::Sarks,
                                 seqs::OrderedDict{<:AbstractString,<:AbstractString};
                                 k::Int64=25,
                                 fill::Symbol=:global
                                 )::Matrix
    out = sarksPositionalScore(this, seqs, k=k)
    if fill == :individual
        return Matrix(c1idx(padAndUnstack(out), collect(keys(seqs)))[:, 2:end])
    else
        return Matrix(c1idx(padAndUnstack(out, fill=mean(this.scores)),
                            collect(keys(seqs)))[:, 2:end])
    end
end

function predictSarks(sarks::Sarks,
                      seqs::OrderedDict{<:AbstractString,<:AbstractString};
                      k::Int64=25,
                      fill::Symbol=:global
                      )::Tuple{Matrix{Float32},Matrix{Float32}}
    posScores = sarksPositionalScore(
        sarks,
        OrderedDict([id => Array{UInt8}(seqs[id]*"\$") for id in keys(seqs)]),
        k = k
    )
    if fill == :individual
        posScoreMat = Matrix{Float32}(c1idx(
            padAndUnstack(posScores),
            collect(keys(seqs))
        )[:, 2:end])
        posGiniMat = Matrix{Float32}(c1idx(
            padAndUnstack(posScores, :id, :bin, :gini),
            collect(keys(seqs))
        )[:, 2:end])
        return (posScoreMat, posGiniMat)
    else
        posScoreMat = Matrix{Float32}(c1idx(
            padAndUnstack(posScores,
                          fill=mean(sarks.scores)),
            collect(keys(seqs))
        )[:, 2:end])
        posGiniMat = Matrix{Float32}(c1idx(
            padAndUnstack(posScores, :id, :bin, :gini,
                          fill=Float64(mean(sarks.windGini))),
            collect(keys(seqs))
        )[:, 2:end])
        return (posScoreMat, posGiniMat)
   end
end



## =============================================================================
mutable struct SarksBoost
    learningRate::Float64
    minGini::Float64
    sarkses::Array{Sarks}
    coef::Array{Float64}
    intercept::Float64
end

function simpleLinReg(x::Vector{Float64}, y::Vector{Float64})::Tuple{Float64, Float64}
    ymu = mean(y)
    dx = x .- mean(x)
    dy = y .- ymu
    slope = sum(dx .* dy) / sum(dx .* dx)
    intercept = ymu - slope * mean(x)
    return (intercept, slope)
end

function boostSarks(iterations::Int64,
                    seqs::AbstractDict{<:AbstractString,<:AbstractString},
                    scoreMap::AbstractDict{<:AbstractString,Float64},
                    halfWindow::Int64,
                    minGini::Float64=-Inf;
                    shrinkage::Float64=0.0,
                    spatialLength::Int64=0,
                    adjustForCharFreqs::Bool=false)::SarksBoost
    lr = 1.0 - shrinkage
    sarks1 = Sarks(seqs, scoreMap, halfWindow, spatialLength,
                   adjustForCharFreqs = adjustForCharFreqs)
    ordScores = [scoreMap[b] for b in sarks1.blocks]
    ordSeqs = [String(seqs[b]) for b in sarks1.blocks]
    resids0 = ordScores .- mean(ordScores)
    spatMeans1 = collect(values(spatialMeans(sarks1, minGini)))
    sarks1Mod = simpleLinReg(spatMeans1, resids0)
    out = SarksBoost(lr, minGini,
                     [sarks1], [lr*sarks1Mod[2]],
                     mean(ordScores) + lr*sarks1Mod[1])
    if iterations > 1
        pseudoResids = resids0 - lr * (sarks1Mod[1] .+ sarks1Mod[2]*spatMeans1)
        for iter = 2:iterations
            sarksIter = rescore(
                sarks1,
                OrderedDict([b => pseudoResids[i]
                             for (i, b) in enumerate(sarks1.blocks)])
            )
            spatMeansIter = collect(values(spatialMeans(sarksIter, minGini)))
            iterMod = simpleLinReg(spatMeansIter, pseudoResids)
            push!(out.sarkses, sarksIter)
            push!(out.coef, lr*iterMod[2])
            out.intercept += lr*iterMod[1]
            pseudoResids -= lr * (iterMod[1] .+ iterMod[2]*spatMeansIter)
        end
    end
    return out
end

function sarksPositionalScoreMat(this::SarksBoost,
                                 seqs::OrderedDict{<:AbstractString,<:AbstractString};
                                 k::Int64)::Matrix
    out = this.coef[1] * sarksPositionalScoreMat(this.sarkses[1], seqs, k=k);
    for i = 2:length(this.sarkses)
        out += this.coef[i] * sarksPositionalScoreMat(this.sarkses[i], seqs, k=k);
    end
    out .+= this.intercept
    return out
end


end # module Sarkses
