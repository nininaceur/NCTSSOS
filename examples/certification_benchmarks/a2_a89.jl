
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using DynamicPolynomials
using NCTSSOS
import Base.Filesystem: basename, splitext
function load_bell_matrices(dir::AbstractString)
    files = sort(filter(f -> occursin(r"^A\d+\.txt$", f), readdir(dir)))
    out   = NamedTuple[]
    for f in files
        lines = filter(!isempty, strip.(readlines(joinpath(dir, f))))
        n = parse(Int, replace(lines[1], r"[^\d]" => ""))
        m = parse(Int, replace(lines[2], r"[^\d]" => ""))
        push!
        M = zeros(Float64, n + 1, m + 1)
        _fill_matrix!(M, lines[3])                         
        push!(out, (name = splitext(basename(f))[1], n = n, m = m, M = M))
    end
    return out
end

function _fill_matrix!(M::AbstractMatrix, expr::AbstractString)
    expr = replace(expr, r"\s+" => "")     
    expr = replace(expr, '*' => "")         
    expr = replace(expr, '-' => "+-")     
    for tok in filter(!isempty, split(expr, '+'))
        coeff, vars = _split_coeff(tok)
        xi, yj = 0, 0                       
        i = 1
        while i ≤ lastindex(vars)
            sym = vars[i]; i += 1
            idx, i = _read_int(vars, i)
            if sym == 'X'; xi = idx
            elseif sym == 'Y'; yj = idx
            end
        end
        M[xi + 1, yj + 1] += coeff       
    end
end

function _split_coeff(tok::AbstractString)
    k = findfirst(c -> c == 'X' || c == 'Y', tok)
    if k === nothing                   
        return (_parse_num(tok), "")
    else
        return (_parse_num(tok[1:k-1]), tok[k:end])
    end
end

function _parse_num(s::AbstractString)
    isempty(s)  && return 1.0
    s == "-"    && return -1.0
    num = match(r"^[+\-]?\d+(?:\.\d+)?(?:/\d+)?", s)
    str = num === nothing ? "1" : num.match
    occursin('/', str) ?
        let p = split(str, '/'); parse(Float64, p[1]) / parse(Float64, p[2]) end :
        parse(Float64, str)
end

function _read_int(s::AbstractString, i::Int)
    j = i
    while j ≤ lastindex(s) && isdigit(s[j]); j += 1; end
    return parse(Int, s[i:j-1]), j
end

function build_pop(rec::NamedTuple)
    n, m, M = rec.n, rec.m, rec.M     
    @ncpolyvar A[1:n] B[1:m]           

    Av = [1.0; A]                     
    Bv = [1.0; B]                     

    pop = sum(M[i,j] * Av[i] * Bv[j] for i in eachindex(Av), j in eachindex(Bv))

    return A, B, pop    
end

function generate_bound_data(A, r, tol; index=1:length(A), filename=joinpath(@__DIR__, "bell_bound_data_r3.txt"))
    names = String[]
    numerical_bounds = BigFloat[]
    correct_bounds = BigFloat[]
    bound_differences = BigFloat[]
    eigs_proj = BigFloat[]
    eigs_raw = BigFloat[]
    diffs_raw = BigFloat[]
    allocs = Float64[]
    glengths = Int[]
    ns = Int[]
    ms = Int[]

    # append mode if file exists, write mode otherwise
    mode = isfile(filename) ? "a" : "w"
    io = open(filename, mode)

    for i in index
        clear_caches!()
        n = A[i].n
        m = A[i].m
        name = A[i].name

        println("\n" * name)

        A_pop = build_pop(A[i])
        X, Y, pop = A_pop
        vars = [X; Y]

        alloc = @allocated begin
            data = rational_certificate(pop, [], [], vars, r;
                partition=length(X), constraint="projector",
                QUIET=true, QUIETTS=true, tol=tol)
        end

        new_bound, old_bound, bdiff, eigproj, eigraw, diffraw, glength = data
        new_bound = -new_bound
        old_bound = -old_bound

        # push to arrays
        push!(names, name)
        push!(numerical_bounds, old_bound)
        push!(correct_bounds, new_bound)
        push!(bound_differences, bdiff)
        push!(eigs_proj, eigproj)
        push!(eigs_raw, eigraw)
        push!(diffs_raw, diffraw)
        push!(glengths, glength)
        push!(ns, n)
        push!(ms, m)
        push!(allocs, alloc)

        println(io,
            name, " ",
            old_bound, " ",
            new_bound, " ",
            bdiff, " ",
            eigproj, " ",
            eigraw, " ",
            diffraw, " ",
            glength, " ",
            n, " ",
            m, " ",
            alloc
        )
        flush(io) 
    end

    close(io)

    return names, numerical_bounds, correct_bounds, bound_differences,
           eigs_proj, eigs_raw, diffs_raw, glengths, ns, ms, allocs
end

function generate_bound_data_sparse(A, r, tol; index=1:length(A), QUIET=true)
    names = Vector{String}([])
    numerical_bounds = Vector{BigFloat}([])
    correct_bounds = Vector{BigFloat}([])
    bound_differences = Vector{BigFloat}([])
    diffs_raw = Vector{BigFloat}([])
    glengths = Vector{Int}([])
    ns = Vector{Int}([])
    ms = Vector{Int}([])
    nCliquesList = Vector{Int}([])

    for i in index
        let      
            n = A[i].n
            m = A[i].m
            name = A[i].name
            println("\n" * name)
            A_pop = build_pop(A[i])
            X, Y, pop = A_pop
            vars = [X;Y]
            data = rational_certificate_sparse(pop, [], [], vars, r; partition=length(X), constraint="projector", QUIET=QUIET, tol=tol)
            new_bound, old_bound, total_shift, nCliques, RHS_per_clique, Gproj_per_clique =
                data.newbound, data.oldbound, data.totalshift,
                data.nCliques, data.RHS_per_clique, data.Gproj_per_clique

            new_bound = -new_bound
            old_bound = -old_bound
            total_shift = -total_shift

            push!(names, name)
            push!(numerical_bounds, old_bound)
            push!(correct_bounds, new_bound)
            push!(bound_differences, total_shift)
            push!(diffs_raw, sum(abs.(values(merged_coeffs(RHS_per_clique[1] - pop)))))
            push!(glengths, length(Gproj_per_clique[1][1,:]))
            push!(ns, n)
            push!(ms, m)
            push!(nCliquesList, nCliques)
        end
    end
    return names, correct_bounds, numerical_bounds,
           bound_differences, diffs_raw,
           glengths, ns, ms, nCliquesList
end

function read_bound_data(filename)
    names             = String[]
    numerical_bounds  = BigFloat[]
    correct_bounds    = BigFloat[]
    bound_differences = BigFloat[]
    eigs_proj         = BigFloat[]
    eigs_raw          = BigFloat[]
    diffs_raw         = BigFloat[]
    glengths          = Int[]
    ns                = Int[]
    ms                = Int[]
    allocs            = Float64[]

    open(filename, "r") do io
        for line in eachline(io)
            startswith(line, "#") && continue   # skip header

            parts = split(line)
            length(parts) == 11 || error("Malformed line: $line")

            name            = parts[1]
            old_bound       = parse(BigFloat, parts[2])
            new_bound       = parse(BigFloat, parts[3])
            bdiff           = parse(BigFloat, parts[4])
            eigproj         = parse(BigFloat, parts[5])
            eigraw          = parse(BigFloat, parts[6])
            diffraw         = parse(BigFloat, parts[7])
            glen            = parse(Int,       parts[8])
            n               = parse(Int,       parts[9])
            m               = parse(Int,       parts[10])
            alloc_bytes     = parse(Float64,   parts[11])

            push!(names,             name)
            push!(numerical_bounds,  old_bound)
            push!(correct_bounds,    new_bound)
            push!(bound_differences, bdiff)
            push!(eigs_proj,         eigproj)
            push!(eigs_raw,          eigraw)
            push!(diffs_raw,         diffraw)
            push!(glengths,          glen)
            push!(ns,                n)
            push!(ms,                m)
            push!(allocs,            alloc_bytes)
        end
    end

    return (
        names,
        numerical_bounds,
        correct_bounds,
        bound_differences,
        eigs_proj,
        eigs_raw,
        diffs_raw,
        glengths,
        ns,
        ms,
        allocs
    )
end

#############################################################

bell_dir = joinpath(@__DIR__, "Bell_A2_A89")

A = load_bell_matrices(bell_dir)
A = sort(A; by = a -> parse(Int, match(r"\d+", a.name).match))

# Simple test
println("********** A[2] **********")

A_pop = build_pop.(A)
A_pop_test = A_pop[1]
X, Y, pop = A_pop_test
vars = [X;Y]
A_pop
data_test = @allocated rational_certificate(pop, [], [], vars, 2; partition=length(X), constraint="projector", QUIET=true, tol=10e-20, eigprec=128)