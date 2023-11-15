using .AirfoilPolars: Polar

"""
    get_filename_ext(fname::String)::Tuple{String, String}

Get the filename and extension from a path.

# Returns 
- `fname_noext`: Filename w/o extension
- `ext`: Extension
"""
function get_filename_ext(fname::String)::Tuple{String, String}
    fname = split(fname, "/")[end]
    ext = split(fname, ".")[end]
    fname_noext = fname[1:end-length(ext)-1]
    return fname_noext, ext
end

function _load_csv(fname)
    df = CSV.read(fname, DataFrame)
    n_data = length(df[:,1])

    re = NaN64
    for n in ["re", "Re", "RE"]
        try 
            re = df[!, n][1]
            break
        catch 
            continue
        end
    end

    mach = NaN64
    for n in ["mach", "Mach", "MACH"]
        try 
            mach = df[!, n][1]
            break
        catch 
            continue
        end
    end

    n_crit = NaN64
    for n in ["n_crit", "N_crit", "N_Crit", "N_Crit"]
        try 
            n_crit = df[!, n][1]
            break
        catch 
            continue
        end
    end

    xtrip = (NaN64, NaN64)
    for n in ["xtrip", "Xtrip", "XTRIP", "xtrip_", "Xtrip_", "XTRIP_"]
        try 
            xtrip1 = df[!, n*"1"][1]
            xtrip2 = df[!, n*"1"][1]
            xtrip = (xtrip1, xtrip2)
            break
        catch 
            continue
        end
    end

    name_airfoil = "NONAME"
    for n in ["name", "Name", "NAME"]
        try 
            name_airfoil = df[!, n][1]
            break
        catch 
            continue
        end
    end

    alpha = zeros(n_data)
    cl = zeros(n_data)
    cd = zeros(n_data)
    cm = zeros(n_data)

    for n in ["alpha", "Alpha", "ALPHA"]
        try 
            alpha = df[!, n]
            break
        catch 
            continue
        end
    end

    for n in ["cl", "CL", "Cl"]
        try 
            cl = df[!, n]
            break
        catch 
            continue
        end
    end

    for n in ["cd", "CD", "Cd"]
        try 
            cd = df[!, n]
            break
        catch 
            continue
        end
    end

    for n in ["cm", "CM", "Cm"]
        try 
            cm = df[!, n]
            break
        catch 
            continue
        end
    end

    return Polar(re, alpha, cl, cd, cm, mach, n_crit, xtrip, name_airfoil)
end

function _load_plr(fname)

    alpha = Vector{Float64}()
    cl = Vector{Float64}()
    cd = Vector{Float64}()
    cm = Vector{Float64}()
    re = NaN64
    mach = NaN64
    n_crit = NaN64
    xtrip = (NaN64, NaN64)
    name_airfoil = "NONAME"

    open(fname) do f
        # -----------------------------
        # Header
        # -----------------------------
        readline(f)
        info = readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        name_airfoil = split(readline(f))[1]
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        re = parse(Float64, split(readline(f))[2])
        readline(f)
        readline(f)
        readline(f)

        # -----------------------------
        # Data
        # -----------------------------
        for line in eachline(f)
            parts = split(line)
            push!(alpha, parse(Float64, parts[1]))
            push!(cl, parse(Float64, parts[2]))
            push!(cd, parse(Float64, parts[3]))
            push!(cm, parse(Float64, parts[4]))
        end
    end
    return Polar(re, alpha, cl, cd, cm, mach, n_crit, xtrip, name_airfoil)
end 

function load(fname)
    ext = get_filename_ext(fname)[2]

    if ext == "csv"
        return _load_csv(fname)
    elseif ext == "plr"
        return _load_plr(fname)
    else 
        error("Unknown file extension: $ext")
    end
end

function save(p::Polar, fname::String)
    n_alpha = length(p.alpha)

    ext = get_filename_ext(fname)[2]

    if ext == "csv"
        df = DataFrame(
            alpha = p.alpha,
            cl = p.cl,
            cd = p.cd,
            cm = p.cm,
            re = [p.Re for _ in n_alpha],
            mach = [p.M for _ in n_alpha],
            n_crit = [p.n_crit for _ in n_alpha],
            xtrip1 = [p.xtrip[1] for _ in n_alpha],
            xtrip2 = [p.xtrip[2] for _ in n_alpha],
        )
        
        CSV.write(fname, df)
    else
        error("Unknown file extension: $ext")
    end
end