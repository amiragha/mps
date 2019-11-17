push!(LOAD_PATH, "./src/")

using DelimitedFiles
using JLD2

using QuantumModels
using MatrixProductStateTools
using ExactDiagonalizationTools
using SymTensors
using GaussianFermions
using GutzwillerMPS

ly = 6
lx = 24
L = lx * ly
t1, t2, t3, mu = 1.0, 1.0, 1.0, 0.0
boundaries = [:OBC]

for boundary in boundaries
    println("Statring boundary condition $boundary ...")
    hamil = generatebdg(
        triangularhopping((ly, lx), t1, t2, t3, mu, boundary=boundary)
    )
    cm = correlationmatrix(hamil, div(L, 2))
    fgs = generate_fishmangates(cm)

    sz, sp, sm = spinoperators(1/2, symmetry=:U1)


    infos = "$boundary,ly=$ly,lx=$lx,t1=$t1,t2=$t2,t3=$t3,mu=$mu"
    for m in [100, 200]

        print("Fermionic mps m=$m ")
        filename = "mps,triangular,hop,$infos,m=$m.jld2"
        if isfile(filename)
            println("already exists! loading...")
            @load filename mps
        else
            println("making ...")
            mps = fishman2mps(fgs, m, symmetry=:U1)
            @save filename mps

            zzdata = real.(measure_2point(mps, sz, sz))
            pmdata = real.(measure_2point(mps, sp, sm))
            vnedata = entanglemententropy(mps)

            println("Measurements for hop m = $m ...")
            writedlm( "zz,triangular,hop,$infos,m=$m.dat", zzdata)
            writedlm( "pm,triangular,hop,$infos,m=$m.dat", pmdata)
            writedlm( "vne,triangular,hop,$infos,m=$m.dat", vnedata)
        end

        for M in [600, 900]
            println("ZipandGutzwiller for M = $M ...")

            filename = "mps,triangular,gutzF23,$infos,m=$m,M=$M.jld2"
            if isfile(filename)
                println("Already exists! dismiss.")
            else
                println("Making ...")
                mpsgutz = zipandgutzwiller!(mps, mps, mode=:F23, maxdim=M)
                normalize!(mpsgutz)
                @save  filename mpsgutz
                println("Measurements for M = $M ...")
                zzdata = real.(measure_2point(mpsgutz, sz, sz))
                pmdata = real.(measure_2point(mpsgutz, sp, sm))
                vnedata = entanglemententropy(mpsgutz)
                writedlm("zz,triangular,gutzF23,$infos,m=$m,M=$M.dat", zzdata)
                writedlm("pm,triangular,gutzF23,$infos,m=$m,M=$M.dat", pmdata)
                writedlm("vne,triangular,gutzF23,$infos,m=$m,M=$M.dat", vnedata)
            end
        end
    end
end
