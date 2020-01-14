push!(LOAD_PATH, "./src/")

using DelimitedFiles
using JLD2

using QuantumModels
using MatrixProductStateTools
using ExactDiagonalizationTools
using SymTensors
using GaussianFermions
using GutzwillerMPS
using TensorNetAlgs

ly = 6
lx = 30
L = lx * ly
t1, t2, t3, mu = 1.0, 1.0, 1.0, 0.0
j, k = 2.0, 0.6
boundary = (:PBC, :OBC)
b_str = "PBCY"

hamil = generatebdg(
    triangularhopping((ly, lx), t1, t2, t3, mu, boundary=boundary)
)
cm = correlationmatrix(hamil, div(L, 2))
fgs = corrmat2gmps(cm)

sz, sp, sm = spinoperators(1/2, symmetry=U1)


infos = "$b_str,ly=$ly,lx=$lx,t1=$t1,t2=$t2,t3=$t3,mu=$mu"
for m in [2000]

    print("Fermionic mps m=$m ")
    filename = "mps,triangular,hop,$infos,m=$m.jld2"
    if isfile(filename)
        println("already exists! loading...")
        @load filename mps
    else
        println("making ...")
        mps = gmps2mps(fgs, m, symmetry=U1)
        #@save filename mps

        # zzdata = measure(mps, sz, sz)
        # pmdata = measure(mps, sp, sm)
        vnedata = entanglemententropy(mps)

        println("Measurements for hop m = $m ...")
        # writedlm( "zz,triangular,hop,$infos,m=$m.dat", zzdata)
        # writedlm( "pm,triangular,hop,$infos,m=$m.dat", pmdata)
        writedlm( "vne,triangular,hop,$infos,m=$m.dat", vnedata)
    end

    for M in [12000]
        println("ZipandGutzwiller for M = $M ...")

        #filename = "mps,triangular,gutzF23,$infos,m=$m,M=$M.jld2"
        if isfile(filename)
            println("Already exists! dismiss.")
        else
            println("Making ...")
            mpsgutz = zipandgutzwiller!(mps, mps, mode=:F23, maxdim=M)
            normalize!(mpsgutz)
            #@save  filename mpsgutz
            println("Measurements for M = $M ...")
            # zzdata = measure(mpsgutz, sz, sz)
            # pmdata = measure(mpsgutz, sp, sm)
            vnedata = entanglemententropy(mpsgutz)
            # writedlm("zz,triangular,gutzF23,$infos,m=$m,M=$M.dat", zzdata)
            # writedlm("pm,triangular,gutzF23,$infos,m=$m,M=$M.dat", pmdata)
            writedlm("vne,triangular,gutzF23,$infos,m=$m,M=$M.dat", vnedata)

            mpo = generatempo(triangularspinmodel((ly, lx), j, j, j, k, k, k,
                                                  boundary=boundary, symmetry=U1))
            reducempo!(mpo)
            center_at!(mpsgutz, 1)
            env = initialenv(mps, mpo)
            dmrg2sitesweep!(mps, mpo, env, maxdim=M, verbose=true)

            dmrginfo = "j=$j,k=$k,M=$M"
            # zzdata = measure(mpsgutz, sz, sz)
            # pmdata = measure(mpsgutz, sp, sm)
            vnedata = entanglemententropy(mpsgutz)
            # writedlm("zz,triangular,1dmrg,$dmrginfo,gutzF23,$infos,m=$m,M=$M.dat", zzdata)
            # writedlm("pm,triangular,1dmrg,$dmrginfo,gutzF23,$infos,m=$m,M=$M.dat", pmdata)
            writedlm("vne,triangular,1dmrg,$dmrginfo,gutzF23,$infos,m=$m,M=$M.dat", vnedata)

        end
    end
end
