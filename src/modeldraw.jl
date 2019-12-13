function tikzlattice(model::UnitCellQModel,
                     filename::String)
    isfile(filename) && error("file already exists: $filename")
    D = dimension(model)
    D == 2 || error("only 2D lattices can be drawn with tikz!")
    lattice = model.lattice

    bondcolors = [
        "blue", "orange", "cyan", "olive", "red", "brown", "green", "yellow" ]
    ballcolors = [
        "green", "blue", "red", "yellow", "brown", "cyan", "orange", "olive"]

    a1 = lattice.unitc.as[1]
    a2 = lattice.unitc.as[2]

    open(filename, "w") do f
        write(f, """\\documentclass[tikz, border=1mm]{standalone}
  \\usepackage{tikz,pgf}
  \\usetikzlibrary{calc,arrows,arrows.meta}
  \\usetikzlibrary{decorations.markings}

  \\begin{document}

  \\begin{tikzpicture}[
        x={$a1},
        y={$a2},
        ]
  """)

        write(f, "  \\tikzset{\n")

        bonds  = [inter for inter in model.inters if support(inter)==2]
        n3bodys = [inter for inter in model.inters if support(inter)==3]
        n4odys = [inter for inter in model.inters if support(inter)==4]

        for i in eachindex(bonds)
            write(f, "    bond$i/.style={thin, double, $(bondcolors[i])!70!black},\n")
        end
        for i in eachindex(lattice.unitc.n)
            write(f, "    ball$i/.style={inner color=blue, ball color=$(ballcolors[i])!20!black},\n")
        end

        write(f, """
    ball/.style={inner color=blue, ball color=green!20!black}
  }

  \\def\\l{1cm}
  \\def\\radius{\\l/10}\n
""")

        pad = 0.3
        ly, lx = model.lattice.sizes
        write(f, "  \\clip (-$pad, -$pad) -- ++(0.0, $(ly-1+2*pad)) -- ++($(lx-1+2*pad), 0.0) -- ++(0.0, -$(ly-1+2*pad)) -- cycle;\n")

        for is in Iterators.product([0:l+1 for l in model.lattice.sizes]...)
            for index in eachindex(bonds)
                interaction = bonds[index]
                #support(interaction) == 2 || continue
                n1, n2 = interaction.ucidxs
                off1, off2 = interaction.offsets
                #if lattice.bc == :OBC
                x_uc1 = [off1[i] .+ is[i] for i=1:D]
                x_uc2 = [off2[i] .+ is[i] for i=1:D]
                if all([0 <= x_uc1[i] <= lattice.sizes[i]+1 for i=1:D]) &
                    all([0 <= x_uc2[i] <= lattice.sizes[i]+1 for i=1:D])
                    one = reduce(coordsum,
                                 #[(x_uc1 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                 [reverse(x_uc1 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                 init=(0,0)).+ lattice.unitc.sites[n1]
                    two = reduce(coordsum,
                                 #[(x_uc2 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                 [reverse(x_uc2 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                 init=(0,0)).+ lattice.unitc.sites[n2]
                    write(f, "  \\draw [bond$index] ($(one[1])*\\l, $(one[2])*\\l) -- ($(two[1])*\\l, $(two[2])*\\l);\n")
                end
                #else
                #    error("boundary not supported yet!")
                #end
            end
        end


        for ls in Iterators.product([1:l for l in model.lattice.sizes]...)
            x_uc = reduce(coordsum,
                          #[(ls .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                          [reverse(ls .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                          init=(0,0))
            for n in 1:lattice.unitc.n
                x_site = x_uc .+ lattice.unitc.sites[n]
                write(f, "  \\shade [ball$n] ($(x_site[1])*\\l, $(x_site[2])*\\l) circle (\\radius);\n")
            end
        end

        for index in eachindex(bonds)
            interaction = bonds[index]
            ns = interaction.ucidxs
            offs = interaction.offsets
            for is in Iterators.product([1:l for l in model.lattice.sizes]...)
                x_ucs = [[offs[n][i] .+ is[i] for i=1:D] for n in 1:length(ns)]
                if all([
                    all([1 <= x_ucs[n][i] <= lattice.sizes[i]
                         for i=1:D])
                    for n in 1:length(ns)])
                    sites = [reduce(coordsum,
                                    #[(x_uc1 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                    [reverse(x_ucs[n] .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                    init=(0,0)) .+ lattice.unitc.sites[ns[n]]
                             for n in 1:length(ns)]
                    #write(f, "  \\draw [->, shorten >=\\radius, shorten <=\\radius] ($(sites[1][1])*\\l, $(sites[1][2])*\\l) to [out=160, in=20, looseness=1] ($(sites[2][1])*\\l, $(sites[2][2])*\\l);\n")
                    mid = (sites[1] .+ sites[2]) ./ 2
                    write(f, "  \\node at ($(mid[1])*\\l, $(mid[2])*\\l) {\$t_$index\$};\n")
                    break
                end
            end
        end
        write(f, "\\end{tikzpicture}\n")
        write(f, "\\end{document}")
    end
end
