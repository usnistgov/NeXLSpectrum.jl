using NeXLMatrixCorrection, NeXLSpectrum
using PrettyTables

det = BasicEDS(4096, 0.0, 5.0, 132.0, 1, Dict(KShell => n"He", LShell => n"K", MShell => n"Cs", NShell => n"Pu"))

open(joinpath(homedir(), "Desktop", "suitable_references.tex"),"w") do io
    for el in elements[5:92]
        println(io,"\\subsection{Element: $(name(el))}")
        pretty_table(io,
            suitability(el, det, minC=0.01, latex=true),
            nosubheader=true, backend=:latex,
            label="The suitability of references for $(name(el))",
            wrap_table=false, table_type=:longtable )
        println(io)
    end
end
