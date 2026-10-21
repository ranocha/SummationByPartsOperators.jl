"""
    VretinarisSchnetter2026()

Source of the fourth- and sixth-order spherical summation-by-parts operators
constructed from Cartesian SBP operators on a mirrored uniform grid.

- Vretinaris, Schnetter (2026)
  High-order summation-by-parts schemes for first-order hyperbolic systems in
  curvilinear coordinates with singularities.
  Physical Review D.
  [doi:10.1103/hg1v-qxsj](https://link.aps.org/doi/10.1103/hg1v-qxsj)
"""
struct VretinarisSchnetter2026 <: SourceOfCoefficients end

function Base.show(io::IO, source::VretinarisSchnetter2026)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
              "Vretinaris, Schnetter (2026)\n",
              "  High-order summation-by-parts schemes for first-order hyperbolic systems\n",
              "  in curvilinear coordinates with singularities.\n",
              "  Physical Review D. doi:10.1103/hg1v-qxsj.")
    end
end
