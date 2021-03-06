export noether_normalization, normalization
export isreduced

##############################################################################
#
# Properties of affine algebras
#
##############################################################################

@doc Markdown.doc"""
    isreduced(A::MPolyQuo)

Return `true` if `A` is reduced, `false` otherwise.
"""
function isreduced(A::MPolyQuo{T}) where T
  I = A.I
  return I == radical(I)
end

##############################################################################
#
# Algebra Containment
#
##############################################################################


##############################################################################
#
# Properties of maps of affine algebras
#
##############################################################################



##############################################################################
#
# Normalization
#
##############################################################################

@doc Markdown.doc"""
    normalize(A::MPolyQuo)

Finds the normalization of a reduced affine algebra over a perfect field $K$:
Given the quotient $A=R/I$ of a multivariate polynomial ring $R$ over $K$
modulo a radical ideal $I$, compute the integral closure $\overline{A}$ 
of $A$ in its total ring of fractions $Q(A)$, together with the embedding 
$f: A \rightarrow \overline{A}$. The function relies on the algorithm 
of Greuel, Laplagne, and Seelisch which proceeds by finding a suitable decomposition 
$I=I_1\cap\dots\cap I_r$ into radical ideals $I_k$, together with
the normalization maps $f_k: R/I_k \rightarrow A_k=\overline{R/I_k}$, such that 

$f=f_1\times \dots\times f_r: A \rightarrow A_1\times \dots\times A_r=\overline{A}$

is the normalization map of $A$. For each $k$, the function specifies two representations
of $A_k$: It returns an array of triples $(A_k, f_k, \mathfrak a_k)$,
where $A_k$ is represented as an affine $K$-algebra, and $f_k$ as a map of affine $K$-algebras.
The third entry $\mathfrak a_k$ is a tuple $(d_k, J_k)$, consisting of an element
$d_k\in A$ and an ideal $J_k\subset A$, such that $\frac{1}{d_k}J_k = A_k$ 
as $A$-submodules of the total ring of fractions of $A$.

By default, as a first step on its way to find the decomposition $I=I_1\cap\dots\cap I_r$, 
the algorithm computes an equidimensional decomposition of the radical ideal $I$.
Alternatively, if specified by `alg=:primeDec`, the algorithm computes $I=I_1\cap\dots\cap I_r$
as the prime decomposition of the radical ideal $I$. 

If `alg=:withDelta` is specified, the algorithm computes additionally the delta invariant of 
$A$, that is, the dimension $\dim_K(\overline{A}/A)$. More precisely, it returns a tuple
consisting of an array containing the delta invariants of the $A_k$ and an integer, the
(total) delta invariant of $A$.  The return value -1 indicates that the delta invariant is infinite.

CAVEAT: The function does not check whether $A$ is reduced. Use `isreduced(A)` in case 
you are unsure (this may take some time).
"""
function normalize(A::MPolyQuo{T}) where T
  I = A.I
  singular_assure(I)
  l = Singular.LibNormal.normal(I.gens.S)
  return [
    begin
      newR = l[1][i][1]
      newA, newAmap = quo(newR, MPolyIdeal(newR, l[1][i][2][:norid]))
      hom = AlgebraHomomorphism(A, newA, map(newAmap, gens(l[1][i][2][:normap])))
      idgens = map(p->_badpolymap(p, A.R), gens(l[2][i]))
      (newA, hom, (A(idgens[end]), ideal(A, idgens)))
    end
    for i in 1:length(l[1])]
end


##############################################################################
#
# Noether Normalization
#
##############################################################################

@doc Markdown.doc"""
    noether_normalization(A::MPolyQuo)

Given an affine algebra $A=R/I$ over a field $K$, the function returns  a triple $(L,f,g)$, such that:
$L$ is an array of $d=\dim A$ elements of $A$, all represented by linear forms in $R$, and
such that $K[L]\rightarrow A$ is a Noether normalization for $A$; $f: A \rightarrow A$ is an
automorphism of $A$, induced by a linear change of coordinates of $R$, and mapping the
$f_i$ to the last $d$ variables of $A$; and $g = f^{-1}$.
"""
function noether_normalization(A::MPolyQuo)
 I = A.I
 R = base_ring(I)
 singular_assure(I)
 l = Singular.LibAlgebra.noetherNormal(I.gens.S)
 i1 = [R(x) for x = gens(l[1])]
 i2 = [R(x) for x = gens(l[2])]
 m = matrix([[coeff(x, y) for y = gens(R)] for x = i1])
 mi = inv(m)
 h1 = AlgebraHomomorphism(A, A, map(A, i1))
 h2 = AlgebraHomomorphism(A, A, map(A, collect(matrix([gens(R)])'*map_entries(R, mi))))
 return map(x->h2(A(x)), i2), h1, h2
end



