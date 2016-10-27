#include "lc.lp".

#const n = 8.
n(1..n).
:- not &distinct { q(X) : n(X) }.
:- &sum { q(X); -q(Y) } = X-Y, n(X), n(Y), X != Y.
:- &sum { q(X); -q(Y) } = Y-X, n(X), n(Y), X != Y.
&assign { q(1) := 1 }    :- not &sum { q(1) } != 1.
&assign { q(X) := 1..n } :- n(X), X > 1.

{ r }.
&assign { q(1) := 4 } :- r.
