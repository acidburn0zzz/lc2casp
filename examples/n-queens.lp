#include "lc.lp".

#const n = 8.
n(1..n).
:- not &distinct { q(X)+0 : n(X) }.
:- not &distinct { q(X)+X : n(X) }.
:- not &distinct { q(X)-X : n(X) }.
&assign { q(X) := 1..n } :- n(X).
