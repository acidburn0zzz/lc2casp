#include "lc.lp".

p(a;b;c).
&show { X/0 : p(X) }.
&assign { X:=1..3 } :- p(X).

:- not &distinct { X : p(X) }.
