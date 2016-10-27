#include "lc.lp".

p(1..3).
&show { x/1 : p(X) }.
&assign { x(X):=2..3 } :- p(X).

:- &sum { x(X) : p(X) } > 7.
