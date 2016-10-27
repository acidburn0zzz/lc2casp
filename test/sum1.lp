#include "lc.lp".

p(a;b;c).
&show { X/0 : p(X) }.
&assign { X:=2..3 } :- p(X).

:- &sum { X : p(X) } > 7.
