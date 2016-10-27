#include "lc.lp".

&assign { x := 1..2 }.

x :- &sum { x } > 1.
y :- &sum { y } > 1.
xy :- &sum { x; y } > 1.
