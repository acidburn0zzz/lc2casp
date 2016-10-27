#include "lc.lp".

2 { a; b; c }.

&assign { x := 1..5 } :- a.
&assign { y := 2..5 } :- b.
&assign { z := 3..5 } :- c.

&minimize { x; y; z }.
