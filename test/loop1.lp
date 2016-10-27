#include "lc.lp".

1 { a; b }.

&assign { x := 1..2; y := 1..2 } :- a.
&assign { y := 1..2 } :- b, &sum { x } <= 2.
&assign { x := 1..2 } :- &sum { y } >= 1.

 :- &sum { x; y } <= 3.
