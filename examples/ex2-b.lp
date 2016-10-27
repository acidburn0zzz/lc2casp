#include "lc.lp".
&assign { x := 1 }.
&assign { y := x - 1 } :- not &sum { 1 } <= z.
