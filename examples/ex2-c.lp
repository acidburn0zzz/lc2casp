#include "lc.lp".
&assign { z := 0..3 }.
&assign { x := 1 }.
&assign { y := x - 1 } :- not &sum { 1 } <= z.
