#include "lc.lp".

&assign { x := y; x := 1 }.
&assign { y := z+3 } :- &sum {z} > 1.
&assign { z := 1..2 }.
