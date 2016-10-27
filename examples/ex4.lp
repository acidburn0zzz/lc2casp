#include "lc.lp".

&assign { z := x; t := y }.
&assign { x := 1 }.
&assign { y := 1 } :- &sum { z } = 1.
&assign { z := 1 } :- &sum { y } = 1.
