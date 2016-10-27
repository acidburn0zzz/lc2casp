#include "lc.lp".

&assign { x := 1..3 }.
&assign { y := 1..3 }.
&assign { z := 1..3 }.

:- not &sum { x; y } >= 4.
:- not &sum { y; z } >= 4.

&minimize { x; y; z }.
