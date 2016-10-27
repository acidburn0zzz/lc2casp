#include "lc.lp".

&assign { a := 1 } :- not &sum { b } = b.
&assign { b := 1 } :- not &sum { a } = a.

&assign { c := 1 } :- not &sum { d } = d.
&assign { d := 1 } :- not &sum { c } = c.

&assign { e := a; f := c }.
