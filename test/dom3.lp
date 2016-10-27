#include "lc.lp".

&show { x/0 }.
{c}.
&assign { x:=1; x:=3 } :- c.
