#include "lc.lp".

#const n = 8.
time(0..n).
step(I,I+1) :- time(I), I < n.

&assign { s(I') := s(I)+D } :- step(I,I'), acc(D,I').
&assign { s(I') := s(I)-D } :- step(I,I'), slow(D,I').

% Inertia for the speed
&assign { s(I') := s(I) } :- step(I,I'), not &sum{ s(I') } != s(I).

def_s(I) :- time(I), &sum{s(I); -s(I)}=0.

% The position keeps moving with the speed, if there is a speed
&assign { p(I') := p(I)+s(I) } :- step(I,I'), def_s(I).

&assign{ rdpos := 400000 }.
%&assign { rdlimit := 90000 }.
%&assign { rdpos := 300000..400000 }.

fine(I') :- step(I,I'), &sum{ p(I) } < rdpos, &sum{ p(I') } >= rdpos, &sum{ s(I') } > rdlimit.

% initial position 0 and initial speed 50 Km/h
&assign { p(0) := 0 }.
&assign { s(0) := 80000 }.
% We accelerate 11.35 Km/h 
acc(11350,4).
% We slow down 2.301 Km/h
slow(2301,6).

#show fine/1.
