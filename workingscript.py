# Currently in progress calculations

from wickcontract import *
from PSfields import *
V1 = V(1)

U2 = U(2)

L21=wick(multiply(V(1),U(2)))
for t in L21.terms:
    for f in t.fields:
        f.arg = 'z_1' if f.arg else '';

UU3 = wick(multiply(U(2),U(3)));
for t in UU3.terms:
    for f in t.fields:
        f.arg = 'z_3' if f.arg else '';


L2331 = wick(V1,divide(UU(2,3,3),'(z_2-z_3)'));
#L2331 = wick(V1,UU3)


ex21 = wick(Q(1),multiply(A(1,'\\alpha',2),W(2,'\\alpha',2)));
for t in ex21.terms:
    for f in t.fields:
        f.arg = 'z_1' if f.arg else '';

ex2131 = wick(ex21,U(3))
print(ex2131.str(False))
