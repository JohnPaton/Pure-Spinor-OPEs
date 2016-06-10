
from copy import deepcopy as dcopy
from basics import *
from operations import *
from itertools import *
from PSfields import *
from PSOPEs import *

# return the powerset of a given set (from itertools documentation)
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# do wick contractions for the fields in a single term
def twick(termin):
    t = dcopy(termin);
    n = len(t.fields);

    opes = [];
    
    for i in range(0,n):
        for j in range(0,n):
            o = ope(t.fields[i],t.fields[j])
            if o.terms and ([j,i] not in opes):
                opes.append([i,j])
                
    opesets = [list(i) for i in powerset(opes)];
    opesets2 = [];

    for oset in opesets:
        bad = False
        osetflat = list(chain(*oset))
        for i in osetflat:
            rest = osetflat[0:osetflat.index(i)]+osetflat[osetflat.index(i)+1:]
            
            if i in rest:
                bad = True;
        if not bad:
            opesets2.append(oset);

    opesets = opesets2;          
    
    e = exp([])
    for oset in opesets:
        enew = exp(dcopy(t));
        if oset:
            for c in oset:
                for tnew in enew.terms:
                    ofields = [t.fields[c[0]],t.fields[c[1]]];
                    operes = ope(ofields[0],ofields[1]);

                    nt = len(operes.terms);
                    newts = [tnew];

                    for n in range(1,nt):
                        newts.append(dcopy(tnew));

                    i,j = -1,-1

                    for ti in range(0,nt):
                        tins = newts[ti]
                        for f in tins.fields:
                            if fieldeq(f,ofields[0]):
                                i = tins.fields.index(f);
                            if fieldeq(f,ofields[1]):
                                j = tins.fields.index(f);
                        f1 = tins.fields[i];
                        f2 = tins.fields[j];
                        tins.shiftto(j,i+1);
                        i = tins.fields.index(f1)
                        j = tins.fields.index(f2)
                        tins.shiftto(i,j-1);
                        i = tins.fields.index(f1)
                        j = tins.fields.index(f2)
                        tins.opeinsert(i,operes.terms[ti]);

                if nt > 1:
                    enew.terms += newts[1:];

                
            e.terms += enew.terms;

    return e

# do wick contractions for each term in e1, or in e1*e2
def wick(e1,e2=exp([])):
    if e2.terms:
        e = multiply(e1,e2);
    else:
        e = e1;
        
    eout = exp([]);
    for t in e.terms:
        enew = twick(t);
        eout.terms += enew.terms

    return eout
