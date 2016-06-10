# multiplication, division, addition, and SYM derivatives
# of terms and expressions, check fields for equivalence

from basics import *
from copy import deepcopy as dcopy

# check for field equality, since fields are often deep copies. Fields equal
# if they have the same indices, argument, and name
def fieldeq(f1,f2):
    if type(f1) is not type(f2):
        return False
    names = f1.name == f2.name
    indups = f1.indup == f2.indup;
    inddns = f1.inddn == f2.inddn;
    args = f1.arg == f2.arg;
    flags = f1.opeflag == f2.opeflag;

    fields = names and indups and inddns and args and flags;

    if type(f1) is field:
        return fields
    elif type(f1) is symfield:
        labels = f1.label == f2.label;
        derivinds= f1.derivind == f2.derivind;

# add two expressions
def add(e1,e2=exp([])):
    e = exp([])
    if e2.terms:
        a = dcopy(e1)
        b = dcopy(e2)
        e.terms = a.terms + b.terms;
    elif type(e1) is list:
        for i in e1:
            a = dcopy(i);
            if type(a) is term:
                a = exp(a)
            e.terms += a.terms;

    return e



# multiply: two terms yeilds another term, expressions yeild expressions
def multiply(a,b):

    if type(a) is field or type(a) is symfield:
        e1 = exp(term(a)); 
    elif type(a) is term:
        e1 = exp(a);
    elif type(a) is exp:
        e1 = dcopy(a);
        
    if type(b) is field or type(b) is symfield:
        e2 = exp(term(b));
    elif type(b) is term:
        e2 = exp(b);
    elif type(b) is exp:
        e2 = dcopy(b);

    terms = [];
    for t1 in e1.terms:
        for t2 in e2.terms:
            
            d2 = t2.dummies();
            i1 = t1.inds();
            for d in d2:
                if d in i1:
                    t2.changedummy(d,excl = i1);
                
            d1 = t1.dummies();
            i2 = t2.inds();
            for d in d1:
                if d in i2:
                    t1.changedummy(d,excl = i2);
            
            fields = t1.fields+t2.fields;
            sign = t1.sign*t2.sign;
            denom = t1.denom + t2.denom;
            pref = t1.pref + t2.pref;
            terms += [term(fields,sign,pref,denom)];

    return exp(terms);

# concatenate the string d onto the denominator of the term n or all the
# denominators of the expression n
def divide(n,d): 
    if type(n) is term:
        e = exp(n);
    else:
        e = dcopy(n);

    for t in e.terms:
        t.denom += d;

    return e;

# Take SYM derivative with index ind of expression ein
def symD(ein,ind):
    e = exp(ein) if type(ein) is term else dcopy(ein);
    eout = exp([]);
    for t in e.terms:
        ds = t.dummies();
        if ind in ds:
            t.changedummy(ind,excl=[ind]);
        for f in t.fields:
            if type(f) is symfield:
                i = t.fields.index(f);
                dt = dcopy(t);
                dt.fields[i].derivind = [ind] + dt.fields[i].derivind;
                eout.terms.append(dt);
    return eout
