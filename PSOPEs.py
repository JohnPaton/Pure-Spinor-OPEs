# the Pure Spinor OPEs

from copy import deepcopy as dcopy
from basics import *
from operations import *
from PSfields import * # Use these now that they have been defined

# return the OPE of f1 and f2 if they are at differnt positions
# Current OPEs are those of the pure spinor formalism but this could be extended
# Worldsheet are handled "by hand" just by increasing the order of the pole
# and potentially changing the sign
# TO DO: Update to use field function definitions instead of defining new
#        fields each time
# See thesis for list of OPEs
def ope(f1,f2):
    t = term([f1,f2]);
    e = exp([]);
    if f1.arg == f2.arg:
        pass
    elif (not f1.arg) or (not f2.arg):
        pass
    
    elif f1.name == 'd': 
        if type(f2) is symfield:
            Df2 = dcopy(f2);
            Df2.derivind = f1.inddn + Df2.derivind;
            Df2.ferm = (Df2.ferm + 1) % 2;
            den = '('+f1.arg+'-'+f2.arg+')'
            e = exp(term(Df2,denom=den))

        elif f2.name == 'd':
            b = t.choosedummy(alphabet('m'))
            gamma = field('\\gamma', inddn = [b]+f1.inddn+f2.inddn)
            Pi = field('\\Pi',indup = b, arg = f2.arg)
            den = '('+f1.arg+'-'+f2.arg+')'
            e = exp(term([gamma,Pi],denom = den))
        
        elif f2.name == '\\Pi':
            b = t.choosedummy(alphabet(f1.inddn[0]));
            gamma = field('\\gamma',f2.indup,f1.inddn+[b])
            dtheta = field('\\partial\\theta',b,ferm=1,arg=f2.arg)
            den = '('+f1.arg+'-'+f2.arg+')'
            e = exp(term([gamma, dtheta],denom=den))

        elif f2.name == '\\theta':
            delta = field('\\delta',f1.inddn,f2.indup)
            den = '('+f1.arg+'-'+f2.arg+')'
            e = exp(term(delta,denom=den))

        elif f2.name == '\\partial\\theta':
            delta = field('\\delta',f1.inddn,f2.indup)
            den = '('+f1.arg+'-'+f2.arg+')^2'
            e = exp(term(delta,denom=den))
            
    elif f1.name == '\\Pi':
        if type(f2) is symfield:
            k = symfield('k',indup = f1.indup, label=f2.label)
            den = '('+f1.arg+'-'+f2.arg+')'
            e = exp(term([k,dcopy(f2)],sign=-1,denom=den));
            
        elif f2.name=='\\Pi':
            eta = field('\\eta',f1.indup+f2.indup);
            den = '('+f1.arg+'-'+f2.arg+')^2'
            e = exp(term(eta,denom=den,sign=-1))
            
    elif f1.name == 'N':
        if f2.name == '\\lambda':
            a = t.choosedummy(alphabet('\\alpha'));
            pre = '\\frac{1}{2}';
            gamma = field('\\gamma',f1.indup+f2.indup,inddn=a)
            lam = field('\\lambda',a,arg = f2.arg)
            den = '('+f1.arg+'-'+f2.arg+')'
            e = exp(term([gamma,lam],sign=-1,pref=pre,denom=den))

        elif f2.name == 'N':
            [m,n] = f1.indup;
            [p,q] = f2.indup;
            den1 = '('+f1.arg+'-'+f2.arg+')'
            den2 = '('+f1.arg+'-'+f2.arg+')^2'
            terms = [];

            s = 1;
            for i in [[m,p,n,q],[m,q,n,p],[n,q,m,p],[n,p,m,q]]:
                N = field('N',indup=i[0:2],arg=f2.arg)
                eta = field('\\eta',indup=i[2:]);
                tnew = term([N,eta],sign=s,denom=den1)
                terms.append(tnew)
                s*=-1
                
            s=-1
            for i in [[n,p,m,q],[n,q,m,p]]:
                pre = '3';
                eta1 = field('\\eta',indup=i[0:2])
                eta2 = field('\\eta',indup=i[2:])
                tnew = term([eta1,eta2],sign=s,denom=den2,pref=pre)
                terms.append(tnew)
                s*= -1

            e = exp(terms)
        

    elif f1.name == 'J' and f2.name == '\\lambda':
        lam = dcopy(f2);
        lam.arg = f1.arg;
        den = '('+f1.arg+'-'+f2.arg+')'
        e = exp(term(lam,denom=den));

    elif f1.name == '\\omega' and f2.name == '\\lambda':
        den = '('+f1.arg+'-'+f2.arg+')'
        e = exp(term(field('\\delta',f2.indup,f1.inddn),denom=den));

    elif f2.name == '\partial\\lambda':
        lam = dcopy(f2);
        lam.name = '\\lambda';
        den = '('+f1.arg+'-'+f2.arg+')'
        e = ope(f1,lam);
        for t in e.terms:
            t.denom += den;
            for f in t.fields:
                if f.name == '\\lambda':
                    f.name = '\partial\\lambda';
        
            
    if e.terms:
        for t in e.terms:
            for f in t.fields:
                f.opeflag = True;

    return e
