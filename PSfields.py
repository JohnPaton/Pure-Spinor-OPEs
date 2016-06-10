# Functions to create PS and SYM (composite) field instances

from basics import *
from operations import *

def Am(lab,ind,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg

    return symfield('A',inddn=ind,label=lab,arg=a)

def A(lab,ind,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg

    f = 1 if '\\alpha' in alphabet(ind) else 0;
    
    return symfield('A',inddn=ind,label=lab,arg=a,ferm=f)

def F(lab,ind,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg
    return symfield('\\mF',inddn=ind,label=lab,arg=a)

def W(lab,ind,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg
    return symfield('W',indup=ind,label=lab,arg=a,ferm=1)

def gamma(indup,inddn):
    return field('\\gamma',indup,inddn)

def k(lab,indup,inddn=''):
    return symfield('k',indup=indup,inddn=inddn,label=lab)

def d(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field('d',inddn = ind,arg=a,ferm=1)

def Pi(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field('\\Pi',ind,arg=a)

def N(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field('N',ind,arg=a);

def theta(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field('\\theta',ind,arg=a,ferm=1);

def dtheta(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field(r'\partial\theta',ind,arg=a,ferm=1)

def lam(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field('\\lambda',ind,arg=a)

def dlam(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field('\\partial\\lambda',ind,arg=a)

def eta(indup='',inddn=''):
    return field('\\eta',indup,inddn);

def delta(indup,inddn):
    return field('\\delta',indup,inddn);

def omega(ind,arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    return field('\\omega',inddn=ind,arg=a)

def T(arg):
    a = 'z_'+str(arg) if type(arg) is int else arg;
    out = exp([]);
    out.add(term([eta(inddn='m n'), Pi('m',a), Pi('n',a)],sign=-1,pref = '\half'))
    out.add(term([d('\\alpha',a),dtheta('\\alpha',a)],sign=-1))
    out.add(term([omega('\\alpha',a),dlam('\\alpha',a)]))

    return out


def V(lab,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg
    return multiply(A(lab,'\\alpha',a),lam('\\alpha',a));
    
def U(lab,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg

    Ui = exp([term([dtheta('\\alpha',a),A(lab,'\\alpha',a)]),
              term([Pi('m',a),Am(lab,'m',a)]),
              term([d('\\alpha',a),W(lab,'\\alpha',a)]),
              term([F(lab,'m n',a), N('m n',a)],pref = '\\half ')])
    return Ui

def Q(arg):
    if type(arg) is int:
        a = 'z_'+str(arg)
    else:
        a = arg

    return multiply(lam('\\alpha',a),d('\\alpha',a));

    

def G(lab,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg

    eout = exp([term([dtheta('\\alpha',a),A(lab,'\\alpha',a)]),
                term([Am(lab,'m',a),Pi('m',a)]),
                term([d('\\alpha',a),W(lab,'\\alpha',a)])]);
    return eout;

# Simplified UU OPE expression
def UU(lab1,lab2,arg=''):
    if type(arg) is int:
        a = 'z_'+str(arg)
    elif not arg:
        a = 'z_'+str(lab)
    else:
        a = arg

    eout = exp([multiply(term([W(lab1,'\\alpha',a),gamma('',r'm \alpha \beta'),
                               W(lab2,'\\beta',a)]),
                         exp([term(Pi('m',a)),
                              term([k(lab1,'','n'),N('n m',a)]),
                              term([k(lab2,'','n'),N('n m',a)])])),
                term([N('m n',a),F(lab1,'m p',a),F(lab2,'n p',a)],sign=-1),
                multiply(W(lab1,'\\alpha',a),symD(G(lab2,a),'\\alpha')),
                term([A(lab2,'m',a),W(lab1,'\\alpha',a),
                      gamma('m',r'\alpha \beta'),dtheta('\\beta',a)]),
                multiply(term([A(lab1,'m',a),k(lab2,indup='m')]),
                         multiply(U(lab2,a), term([],sign=-1))),
                multiply(W(lab2,'\\alpha',a),symD(G(lab2,a),'\\alpha')),
                term([A(lab1,'m',a),W(lab2,'\\alpha',a),
                      gamma('m',r'\alpha \beta'),dtheta('\\beta',a)],sign=-1),
                multiply(term([A(lab2,'m',a),k(lab1,indup='m')],sign=-1),
                         U(lab1,a))])

    return eout
