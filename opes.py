from copy import deepcopy as dcopy
from itertools import *

# worldsheet field class
# indices can be a list of strings or a single string separated by spaces
class field:
    def __init__(self,name,indup=[],inddn=[],arg='',ferm=0,opeflag=False):
        self.name = name;
        self.indup = indup.split() if type(indup) is str else indup;
        self.inddn = inddn.split() if type(inddn) is str else inddn;
        self.arg = arg;
        self.ferm = ferm;
        self.opeflag = opeflag;

    # LaTeX output of this field, optionally display arguments
    def str(self,args=True):
        indups = ' '.join(self.indup);
        inddns = ' '.join(self.inddn);

        arg = ''
        if self.arg and args:
            arg = '('+self.arg+')';
            
        s = self.name;
        if indups:
            s+= '^{'+indups+'}'
        if inddns:
            s+='_{'+inddns+'}'

        arg = ''
        if self.arg and args:
            arg = '('+self.arg+')';

        s+= arg;
        return s

# super Yang-Mills fields additionally have labels, and can be acted on by
# the supersymmetric derivative
class symfield(field):
    def __init__(self,name,label,indup=[],inddn=[],arg='',ferm=0,opeflag=False,
                 derivind = []):
        field.__init__(self,name,indup,inddn,arg,ferm,opeflag)
        self.label = label;
        self.derivind = derivind.split() if type(derivind) is str else derivind;

    def str(self,args=True):
        indups = ' '.join([str(self.label)+'\\,']+self.indup);
        inddns = ' '.join(self.inddn);

        arg = ''
        if self.arg and args:
            arg = '('+self.arg+')';

        ds = '';

        for dind in self.derivind:
            ds += 'D_{'+dind+'}'
            
        s = ds + self.name+'^{'+indups+'}'
        if inddns:
            s += '_{'+inddns+'}'
        s += arg;
        return s



# a term is a product of (SYM) fields, given in a list.
# Potentially with a sign, prefactor,
# and denominator. Note that no "math" is performed on the pref or denom;
# they are just handled as strings. 
class term:
    def __init__(self,fields,sign=1,pref='',denom=''):
        self.sign = sign;
        self.fields = self.__setfields(fields);
        self.denom = denom; # denominato string
        self.pref = pref # prefactor string

    # push the field at index ind back by one position, change sign if
    # exchanging fermions
    def shiftback(self,ind):
        if ind >= len(self.fields)-1:
            pass
        else:
            f1 = self.fields[ind];
            f2 = self.fields[ind+1];

            self.sign *= (-1)**(f1.ferm * f2.ferm);

            self.fields[ind], self.fields[ind+1] = f2, f1;

    
    # push the field at index ind back by one position, change sign if
    # exchanging fermions
    def shiftforward(self,ind):
        if ind <= 0:
            pass
        else:
            f2 = self.fields[ind];
            f1 = self.fields[ind-1];

            self.sign *= (-1)**(f2.ferm * f1.ferm);

            self.fields[ind], self.fields[ind-1] = f1, f2;

    # shift the field at index ind n places left (negative) or right (positive)
    def shift(self,ind,n):
        if n<0:
            for i in range(ind,ind+n,-1):
                self.shiftforward(i);
        elif n>0:
            for i in range(ind,ind+n):
                self.shiftback(i);

    # shift the field at position ind to destination dest
    def shiftto(self,ind,dest):
        if dest <0:
            dest = len(self.fields);
            
        n = dest - ind;
        self.shift(ind,n);

    def __setfields(self,fields):
        fins = fields if type(fields) is list else [fields];
        fnews = [];

        # make copies of all fields so we don't change the field objects
        for fin in fins:
            fnew = dcopy(fin);
            fnews.append(fnew);
            
        return fnews;

    # return list of all field indices in the term
    def inds(self):
        indups = [];
        inddns = [];

        for f in self.fields:
            indups += f.indup;
            inddns += f.inddn;
            if type(f) is symfield and f.derivind:
                inddns += f.derivind
                
        return indups + inddns;

    # return list of all dummy (repeated once up and once down) indices
    def dummies(self):
        indups = [];
        inddns = [];
        dummies = [];
        for f in self.fields:
            indups += f.indup;
            inddns += f.inddn;
            if type(f) is symfield:
                inddns += f.derivind;

        for i in indups:
            if i in inddns:
                dummies.append(i);

        return dummies

    # change a dummy index to the new one if provided, else the next one in the
    # appropriate alphabet
    def changedummy(self,old,new='',excl =['O','o','\\omicron']): # exclude some indices
        ds = self.dummies();
        if new:
            for f in self.fields:
                if old in f.indup:
                    f.indup[f.indup.index(old)] = new;
                if old in f.inddn:
                    f.inddn[f.inddn.index(old)] = new;
        else:
            alph = alphabet(old);
            i = alph.index(old);
            
            for f in self.fields:
                swapped = False;
                if old in f.indup:
                    inew = i+1 if i+1 != len(alph) else 0;

                    while not swapped:
                        n = alph[inew];
                        if (n not in ds) and (n not in excl):
                            f.indup[f.indup.index(old)] = n;
                            swapped = True;
                        else:
                            inew += 1;
                            
            for f in self.fields:
                if old in f.inddn:
                    f.inddn[f.inddn.index(old)] = n;

    # provide a new dummy from the appropriate alphabet that has not been used
    # yet in the term
    def choosedummy(self,alph,n=1):
        inds = self.inds();
        ds = [];
        for l in alph:
            if (l not in inds) and (len(ds) < n):
                ds.append(l);

        return ds[0] if n == 1 else ds;

    # replace the fields at positions ind1 and ind1+1 with the fields in
    # the opeterm, and combine denominators, signs, and prefactors
    def opeinsert(self,ind1,opeterm):
        self.sign *= opeterm.sign;
        self.pref += opeterm.pref;
        self.denom += opeterm.denom;

        ds = opeterm.dummies();
        
        for d in ds:
            if d in self.inds():
                opeterm.changedummy(d,excl = self.inds() +['o','\\ommicron','O'])
      
        nfields = self.fields[0:ind1] + opeterm.fields + self.fields[ind1+2:]
        self.fields = nfields;

    # provide LaTeX output string for the entire term
    def str(self,args=True):
        s = ''
        if self.sign == 1:
            s = '+';
        elif self.sign == -1:
            s = '-';
            
        if self.pref:
            s += self.pref;
            
        if not self.denom:
            for f in self.fields:
                s += f.str(args);
        else:
            s += '\\frac{'
            for f in self.fields:
                s += f.str(args);
            s +='}{'+self.denom+'}'
        return s

# an expression is a list of terms with nice LaTeX output
class exp:
    def __init__(self,terms):
        self.terms = self.__setterms(terms);

    # make copies of all terms so we don't change the original term objects
    def __setterms(self,terms):
        if type(terms) is term:
            ts = [dcopy(terms)];
        elif type(terms) is exp:
            ts = [dcopy(t) for t in terms.terms];
        elif type(terms) is list:
            ts = [];
            for t in terms:
                if type(t) is term:
                    ts += [dcopy(t)];
                elif type(t) is exp:
                    ts += [dcopy(t) for t in t.terms];
        return ts;
                    
    # provide LaTeX output. Optionally display field arguments, start a new line
    # every splitevery terms, and align all the lines using &
    def str(self,args=True,splitevery=3,align=True):
        n = 0
        s = '';
        if not self.terms:
            return s

        al = r'&\ ' if align else '';
        for t in self.terms:
            n +=1;
            nl = '';
            if splitevery:
                nl = r'\\'+'\n'+al if n%splitevery == 0 else '';
                
                
            s += t.str(args)+nl;

        if s[0] == '+':
            s = s[1:]; # remove leading +

        if align:
            s = al + s;

        return s

    # add more terms to the expression
    def add(self,e):
        if type(e) is term:
            self.terms += [t];
        elif type(e) is exp:
            self.terms += e.terms;

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

# returns the alphabet to which l belongs as a list, e.g. for choosing dummies
def alphabet(l):
    latin = list('abcdefghijklmnpqrstuvwxyz') 
    caplat = list('ABCDEFGHIJKLMNPQRSTUVWXYZ')
    greek = ['\\alpha','\\beta','\\gamma','\\delta','\\epsilon','\\zeta',
             '\\eta','\\theta','\\iota','\\kappa','\\lambda','\\mu','\\nu',
             '\\xi','\\pi','\\rho','\\sigma','\\tau','\\upsilon','\\phi',
             '\\chi','\\psi','\\omega']

    if l in caplat:
        alph =  caplat
    elif l in greek:
        alph =  greek
    else:
        alph = latin

    i = alph.index(l);
    out = alph[alph.index(l):] + alph[0:alph.index(l)]
    return out

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

        return fields and labels and derivinds

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
        
    

##def diff(first, second):
##        second = set(second)
##        return [item for item in first if item not in second]

# return the powerset of a given set (from itertools documentation)
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

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


    
    
