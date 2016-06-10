from copy import deepcopy as dcopy

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
