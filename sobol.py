#!/usr/bin/env python

import math, re, csv
from scipy.stats import qmc
from itertools import product

def parseBatchParams (b):
    ''' read a batch.py file for NetPyNE param search; returning list of (name, valueList, [indexed]) where optional indexed means to use all the values '''
    with open(b,'r') as fb: lines = fb.readlines()
    p = re.compile(r'''\s +params[^a-z]+([^]']+)'\]\s *=\s *(\[[^]]+\])\s *#*\s *(indexed)*''') # 'indexed' is the keyword to NOT sobolize eg amp or cellnum
    bl, pl = [(i, p.match(l)) for i,l in enumerate(lines)], [] # bl: lines that match regexp
    for i,m in bl:
        if m:
            try: 
                pl.append((m.group(1), eval(m.group(2)), m.group(3))) # strings: name, valueList, [indexed]
            except Exception as e:
                print(f"ERROR >>>{e}<<<\n\tunable to evaluate '{m.group(2)}':\n\tline {i}: {m.string.strip()}")
    return pl

def sobcall (pl, num, seed=33):
    ''' determine the min, max of sobolized params and do the combos with indexed params '''
    labels, Mins, Maxs, ilabels, ivals = [],[],[],[],[]
    for x in pl:
        if (not x[2]):
            labels.append(x[0]); Mins.append(min(x[1])); Maxs.append(max(x[1]))
        else:
            ilabels.append(x[0]); ivals.append(x[1])
    print(f'Mins/Maxs for {labels}: {Mins}, {Maxs}')
    sobolVals = sob(len(Mins), num, seed=seed)
    scaledVals = qmc.scale(sobolVals, Mins, Maxs) # only for those that are not 'indexed'
    allVals = list(zip(*scaledVals)) + ivals
    combos = list(product(*allVals))
    combos.insert(0, labels+ilabels)
    return combos

def sob (dim=4, num=4096, seed=1234):
    sm = qmc.Sobol(d=dim, scramble=True, seed=seed)
    m = math.floor(math.log(num)/math.log(2) + 0.99) # round up to nearest power of 2
    if 2**m != num: print(f'{2**m} samples (2^{m} ; {num} requested)')
    vals = sm.random_base2(m=m) # 2^m points
    return vals

def getArgs ():
    import argparse
    msg = '''Generate at least cnt lines of sobol samples based on params definitions in batch.py.
             Params may be indicated as '# indexed' in which case the batch values given will be used.
             Other params will be scaled using a sobol quasi-monte carlo distribution'''
    parser = argparse.ArgumentParser(description=msg, )
    # parser.add_argument('--dim', nargs='?', type=int, default=4, help='dim of space to be sampled')
    parser.add_argument('--cnt', nargs='?', type=int, default=10, help='num of samples: will be rounded up to a power of 2')
    # parser.add_argument("-r", default='sobol.csv', help='raw output from sobol call (default ./sobol.csv)')
    parser.add_argument("-f", default='params.csv', help='file for saving param lists (default ./params.csv)')
    parser.add_argument("-s", default=1234, type=int, help='seed')
    # parser.add_argument("-v", action='store_true', default=False, help='output to terminal')
    parser.add_argument("-b", default='batch.py', help='name of batchfile with "params" ranges (default ./batch.py)')
    return parser.parse_args()

    np.savetxt(ag.f, allvals, delimiter=',', fmt='%.10f')

if __name__ == '__main__':
    ag = getArgs()
    allvals = sobcall(parseBatchParams(ag.b), ag.cnt)
    with open(ag.f, 'w') as f:
        wr = csv.writer(f)
        wr.writerows(allvals)
