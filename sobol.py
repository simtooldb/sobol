#!/usr/bin/env python

import math, numpy as np, re
from scipy.stats.qmc import Sobol as sbl
from scipy.stats import qmc
from itertools import product

def parseBatchParams (b):
    ''' read a batch.py file for NetPyNE param search; returning list of (name, valueList, [indexed]) where optional indexed means to use all the values '''
    with open(b,'r') as fb: lines = fb.readlines()
    p = re.compile('''\s +params[^a-z]+([^]']+)'\]\s *=\s *(\[[^]]+\])\s *#*\s *(indexed)*''') # 'indexed' is the keyword to NOT sobolize eg amp or cellnum
    bl, pl = [(i, p.match(l)) for i,l in enumerate(lines)], [] # bl: lines that match regexp
    for i,m in bl:
        if m:
            try: 
                pl.append((m.group(1), eval(m.group(2)), m.group(3))) # strings: name, valueList, [indexed]
            except Exception as e:
                print("ERROR >>>%s<<<\n\tunable to evaluate '%s':\n\tline %d:'%s'"%(e, m.group(2), i, m.string.strip()))
    return pl
        
def sobcall (pl, num, seed=1234):
    ''' determine the min, max of sobolized params and do the combos with indexed params '''
    labels, mins, maxs = [], [], []
    for x in pl:
        if not x[2]:
            labels.append(x[0]); mins.append(min(x[1])); maxs.append(max(x[1]))
    print('mins/maxs for %s: '%(labels), mins, maxs)
    vals = sob(len(mins), num, seed=seed)
    svals = qmc.scale(vals, mins, maxs) # only for those that are not 'indexed'
    
    return svals
    
def sob (dim=4, num=4096, f=None, seed=1234):
    sm = sbl(d=dim, scramble=True, seed=seed)
    m = math.floor(math.log(num)/math.log(2) + 0.99) # round up to nearest power of 2
    if 2**m != num: print('%d samples (2^%d ; %d requested)'%(2**m, m, num))
    vals = sm.random_base2(m=m) # 2^m points
    if f: 
        print('saving raw sobol values to %s'%(f))
        np.savetxt(f, vals, delimiter=',', fmt='%.10f')
    return vals

def getArgs ():
    import argparse
    msg = "Generate at least cnt lines of sobol samples of size dim and place in file sobol.out"
    parser = argparse.ArgumentParser(description=msg, )    
    parser.add_argument('dim', nargs='?', type=int, default=4, help='dim of space to be sampled')
    parser.add_argument('cnt', nargs='?', type=int, default=10, help='num of samples: will be rounded up to a power of 2')
    parser.add_argument("-r", default='sobol.csv', help='raw output from sobol call (default ./sobol.csv)')
    parser.add_argument("-f", default='params.csv', help='file for saving param lists (default ./params.csv)')
    parser.add_argument("-s", default=1234, type=int, help='seed')
    parser.add_argument("-v", action='store_true', default=False, help='output to terminal')
    parser.add_argument("-b", default='batch.py', help='name of batchfile with "params" ranges (default ./batch.py)')
    return parser.parse_args()

if __name__ == '__main__':
    ag = getArgs()
    svals = sobcall(parseBatchParams(ag.b), ag.cnt)
