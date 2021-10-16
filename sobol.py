#!/usr/bin/env python

import math, numpy as np, re
from scipy.stats.qmc import Sobol as sbl
from scipy.stats import qmc
 
def parseBatchParams(b):
    with open(b,'r') as fb: lines=fb.readlines()
    params = [x for x in lines if 
    return lines
        
def sob (dim=4, num=4096, f=None, seed=1234):
    sm = sbl(d=dim, scramble=True, seed=seed)
    m = math.floor(math.log(num)/math.log(2) + 0.99) # round up to nearest power of 2
    if 2**m != num: print('Saving %d samples (2^%d ; %d requested)'%(2**m, m, num))
    vals = sm.random_base2(m=m) # 2^m points
    if f: np.savetxt(f, vals, delimiter=',', fmt='%.10f')
    return vals

def getArgs ():
    import argparse
    msg = "Generate at least cnt lines of sobol samples of size dim and place in file sobol.out"
    parser = argparse.ArgumentParser(description=msg, )    
    parser.add_argument('dim', nargs='?', type=int, default=4, help='dim of space to be sampled')
    parser.add_argument('cnt', nargs='?', type=int, default=16, help='num of samples: will be rounded up to a power of 2')
    parser.add_argument("-o", default='sobol.csv', help='output filename (default ./sobol.csv)')
    parser.add_argument("-s", default=1234, type=int, help='seed')
    parser.add_argument("-v", action='store_true', default=False, help='output to terminal')
    parser.add_argument("-b", default='batch.py', help='name of batchfile (default ./batch.py')
    return parser.parse_args()

if __name__ == '__main__':
    ag = getArgs()
    parseBatchParams(ag.b)
    vals = sob(dim=ag.dim, num=ag.cnt, seed=ag.s, f=ag.o if ag.o else None)
    if ag.v: print(vals)
