import csv
from scipy.stats.qmc import Sobol as sbl
import scipy.stats.qmc as qmc

sm = sbl(d=5, scramble=True) # should be power of 2
vals = sm.random_base2(m=4) # 2^m points

