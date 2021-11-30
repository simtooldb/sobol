# sobol

to be used with netpyne batch 'list' functionality -- generates a list based on 'params[]' identified in '-b' file with additional key word comment '# indexed'
to indicate use of discrete values (from 'grid') rather than continuous values between min,max of that param (e.g., IClamp.amp would typically be indexed)

usage: sobol.py [-h] [--cnt [CNT]] [-f F] [-s S] [-b B]

Generate at least cnt lines of sobol samples based on params definitions in batch.py. Params may be indicated as '# indexed' in which case the batch values given will be used. Other
params will be scaled using a sobol quasi-monte carlo distribution

optional arguments:
  -h, --help   show this help message and exit
  --cnt [CNT]  num of samples: will be rounded up to a power of 2
  -f F         file for saving param lists (default ./params.csv)
  -s S         seed
  -b B         name of batchfile with "params" ranges (default ./batch.py)
