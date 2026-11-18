"""Study 1 oracle-knot power experiment. Requires numpy, pandas, scipy.
Run: python run.py
"""
import numpy as np
import pandas as pd
from scipy.stats import t as tdist
from pathlib import Path
out=Path(__file__).resolve().parent
rng=np.random.default_rng(20261118)
N=474; B=5000; knot=4.75; falls=[0,.0184,.08,.20,.50,1.00]
rows=[]
for fall in falls:
    rejected=negative=0; n_tail=[]; slopes=[]; corrs=[]
    for j in range(B):
        z=rng.standard_normal(N)
        x=np.clip(np.round((3.335+.751*z)*7)/7,1,5)
        tenure=-.1724*z+np.sqrt(1-.1724**2)*rng.standard_normal(N)
        mu=.27*(np.minimum(x,knot)-3.335)-(fall/(5-knot))*np.maximum(x-knot,0)
        y=np.clip(np.round((3.695+mu+.72*rng.standard_normal(N))*4)/4,1,5)
        X=np.column_stack([np.ones(N),x,np.maximum(x-knot,0),tenure])
        inv=np.linalg.inv(X.T@X); beta=inv@(X.T@y)
        residual=y-X@beta; s2=(residual@residual)/(N-4)
        c=np.array([0,1,1,0]); se=np.sqrt(s2*(c@inv@c))
        slope=c@beta; p=tdist.cdf(slope/se,df=N-4)
        rejected+=p<.05; negative+=slope<0
        n_tail.append((x>knot).sum()); slopes.append(slope); corrs.append(np.corrcoef(x,y)[0,1])
    rows.append(dict(true_fall_points=fall,replicates=B,significantly_negative_pct=100*rejected/B,
       negative_estimate_pct=100*negative/B,median_estimated_tail_slope=np.median(slopes),
       median_n_above_4_75=np.median(n_tail),mean_n_above_4_75=np.mean(n_tail),mean_xy_correlation=np.mean(corrs)))
    print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(out/'oracle_power_results.csv',index=False)
