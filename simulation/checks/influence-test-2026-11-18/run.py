import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t as student_t
out = Path(__file__).resolve().parent
rng = np.random.default_rng(20261118)
N, B, knot = 474, 1200, 4.75

def fit(x, y, tenure):
    X = np.column_stack([np.ones(len(x)), x, np.maximum(x-knot, 0), tenure])
    if np.linalg.matrix_rank(X) < 4:
        return np.nan, np.nan
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    c = np.array([0, 1, 1, 0])
    slope = float(c @ beta)
    inv = np.linalg.inv(X.T @ X)
    resid = y - X @ beta
    se = np.sqrt((resid @ resid)/(len(y)-4) * (c @ inv @ c))
    return slope, bool(student_t.cdf(slope/se, df=len(y)-4) < .05)

records = []
for truth, fall in [('plateau', 0), ('published_size_decline', .0184)]:
    for rep in range(B):
        z = rng.normal(size=N)
        x = np.clip(np.round((3.335+.751*z)*7)/7, 1, 5)
        tenure = -.1724*z + np.sqrt(1-.1724**2)*rng.normal(size=N)
        mu = .27*(np.minimum(x,knot)-3.335) - (fall/.25)*np.maximum(x-knot,0)
        y = np.clip(np.round((3.695+mu+.72*rng.normal(size=N))*4)/4,1,5)
        base, sig = fit(x,y,tenure)
        tail = np.where(x > knot)[0]
        if len(tail)<5 or not np.isfinite(base):
            continue
        for k in [1,3,5]:
            idx = rng.choice(tail,k,replace=False)
            for direction in [-1,1]:
                y2 = y.copy()
                y2[idx] = np.clip(y2[idx] + direction*.25,1,5)
                changed, sg = fit(x,y2,tenure)
                records.append((truth,rep,len(tail),k,f'shift_{direction:+d}_0.25',base,changed,sig,sg))
            for label,idx2 in [('delete_random',idx),('delete_highest_x',tail[np.argsort(x[tail])[-k:]])]:
                keep = np.ones(N,bool)
                keep[idx2] = False
                changed, sg = fit(x[keep],y[keep],tenure[keep])
                records.append((truth,rep,len(tail),k,label,base,changed,sig,sg))

df = pd.DataFrame(records,columns=['truth','rep','n_tail','n_modified','intervention','base_slope','new_slope','base_sig','new_sig'])
df['unestimable'] = ~np.isfinite(df.new_slope)
df['sign_flip'] = np.where(df.unestimable,np.nan,(df.base_slope<0)!=(df.new_slope<0))
df['sig_flip'] = np.where(df.unestimable,np.nan,df.base_sig!=df.new_sig)
df['abs_slope_change'] = abs(df.new_slope-df.base_slope)
summ = df.groupby(['truth','n_modified','intervention']).agg(n=('rep','size'),unestimable_pct=('unestimable',lambda s:100*s.mean()),sign_flip_pct=('sign_flip',lambda s:100*s.mean()),sig_flip_pct=('sig_flip',lambda s:100*s.mean()),median_abs_slope_change=('abs_slope_change','median')).reset_index()
df.to_csv(out/'per_run.csv',index=False)
summ.to_csv(out/'summary.csv',index=False)
print(summ.to_string(index=False))
