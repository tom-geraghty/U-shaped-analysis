import pandas as pd
from math import sqrt
from pathlib import Path
source=Path(__file__).resolve().parent.parent/'eldor_python_extended'/'per_run.csv'
d=pd.read_csv(source)
def wilson(k,n,z=1.95996398454):
 p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return 100*(c-h),100*(c+h)
rows=[]
for (sc,w),g in d.groupby(['scenario','width']):
 for threshold in [.0092,.0184,.0368]:
  test=(g.slope_at_4>0)&(g.endpoint_slope<0)&(g.fall>=threshold)
  k=int(test.sum());n=len(g);l,u=wilson(k,n)
  rows.append(dict(scenario=sc,width=w,threshold=threshold,k=k,n=n,rate_pct=round(100*k/n,1),ci95_low_pct=round(l,1),ci95_high_pct=round(u,1)))
out=Path(__file__).resolve().parent
pd.DataFrame(rows).to_csv(out/'joint_claim_rates.csv',index=False)
print(pd.DataFrame(rows).query('threshold==0.0184').to_string(index=False))
