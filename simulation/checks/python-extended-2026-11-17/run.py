import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.integrate import trapezoid
from scipy.stats import beta
from statsmodels.nonparametric.smoothers_lowess import lowess
from pathlib import Path
P=Path(__file__).resolve().parent
D=np.load(P/'digitised_all.npy',allow_pickle=True).item()['Fig 1 (Study 1, in-role)']
TX=np.linspace(1,5,101); TY=np.interp(TX,D['lx'],D['ly']); TRANGE=TY.max()-TY.min()
rng=np.random.default_rng(314159)
N=474; reps=int(os.environ.get('REPS','100'))
rows=[];curves=[]
for scenario in ['linear','diminishing','plateau','downturn']:
 for i in range(reps):
  z=rng.normal(size=N); x=np.clip(np.round((3.335+.751*z)*7)/7,1,5)
  t=-.1724*z+np.sqrt(1-.1724**2)*rng.normal(size=N)
  if scenario=='linear': g=x
  elif scenario=='diminishing':g=np.log1p(2*np.maximum(x-1,0))
  elif scenario=='plateau':g=np.minimum(x,3.335)
  else:g=-(x-4.75)**2
  g=(g-g.mean())/g.std(); noise=rng.normal(size=N)
  # Calibrate observed correlation to published 0.2595 using same fixed noise draw.
  def response(a):
   lat=a*g+noise
   lat=(lat-lat.mean())/lat.std()
   return np.clip(np.round((3.695+.752*lat)*4)/4,1,5)
  lo,hi=0.,4.
  for _ in range(17):
   mid=(lo+hi)/2
   if np.corrcoef(x,response(mid))[0,1]<.2595:lo=mid
   else:hi=mid
  y=response((lo+hi)/2)
  for width in [1.,2.]:
   # Gaussian kernel on original units; regularization fixed at 0.05, not official KRLS.
   X=np.column_stack([x,t]); dist=((X[:,None,:]-X[None,:,:])**2).sum(axis=2)
   K=np.exp(-dist/(2*width**2));K.flat[::N+1]+=.05
   try: alpha=cho_solve(cho_factor(K,check_finite=False),y,check_finite=False)
   except Exception:continue
   K.flat[::N+1]-=.05
   deriv=np.sum(-(x[:,None]-x[None,:])/(width**2)*K*alpha[None,:],axis=1)
   # Official paper uses lowess; f=.8 and robust it=0 approximation
   sm=lowess(deriv,x,frac=.8,it=0,return_sorted=True)
   sx,ix=np.unique(sm[:,0],return_index=True);sy=sm[ix,1]
   # Match only observed support, avoid unobserved extrapolation.
   mask=(TX>=sx[0])&(TX<=sx[-1]);tx=TX[mask]; target=TY[mask]
   pred=np.interp(tx,sx,sy)
   rmse=np.sqrt(np.mean((pred-target)**2))
   # Also scale-free shape comparison (allows multiplicative y slope)
   scale=(pred@target)/(pred@pred) if pred@pred>1e-10 else 0
   shape_rmse=np.sqrt(np.mean((scale*pred-target)**2))
   # upper-end fall from last zero crossing to observed max
   fall=0.
   if sy[-1]<0:
    k=np.flatnonzero(sy>=0)
    start=int(k[-1]) if len(k) else 0
    fall=-trapezoid(np.minimum(sy[start:],0),sx[start:])
   # How close: slope + fall + peak and RMSE jointly, no post hoc hypothesis probabilities
   row=dict(scenario=scenario,rep=i,width=width,cor=np.corrcoef(x,y)[0,1],xmean=x.mean(),xsd=x.std(ddof=1),ymean=y.mean(),ysd=y.std(ddof=1),xmin=x.min(),xmax=x.max(),rmse=rmse,shape_rmse=shape_rmse,fall=fall,negative_endpoint=int(sy[-1]<0),endpoint_slope=sy[-1],slope_at_4=np.interp(4,sx,sy),peak_x=sx[np.argmax(sy)],scale=scale)
   rows.append(row)
   if i<3:curves.extend([dict(scenario=scenario,rep=i,width=width,x=a,estimated_slope=b,target_slope=c) for a,b,c in zip(tx,pred,target)])
 if (len(rows)>0):print(scenario,'complete',flush=True)
pd.DataFrame(rows).to_csv(P/'per_run.csv',index=False)
pd.DataFrame(curves).to_csv(P/'example_curves.csv',index=False)
df=pd.DataFrame(rows)
s=[]
for (sc,w),g in df.groupby(['scenario','width']):
 s.append(dict(scenario=sc,width=w,n=len(g),corr_mean=g.cor.mean(),rmse_median=g.rmse.median(),rmse_q25=g.rmse.quantile(.25),shape_rmse_median=g.shape_rmse.median(),negative_endpoint_pct=100*g.negative_endpoint.mean(),fall_at_least_published_pct=100*(g.fall>=.0184).mean(),joint_negative_and_fall_pct=100*((g.negative_endpoint==1)&(g.fall>=.0184)).mean(),rmse_below_0_15_pct=100*(g.rmse<.15).mean(),rmse_below_0_2_pct=100*(g.rmse<.2).mean()))
pd.DataFrame(s).to_csv(P/'summary.csv',index=False)
print(pd.DataFrame(s).to_string(index=False))
