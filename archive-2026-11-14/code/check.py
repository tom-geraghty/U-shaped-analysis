# Reproduce the paper's simple slopes from its own coefficients, then derive what else they imply.
def slopes(b1,b2,b4,b5,psm,pss,cam,cas,label):
    print(label)
    for cl,ca in (("lowCA",cam-cas),("highCA",cam+cas)):
        lin=b1+b4*ca; quad=b2+b5*ca
        tp = -lin/(2*quad)
        s=[round(lin+2*quad*x,3) for x in (psm-pss,psm,psm+pss)]
        print(f"  {cl}: slopes at -1SD/mean/+1SD PS = {s}; turning point PS={tp:.2f} ({'max' if quad<0 else 'min'})")
# Study 4 Table 7 model 3; Table 6 means
slopes(13.73,-1.143,-2.978,0.252,5.286,0.824,3.827,0.900,"Study 4 unit performance (paper: low 1.384/0.714/0.045; high 0.07/0.147/0.224)")
# Study 5 Table 10 model 3 (service) and 6 (productivity); Table 9 means
slopes(10.39,-1.771,-2.390,0.403,3.390,0.425,4.536,0.409,"Study 5 service (paper: low -0.12/-0.212/-0.305; high -0.121/0.066/0.254)")
slopes(703.2,-116.0,-159.8,25.84,3.390,0.425,4.536,0.409,"Study 5 productivity (paper: low -11.62/-19.55/-27.49; high -17.0/-6.97/3.06)")
# Study 4: predicted gap between high-CA and low-CA units across PS
b3,b4,b5=8.508,-2.978,0.252; d=1.8
print("Study 4: predicted (highCA - lowCA) unit performance by PS level:")
for x in (2.5,3.5,4.462,5.286,6.11,7.0):
    print(f"  PS={x}: {d*(b3+b4*x+b5*x*x):+.2f}")
print("Quadratic checks: turning points (min unless noted)")
for lab,b1,b2 in (("S1",-0.824,0.166),("S2 service",-0.597,0.0858),("S2 in-role",-0.212,0.0491),("S3",-0.803,0.129),("S4 (max)",1.780,-0.133),("S5 service",-0.632,0.0937)):
    print(f"  {lab}: {-b1/(2*b2):.2f}")
# effect sizes Study 5 low CA: 1 SD of PS
print("S5 lowCA at mean: 1SD PS ->", round(-0.212*0.425,3),"service pts =",round(-0.212*0.425/0.487,2),"SD;", round(-19.553*0.425,1),"prod =",round(-19.553*0.425/93.778,2),"SD")
import statistics as st
from math import erf,sqrt
print("S4 peak 6.69 is z=",round((6.69-5.286)/0.824,2))
