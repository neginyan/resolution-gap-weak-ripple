import numpy as np, math, random
from decimal import Decimal, getcontext
getcontext().prec=60
rng=np.random.default_rng(1)
def ulp_err(approx, exact):
    return abs(Decimal(float(approx))-exact)/Decimal(np.spacing(abs(float(approx))))
res={}
N=20000
xs=np.concatenate([rng.uniform(-700,700,N), rng.uniform(-1,1,N)*rng.choice([1e-8,1e-3,1,10],N)])
m=0
for x in xs:
    m=max(m,ulp_err(np.exp(np.float64(x)), Decimal(float(x)).exp()))
res['exp']=m
pos=np.exp(rng.uniform(-700,700,2*N))
m=0
for x in pos: m=max(m,ulp_err(np.log(x), Decimal(float(x)).ln()))
res['log']=m
sm=rng.uniform(-0.9,5,N)*rng.choice([1e-12,1e-6,1e-2,1],N)
m=0;m2=0
for x in sm:
    d=Decimal(float(x))
    m=max(m,ulp_err(np.log1p(x),(1+d).ln()))
    m2=max(m2,ulp_err(np.expm1(x),d.exp()-1))
res['log1p']=m;res['expm1']=m2
m=0
for x in pos: m=max(m,ulp_err(np.sqrt(x),Decimal(float(x)).sqrt()))
res['sqrt']=m
# also math module
m=0
for x in xs[:5000]: m=max(m,ulp_err(math.exp(float(x)),Decimal(float(x)).exp()))
res['math.exp']=m
print(np.__version__, {k:round(float(v),3) for k,v in res.items()})
np.show_config() if False else None
