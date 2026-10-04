"""Saltelli/Jansen Sobol indices (first-order and total) and local elasticities for the marked process."""
import numpy as np, json, time
from scipy.stats import qmc
import marked_process as mp

# inputs: name, type, range (log10 bounds for 'log', linear for 'lin', categorical list)
INPUTS=[('N_C','lin',(4,20)),('lam_sig','lin',(0.2,1.2)),('beta0','log',(-2,1)),('r0','log',(1.5,3.0)),
        ('L_C_med','log',(-2,0.5)),('L_info','log',(-2,0.7)),('p_M','lin',(0,1)),('lam_mult_M','lin',(1,5)),('fam','cat',['lognormal','hardstep','triggered'])]
OUTS=['mutual_now','oneway_frac','scc_max','isolated']
def decode(u):
    kw={}
    for (n,t,r),x in zip(INPUTS,u):
        if t=='lin': kw[n]=r[0]+x*(r[1]-r[0])
        elif t=='log': kw[n]=10**(r[0]+x*(r[1]-r[0]))
        else: kw[n]=r[int(min(x*len(r),len(r)-1))]
    kw['N_C']=int(round(kw['N_C'])); return kw
def model(u, runs=40, seed=None):
    kw=decode(u); fam=kw.pop('fam')
    if seed is not None: mp.rng=np.random.default_rng(seed)
    s,_=mp.run(fam,runs=runs,**kw)
    return np.array([np.nanmean(s[q]) if np.isfinite(np.nanmean(s[q])) else 0.0 for q in OUTS])

def sobol(Ns=256, runs=40):
    k=len(INPUTS)
    sampler=qmc.Sobol(d=2*k,scramble=True,seed=3); X=sampler.random(Ns)
    A=X[:,:k]; B=X[:,k:]
    fA=np.array([model(a,runs) for a in A]); fB=np.array([model(b,runs) for b in B])
    S1=np.zeros((k,len(OUTS))); ST=np.zeros((k,len(OUTS)))
    V=np.var(np.vstack([fA,fB]),axis=0)
    for i in range(k):
        ABi=A.copy(); ABi[:,i]=B[:,i]
        fABi=np.array([model(x,runs) for x in ABi])
        S1[i]=np.mean(fB*(fABi-fA),axis=0)/V                    # Saltelli 2010 estimator
        ST[i]=0.5*np.mean((fA-fABi)**2,axis=0)/V                 # Jansen estimator
        print('input',INPUTS[i][0],'done',flush=True)
    return S1,ST,V
def elasticities(runs=400, h=0.2):
    base={'N_C':10,'lam_sig':0.5,'beta0':0.4,'r0':150.0,'L_C_med':0.05,'L_info':0.5,'p_M':0.5,'lam_mult_M':3.0}
    def f(kw,seed=11):
        mp.rng=np.random.default_rng(seed); s,_=mp.run('lognormal',runs=runs,**kw)
        return np.array([np.nanmean(s[q]) if np.isfinite(np.nanmean(s[q])) else 0.0 for q in OUTS])
    y0=f(base); E={}
    for n in base:
        up=dict(base); dn=dict(base)
        if n=='N_C': up[n]=12; dn[n]=8; dl=np.log(12/8)
        elif n=='p_M': up[n]=0.6; dn[n]=0.4; dl=np.log(0.6/0.4)
        else: up[n]=base[n]*(1+h); dn[n]=base[n]*(1-h); dl=np.log((1+h)/(1-h))
        yu=f(up); yd=f(dn)
        E[n]=((np.log(np.maximum(yu,1e-6))-np.log(np.maximum(yd,1e-6)))/dl).tolist()   # d ln Y / d ln x
    return y0.tolist(),E

if __name__=='__main__':
    t=time.time()
    S1,ST,V=sobol()
    print('sobol time',time.time()-t)
    res={'inputs':[i[0] for i in INPUTS],'outputs':OUTS,'S1':S1.tolist(),'ST':ST.tolist(),'V':V.tolist()}
    for j,q in enumerate(OUTS):
        order=np.argsort(-ST[:,j])
        print(q,' | '.join(f"{INPUTS[i][0]}: S1={S1[i,j]:.2f} ST={ST[i,j]:.2f}" for i in order))
    y0,E=elasticities(); res['baseline']=y0; res['elasticities']=E
    print('baseline',dict(zip(OUTS,np.round(y0,3))))
    for n,e in E.items(): print(f"{n:11s} dlnY/dlnx: "+'  '.join(f"{q}={v:+.2f}" for q,v in zip(OUTS,e)))
    json.dump(res,open('sobol.json','w'),indent=1); print('total',time.time()-t)
