"""Within-family Sobol (first/total, bootstrap CIs), group indices, stochastic variance fraction."""
import numpy as np, json, time
from scipy.stats import qmc
import marked_process as mp
INPUTS=[('N_C','lin',(4,20)),('lam_sig','lin',(0.2,1.2)),('beta0','log',(-2,1)),('r0','log',(1.5,3.0)),
        ('L_C_med','log',(-2,0.5)),('L_info','log',(-2,0.7)),('p_M','lin',(0,1)),('lam_mult_M','lin',(1,5))]
GROUPS={'abundance':['N_C'],'persistence':['L_C_med','L_info','p_M'],'rates':['lam_sig','lam_mult_M'],'communication':['beta0','r0']}
OUTS=['p_contact','mutual_now','oneway_frac','gscc','isolated']
def decode(u):
    kw={}
    for (n,t,r),x in zip(INPUTS,u):
        kw[n]= r[0]+x*(r[1]-r[0]) if t=='lin' else 10**(r[0]+x*(r[1]-r[0]))
    kw['N_C']=int(round(kw['N_C'])); return kw
def model(fam,u,runs=40):
    s,_=mp.run(fam,runs=runs,**decode(u))
    with np.errstate(all='ignore'):
        y=np.array([np.nanmean(s[q]) for q in OUTS])      # NaN where undefined at every run (hurdle: f_one^+ with no contacted pair)
        v=np.array([np.nanvar(s[q]) for q in OUTS])
    return y,v
def analyse(fam,Ns=192,runs=40,boot=300):
    k=len(INPUTS); X=qmc.Sobol(d=2*k,scramble=True,seed=5).random(Ns); A=X[:,:k]; B=X[:,k:]
    fA=[];vA=[]
    for a in A: y,v=model(fam,a,runs); fA.append(y); vA.append(v)
    fA=np.array(fA); vA=np.array(vA); fB=np.array([model(fam,b,runs)[0] for b in B])
    n_undef=int(np.isnan(fA).sum(0)[2]+np.isnan(fB).sum(0)[2])
    names=[i[0] for i in INPUTS]; cols={n:i for i,n in enumerate(names)}
    sets={**{n:[cols[n]] for n in names}, **{'G:'+g:[cols[n] for n in m] for g,m in GROUPS.items()}}
    fAB={}
    for key,idx in sets.items():
        ABi=A.copy(); ABi[:,idx]=B[:,idx]; fAB[key]=np.array([model(fam,x,runs)[0] for x in ABi])
    # impute undefined f_one^+ points by the mean of the defined ones (removes, rather than adds, variance); counts reported
    mu=np.nanmean(np.vstack([fA,fB]),axis=0)
    fA=np.where(np.isnan(fA),mu,fA); fB=np.where(np.isnan(fB),mu,fB); vA=np.where(np.isnan(vA),0.0,vA)
    for k in fAB: fAB[k]=np.where(np.isnan(fAB[k]),mu,fAB[k])
    def indices(rows):
        fa,fb=fA[rows],fB[rows]; V=np.var(np.vstack([fa,fb]),axis=0)
        S1={k:np.mean(fb*(fAB[k][rows]-fa),axis=0)/V for k in fAB}
        ST={k:0.5*np.mean((fa-fAB[k][rows])**2,axis=0)/V for k in fAB}
        return S1,ST,V
    S1,ST,V=indices(np.arange(Ns))
    rng=np.random.default_rng(1); bs1={k:[] for k in fAB}; bst={k:[] for k in fAB}
    for _ in range(boot):
        rows=rng.integers(0,Ns,Ns); s1,st,_=indices(rows)
        for k in fAB: bs1[k].append(s1[k]); bst[k].append(st[k])
    ci=lambda d: {k:np.percentile(np.array(v),[2.5,97.5],axis=0).tolist() for k,v in d.items()}
    # stochastic fraction: E_X[Var(Y|X)] / Var(Y)  (Var(Y|X) from the 40 replications at A points)
    R_stoch=np.mean(vA,axis=0)/(np.var(fA,axis=0)+np.mean(vA,axis=0))
    return {'S1':{k:v.tolist() for k,v in S1.items()},'ST':{k:v.tolist() for k,v in ST.items()},'S1_ci':ci(bs1),'ST_ci':ci(bst),'V':V.tolist(),'R_stoch':R_stoch.tolist(),'R_mean':((R_stoch/runs)/((1-R_stoch)+R_stoch/runs)).tolist(),'n_undef_oneway':n_undef,'Ns':Ns}
if __name__=='__main__':
    import sys
    fams=sys.argv[1:] or ['lognormal','hardstep','triggered']
    t=time.time(); res={'outputs':OUTS,'inputs':[i[0] for i in INPUTS],'groups':GROUPS}
    for fam in fams:
        res[fam]=analyse(fam,Ns=128,runs=30,boot=300); print(fam,'done',round(time.time()-t))
        for j,q in enumerate(OUTS):
            order=sorted(res[fam]['ST'],key=lambda k:-res[fam]['ST'][k][j])
            print(f"  {q:12s} R_stoch={res[fam]['R_stoch'][j]:.2f} R_mean={res[fam]['R_mean'][j]:.3f} | "+' '.join(f"{k}:{res[fam]['ST'][k][j]:.2f}[{res[fam]['ST_ci'][k][0][j]:.2f},{res[fam]['ST_ci'][k][1][j]:.2f}]" for k in order[:7]))
    json.dump(res,open('sobol3_'+'_'.join(fams)+'.json','w'),indent=1); print('total',round(time.time()-t))
