import numpy as np, itertools, json, time
import marked_process as mp
np.set_printoptions(precision=3)
out={}
t=time.time()
# 1. kernel robustness (baseline: lognormal family, short carbon, machines persist)
out['kernel']={}
for k in ['exp','inv2','step']:
    s,_=mp.run('lognormal',runs=250,kernel=k)
    out['kernel'][k]={q:float(np.nanmean(s[q])) for q in ['mutual_now','oneway_frac','scc_max','isolated','living_frac']}
    print('kernel',k,out['kernel'][k])
# 2. two-level factorial, main-effect variance decomposition
factors={'N_C':[5,15],'lam_sig':[0.25,1.0],'beta0':[0.1,1.0],'r0':[50.0,500.0],'fam':['lognormal','hardstep'],'L_C_med':[0.05,1.0],'L_info':[0.05,2.0],'p_M':[0.0,0.5]}
names=list(factors); combos=list(itertools.product(*[factors[n] for n in names]))
Y={'mutual_now':[],'oneway_frac':[],'scc_max':[]}; X=[]
for c in combos:
    kw=dict(zip(names,c)); fam=kw.pop('fam')
    s,_=mp.run(fam,runs=40,**kw)
    X.append(c)
    for q in Y: Y[q].append(float(np.nanmean(s[q])) if np.isfinite(np.nanmean(s[q])) else 0.0)
print('factorial done',time.time()-t)
out['sobol']={}
for q in Y:
    y=np.array(Y[q]); V=y.var(); S={}
    for i,n in enumerate(names):
        lv=[c[i] for c in X]; m=[y[[j for j in range(len(X)) if lv[j]==v]].mean() for v in factors[n]]
        S[n]=float(np.var(m)/V) if V>0 else 0.0
    out['sobol'][q]=dict(sorted(S.items(),key=lambda kv:-kv[1]))
    print(q, out['sobol'][q])
# 3. phase diagram: persistence scale P (L_C=P/10, L_M=P, L_info=P) vs bandwidth beta0; outputs
Ps=np.logspace(-2,1,7); Bs=np.logspace(-2,1,7)
grid={q:np.zeros((len(Ps),len(Bs))) for q in ['mutual_now','oneway_frac','scc_max','isolated']}
for i,P in enumerate(Ps):
    for j,B in enumerate(Bs):
        s,_=mp.run('lognormal',runs=80,beta0=B,L_C_med=P/10,L_C_sig=0.5,L_M=P,L_info=P,p_M=0.5)
        for q in grid: grid[q][i,j]=np.nanmean(s[q]) if np.isfinite(np.nanmean(s[q])) else 0.0
    print('phase row',i,time.time()-t)
out['phase']={'P':Ps.tolist(),'beta0':Bs.tolist(),**{q:grid[q].tolist() for q in grid}}
json.dump(out,open('experiments.json','w'),indent=1)
print('all done',time.time()-t)
