import numpy as np, json, marked_process as mp
Q=['p_contact','mutual_now','oneway_frac','gscc','isolated']
out={}
for reg,base in [('baseline',{}),('highpersist',dict(L_info=5.0,L_M=5.0))]:
    for lab,kw in [('A+P',{}),('P only',dict(active=False,beta0=0.0)),('A only',dict(kappa0=0.0)),('none',dict(active=False,beta0=0.0,kappa0=0.0)),('A + P(kappa0=1)',dict(kappa0=1.0)),('A + P(kappa0=0.3)',dict(kappa0=0.3))]:
        mp.rng=np.random.default_rng(7); s,_=mp.run('lognormal',runs=400,**base,**kw)
        out[reg+'/'+lab]={q:float(np.nanmean(s[q])) if np.isfinite(np.nanmean(s[q])) else None for q in Q}
        print(reg,lab,{q:(round(v,3) if v is not None else None) for q,v in out[reg+'/'+lab].items()})
json.dump(out,open('channels.json','w'),indent=1)
