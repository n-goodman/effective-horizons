import numpy as np, json, time, marked_process as mp
a_med=2.426; ys=np.logspace(np.log10(0.15),np.log10(1.5),11)
out={'y':ys.tolist(),'N':{}}; t=time.time()
for N,runs in [(5,400),(10,400),(20,300),(40,150),(80,60)]:
    box=500.0*(N/10)**(1/3); res={q:[] for q in ['gscc','isolated']}
    for y in ys:
        mp.rng=np.random.default_rng(int(1000*y)+N); s,_=mp.run('lognormal',N_C=N,box=box,runs=runs,L_info=y*a_med,L_M=1.0,L_C_med=0.05,L_C_sig=1.0,p_M=0.5)
        for q in res: res[q].append(float(np.nanmean(s[q])))
    out['N'][str(N)]={'box':box,'runs':runs,**res}; print(N,round(time.time()-t),np.round(res['gscc'],2).tolist(),flush=True)
json.dump(out,open('nscale2.json','w'),indent=1)
