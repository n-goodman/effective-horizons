import numpy as np, json, time, marked_process as mp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
fam='lognormal'; sig=0.5; box=500.0; r0=150.0
d=json.load(open('phase3.json')); a_med=d['a_med']; lam_std=d['lam_std']; kfac=d['kfac']
xs=np.array(d['x']); ys=np.array(d['y'])
out={'x':xs.tolist(),'y':ys.tolist(),'a_med':a_med,'lam_std':lam_std,'kfac':kfac}
t=time.time()
for vary in ['info','synth']:
    res={q:np.zeros((len(ys),len(xs))) for q in ['isolated','oneway_frac','gscc','p_contact','mutual_now','scc_max']}
    for i,y in enumerate(ys):
        for j,x in enumerate(xs):
            mp.rng=np.random.default_rng(200+i*20+j); beta0=x*lam_std/kfac
            kw=dict(L_info=y*a_med,L_M=1.0) if vary=='info' else dict(L_info=0.5,L_M=y*a_med)
            s,_=mp.run(fam,runs=100,beta0=beta0,lam_sig=sig,L_C_med=0.05,L_C_sig=1.0,p_M=0.5,box=box,r0=r0,**kw)
            for q in res: res[q][i,j]=np.nanmean(s[q])
    out[vary]={q:v.tolist() for q,v in res.items()}; print(vary,round(time.time()-t),flush=True)
json.dump(out,open('phase4.json','w'),indent=1)
navy='#0D1F2D'; gold='#B8962E'; slate='#4A5568'; cream='#F3F1EA'; cmap=LinearSegmentedColormap.from_list('house',[cream,gold,navy])
X,Y=np.meshgrid(np.where(xs>0,xs,xs[1]/3),np.where(ys>0,ys,ys[1]/3))
for vary,ylab,fname in [('info',r'information persistence / median head start  $L^{\rm info}/\tilde a$','phase4'),('synth',r'synthetic lifetime / median head start  $L_M/\tilde a$','phase4b')]:
    fig,axes=plt.subplots(1,3,figsize=(13,3.9))
    for ax,(q,title) in zip(axes,[('isolated','isolated fraction $I_F$ of existing lineages'),('oneway_frac',r'one-way share of contacted pairs $f^{+}_{\rm one}$'),('gscc',r'cluster order parameter $G_{\rm SCC}$')]):
        Z=np.array(out[vary][q]); cf=ax.contourf(X,Y,Z,levels=np.linspace(0,1,11),cmap=cmap); cs=ax.contour(X,Y,Z,levels=[0.1,0.25,0.5,0.75],colors='white',linewidths=0.7); ax.clabel(cs,fontsize=7,fmt='%.2f')
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_title(title,color=navy,fontsize=10); ax.set_xlabel(r'typical-pair bandwidth / rate dispersion  $\tilde\beta/\sigma_\lambda$',color=navy,fontsize=9); ax.tick_params(colors=navy)
        cb=plt.colorbar(cf,ax=ax,fraction=0.046); cb.ax.tick_params(colors=navy,labelsize=8)
        for sp in ax.spines.values(): sp.set_color(slate)
    axes[0].set_ylabel(ylab,color=navy,fontsize=9)
    plt.tight_layout(); plt.savefig(f'../{fname}_fig.pdf',bbox_inches='tight'); plt.savefig(f'fig_{fname}.png',dpi=160,bbox_inches='tight')
print('done')
