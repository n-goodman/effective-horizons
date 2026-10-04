import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
navy='#0D1F2D'; gold='#B8962E'; slate='#4A5568'; cream='#F3F1EA'
res={}
for f in ['lognormal','hardstep','triggered']: res.update(json.load(open(f'sobol3_{f}.json')))
inputs=res['inputs']; outs=res['outputs']; groups=['G:'+g for g in res['groups']]
labels={'N_C':'$N$','lam_sig':r'$\sigma_{\ln\lambda}$','beta0':r'$\beta_0$','r0':'$r_0$','L_C_med':'$L_C$','L_info':r'$L^{\rm info}$','p_M':'$p_M$','lam_mult_M':r'$\lambda_M/\lambda_C$',
        'G:abundance':'abund.','G:persistence':'persist.','G:rates':'rates','G:communication':'comm.'}
onames={'p_contact':r'$P_{\rm contact}$','mutual_now':'P(mutual contact by now)','oneway_frac':r'$f^{+}_{\rm one}$','gscc':r'$G_{\rm SCC}$','isolated':'isolated fraction'}
fams=['lognormal','hardstep','triggered']; fnames={'lognormal':'fixed delay','hardstep':'hard step','triggered':'triggered'}
fig,axes=plt.subplots(3,5,figsize=(15,7.2),sharey=True)
keys=inputs+groups; xpos=np.arange(len(keys)); xpos[len(inputs):]+=1
for i,f in enumerate(fams):
    r=res[f]
    for j,q in enumerate(outs):
        ax=axes[i,j]
        st=np.array([r['ST'][k][j] for k in keys]); lo=np.array([r['ST_ci'][k][0][j] for k in keys]); hi=np.array([r['ST_ci'][k][1][j] for k in keys])
        s1=np.array([r['S1'][k][j] for k in keys])
        cols=[gold]*len(inputs)+[navy]*len(groups)
        ax.bar(xpos,np.clip(st,0,None),color=cols,alpha=0.85,width=0.7)
        ax.errorbar(xpos,st,yerr=[np.clip(st-lo,0,None),np.clip(hi-st,0,None)],fmt='none',ecolor=slate,elinewidth=0.8,capsize=2)
        ax.plot(xpos,np.clip(s1,0,None),'o',color='white',mec=navy,ms=3.5)
        ax.set_xticks(xpos); ax.set_xticklabels([labels[k] for k in keys],rotation=70,fontsize=7,color=navy)
        ax.tick_params(colors=navy,labelsize=8); ax.set_ylim(0,1.3); ax.axhline(1,color=slate,lw=0.5,ls=':')
        if i==0: ax.set_title(onames[q],color=navy,fontsize=10)
        if j==0: ax.set_ylabel(fnames[f]+'\n$S_T$ (bars), $S_1$ (dots)',color=navy,fontsize=9)
        ax.text(0.98,0.95,f"$R_{{\\rm stoch}}={r['R_stoch'][j]:.2f}$, $R_{{\\rm mean}}={r['R_mean'][j]:.2f}$",transform=ax.transAxes,ha='right',va='top',fontsize=8,color=slate)
        for s in ax.spines.values(): s.set_color(slate)
plt.tight_layout(); plt.savefig('../sobol3_fig.pdf',bbox_inches='tight'); plt.savefig('fig_sobol3.png',dpi=150,bbox_inches='tight'); print('ok')
