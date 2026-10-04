import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sobol; sobol.OUTS=['p_contact','mutual_now','oneway_frac','gscc','isolated']; OUTS=sobol.OUTS; elasticities=sobol.elasticities
y0,E=elasticities()
print('baseline',dict(zip(OUTS,np.round(y0,3))))
for n,e in E.items(): print(f"{n:11s} "+'  '.join(f"{q}={v:+.2f}" for q,v in zip(OUTS,e)))
json.dump({'baseline':y0,'elasticities':E,'outputs':OUTS},open('elast06.json','w'),indent=1)
navy='#0D1F2D'; gold='#B8962E'; slate='#4A5568'
labels={'N_C':'$N$','lam_sig':r'$\sigma_{\ln\lambda}$','beta0':r'$\beta_0$','r0':'$r_0$','L_C_med':'$L_C$','L_info':r'$L^{\rm info}$','p_M':'$p_M$','lam_mult_M':r'$\lambda_M/\lambda_C$'}
onames=[r'$P_{\rm contact}$','P(mutual by now)',r'$f^{+}_{\rm one}$',r'$G_{\rm SCC}$','isolated']
names=list(E); M=np.array([E[n] for n in names])
fig,ax=plt.subplots(figsize=(7,3.2)); w=0.16; x=np.arange(len(names)); cols=['#C9B37A',navy,gold,slate,'#9AA5B1']
for j in range(5): ax.bar(x+(j-2)*w,M[:,j],w,color=cols[j],label=onames[j])
ax.axhline(0,color=slate,lw=0.6); ax.set_xticks(x); ax.set_xticklabels([labels[n] for n in names],color=navy); ax.set_ylabel(r'$\partial\ln Y/\partial\ln x$ at baseline',color=navy,fontsize=9)
ax.tick_params(colors=navy,labelsize=8); ax.legend(fontsize=7,frameon=False,ncol=2)
for s in ax.spines.values(): s.set_color(slate)
plt.tight_layout(); plt.savefig('../elast_fig.pdf',bbox_inches='tight'); plt.savefig('fig_elast.png',dpi=150,bbox_inches='tight')
