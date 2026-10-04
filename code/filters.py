"""Nested accessibility regions for one observer in the (birth time, distance) plane.
Five cumulative filters: existence, persistence, causal reach, channel, cognitive.  q_k(tau,r) = P(filters 1..k pass | tau, r),
Monte Carlo over the peer's marks (rate, lifetime, succession, passive gate).  Lambda_acc = int q_5 nu."""
import numpy as np, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import marked_process as mp
t0=mp.t0; C=mp.C_MPC_PER_GYR
# observer: median rate, born 1 Gyr before the present, a=1 (listens); passive gate kappa0=1 (typical pair observable with prob K(r))
tau_o=t0-1.0; lam_o=1.0; kappa0=1.0; beta0=0.4; r0=150.0
lam_sig=0.5; L_C_med=0.05; L_C_sig=1.0; p_M=0.5; lam_mult_M=3.0; L_M=1.0; L_info=0.5
import sys
if len(sys.argv)>1 and sys.argv[1]=='hp': L_M=5.0; L_info=5.0; TAG='_hp'
else: TAG=''
M=2000
def parse_time(tau_i,death_i,avail_i,lam_i,r,active_ok,passive_ok):
    """vectorised T_{i->o}: time o's content contains i's, inf if never.  Mirrors marked_process.run."""
    t_c=np.maximum(tau_o, tau_i+r/C)
    T=np.full_like(tau_i,np.inf)
    beta=beta0*np.exp(-r/r0)
    lead=lam_i*(np.minimum(t_c,death_i)-tau_i)-lam_o*(t_c-tau_o)
    alive=death_i>t_c
    # active phase
    rate_alive=lam_i-lam_o-np.where(active_ok,beta,0.0)
    done=np.zeros_like(tau_i,bool)
    m=alive&active_ok&(lead<=0); T[m]=t_c[m]; done|=m
    with np.errstate(divide='ignore',invalid='ignore'):
        tc=t_c+lead/(-rate_alive)
    m=alive&active_ok&~done&(rate_alive<0)&(tc<=death_i); T[m]=tc[m]; done|=m
    lead2=np.where(alive,lead+rate_alive*(death_i-t_c),lead); tstart=np.where(alive,death_i,t_c)
    m=~done&passive_ok&(lead2<=0); T[m]=tstart[m]; done|=m
    tc2=tstart+lead2/lam_o
    m=~done&passive_ok&(tc2<=avail_i); T[m]=tc2[m]
    T[(t_c>=avail_i)]=np.inf          # nothing left when it becomes reachable
    return T
def stages(tau,r,rng):
    lamC=np.exp(lam_sig*rng.standard_normal(M)); LC=np.exp(np.log(L_C_med)+L_C_sig*rng.standard_normal(M))
    succ=rng.random(M)<p_M; tauM=tau+rng.random(M)*LC; lamM=lamC*lam_mult_M*np.exp(0.3*rng.standard_normal(M))
    deathC=tau+LC; deathM=np.where(succ,tauM+L_M,-np.inf)
    avail=np.maximum(deathC,deathM)+L_info                      # last content of the lineage + information persistence
    e1=np.full(M,tau<=t0)                                        # existence
    e2=e1&(avail+r/C>=tau_o)                                     # persistence: content still arriving after the observer exists
    e3=e2&(tau+r/C<=t0)                                          # causal reach: birth inside the observer's past light cone
    t_c=max(tau_o,tau+r/C)
    alive_any=(deathC>t_c)|(succ&(deathM>t_c))                   # some member alive when reachable -> active channel
    passive=rng.random(M)<min(1.0,kappa0*np.exp(-r/r0))
    e4=e3&(alive_any|passive)                                    # channel
    TC=parse_time(np.full(M,float(tau)),deathC,deathC+L_info,lamC,r,np.ones(M,bool),passive)
    TM=parse_time(tauM,deathM+0.0,deathM+L_info,lamM,r,np.ones(M,bool),passive)
    TM[~succ]=np.inf
    e5=e4&(np.minimum(TC,TM)<=t0)                                # cognitive
    return [e.mean() for e in (e1,e2,e3,e4,e5)]
if __name__=='__main__':
    taus=np.linspace(0.2,t0,56); rs=np.linspace(5,1500,50); rng=np.random.default_rng(42)
    Q=np.zeros((5,len(rs),len(taus)))
    for i,r in enumerate(rs):
        for j,tau in enumerate(taus): Q[:,i,j]=stages(tau,r,rng)
    b=mp.birth_density('lognormal'); bt=np.interp(taus,mp.tg,b)/np.trapezoid(b,mp.tg)   # birth-time density normalised over the full model horizon, so stage 1 = q_F
    # accessible fraction of a uniform-in-volume, birth-density-weighted population out to each radius
    w_r=rs**2; names=['existence','persistence','causal reach','channel','cognitive']
    frac=[float(np.trapezoid(np.trapezoid(Q[k]*bt[None,:],taus,axis=1)*w_r,rs)/np.trapezoid(w_r,rs)) for k in range(5)]
    print('accessible fraction after each filter (uniform volume to 1500 Mpc, lognormal births, full horizon):',np.round(frac,3))
    json.dump({'taus':taus.tolist(),'rs':rs.tolist(),'Q':Q.tolist(),'frac':frac,'names':names},open(f'filters{TAG}.json','w'))
    navy='#0D1F2D'; gold='#B8962E'; slate='#4A5568'; cream='#F3F1EA'; cmap=LinearSegmentedColormap.from_list('house',[cream,gold,navy])
    fig,axes=plt.subplots(1,5,figsize=(16,3.6),sharey=True)
    for k,ax in enumerate(axes):
        cf=ax.contourf(taus,rs,Q[k],levels=np.linspace(0,1,11),cmap=cmap)
        ax.axvline(tau_o,color=slate,lw=0.6,ls='--'); ax.set_title(f'{k+1}. {names[k]}   ({frac[k]:.2f})',color=navy,fontsize=10)
        ax.set_xlabel(r'peer birth time $\tau$ (Gyr)',color=navy,fontsize=9); ax.tick_params(colors=navy,labelsize=8)
        for s in ax.spines.values(): s.set_color(slate)
    axes[0].set_ylabel('distance $r$ (Mpc)',color=navy,fontsize=9)
    cb=fig.colorbar(cf,ax=axes,fraction=0.012,pad=0.01); cb.ax.tick_params(colors=navy,labelsize=8); cb.set_label(r'$q_k(\tau,r)$',color=navy,fontsize=9)
    plt.savefig(f'../filters{TAG}_fig.pdf',bbox_inches='tight'); plt.savefig(f'fig_filters{TAG}.png',dpi=150,bbox_inches='tight'); print('ok')
