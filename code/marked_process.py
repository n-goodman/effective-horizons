"""v2: marked spacetime birth-and-contact model.
Opportunity process (SFR x metallicity) -> emergence hazard families -> reset (sterilisation) hazard
-> carbon lineage births; machine lineage as a transition from carbon; lifetimes and information
persistence; spatial positions with light travel and distance-dependent bandwidth; directed meeting
graph with first-edge, first-mutual-edge and component statistics. All inputs are scenario
parametrisations, not fitted curves."""
import numpy as np
from scipy.optimize import brentq
rng=np.random.default_rng(7)

# cosmology / opportunity process ---------------------------------------------------
H0=67.7; Om=0.31; OL=1-Om; tH=977.8/H0
def t_of_z(z): return (2/(3*np.sqrt(OL)))*tH*np.arcsinh(np.sqrt(OL/Om)*(1+z)**-1.5)
t0=t_of_z(0.0)
def z_of_t(t): return brentq(lambda z: t_of_z(z)-t, 0, 60)
def sfr_z(z): return 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)
TMAX=25.0
tg=np.linspace(0.2,TMAX,1000); dt=tg[1]-tg[0]
sfr=np.array([sfr_z(z_of_t(t)) if t<=t0 else sfr_z(0.0)*np.exp(-(t-t0)/3.5) for t in tg])
csf=np.cumsum(sfr)*dt; Zrel=np.minimum(1.0, csf/np.interp(8.0,tg,csf))
R_planet=sfr*Zrel                                    # opportunity density in time (units arbitrary)

# reset (sterilisation) hazard, per Gyr, decreasing with cosmic time -------------------
def gamma_reset(t, g0=0.6, tref=5.0, tau=3.0): return g0*np.exp(-(t-tref)/tau)
# emergence hazard families: h(u, t) = hazard at planet age u, cosmic time t ----------
def family(name):
    if name=='lognormal':      # fixed delay distribution, median 4.5 Gyr
        med,sig=4.5,0.25
        pdf=lambda u: np.where(u>0,np.exp(-(np.log(np.maximum(u,1e-9))-np.log(med))**2/(2*sig**2))/(np.maximum(u,1e-9)*sig*np.sqrt(2*np.pi)),0)
        return lambda u,t: pdf(u)                     # treated directly as completion density
    if name=='hardstep':       # Carter n hard steps, completion density ~ u^(n-1) on [0,Th], total prob p_e
        n,Th=5,10.0
        return lambda u,t: np.where((u>0)&(u<Th), n*u**(n-1)/Th**n, 0.0)
    if name=='triggered':      # hazard switched on by a cosmic-time trigger (e.g. metallicity/oxygen analogue)
        h0=0.3; S=lambda t: 1/(1+np.exp(-(t-9.0)/1.0))
        # completion density at age u for a planet formed at t-u: h(t) * exp(-int h)
        def dens(u,t):
            tp=t-u; out=np.zeros_like(u)
            ss=np.linspace(0,1,40)
            for k in np.where(u>0)[0]:
                s=tp[k]+ss*u[k]; H=h0*np.trapezoid(S(s),s)
                out[k]=h0*S(t)*np.exp(-H)
            return out
        return dens
def birth_density(fam):
    """carbon birth density Lambda_C(t) = int R_planet(tp) * surv_reset(tp->t) * dens(t-tp, t) dtp"""
    dens=family(fam); out=np.zeros_like(tg)
    # survival against resets from tp to t
    G=np.cumsum(gamma_reset(tg))*dt
    for j,t in enumerate(tg):
        tp=tg[:j+1]; u=t-tp
        surv=np.exp(-(G[j]-G[:j+1]))
        d=dens(u,t)
        out[j]=np.sum(R_planet[:j+1]*surv*d)*dt
    return out
def sample_times(b,n):
    cdf=np.cumsum(b); cdf/=cdf[-1]; return np.interp(rng.random(n),cdf,tg)

# marks ---------------------------------------------------------------------------------
C_MPC_PER_GYR=306.6
KERNELS={'exp':lambda r,r0: np.exp(-r/r0),'inv2':lambda r,r0: (1+r/r0)**-2,'step':lambda r,r0: float(r<r0)}
_BCACHE={}
def run(fam, N_C=10, box=500.0, beta0=0.4, r0=150.0, lam_med_C=1.0, lam_mult_M=3.0, lam_sig=0.5,
        L_C_med=0.05, L_C_sig=1.0, p_M=0.5, L_M=np.inf, L_info=0.5, runs=400, t_eval=None, kernel='exp', beta_passive=None, kappa0=None, active=True):
    t_eval=t0 if t_eval is None else t_eval
    if fam not in _BCACHE: _BCACHE[fam]=birth_density(fam)
    b=_BCACHE[fam]; K=KERNELS[kernel]
    stats=dict(first_edge=[],first_mutual=[],oneway_frac=[],scc_max=[],isolated=[],n_alive=[],living_frac=[],mutual_now=[],born_frac=[],p_contact=[],gscc=[])
    bp=beta0 if beta_passive is None else beta_passive
    for _ in range(runs):
        tau=sample_times(b,N_C); pos=rng.random((N_C,3))*box
        lam=np.exp(np.log(lam_med_C)+lam_sig*rng.standard_normal(N_C))
        L=np.exp(np.log(L_C_med)+L_C_sig*rng.standard_normal(N_C))      # carbon lifetimes (Gyr)
        kind=['C']*N_C
        # machine descendants: transition with prob p_M at a fraction of the carbon lifetime
        for i in range(N_C):
            if rng.random()<p_M:
                tau=np.append(tau,tau[i]+rng.random()*L[i]); pos=np.vstack([pos,pos[i]+rng.normal(0,1,3)])
                lam=np.append(lam,lam[i]*lam_mult_M*np.exp(0.3*rng.standard_normal())); L=np.append(L,L_M); kind.append('M')
        n=len(tau); death=tau+L; avail=death+L_info
        T=np.full((n,n),np.inf)                      # T[i,j]: time j parses i
        for i in range(n):
            for j in range(n):
                if i==j: continue
                r=np.linalg.norm(pos[i]-pos[j]); t_causal=max(tau[j], tau[i]+r/C_MPC_PER_GYR)
                if t_causal>=avail[i]: continue
                beta=beta0*K(r,r0); beta_p=bp*K(r,r0)
                # access gates: active channel (testimony while i alive) and passive observability of i's record
                passive_open = True if kappa0 is None else (rng.random()<min(1.0,kappa0*K(r,r0)))
                if (not passive_open) and t_causal>=death[i]: continue
                # lead of i over j at t_causal: i's content minus j's content (rungs), both from rung 0
                lead=lam[i]*(min(t_causal,death[i])-tau[i]) - lam[j]*(t_causal-tau[j])
                if death[i]>t_causal:                       # i alive at causal availability
                    rate_alive=lam[i]-lam[j]-(beta if active else 0.0)
                    if active:
                        if lead<=0: T[i,j]=t_causal; continue
                        if rate_alive<0:
                            tc=t_causal+lead/(-rate_alive)
                            if tc<=death[i]: T[i,j]=tc; continue
                    lead=lead+rate_alive*(death[i]-t_causal); tstart=death[i]
                else: tstart=t_causal
                if not passive_open: continue
                if lead<=0: T[i,j]=tstart; continue
                rate_dead=-lam[j]-beta_p
                tc=tstart+lead/(-rate_dead)
                if tc<=avail[i]: T[i,j]=tc
        finite=np.isfinite(T)
        stats['first_edge'].append(T.min() if finite.any() else np.inf)
        M=np.maximum(T,T.T); mut=M[np.triu_indices(n,1)]
        stats['first_mutual'].append(mut.min() if np.isfinite(mut).any() else np.inf)
        born=tau<=t_eval; nb=int(born.sum()); stats['born_frac'].append(nb/n)
        E=(T<=t_eval); E=E[np.ix_(born,born)]
        if nb<2:
            for q in ['oneway_frac','living_frac','scc_max','isolated','mutual_now','n_alive','p_contact','gscc']: stats[q].append(np.nan if q in ('oneway_frac','living_frac') else 0.0)
            continue
        iu=np.triu_indices(nb,1); any_pair=(E|E.T)[iu]; both=(E&E.T)[iu]
        stats['oneway_frac'].append((any_pair.sum()-both.sum())/any_pair.sum() if any_pair.sum() else np.nan)
        stats['mutual_now'].append(float(both.any())); stats['p_contact'].append(float(any_pair.any()))
        alive=(death[born]>t_eval)
        stats['living_frac'].append((E[alive,:].sum()/E.sum()) if E.sum() else np.nan)
        R=E.copy()|np.eye(nb,dtype=bool)
        for k in range(nb): R=R|(R[:,[k]]&R[[k],:])
        scc=R&R.T; seen=np.zeros(nb,bool); sizes=[]
        for i in range(nb):
            if not seen[i]: comp=scc[i]; sizes.append(comp.sum()); seen|=comp
        stats['scc_max'].append(max(sizes)/nb)                 # fraction of BORN lineages
        stats['gscc'].append((max(sizes)-1)/(nb-1))            # order parameter: 0 = no nontrivial SCC, 1 = all existing lineages
        stats['isolated'].append(np.mean(~(E.any(0)|E.any(1))))  # among BORN lineages
        stats['n_alive'].append(np.mean(alive))
    return {k:np.array(v) for k,v in stats.items()}, b

if __name__=='__main__':
    for fam in ['lognormal','hardstep','triggered']:
        b=birth_density(fam); m=tg<=t0
        peak=tg[np.argmax(b)]; frac=np.trapezoid(b[m],tg[m])/np.trapezoid(b,tg)
        w=b[m]/b[m].sum(); a=t0-tg[m]; med=a[np.searchsorted(np.cumsum(w),0.5)]
        print(f"[{fam:9s}] carbon birth peak {peak:.1f} Gyr; fraction of births within model horizon ({TMAX} Gyr) before now {frac:.2f}; median head start of existing peers {med:.2f} Gyr")
    print()
    for fam in ['lognormal','hardstep','triggered']:
        for lab,kw in [('persistent carbon (L_C=1 Gyr)',dict(L_C_med=1.0,L_C_sig=0.5)),('short carbon (L_C=50 Myr), machines persist',dict(L_C_med=0.05,L_C_sig=1.0)),('short carbon, no machines',dict(L_C_med=0.05,L_C_sig=1.0,p_M=0.0))]:
            s,_=run(fam,runs=300,**kw)
            fe=s['first_edge']; fm=s['first_mutual']
            print(f"[{fam:9s}] {lab:44s} first directed edge median {np.nanmedian(np.where(np.isfinite(fe),fe,np.nan)):5.2f} Gyr (P<now {np.mean(fe<=t0):.2f}); first mutual median {np.nanmedian(np.where(np.isfinite(fm),fm,np.nan)):5.2f} (P<now {np.mean(fm<=t0):.2f}); one-way edge fraction now {np.nanmean(s['oneway_frac']):.2f}; largest SCC {np.mean(s['scc_max']):.2f}; isolated {np.mean(s['isolated']):.2f}; alive now {np.mean(s['n_alive']):.2f}")
