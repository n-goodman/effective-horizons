"""Start-time model for cognitive growth: cosmic star formation -> planets -> intelligence,
substrate-dependent stability windows, promotion rates, and meeting times. First pass, all
astrophysical inputs are crude parametrisations flagged in the notes."""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm
rng=np.random.default_rng(1)

# ---- cosmology (flat LCDM, Planck-ish) ----
H0=67.7; Om=0.31; OL=1-Om; tH=977.8/H0   # Gyr
def t_of_z(z): return (2/(3*np.sqrt(OL)))*tH*np.arcsinh(np.sqrt(OL/Om)*(1+z)**-1.5)
t0=t_of_z(0.0)
def z_of_t(t): return brentq(lambda z: t_of_z(z)-t, 0, 60)
# ---- Madau-Dickinson cosmic SFR density ----
def sfr_z(z): return 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)
tgrid=np.linspace(0.2,40.0,1200)
sfr_t=np.array([sfr_z(z_of_t(t)) if t<=t0 else sfr_z(0.0)*np.exp(-(t-t0)/3.5) for t in tgrid])  # future SFR: exponential decline, tau=3.5 Gyr (assumption)
csf=np.cumsum(sfr_t)*(tgrid[1]-tgrid[0]); csf/=csf[-1]
# metallicity proxy: cumulative star formation, normalised to 1 at t=8 Gyr (solar-ish)
Zrel=np.minimum(1.0, csf/np.interp(8.0,tgrid,csf))
planet_rate=sfr_t*Zrel                           # terrestrial-planet formation history (crude Lineweaver-like)
# ---- substrate stability measures S_s(t) (fraction of universe stable enough) ----
S_c=lambda t: 1/(1+np.exp(-(t-7.0)/1.5))         # carbon: GRB/metallicity-limited, mid 7 Gyr
S_m=lambda t: np.ones_like(t)                    # machine: not star-bound
S_x=lambda t: 1/(1+np.exp(-(t-12.0)/1.0))        # exotic late substrate: cold universe
# ---- delay from planet formation to agentic intelligence ----
def birth_density(S, med_delay=4.5, sig=0.25):
    b=np.zeros_like(tgrid)
    dt=tgrid[1]-tgrid[0]
    for i,tf in enumerate(tgrid):
        tau=tgrid; d=tau-tf
        pdf=np.where(d>0, np.exp(-(np.log(np.maximum(d,1e-9))-np.log(med_delay))**2/(2*sig**2))/(np.maximum(d,1e-9)*sig*np.sqrt(2*np.pi)),0)
        b+=planet_rate[i]*pdf*dt
    return b*S(tgrid)
b_c=birth_density(S_c)
b_m=b_c.copy()                                   # post-biological: same origin times (delta << Gyr)
b_x=birth_density(S_x, med_delay=2.0, sig=0.3)   # exotic: faster to intelligence once stable
def sample_tau(b, n):
    cdf=np.cumsum(b); cdf/=cdf[-1]
    return np.interp(rng.random(n), cdf, tgrid)
def rates(n, med, sig=0.5): return np.exp(np.log(med)+sig*rng.standard_normal(n))

# ---- 1. head starts of peers born before now ----
def headstart_stats(b):
    mask=tgrid<=t0
    w=b[mask]/b[mask].sum(); a=t0-tgrid[mask]
    mean=(w*a).sum(); 
    cdf=np.cumsum(w); med=a[np.searchsorted(cdf,0.5)]
    return mean, med, (w*(a<1)).sum(), (w*(a>5)).sum()
for name,b in [('carbon',b_c),('machine',b_m),('exotic',b_x)]:
    print(f"{name:8s} head start: mean {headstart_stats(b)[0]:.2f} Gyr, median {headstart_stats(b)[1]:.2f}, P(a<1Gyr)={headstart_stats(b)[2]:.3f}, P(a>5Gyr)={headstart_stats(b)[3]:.3f}, fraction of all-time births before now {np.trapezoid(b[tgrid<=t0],tgrid[tgrid<=t0])/np.trapezoid(b,tgrid):.2f}")

# ---- 2. our silence probability vs budget B, N peers of a substrate ----
beta=0.4; lam0=1.0   # our rate = carbon median
def silence(b, lam_med, N, B, runs=4000):
    cnt=0
    for _ in range(runs):
        tau=sample_tau(b*(tgrid<=t0), N); a=t0-tau; lam=rates(N, lam_med)
        tc=np.where(lam<lam0+beta, lam*a/np.maximum(lam0+beta-lam,1e-12), np.inf)
        cnt+= (tc>B).all()
    return cnt/runs
print("\nsilence probability (our rate = carbon median, beta=0.4):")
for B in [0.001,0.01,0.1,1.0]:
    row=[]
    for name,b,lm in [('carbon',b_c,1.0),('machine',b_m,3.0),('exotic',b_x,3.0)]:
        row.append(f"{name} N=10: {silence(b,lm,10,B):.2f}  N=30: {silence(b,lm,30,B):.2f}")
    print(f"B={B:>6} Gyr | "+" | ".join(row))

# ---- 3. cosmic meeting epochs among a population ----
def meetings(Ns, runs=2000, beta=0.4):
    """Ns: dict substrate->(density, lam_med, N). Returns first-meeting epochs and counts by now."""
    firsts=[]; counts=[]; pairtypes={}
    for _ in range(runs):
        taus=[];lams=[];kinds=[]
        for k,(b,lm,N) in Ns.items():
            taus+=list(sample_tau(b,N)); lams+=list(rates(N,lm)); kinds+=[k]*N
        taus=np.array(taus); lams=np.array(lams); kinds=np.array(kinds)
        order=np.argsort(taus)
        T=[]; types=[]
        n=len(taus)
        for ii in range(n):
            for jj in range(ii+1,n):
                i,j=order[ii],order[jj]   # i earlier
                a=taus[j]-taus[i]
                if lams[j]+beta>lams[i]:
                    tc=lams[i]*a/(lams[j]+beta-lams[i]); T.append(taus[j]+tc); types.append(kinds[i]+'-'+kinds[j])
        T=np.array(T)
        if len(T): 
            firsts.append(T.min()); counts.append((T<=t0).sum())
            k=types[int(T.argmin())]; pairtypes[k]=pairtypes.get(k,0)+1
        else: firsts.append(np.inf); counts.append(0)
    return np.array(firsts), np.array(counts), pairtypes
for N in [3,10,30]:
    f,c,pt=meetings({'carbon':(b_c,1.0,N),'machine':(b_m,3.0,N//3 if N>=3 else 1),'exotic':(b_x,3.0,max(1,N//10))})
    print(f"\nN_carbon={N}: first meeting epoch median {np.median(f):.2f} Gyr (10-90%: {np.percentile(f,10):.2f}-{np.percentile(f,90):.2f}); P(first meeting before now)={np.mean(f<=t0):.2f}; mean meetings by now {c.mean():.1f}; first pair types {pt}")
np.save('grid.npy', np.vstack([tgrid,sfr_t,planet_rate,b_c,b_x]))
