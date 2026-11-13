# analysis_scripts.py - drivers used for Sections 10-12 (Python 3 + numpy; the engine itself stays dependency-free)
# ---- exp.py: thin wrapper around the compiled engine (rpt_adapt) and a bisection for the persistence threshold ----
import subprocess, re, json, numpy as np, os
EXE=os.environ.get('RPT_ADAPT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rpt_adapt'))   # build: g++ -std=c++20 -O2 rpt_adapt.cpp -o rpt_adapt
BEST=dict(beta_T=8.129e-7,alpha=0.05018,p_T=2877,tip_dose=3.794e7,tip_time=0.2602)   # screened-best sampled profile
REFINED=dict(beta_T=2.0e-6,alpha=0.05,p_T=1779,tip_dose=2.24e6,tip_time=0.02865)    # Nelder-Mead refined profile
def run(over, traj=None, adaptive=0):
    args=[f"{k}={v}" for k,v in over.items()]+[f"adaptive={adaptive}"] if 'adaptive' not in over else [f"{k}={v}" for k,v in over.items()]
    if traj: args.append(f"traj={traj}")
    o=subprocess.run([EXE]+args,capture_output=True,text=True).stdout
    d={k:float(v) for k,v in re.findall(r'(\w+)=([\d.eE+-]+|nan|inf)',o)}
    return d
def ext(over,**kw): return run(over,**kw)['extinct']==1
def cT_star(prof, extra=None, lo=0.01, hi=3.0, iters=22, adaptive=0):
    extra=extra or {}
    f=lambda c: ext({**prof,**extra,'c_T':c},adaptive=adaptive)
    if not f(lo): return np.nan
    if f(hi): return hi
    for _ in range(iters):
        mid=(lo*hi)**0.5
        if f(mid): lo=mid
        else: hi=mid
    return lo


# ---- tail-stochastic extinction check (Section 12.2): Gillespie simulation of the low-count WT subsystem in an ODE environment ----
import numpy as np
rng=np.random.default_rng(20261109)
VOL=5000.0   # host blood volume, mL
def tail_sim(env, pars, N_start, reps, t_end=30.0, esc=1500):
    t=env[:,0]; T=env[:,1]; EW=env[:,2]; IW=env[:,3]; EC=env[:,4]; IC=env[:,5]; P=env[:,6]; VW=env[:,7]; VT=env[:,8]; C=env[:,9]
    tot=VOL*(VW+EW+IW+EC+IC); ipk=int(np.argmax(t>=0.3)); tmin=float(tot[ipk:].min())
    N_start=max(N_start,1.5*tmin)
    k0=next(k for k in range(ipk,len(t)) if tot[k]<N_start)
    bW,bT,kW,kCo,dW,dCo,pW,al,cW,dN,KF=[pars[x] for x in ('bW','bT','kW','kCo','dW','dCo','pW','alpha','cW','dN','KF')]
    ext=0; esc_n=0
    for _ in range(reps):
        X=[rng.poisson(VOL*EW[k0]),rng.poisson(VOL*IW[k0]),rng.poisson(VOL*EC[k0]),rng.poisson(VOL*IC[k0]),rng.poisson(VOL*VW[k0])]  # E,I,Ec,Ic,V
        k=k0; tc=t[k0]; done=False
        while k<len(t)-1 and not done:
            dt=t[k+1]-t[k]; cyt=dN*C[k]/(KF+C[k]); tt=0.0
            while True:
                E,I,Ec,Ic,V=X
                r=[bW*T[k]*V, kW*E, bT*VT[k]*E, (dW+cyt)*I, bT*VT[k]*I, bW*P[k]*V, kCo*Ec, (dCo+cyt)*Ic, pW*I+al*pW*Ic, cW*V]
                R=sum(r)
                if R<=0: tt=dt; break
                tau=rng.exponential(1.0/R)
                if tt+tau>dt: break
                tt+=tau
                u=rng.random()*R; acc=0
                for j,rj in enumerate(r):
                    acc+=rj
                    if u<acc: break
                if j==0: X[0]+=1
                elif j==1: X[0]-=1; X[1]+=1
                elif j==2: X[0]-=1; X[2]+=1
                elif j==3: X[1]-=1
                elif j==4: X[1]-=1; X[3]+=1
                elif j==5: X[2]+=1
                elif j==6: X[2]-=1; X[3]+=1
                elif j==7: X[3]-=1
                elif j==8: X[4]+=1
                else: X[4]-=1
                if X[4]+X[0]+X[1]+X[2]+X[3] > esc: done=True; break
            if sum(X)==0: done=True
            k+=1
        if sum(X)==0: ext+=1
        elif sum(X)>esc: esc_n+=1
    return ext, esc_n, float(t[k0]), int(tot[k0])

# ---- branching-process establishment probability and its Monte Carlo validation (Section 12.1) ----
def q_theory(p,dl,R): return min(1.0,(p+dl)/(p+R*dl))
def mc(p,dl,c,R,n,thr=400):
    bTs=R*c*dl/p; est=0
    for _ in range(n):
        V,I=1,0
        while True:
            r1=c*V; r2=bTs*V; r3=dl*I; r4=p*I; R_=r1+r2+r3+r4
            if R_==0: break
            u=rng.random()*R_
            if u<r1: V-=1
            elif u<r1+r2: I+=1
            elif u<r1+r2+r3: I-=1
            else: V+=1
            if V>=thr: est+=1; break
    return est/n
