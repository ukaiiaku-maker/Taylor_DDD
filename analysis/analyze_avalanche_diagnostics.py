#!/usr/bin/env python3
"""Postprocess DDD histories for avalanche-like depinning/plastic-burst diagnostics."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path


def sf(x, default=float('nan')):
    try:
        if x is None: return default
        y = float(x)
        return y if math.isfinite(y) else default
    except Exception:
        return default

def si(x, default=0):
    try: return int(float(x))
    except Exception: return default

def finite(vs): return [v for v in vs if math.isfinite(v)]

def median(vs):
    vs = sorted(finite(vs))
    if not vs: return float('nan')
    n=len(vs); m=n//2
    return vs[m] if n%2 else 0.5*(vs[m-1]+vs[m])

def gini(vs):
    vs = sorted([max(0.0, v) for v in vs if math.isfinite(v)])
    if not vs: return float('nan')
    s=sum(vs)
    if s <= 0: return 0.0
    n=len(vs)
    return (2.0*sum((i+1)*v for i,v in enumerate(vs))/(n*s)) - (n+1.0)/n

def analyze_run(run_dir: Path, quiet_steps=2, window_steps=100, min_events=2, plastic_factor=1.0, write=True):
    hist = run_dir / 'single_glider_history.csv'
    if not hist.exists(): raise FileNotFoundError(hist)
    params = {}
    pp = run_dir / 'clean_arrhenius_params.json'
    if pp.exists():
        try: params = json.loads(pp.read_text())
        except Exception: params = {}
    target = sf(params.get('target_strain'), float('inf'))
    rows=[]
    with open(hist, newline='') as f:
        for r in csv.DictReader(f):
            if None in r: continue
            e=sf(r.get('eps_total'))
            if not math.isfinite(e): continue
            if math.isfinite(target) and e > 1.05*max(target,1e-300): continue
            rows.append(r)
    if not rows: raise RuntimeError(f'no valid rows in {hist}')
    steps=[si(r.get('step'), i+1) for i,r in enumerate(rows)]
    times=[sf(r.get('time_s'),0.0) for r in rows]
    eps=[sf(r.get('eps_total'),0.0) for r in rows]
    tau=[sf(r.get('tau_MPa')) for r in rows]
    dtau=[sf(r.get('d_tau_step_MPa'),0.0) for r in rows]
    depp=[sf(r.get('d_eps_p'),0.0) for r in rows]
    de=[sf(r.get('d_eps_total'),0.0) for r in rows]
    nd=[sf(r.get('n_depin'),0.0) for r in rows]
    nc=[sf(r.get('n_capture'),0.0) for r in rows]
    live=[sf(r.get('n_live_pins'),0.0) for r in rows]
    dt=sf(params.get('dt'))
    rdt=[]
    for r in rows:
        rr=sf(r.get('crossing_rate_max_s'))
        rdt.append(rr*dt if math.isfinite(rr) and math.isfinite(dt) else float('nan'))
    med_de=median([x for x in de if x>0])
    if not math.isfinite(med_de) or med_de <= 0:
        med_de=abs(sf(params.get('strain_rate'),0.0)*dt) if math.isfinite(dt) else 0.0
    active=[(a>0) or (dp > max(plastic_factor*max(med_de,0.0), dd)) for a,dp,dd in zip(nd,depp,de)]
    aval=[]; start=None; last=None
    for i,a in enumerate(active):
        if a:
            if start is None: start=i
            last=i
        elif start is not None and (i-last)>quiet_steps:
            aval.append((start,last)); start=last=None
    if start is not None: aval.append((start,last))
    avrows=[]
    for k,(i0,i1) in enumerate(aval,1):
        idx=range(i0,i1+1)
        size=sum(max(0.0,nd[i]) for i in idx)
        psize=sum(max(0.0,depp[i]) for i in idx)
        isize=sum(max(0.0,de[i]) for i in idx)
        drop=-sum(min(0.0,dtau[i]) for i in idx)
        rise=sum(max(0.0,dtau[i]) for i in idx)
        tseg=[tau[i] for i in idx if math.isfinite(tau[i])]
        ptv=0.0
        if len(tseg)>=2:
            pk=tseg[0]
            for tv in tseg:
                pk=max(pk,tv); ptv=max(ptv, pk-tv)
        rd=[rdt[i] for i in idx if math.isfinite(rdt[i])]
        avrows.append({
            'avalanche_id':k, 'start_step':steps[i0], 'end_step':steps[i1],
            'start_time_s':times[i0], 'end_time_s':times[i1],
            'start_eps_total':eps[i0], 'end_eps_total':eps[i1],
            'duration_steps':steps[i1]-steps[i0]+1, 'duration_s':times[i1]-times[i0],
            'active_steps':sum(1 for i in idx if active[i]), 'event_size_depin':size,
            'capture_size':sum(max(0.0,nc[i]) for i in idx), 'plastic_size':psize,
            'imposed_size':isize, 'plastic_over_imposed':psize/isize if isize>0 else float('nan'),
            'stress_drop_negsum_MPa':drop, 'stress_rise_possum_MPa':rise,
            'net_d_tau_MPa':sum(dtau[i] for i in idx if math.isfinite(dtau[i])),
            'peak_to_valley_stress_drop_MPa':ptv,
            'max_n_depin_step':max([max(0.0,nd[i]) for i in idx] or [0.0]),
            'max_d_eps_p_step':max([max(0.0,depp[i]) for i in idx] or [0.0]),
            'min_d_tau_step_MPa':min([dtau[i] for i in idx if math.isfinite(dtau[i])] or [float('nan')]),
            'max_rate_dt':max(rd) if rd else float('nan'), 'n_live_pins_median':median([live[i] for i in idx])
        })
    event_av=[a for a in avrows if a['event_size_depin'] >= min_events]
    total=sum(max(0.0,x) for x in nd)
    esizes=[a['event_size_depin'] for a in event_av]
    psizes=[a['plastic_size'] for a in avrows]
    drops=[a['stress_drop_negsum_MPa'] for a in avrows]
    depvals=[max(0.0,x) for x in nd]
    topn=max(1,int(math.ceil(0.01*len(depvals))))
    topfrac=sum(sorted(depvals, reverse=True)[:topn])/total if total>0 else float('nan')
    wins=[sum(depvals[i:i+window_steps]) for i in range(0,len(depvals),window_steps)]
    mw=sum(wins)/len(wins) if wins else 0.0
    vw=sum((x-mw)**2 for x in wins)/len(wins) if wins else 0.0
    fano=vw/mw if mw>0 else float('nan')
    depidx=[i for i,x in enumerate(nd) if x>0]
    waits=[steps[j]-steps[i] for i,j in zip(depidx[:-1], depidx[1:])]
    wcv=float('nan')
    if len(waits)>=2:
        mwait=sum(waits)/len(waits); vwait=sum((w-mwait)**2 for w in waits)/len(waits)
        wcv=math.sqrt(vwait)/mwait if mwait>0 else float('nan')
    largest=max(esizes) if esizes else 0.0
    lfrac=largest/total if total>0 else float('nan')
    lplast=max(psizes) if psizes else 0.0
    ldrop=max(drops) if drops else 0.0
    reasons=[]
    if total>=20 and largest>=max(5,0.05*total): reasons.append('one avalanche contains >=5% of all depinning events')
    if math.isfinite(fano) and fano>5.0 and total>=20: reasons.append('depin-count Fano factor > 5')
    if math.isfinite(topfrac) and topfrac>0.25 and total>=20: reasons.append('top 1% of steps contain >25% of depinning events')
    if ldrop>1.0 and lplast>5.0*max(med_de,1e-300): reasons.append('stress-drop/plastic-strain burst detected')
    rdts=finite(rdt)
    summary={'run_dir':str(run_dir), 'T_K':params.get('temperature_K'), 'rho_m2':params.get('forest_rho_m2'),
             'strain_rate_s':params.get('strain_rate'), 'dt_s':params.get('dt'), 'quiet_steps':quiet_steps,
             'window_steps':window_steps, 'min_events':min_events, 'valid_history_rows':len(rows),
             'total_depin':total, 'positive_depin_steps':sum(1 for x in nd if x>0),
             'active_steps':sum(1 for x in active if x), 'n_avalanches_all_activity':len(avrows),
             'n_avalanches_event_min':len(event_av), 'largest_event_avalanche_depin':largest,
             'largest_event_avalanche_fraction':lfrac, 'largest_plastic_avalanche':lplast,
             'largest_stress_drop_negsum_MPa':ldrop, 'event_size_mean':sum(esizes)/len(esizes) if esizes else 0.0,
             'event_size_median':median(esizes), 'event_size_gini':gini(esizes), 'depin_step_gini':gini(depvals),
             'event_count_window_fano':fano, 'top_1pct_steps_event_fraction':topfrac,
             'waiting_time_steps_cv':wcv, 'median_d_eps_total':med_de, 'max_rate_dt':max(rdts) if rdts else float('nan'),
             'avalanche_like':bool(reasons), 'avalanche_like_reasons':reasons}
    if write:
        fields=list(avrows[0].keys()) if avrows else ['avalanche_id']
        with open(run_dir/'avalanche_events.csv','w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(avrows)
        (run_dir/'avalanche_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True))
        with open(run_dir/'avalanche_summary.csv','w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(summary.keys())); w.writeheader(); w.writerow(summary)
        lines=['Avalanche diagnostic','====================']
        for k in ['total_depin','positive_depin_steps','active_steps','n_avalanches_all_activity','n_avalanches_event_min','largest_event_avalanche_depin','largest_event_avalanche_fraction','largest_plastic_avalanche','largest_stress_drop_negsum_MPa','event_count_window_fano','top_1pct_steps_event_fraction','waiting_time_steps_cv','event_size_gini','depin_step_gini','max_rate_dt','avalanche_like']:
            lines.append(f'{k}: {summary.get(k)}')
        if reasons:
            lines.append('avalanche_like_reasons:'); lines += [f'  - {r}' for r in reasons]
        (run_dir/'avalanche_summary.txt').write_text('\n'.join(lines)+'\n')
    return summary

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run-dir',type=Path)
    ap.add_argument('--root',type=Path)
    ap.add_argument('--quiet-steps',type=int,default=2)
    ap.add_argument('--window-steps',type=int,default=100)
    ap.add_argument('--min-events',type=int,default=2)
    ap.add_argument('--plastic-factor',type=float,default=1.0)
    ap.add_argument('--show-table',action='store_true')
    args=ap.parse_args()
    if args.run_dir: dirs=[args.run_dir]
    elif args.root: dirs=sorted({p.parent for p in args.root.rglob('single_glider_history.csv')})
    else: raise SystemExit('Use --run-dir or --root')
    sums=[]
    for d in dirs:
        try: sums.append(analyze_run(d,args.quiet_steps,args.window_steps,args.min_events,args.plastic_factor,True))
        except Exception as e: print(f'WARN {d}: {e}')
    if args.root and sums:
        out=args.root/'avalanche_summary_all.csv'
        keys=list(sums[0].keys())
        with open(out,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(sums)
        print(f'Wrote: {out}')
    if args.show_table:
        cols=['rho_m2','total_depin','n_avalanches_event_min','largest_event_avalanche_depin','largest_event_avalanche_fraction','largest_stress_drop_negsum_MPa','event_count_window_fano','top_1pct_steps_event_fraction','avalanche_like']
        print('\t'.join(cols))
        for s in sums: print('\t'.join(str(s.get(c,'')) for c in cols))
if __name__=='__main__': main()
