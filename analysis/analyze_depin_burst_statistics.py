#!/usr/bin/env python3
"""
Depinning burst / avalanche statistics for OpenDiS histories.

This postprocessor complements the simple avalanche flag by analyzing the
distribution of burst sizes and comparing temporal clustering to a randomized
null model.

Outputs per run:
  depin_burst_events.csv
  depin_burst_event_size_ccdf.csv
  depin_burst_summary.txt/json/csv

Outputs for --root:
  depin_burst_summary_all.csv
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd


def to_num(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def sf(x, default=float("nan")):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def gini(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    x = np.abs(x)
    if x.sum() <= 0:
        return 0.0
    x.sort()
    n = len(x)
    return float((2*np.sum(np.arange(1, n+1)*x)/(n*x.sum())) - (n+1)/n)


def clusters(active_steps, gap):
    s = np.unique(np.asarray(active_steps, int))
    if len(s) == 0:
        return []
    out = []
    a = p = int(s[0])
    for v in s[1:]:
        v = int(v)
        if v - p <= gap:
            p = v
        else:
            out.append((a, p))
            a = p = v
    out.append((a, p))
    return out


def ccdf(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v) & (v > 0)]
    if len(v) == 0:
        return pd.DataFrame(columns=["size", "ccdf", "count_ge_size"])
    xs = np.unique(np.sort(v))
    n = len(v)
    return pd.DataFrame([{"size": float(x), "ccdf": float(np.sum(v >= x)/n), "count_ge_size": int(np.sum(v >= x))} for x in xs])


def fit_power_exp(sizes, min_tail_count=12):
    x = np.asarray(sizes, float)
    x = x[np.isfinite(x) & (x > 0)]
    if len(x) < min_tail_count:
        return {"fit_status": "too_few_bursts", "n_bursts_for_fit": int(len(x))}
    best = None
    for xmin in np.unique(x):
        xt = x[x >= xmin]
        n = len(xt)
        if n < min_tail_count:
            continue
        denom = np.sum(np.log(xt/xmin))
        if denom <= 0:
            continue
        alpha = 1 + n/denom
        xs = np.sort(xt)
        emp = np.arange(1, n+1)/n
        model = 1 - (xs/xmin)**(1-alpha)
        ks = float(np.max(np.abs(emp-model)))
        if best is None or ks < best["power_ks"]:
            best = {"xmin": float(xmin), "tail_n": int(n), "power_alpha": float(alpha), "power_ks": ks}
    if best is None:
        return {"fit_status": "fit_failed", "n_bursts_for_fit": int(len(x))}
    xmin = best["xmin"]
    xt = x[x >= xmin]
    y = xt - xmin
    lam = 1/np.mean(y) if np.mean(y) > 0 else float("inf")
    xs = np.sort(xt)
    emp = np.arange(1, len(xt)+1)/len(xt)
    exp_model = 1 - np.exp(-lam*(xs-xmin)) if np.isfinite(lam) else np.ones_like(xs)
    exp_ks = float(np.max(np.abs(emp-exp_model)))
    alpha = best["power_alpha"]
    ll_power = len(xt)*math.log(alpha-1) - len(xt)*math.log(xmin) - alpha*np.sum(np.log(xt/xmin))
    ll_exp = len(xt)*math.log(lam) - lam*np.sum(y) if np.isfinite(lam) else -float("inf")
    return {"fit_status": "ok", **best, "exp_lambda": float(lam), "exp_ks": exp_ks,
            "loglik_power_minus_exp": float(ll_power-ll_exp), "n_bursts_for_fit": int(len(x))}


def load_run(run_dir):
    hp = run_dir/"single_glider_history.csv"
    ep = run_dir/"single_glider_crossing_events.csv"
    if not hp.exists():
        raise FileNotFoundError(hp)
    h = pd.read_csv(hp, low_memory=False)
    h = to_num(h, ["step","time_s","eps_total","eps_plastic","d_eps_total","d_eps_p_actual",
                   "d_tau_step_MPa","tau_MPa","tau_after_step_MPa","n_depin",
                   "frac_tau_local_capped","tau_local_median_MPa","n_live_pins"])
    h = h.dropna(subset=["step"]).copy()
    h["step"] = h["step"].astype(int)
    if "d_eps_total" not in h.columns:
        h["d_eps_total"] = h["eps_total"].diff().fillna(h["eps_total"]) if "eps_total" in h.columns else np.nan
    if "d_eps_p_actual" not in h.columns:
        h["d_eps_p_actual"] = h["eps_plastic"].diff().fillna(h["eps_plastic"]) if "eps_plastic" in h.columns else np.nan
    if "d_tau_step_MPa" not in h.columns:
        tc = "tau_after_step_MPa" if "tau_after_step_MPa" in h.columns else "tau_MPa"
        h["d_tau_step_MPa"] = h[tc].diff().fillna(0.0)
    ev = None
    if ep.exists():
        ev0 = pd.read_csv(ep, low_memory=False)
        ev0 = to_num(ev0, ["step","tau_local_MPa","rate_s","barrier_eV","work_eV","cross_force_ratio"])
        ev = ev0[ev0["event"].astype(str).eq("depin_cross")].copy() if "event" in ev0.columns else ev0.copy()
        by = ev.groupby("step").size().rename("n_depin_event").reset_index()
        h = h.merge(by, on="step", how="left")
        h["n_depin_event"] = h["n_depin_event"].fillna(0.0)
        if "tau_local_MPa" in ev.columns:
            cap = ev.groupby("step")["tau_local_MPa"].agg(depin_tau_local_max_MPa="max",
                                                          depin_tau_local_median_MPa="median").reset_index()
            h = h.merge(cap, on="step", how="left")
    else:
        h["n_depin_event"] = h["n_depin"].fillna(0.0) if "n_depin" in h.columns else 0.0
    return h, ev


def build_bursts(h, gap=1, active_plastic_ratio=1.0, stress_drop_threshold_MPa=0.0):
    dEt = h["d_eps_total"].to_numpy(float)
    dEp = h["d_eps_p_actual"].to_numpy(float)
    pr = np.divide(dEp, dEt, out=np.zeros_like(dEp), where=np.abs(dEt)>1e-300)
    drop = np.maximum(-h["d_tau_step_MPa"].to_numpy(float), 0.0)
    ndep = h["n_depin_event"].to_numpy(float)
    active = (ndep > 0) | (pr > active_plastic_ratio) | (drop > stress_drop_threshold_MPa)
    out = []
    tc = "tau_after_step_MPa" if "tau_after_step_MPa" in h.columns else "tau_MPa"
    for a,b in clusters(h.loc[active, "step"].to_numpy(int), gap):
        s = h[(h["step"] >= a) & (h["step"] <= b)].copy()
        dep = s["n_depin_event"].to_numpy(float)
        dtau = s["d_tau_step_MPa"].to_numpy(float)
        tau = s[tc].to_numpy(float) if tc in s.columns else np.array([np.nan])
        imposed = float(np.nansum(s["d_eps_total"]))
        plastic = float(np.nansum(s["d_eps_p_actual"]))
        out.append({
            "start_step": int(a), "end_step": int(b), "duration_steps": int(b-a+1),
            "event_size_depin": float(np.nansum(dep)),
            "max_n_depin_step": float(np.nanmax(dep)) if len(dep) else 0.0,
            "plastic_size": plastic, "imposed_size": imposed,
            "plastic_over_imposed": plastic/imposed if imposed > 0 else np.nan,
            "stress_drop_negsum_MPa": float(np.nansum(np.maximum(-dtau,0.0))),
            "stress_rise_possum_MPa": float(np.nansum(np.maximum(dtau,0.0))),
            "peak_to_valley_stress_MPa": float(np.nanmax(tau)-np.nanmin(tau)) if np.isfinite(tau).any() else np.nan,
            "max_depin_tau_local_MPa": sf(s["depin_tau_local_max_MPa"].max()) if "depin_tau_local_max_MPa" in s.columns else np.nan,
            "mean_frac_tau_local_capped": sf(s["frac_tau_local_capped"].mean()) if "frac_tau_local_capped" in s.columns else np.nan,
        })
    return pd.DataFrame(out)


def null_test(h, args, obs_events):
    rng = np.random.default_rng(args.seed)
    counts0 = h["n_depin_event"].to_numpy(float)
    active_counts = counts0[counts0 > 0]
    total = counts0.sum()
    if len(active_counts) == 0 or total <= 0 or args.n_boot <= 0:
        return {}
    obs_largest = float(obs_events["event_size_depin"].max()/total) if len(obs_events) else 0.0
    def metrics(counts):
        hh = h.copy()
        hh["n_depin_event"] = counts
        ev = build_bursts(hh, args.cluster_gap_steps, args.active_plastic_ratio, args.stress_drop_threshold_MPa)
        largest = float(ev["event_size_depin"].max()/counts.sum()) if len(ev) and counts.sum()>0 else 0.0
        sorted_counts = np.sort(counts)[::-1]
        topn = max(1, int(math.ceil(0.01*len(sorted_counts))))
        top = float(sorted_counts[:topn].sum()/counts.sum()) if counts.sum()>0 else 0.0
        steps = hh["step"].to_numpy(int)
        bins = np.arange(steps.min(), steps.max()+args.window_steps+1, args.window_steps)
        win,_ = np.histogram(steps, bins=bins, weights=counts)
        fano = float(np.var(win)/np.mean(win)) if np.mean(win)>0 else 0.0
        return largest, fano, top
    obs_largest, obs_fano, obs_top = metrics(counts0)
    L,F,T = [],[],[]
    n = len(counts0)
    for _ in range(args.n_boot):
        c = np.zeros(n)
        pos = rng.choice(n, size=len(active_counts), replace=False)
        c[pos] = rng.permutation(active_counts)
        a,b,c0 = metrics(c)
        L.append(a); F.append(b); T.append(c0)
    L=np.asarray(L); F=np.asarray(F); T=np.asarray(T)
    return {
        "obs_largest_burst_fraction": obs_largest,
        "obs_event_count_window_fano": obs_fano,
        "obs_top_1pct_steps_event_fraction": obs_top,
        "null_largest_fraction_mean": float(L.mean()),
        "null_largest_fraction_p_ge_obs": float(np.mean(L >= obs_largest)),
        "null_fano_mean": float(F.mean()),
        "null_fano_p_ge_obs": float(np.mean(F >= obs_fano)),
        "null_top1pct_mean": float(T.mean()),
        "null_top1pct_p_ge_obs": float(np.mean(T >= obs_top)),
    }


def meta(run_dir):
    d = {}
    try:
        nm = run_dir.name
        if nm.startswith("T") and "_rho" in nm:
            d["T_K"] = float(nm.split("_rho")[0][1:].replace("p","."))
            d["rho_m2"] = float(nm.split("_rho",1)[1].replace("p","."))
    except Exception:
        pass
    js = run_dir/"run_summary.json"
    if js.exists():
        try:
            j = json.loads(js.read_text())
            for k in ["T_K","rho_m2","strain_rate_s","dt_s","epsp_over_epstotal_final","tau_tail_median_MPa"]:
                if k in j:
                    d[k] = sf(j[k])
        except Exception:
            pass
    return d


def analyze(run_dir, args):
    h, ev = load_run(run_dir)
    bursts = build_bursts(h, args.cluster_gap_steps, args.active_plastic_ratio, args.stress_drop_threshold_MPa)
    bursts.to_csv(run_dir/"depin_burst_events.csv", index=False)
    sizes = bursts["event_size_depin"].to_numpy(float) if len(bursts) else np.array([])
    ccdf(sizes).to_csv(run_dir/"depin_burst_event_size_ccdf.csv", index=False)
    total = float(np.nansum(h["n_depin_event"].to_numpy(float)))
    if ev is not None and len(ev) and "tau_local_MPa" in ev.columns:
        tl = ev["tau_local_MPa"].to_numpy(float)
        frac_cap = float(np.mean(tl >= args.tau_cap_MPa*(1-args.cap_tol)))
        med_tl = sf(np.nanmedian(tl)); p90_tl = sf(np.nanquantile(tl,0.9))
    else:
        frac_cap = med_tl = p90_tl = np.nan
    fit = fit_power_exp(sizes, args.min_tail_count)
    nt = null_test(h, args, bursts)
    summary = {
        "run_dir": str(run_dir), **meta(run_dir),
        "n_steps": int(len(h)),
        "active_step_count": int(np.sum(h["n_depin_event"].to_numpy(float)>0)),
        "total_depin": total,
        "n_bursts": int(len(bursts)),
        "largest_burst_depin": float(np.nanmax(sizes)) if len(sizes) else 0.0,
        "largest_burst_fraction": float(np.nanmax(sizes)/total) if len(sizes) and total>0 else 0.0,
        "median_burst_depin": sf(np.nanmedian(sizes)) if len(sizes) else np.nan,
        "p90_burst_depin": sf(np.nanquantile(sizes,0.9)) if len(sizes) else np.nan,
        "p99_burst_depin": sf(np.nanquantile(sizes,0.99)) if len(sizes) else np.nan,
        "max_stress_drop_negsum_MPa": sf(bursts["stress_drop_negsum_MPa"].max()) if len(bursts) else 0.0,
        "depin_step_gini": gini(h["n_depin_event"].to_numpy(float)),
        "frac_depin_events_at_tau_cap": frac_cap,
        "depin_tau_local_median_MPa": med_tl,
        "depin_tau_local_p90_MPa": p90_tl,
        **fit, **nt
    }
    bursty_vs_null = (sf(summary.get("null_largest_fraction_p_ge_obs"),1) < 0.05 or
                      sf(summary.get("null_fano_p_ge_obs"),1) < 0.05)
    heavy_tail = (summary.get("fit_status") == "ok" and
                  sf(summary.get("loglik_power_minus_exp")) > 0 and
                  sf(summary.get("tail_n")) >= args.min_tail_count)
    summary["burst_distribution_suggestive"] = bool(bursty_vs_null or heavy_tail)
    summary["interpretation"] = (
        "cap_dominated_depin_ensemble" if sf(frac_cap) > 0.5 else
        "bursty_distribution_candidate" if summary["burst_distribution_suggestive"] else
        "smooth_or_finite_small_bursts"
    )
    (run_dir/"depin_burst_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    pd.DataFrame([summary]).to_csv(run_dir/"depin_burst_summary.csv", index=False)
    (run_dir/"depin_burst_summary.txt").write_text("\n".join(["Depinning burst statistics","===========================",""]+[f"{k}: {v}" for k,v in summary.items()])+"\n")
    return summary


def discover(root):
    return sorted(p.parent for p in Path(root).rglob("single_glider_history.csv"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path)
    ap.add_argument("--run-dir", type=Path)
    ap.add_argument("--cluster-gap-steps", type=int, default=1)
    ap.add_argument("--active-plastic-ratio", type=float, default=1.0)
    ap.add_argument("--stress-drop-threshold-MPa", type=float, default=0.0)
    ap.add_argument("--min-tail-count", type=int, default=12)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--window-steps", type=int, default=500)
    ap.add_argument("--tau-cap-MPa", type=float, default=8000.0)
    ap.add_argument("--cap-tol", type=float, default=1e-6)
    ap.add_argument("--show-table", action="store_true")
    args = ap.parse_args()
    if args.root is None and args.run_dir is None:
        ap.error("Provide --root or --run-dir")
    runs = [args.run_dir] if args.run_dir else discover(args.root)
    rows = []
    for rd in runs:
        try:
            rows.append(analyze(rd, args))
        except Exception as e:
            print(f"ERROR {rd}: {e}")
            rows.append({"run_dir": str(rd), "error": str(e)})
    df = pd.DataFrame(rows)
    if args.root:
        out = Path(args.root)/"depin_burst_summary_all.csv"
        df.to_csv(out, index=False)
        print(f"Wrote: {out}")
    if args.show_table and len(df):
        cols = ["rho_m2","total_depin","n_bursts","largest_burst_depin","largest_burst_fraction",
                "p90_burst_depin","p99_burst_depin","max_stress_drop_negsum_MPa",
                "frac_depin_events_at_tau_cap","fit_status","xmin","tail_n","power_alpha",
                "loglik_power_minus_exp","null_largest_fraction_p_ge_obs","null_fano_p_ge_obs",
                "burst_distribution_suggestive","interpretation"]
        cols = [c for c in cols if c in df.columns]
        print(df[cols].to_string(index=False))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
