#!/usr/bin/env python3
"""
MetaDrugSig leave-one-study-out (LOSO) validation of tau^2  --  streaming Python port of loso_functions.R.

Reads DEG_results_v3.tar.gz in ONE sequential pass (no unpacking, no extra disk) and, for every eligible
drug (Tier 1, not stratum C, >=3 GSEs), runs the LOSO folds and writes small per-drug CSV tables that
loso_summarize.R turns into the report.

  python3 loso_stream.py --tar DEG_results_v3.tar.gz --meta drug_metadata_v3.csv --out loso_out
         [--max_folds 10] [--n_sub 300] [--seed 42] [--min_gse 3] [--max_drugs 20] [--selftest]

Same model and metrics as the R version:  REML tau^2 by vectorised Fisher scoring (== metafor::rma REML),
SE = |logFC|/sqrt(F), held-out GSE pooled by fixed-effect IVW, PI = theta +/- 1.96*sqrt(se_h^2+se^2+tau^2).
Restartable: drugs whose outputs exist are skipped.
"""
import argparse, io, os, re, sys, tarfile, time, csv, math, json
import numpy as np
import pandas as pd

try:
    from scipy.special import ndtr as _ndtr
    def sf2(z):                      # two-sided normal p-value
        return 2.0 * _ndtr(-np.abs(z))
except Exception:                    # scipy not available -> math.erfc
    _erfc = np.frompyfunc(math.erfc, 1, 1)
    def sf2(z):
        return np.asarray(_erfc(np.abs(z) / math.sqrt(2.0)), dtype=float)

TAU_LEVELS = ["0", "(0,.01]", "(.01,.03]", "(.03,.1]", "(.1,.3]", "(.3,1]", ">1"]
TH_LEVELS = ["[0,.25)", "[.25,.5)", "[.5,1)", "[1,Inf)"]
K_LEVELS = ["2", "3-4", "5-9", "10+"]
TOPN = (100, 500, 2000)


# ---------------------------------------------------------------- REML (vectorised over genes)
def reml_core(Y, V, maxit=200, tol=1e-9):
    M = np.isfinite(Y) & np.isfinite(V) & (V > 0)
    Y = np.where(M, Y, 0.0); V = np.where(M, V, 1.0)
    k = M.sum(1)
    ok = k >= 2
    wf = M / V
    S1f = wf.sum(1); S2f = (wf ** 2).sum(1)
    with np.errstate(all="ignore"):
        thf = (wf * Y).sum(1) / S1f
        Q = (wf * (Y - thf[:, None]) ** 2).sum(1)
        tau2 = np.maximum(0.0, (Q - (k - 1)) / (S1f - S2f / S1f))
    tau2[~np.isfinite(tau2)] = 0.0
    active = ok.copy()
    with np.errstate(all="ignore"):
        for _ in range(maxit):
            W = M / (V + tau2[:, None])
            S1 = W.sum(1); S2 = (W ** 2).sum(1); S3 = (W ** 3).sum(1)
            th = (W * Y).sum(1) / S1
            R2 = (W ** 2 * (Y - th[:, None]) ** 2).sum(1)
            score = 0.5 * (R2 - (S1 - S2 / S1))
            info = 0.5 * (S2 - 2 * S3 / S1 + S2 ** 2 / S1 ** 2)
            step = score / info
            step[~np.isfinite(step)] = 0.0
            new = np.minimum(np.maximum(0.0, tau2 + step), 1e3)
            delta = np.abs(new - tau2)
            tau2 = np.where(active, new, tau2)
            active = active & (delta > tol)
            if not active.any():
                break
        W = M / (V + tau2[:, None])
        S1 = W.sum(1)
        th = (W * Y).sum(1) / S1
        s2 = (k - 1) * S1f / (S1f ** 2 - S2f)
        I2 = 100 * tau2 / (tau2 + s2)
        se = np.where(ok, 1 / np.sqrt(S1), np.nan)
        se_fe = np.where(ok, 1 / np.sqrt(S1f), np.nan)
    thf = np.where(ok, thf, np.nan); th = np.where(ok, th, np.nan)
    tau2 = np.where(ok, tau2, np.nan); I2 = np.where(ok, I2, np.nan)
    return dict(theta=th, se=se, tau2=tau2, I2=I2, k=k, theta_fe=thf, se_fe=se_fe)


def fast_reml(Y, V, chunk=8000):
    parts = [reml_core(Y[s:s + chunk], V[s:s + chunk]) for s in range(0, Y.shape[0], chunk)]
    return {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}


def bh(p):
    n = len(p); o = np.argsort(p); r = np.empty(n)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r[o] = np.minimum(q, 1.0)
    return r


def cut_tau(x):   # (-Inf,0],(0,.01],...,(1,Inf]   right-closed
    return np.searchsorted(np.array([0, .01, .03, .1, .3, 1.0]), x, side="left")


def cut_th(x):    # [-Inf,.25),[.25,.5),[.5,1),[1,Inf)  left-closed
    return np.searchsorted(np.array([.25, .5, 1.0]), x, side="right")


def cut_k(x):     # (1,2],(2,4],(4,9],(9,Inf]
    return np.searchsorted(np.array([2, 4, 9]), x, side="left")


# ---------------------------------------------------------------- one LOSO fold
def one_fold(Y, V, gse, genes, g, rng, n_sub):
    held = gse == g
    Yt, Vt, gt = Y[:, ~held], V[:, ~held], gse[~held]
    Yh, Vh = Y[:, held], V[:, held]
    okh = np.isfinite(Yh) & np.isfinite(Vh) & (Vh > 0)
    Wh = np.where(okh, 1.0 / np.where(okh, Vh, 1.0), 0.0); Yh = np.where(okh, Yh, 0.0)
    sw = Wh.sum(1)
    with np.errstate(all="ignore"):
        y_h = (Wh * Yh).sum(1) / sw
        se_h = 1 / np.sqrt(sw)
    fit = fast_reml(Yt, Vt)
    ng = np.zeros(len(genes), int)
    for u in np.unique(gt):
        ng += np.isfinite(Yt[:, gt == u]).any(1)
    ok = np.isfinite(fit["tau2"]) & (ng >= 2) & (sw > 0) & np.isfinite(y_h)
    if ok.sum() < 100:
        return None
    j = int(rng.integers(Yt.shape[1]))
    with np.errstate(all="ignore"):
        z1 = np.abs(Yt[:, j]) / np.sqrt(Vt[:, j])
    z1[~np.isfinite(z1)] = np.nan

    theta, se, tau2, I2, k = (fit[x][ok] for x in ("theta", "se", "tau2", "I2", "k"))
    theta_fe, se_fe, z1, y, se_h = fit["theta_fe"][ok], fit["se_fe"][ok], z1[ok], y_h[ok], se_h[ok]
    gn = np.asarray(genes)[ok]
    z = theta / se
    z_fe = theta_fe / se_fe
    fdr = bh(sf2(z))
    broad = (fdr < 0.05) & (np.abs(theta) > 0.5)
    robust = np.abs(theta) > 1.96 * np.sqrt(tau2 + se ** 2)
    cover_pi = np.abs(y - theta) <= 1.96 * np.sqrt(se_h ** 2 + se ** 2 + tau2)
    cover_ci = np.abs(y - theta) <= 1.96 * np.sqrt(se_h ** 2 + se ** 2)
    same = np.sign(y) == np.sign(theta)
    sig_h = np.abs(y / se_h) > 1.96
    sig_same = sig_h & same
    aligned = y * np.sign(theta)
    n = len(theta)

    df = pd.DataFrame(dict(tb=cut_tau(tau2), hb=cut_th(np.abs(theta)), kb=cut_k(k),
                           cover_pi=cover_pi, cover_ci=cover_ci, same=same, sig_h=sig_h, sig_same=sig_same,
                           aligned=aligned, abs_theta=np.abs(theta)))
    b = df.groupby(["tb", "hb", "kb"]).agg(n=("same", "size"), cover_pi=("cover_pi", "sum"), cover_ci=("cover_ci", "sum"),
                                           same=("same", "sum"), sig_h=("sig_h", "sum"), sig_same=("sig_same", "sum"),
                                           aligned=("aligned", "sum"), abs_theta=("abs_theta", "sum")).reset_index()
    b["tau_bin"] = [TAU_LEVELS[i] for i in b.pop("tb")]
    b["th_bin"] = [TH_LEVELS[i] for i in b.pop("hb")]
    b["k_bin"] = [K_LEVELS[i] for i in b.pop("kb")]

    def sumset(m, nm):
        if not m.any():
            return None
        return dict(set=nm, n=int(m.sum()), same=int(same[m].sum()), sig_h=int(sig_h[m].sum()),
                    sig_same=int(sig_same[m].sum()), aligned=float(aligned[m].sum()), sx=float(theta[m].sum()),
                    sy=float(y[m].sum()), sxx=float((theta[m] ** 2).sum()), syy=float((y[m] ** 2).sum()),
                    sxy=float((theta[m] * y[m]).sum()))
    allm = np.ones(n, bool)
    sets = [sumset(m, nm) for m, nm in [(allm, "all_evaluable"), (broad, "Broad"), (robust, "Robust"),
                                        (broad & ~robust, "Broad_not_Robust"), (robust & ~broad, "Robust_not_Broad"),
                                        (broad & robust, "Broad_and_Robust"), (~broad, "not_Broad")]]
    sets = pd.DataFrame([s for s in sets if s])

    with np.errstate(all="ignore"):
        rob = np.abs(theta) / np.sqrt(tau2 + se ** 2)
    sc = dict(single_study_z=(np.isfinite(z1), z1),
              fixed_effect_z=(allm, np.abs(z_fe)),
              random_effects_z=(allm, np.abs(z)),
              effect_size=(allm, np.abs(theta)),
              robust_score=(allm, rob),
              broad_low_tau2=(broad, -tau2),
              broad_low_I2=(broad, -I2))
    az = np.abs(z)
    rows = []
    for nm, (pool, s) in sc.items():
        idx = np.where(pool & np.isfinite(s))[0]
        o = idx[np.lexsort((-az[idx], -s[idx]))]
        for N in TOPN:
            sel = o[:N]
            if len(sel):
                rows.append(dict(method=nm, N=N, n_sel=len(sel), same=same[sel].mean(),
                                 sig_same=sig_same[sel].mean(), aligned=aligned[sel].mean()))
    topn = pd.DataFrame(rows)

    sub = None
    if n_sub > 0:
        ii = rng.choice(n, size=min(n, n_sub), replace=False)
        sub = pd.DataFrame(dict(gene=gn[ii], theta=theta[ii], se=se[ii], tau2=tau2[ii], I2=I2[ii], k=k[ii],
                                theta_fe=theta_fe[ii], se_fe=se_fe[ii], y=y[ii], se_h=se_h[ii], z=z[ii],
                                broad=broad[ii], robust=robust[ii]))
    fold = pd.DataFrame([dict(n_train_comp=int((~held).sum()), n_heldout_comp=int(held.sum()), n_eval=n,
                              n_broad=int(broad.sum()), n_robust=int(robust.sum()), median_tau2=float(np.median(tau2)))])
    for t in (b, sets, topn, sub, fold):
        if t is not None:
            t.insert(0, "fold", g)
    return dict(bin=b, sets=sets, topn=topn, sub=sub, fold=fold)


def loso_one_drug(Y, V, gse, genes, rng, max_folds=10, n_sub=300):
    ug = np.unique(gse)
    if len(ug) < 3:
        return None
    folds = rng.permutation(ug)[:max_folds]
    out = [one_fold(Y, V, gse, genes, g, rng, n_sub) for g in folds]
    out = [o for o in out if o is not None]
    if not out:
        return None
    return {nm: pd.concat([o[nm] for o in out if o[nm] is not None], ignore_index=True)
            for nm in ("bin", "sets", "topn", "sub", "fold") if any(o[nm] is not None for o in out)}


# ---------------------------------------------------------------- assemble matrices from per-comparison tables
def assemble(comps, max_se=10.0):
    """comps: list of (gse, DataFrame[gene, logFC, F]) -> Y, V (genes x comps), gse array, genes"""
    one = []
    for gse, d in comps:
        d = d.dropna(subset=["gene"]).drop_duplicates("gene")
        with np.errstate(all="ignore"):
            se = d["logFC"].abs() / np.sqrt(d["F"])
        m = np.isfinite(d["logFC"]) & np.isfinite(d["F"]) & (d["F"] > 0) & np.isfinite(se) & (se > 0) & (se <= max_se)
        d = d[m]
        if len(d):
            one.append((gse, d["gene"].astype(str).values, d["logFC"].values.astype(float), (se[m].values.astype(float)) ** 2))
    if not one:
        return None
    genes = np.unique(np.concatenate([o[1] for o in one]))
    Y = np.full((len(genes), len(one)), np.nan); V = Y.copy()
    for j, (_, g, y, v) in enumerate(one):
        i = np.searchsorted(genes, g)
        Y[i, j] = y; V[i, j] = v
    return Y, V, np.array([o[0] for o in one]), genes


def write_drug(outdir, drug, res):
    for nm, t in res.items():
        t.insert(0, "drug", drug)
        t.to_csv(os.path.join(outdir, f"{drug}__{nm}.csv" + (".gz" if nm == "sub" else "")), index=False)


# ---------------------------------------------------------------- streaming driver
def stream(args):
    meta = pd.read_csv(args.meta, dtype=str)
    elig = meta[(meta["tier"] == "1") & (meta["stratum_C"] != "TRUE") & (meta["n_gses"].astype(float) >= args.min_gse)]
    elig_ids = set(elig["DrugBank_ID"])
    print(f"eligible by metadata: {len(elig_ids)} drugs", flush=True)
    pd_dir = os.path.join(args.out, "per_drug"); os.makedirs(pd_dir, exist_ok=True)
    done = {f.split("__")[0] for f in os.listdir(pd_dir) if f.endswith("__fold.csv")}
    skipped = []
    stats = dict(drugs_done=0, folders=0, files_read=0)
    t0 = time.time()
    cur = None; cur_id = None; cur_name = None; bucket = []
    seen_ids = set()

    def flush():
        nonlocal bucket, cur_id, cur_name
        if cur_id is None or not bucket:
            bucket = []; return
        asm = assemble(bucket)
        if asm is None:
            skipped.append((cur_id, cur_name, "no usable comparisons")); bucket = []; return
        Y, V, gse, genes = asm
        if len(np.unique(gse)) < args.min_gse:
            skipped.append((cur_id, cur_name, f"only {len(np.unique(gse))} GSEs after QC")); bucket = []; return
        rng = np.random.default_rng([args.seed, int(re.sub(r"\D", "", cur_id) or 0)])
        res = loso_one_drug(Y, V, gse, genes, rng, args.max_folds, args.n_sub)
        if res is None:
            skipped.append((cur_id, cur_name, "no evaluable folds"))
        else:
            write_drug(pd_dir, cur_id, res)
            stats["drugs_done"] += 1
            print(f"[{time.time()-t0:7.0f}s] {cur_id} {cur_name}: {len(np.unique(gse))} GSEs, {Y.shape[1]} comps, "
                  f"{len(res['fold'])} folds  (done {stats['drugs_done']})", flush=True)
        bucket = []

    with tarfile.open(args.tar, mode="r|gz") as tf:
        for m in tf:
            if not m.isfile() or not m.name.endswith(".csv"):
                continue
            parts = m.name.split("/")
            folder = parts[-2]
            if folder != cur:
                flush()
                cur = folder; cur_id = None; cur_name = folder; stats["folders"] += 1
                if stats["folders"] % 100 == 0:
                    print(f"[{time.time()-t0:7.0f}s] scanned {stats['folders']} folders", flush=True)
                if args.max_drugs and stats["drugs_done"] >= args.max_drugs:
                    break
                if args.budget and time.time() - t0 > args.budget:
                    print("time budget reached - rerun the same command to resume", flush=True)
                    break
            f = tf.extractfile(m)
            if cur_id is None or cur_id == "__skip__":
                head = f.read(1500).decode("utf-8", "ignore").split("\n")
                try:
                    hdr = next(csv.reader([head[0]])); row = next(csv.reader([head[1]]))
                    did = row[hdr.index("drug_id")]
                except Exception:
                    did = "__skip__"
                if cur_id is None:
                    if did in elig_ids and did not in done and did not in seen_ids:
                        cur_id = did; seen_ids.add(did)
                    else:
                        cur_id = "__skip__"
                if cur_id == "__skip__":
                    continue
                data = (head and "\n".join(head)).encode() + f.read()
            else:
                data = f.read()
            gse = re.sub(r"^DEGs_", "", re.sub(r"_Cell_Line.*\.csv$", "", os.path.basename(m.name)))
            try:
                d = pd.read_csv(io.BytesIO(data), usecols=["gene", "logFC", "F"])
            except Exception as e:
                print("read error", m.name, e, flush=True); continue
            bucket.append((gse, d)); stats["files_read"] += 1
        else:
            flush()
    if skipped:
        pd.DataFrame(skipped, columns=["drug_id", "drug_name", "why"]).to_csv(os.path.join(args.out, "drugs_skipped.csv"), index=False)
    print(f"finished: {stats}  in {time.time()-t0:.0f}s", flush=True)


# ---------------------------------------------------------------- self-test data (simulated, known truth)
def simulate(G=3000, K=8, tau_sd=0.3, seed=0):
    rng = np.random.default_rng(seed)
    mu = rng.normal(0, 0.4, G) * (rng.random(G) < 0.3)
    tau = np.abs(rng.normal(0.15, tau_sd, G))
    se = rng.uniform(0.1, 0.5, (G, K))
    Y = mu[:, None] + rng.normal(0, 1, (G, K)) * np.sqrt(tau[:, None] ** 2 + se ** 2)
    Y[rng.random((G, K)) < 0.1] = np.nan
    V = se ** 2
    gse = np.array([f"GSE{i//2}" for i in range(K)])
    return Y, V, gse, np.array([f"g{i}" for i in range(G)])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tar"); ap.add_argument("--meta"); ap.add_argument("--out", default="loso_out")
    ap.add_argument("--max_folds", type=int, default=10); ap.add_argument("--n_sub", type=int, default=300)
    ap.add_argument("--seed", type=int, default=42); ap.add_argument("--min_gse", type=int, default=3)
    ap.add_argument("--max_drugs", type=int, default=0)
    ap.add_argument("--budget", type=float, default=0, help="stop after this many seconds (resume by rerunning)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        Y, V, gse, genes = simulate()
        res = loso_one_drug(Y, V, gse, genes, np.random.default_rng(1), 4, 100)
        b = res["bin"]
        print("PI coverage:", round(b.cover_pi.sum() / b.n.sum(), 3), " CI coverage:", round(b.cover_ci.sum() / b.n.sum(), 3))
        print(res["topn"].query("N==500")[["method", "same"]].round(3).to_string(index=False))
        print("SELF-TEST COMPLETE")
    else:
        os.makedirs(a.out, exist_ok=True)
        stream(a)
