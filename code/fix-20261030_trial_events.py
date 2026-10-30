"""fix-20261030 (Oct 2026): task events for nm000429 from the behaviour .mat files, aligned to the Alpha Omega port words.

Sources (all in this dataset):
  sourcedata/dryad-2jm63xsq2/SubXX/BlockYY/<behaviour>.mat  trial_data{i}: start_t, end_t, event_name, event_happened,
      event_time (s from start_t), trial_type (index into trial_type_names), result, rt (behaviour PC clock, s)
  sub-XX/ieeg/*_events.tsv rows with trial_type == 'port' (CPORT__1 words on the RAW time axis)
  code/conversion_report.json  per run: n_samples, sfreq, fs_from_times (begin/end-time check), behaviour file name

Method (no author code exists for this step; every number below is measured from the files, see REPORT):
  1. Trial start = port word 21002. The sequence of 21002 onsets is matched to the behaviour start_t sequence by
     sliding the trial index offset k and taking the median absolute residual after removing the median shift.
     Candidate offsets need >= 10 paired trials. Accept the block only if the best offset pairs every start word or
     every behaviour trial (paired count = min of the two counts, one contiguous run), is unique (median residual
     < 10 ms; next candidate > 50 ms), and every pair lies within 20 ms of the least-squares clock line of step 2
     (no missing or extra start word inside the overlap).
  2. Clock map: least-squares line onset = a * t_behaviour + b over the paired starts (a ~ 0.99995: the two clocks
     drift by ~0.05 ms/s).
  3. Every happened behaviour event (event_time not NaN) of a paired trial becomes a row at
     a * (start_t + event_time) + b, plus one 'trial' row (onset start_t, duration end_t - start_t).
     port_check_ms = distance to the nearest port word of the code that marks the event in every block
     (fp_on 21004, targ_on 21006, targ_chg 21300, in_chg_window/in_knot1_window/response 21101, feedback 21015,
     trial 21002); fp_off has no word within 20 ms (21005 follows it by about 25 ms) -> n/a.
  4. Runs whose stored RAW begin/end times span less time than samples / rate (README "Timing caveat"; the RAW
     axis is exact only up to the first possible discontinuity at TimeEnd) keep only trials that end before
     (n_samples - 1) / fs_from_times.
  5. Runs listed in --skip are not changed (default: sub-02 runs; branch fix-20261018 edits those events.tsv).
Original rows (port, ttl_up) are kept unchanged; rerunning the script on its own output gives the same files.
Usage: python fix-20261030_trial_events.py <bids_root> <dir holding sourcedata/dryad-2jm63xsq2 with the .mat content> [--skip sub-02]
"""
import csv, json, os, re, sys
import numpy as np
import scipy.io as sio

ROOT, MATROOT = sys.argv[1], sys.argv[2]
SKIP = sys.argv[sys.argv.index("--skip") + 1].split(",") if "--skip" in sys.argv else ["sub-02"]
ONLY = sys.argv[sys.argv.index("--only") + 1].split(",") if "--only" in sys.argv else None
CODE = {"trial": 21002, "fp_on": 21004, "targ_on": 21006, "targ_chg": 21300, "in_chg_window": 21101,
        "in_knot1_window": 21101, "response": 21101, "feedback": 21015}
COLS = ["onset", "duration", "trial_type", "value", "sample", "source_channel", "trial", "condition", "result", "rt", "port_check_ms"]


def py(x):
    if isinstance(x, np.ndarray):
        if x.dtype == object:
            return [py(v) for v in x.ravel()]
        if x.dtype.kind == "U":
            return "".join(x.ravel().tolist()) if x.size else ""
        v = x.ravel().tolist()
        return v[0] if len(v) == 1 else v
    return x


def load_trials(p):
    out = []
    for c in sio.loadmat(p)["trial_data"].ravel():
        s = c.flat[0] if (isinstance(c, np.ndarray) and c.dtype == object) else c
        s = s.flat[0]
        d = {n: py(s[n]) for n in ("start_t", "end_t", "event_name", "event_happened", "event_time", "trial_type",
                                   "trial_type_names", "result", "rt")}
        for k in ("event_name", "event_happened", "event_time", "trial_type_names", "result"):
            if not isinstance(d[k], list):
                d[k] = [d[k]]
        out.append(d)
    return out


def match(on, st):
    res = []
    for k in range(-len(on) + 5, len(st) - 4):
        i0, i1 = max(0, -k), min(len(on), len(st) - k)
        if i1 - i0 < 10:
            continue
        a, b = on[i0:i1], st[i0 + k:i1 + k]
        r = np.abs(a - b - np.median(a - b))
        res.append((float(np.median(r)), k, i0, i1, float(r.max())))
    res.sort()
    return res


def fmt(v):
    if isinstance(v, (float, np.floating)):
        v = float(v)
        return "n/a" if v != v else ("%.6f" % v).rstrip("0").rstrip(".")
    return str(v)


rep = json.load(open(os.path.join(ROOT, "code", "conversion_report.json")))
summary = []
for run in rep["runs"]:
    sub, blk, stem = run["sub"], run["block"], run["stem"]
    evp = os.path.join(ROOT, sub, "ieeg", stem + "_events.tsv")
    s = {"run": stem, "behaviour": run["behaviour"]}
    if ONLY and sub not in ONLY:
        continue
    if sub in SKIP:
        s["status"] = "skipped (--skip)"; summary.append(s); continue
    mat = os.path.join(MATROOT, "sourcedata", "dryad-2jm63xsq2", f"Sub{int(sub[4:]):02d}", f"Block{blk:02d}", run["behaviour"])
    T = load_trials(mat)
    rows = [r for r in csv.DictReader(open(evp), delimiter="\t") if r["trial_type"] in ("port", "ttl_up")]
    P = [(float(r["onset"]), int(r["value"])) for r in rows if r["trial_type"] == "port"]
    on = np.array([o for o, v in P if v == 21002]); st = np.array([t["start_t"] for t in T])
    res = match(on, st)
    best, second = res[0], res[1]
    _, k, i0, i1, rmax = best
    s.update(n_behaviour_trials=len(T), n_start_words=len(on), offset_k=k, paired=i1 - i0,
             median_resid_ms=best[0] * 1e3, max_resid_ms=rmax * 1e3, second_best_median_ms=second[0] * 1e3)
    a_, b_ = on[i0:i1], st[i0 + k:i1 + k]
    sl, ic = np.linalg.lstsq(np.vstack([b_, np.ones_like(b_)]).T, a_, rcond=None)[0]
    lin = float(np.abs(a_ - (sl * b_ + ic)).max())
    s.update(slope=float(sl), intercept_s=float(ic), lin_max_resid_ms=lin * 1e3)
    if not (i1 - i0 == min(len(on), len(st)) and best[0] < 0.010 and second[0] > 0.050 and lin < 0.020):
        s["status"] = "not written: alignment not unambiguous"; summary.append(s); continue
    fs = float(run["sfreq"]); n = int(run["n_samples"])
    tmax = n / fs if run["span_ok"] else (n - 1) / float(run["fs_from_times"])
    s["exact_until_s"] = tmax
    codes = {}
    for o, v in P:
        codes.setdefault(v, []).append(o)
    codes = {c: np.array(x) for c, x in codes.items()}
    new = []; ntr = 0; ncut = 0; checks = {}
    for j in range(i0, i1):
        t = T[j + k]; ti = j + k + 1
        cond = t["trial_type_names"][int(t["trial_type"]) - 1]
        res_ = t["result"][0] if t["result"] and t["result"][0] else "n/a"
        rt = float(t["rt"]) if isinstance(t["rt"], (int, float)) else float("nan")
        evs = [("trial", float(t["start_t"]), float(t["end_t"]) - float(t["start_t"]))]
        for nm, h, et in zip(t["event_name"], t["event_happened"], t["event_time"]):
            if h and et == et:
                evs.append((nm, float(t["start_t"]) + float(et), 0.0))
        tr_end = sl * float(t["end_t"]) + ic
        if tr_end > tmax:
            ncut += 1; continue
        ntr += 1
        for nm, tb, du in evs:
            o = sl * tb + ic
            chk = float("nan")
            c = CODE.get(nm)
            if c in codes:
                chk = float(np.min(np.abs(codes[c] - o)) * 1e3)
                checks.setdefault(nm, []).append(chk)
            new.append([round(o, 6), round(du * sl, 6), nm, "n/a", int(round(o * fs)), "behaviour", ti, cond, res_,
                        rt, round(chk, 3) if chk == chk else float("nan")])
    s.update(trials_written=ntr, trials_cut_after_exact_limit=ncut, rows_written=len(new),
             port_check_ms={nm: {"n": len(v), "median": float(np.median(v)), "max": float(np.max(v))} for nm, v in checks.items()})
    if not new:
        s["status"] = "no rows (all paired trials after the exact limit)"; summary.append(s); continue
    old = [[r["onset"], r["duration"], r["trial_type"], r["value"], r["sample"], r["source_channel"]] + ["n/a"] * 5 for r in rows]
    allr = sorted(old + new, key=lambda r: (float(r[0]), 0 if r[5] != "behaviour" else 1))  # stable: original order kept
    with open(evp, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        for r in allr:
            w.writerow([fmt(x) if not isinstance(x, str) else x for x in r])
    s["status"] = "written"; summary.append(s)
    print(stem, json.dumps({k_: s[k_] for k_ in ("offset_k", "paired", "trials_written", "trials_cut_after_exact_limit", "rows_written", "lin_max_resid_ms")}), flush=True)
json.dump(summary, open(os.path.join(ROOT, "code", "fix-20261030_trial_events_report.json"), "w"), indent=1)
print("DONE")
