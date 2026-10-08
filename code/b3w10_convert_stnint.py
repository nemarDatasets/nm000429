"""d3-dryad-stn-interrupt (Dryad doi:10.5061/dryad.2jm63xsq2; London et al., eLife 2021, doi:10.7554/eLife.64893) -> BIDS-iEEG.
Source: Sub01..Sub13.zip (extracted by b3w10_inspect.sh into work/<id>/extract/). Each SubXX/BlockYY holds the AlphaOmega
Neuro Omega MATLAB export of one task block (RT1D*/LT1D*.mat or ao_data.mat), the task/behaviour .mat and Plexon spike
times (.txt). Converted: the wideband RAW streams (CRAW_xx = microelectrode, CMacro_RAW_xx = macro contact) at the raw
rate (44 kHz; 22 kHz for Sub04), int16 x BitResolution / Gain -> µV, float32. The online LFP streams (CLFP, CMacro_LFP,
1.375 kHz), spike segments (CSEG), the behaviour files and spike times stay in sourcedata (renamed: initials/dates removed;
MAT header dates set to day 01).
events.tsv: CPORT__1 (digital port words; row 0 = timestamp in samples of the port clock, row 1 = value) and CDIG_IN_x
rising edges (_Up) on the RAW time axis.
Usage: python b3w10_convert_stnint.py <sourcedata_dl> <bids root>"""
import glob, json, os, re, shutil, sys
import numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b3w10_common import write_run, wtsv, wjson, sha256_file, scrub_mat_header, mat_strings, flag_strings

SRC, OUT = sys.argv[1], sys.argv[2]
X_ = os.path.join(os.path.dirname(SRC.rstrip("/")), "extract")
report = {"runs": [], "privacy_flags": {}, "renamed": {}, "sourcedata": []}
SD = os.path.join(OUT, "sourcedata", "dryad-2jm63xsq2")


def clean_name(fn):
    n = re.sub(r"_(AK|PH)_\d{8}_", "_", fn)       # Sub02_AK_04262017_OR_... -> Sub02_OR_...
    n = re.sub(r"_(AK|PH)_", "_", n)                 # Sub09_PH_HandR_...    -> Sub09_HandR_...
    return n


scans = {}
for bd in sorted(glob.glob(os.path.join(X_, "Sub*", "Block*"))):
    subn = int(re.search(r"Sub(\d+)", bd).group(1)); blk = int(re.search(r"Block(\d+)", bd).group(1))
    sub = f"sub-{subn:02d}"
    files = sorted(os.listdir(bd))
    beh = [f for f in files if re.match(r"(?i)sub\d+_.*\.mat$", f)]
    ao = [f for f in files if f.endswith(".mat") and f not in beh]
    assert len(ao) == 1 and len(beh) == 1, (bd, files)
    side = re.search(r"STN([LR])", beh[0]).group(1)
    hemi = {"L": "left", "R": "right"}[side]
    hand = re.search(r"Hand([LR])", beh[0]).group(1)
    aop = os.path.join(bd, ao[0])
    try:
        who = {n: (shape, cls) for n, shape, cls in sio.whosmat(aop)}
        full = None
    except TypeError:  # scipy whosmat fails on some Alpha Omega exports; load everything instead
        full = sio.loadmat(aop, squeeze_me=True)
        who = {n: (np.shape(v), "") for n, v in full.items() if not n.startswith("__")}
    raw = sorted([n for n in who if re.fullmatch(r"C(RAW|Macro_RAW)_\d+", n)], key=lambda n: (n.startswith("CMacro"), n))
    meta_names = [f"{n}_{k}" for n in raw for k in ("KHz", "KHz_Orig", "BitResolution", "Gain", "TimeBegin", "TimeEnd") if f"{n}_{k}" in who]
    dig = sorted({m.group(1) for n in who for m in [re.fullmatch(r"(CDIG_IN_\d+)_Up", n)] if m})
    dig_names = [f"{d}_{k}" for d in dig for k in ("Up", "Down", "KHz", "TimeBegin") if f"{d}_{k}" in who]
    port = [n for n in ("CPORT__1", "CPORT__1_KHz") if n in who]
    m = full if full is not None else sio.loadmat(aop, variable_names=raw + meta_names + dig_names + port, squeeze_me=True)
    lens = {n: int(np.asarray(m[n]).size) for n in raw}
    t0s = {n: float(m[f"{n}_TimeBegin"]) for n in raw}
    groups = {}
    for n in raw:
        groups.setdefault((lens[n], round(t0s[n], 6)), []).append(n)
    for gi, ((n_s, t0), names) in enumerate(sorted(groups.items()), 1):
        t1 = float(m[f"{names[0]}_TimeEnd"])
        fs_t = (n_s - 1) / (t1 - t0)
        khz = float(m[f"{names[0]}_KHz"])
        fs = float(khz * 1000)  # nominal rate stored by the system
        span_ok = abs(fs_t - fs) / fs < 0.001
        if not span_ok:
            report.setdefault("span_mismatch", []).append(f"Sub{subn:02d}/Block{blk:02d}: {n_s} samples over TimeBegin-TimeEnd {t1 - t0:.1f} s ({fs_t:.1f} Hz vs nominal {fs:.0f} Hz)")
        scale = [float(m[f"{n}_BitResolution"]) / float(m[f"{n}_Gain"]) for n in names]
        X = np.stack([np.asarray(m[n]).ravel() for n in names], axis=1)
        stem = f"{sub}_task-interruptedreach_run-{blk:02d}" + (f"_split-{gi}" if len(groups) > 1 else "")
        D = os.path.join(OUT, sub, "ieeg")
        # events on this RAW time axis
        ev = []
        if "CPORT__1" in m:
            P = np.asarray(m["CPORT__1"]); pk = float(m["CPORT__1_KHz"]) * 1000
            for ts, val in zip(P[0], P[1]):
                on = ts / pk - t0
                if -1 <= on <= n_s / fs + 1:
                    ev.append([round(on, 6), 0, "port", int(val), int(round(on * fs)), "CPORT__1"])
        for d in dig:
            up = np.atleast_1d(m.get(f"{d}_Up", [])); dk = float(m[f"{d}_KHz"]) * 1000; db = float(m[f"{d}_TimeBegin"])
            for u in up:
                on = db - t0 + (float(u) - 1) / dk
                ev.append([round(on, 6), 0, "ttl_up", int(u), int(round(on * fs)), d])
        ev.sort(key=lambda r: (r[0], r[2]))
        types = ["SEEG"] * len(names)
        r = write_run(D, stem, X, fs, names, ["µV"] * len(names), types, {
            "TaskName": "interruptedreach",
            "TaskDescription": "Interrupted reach task on a touchscreen: reach to a target; on some trials the target jumps mid-movement and the participant must change the ongoing movement (article).",
            "PowerLineFrequency": 60,
            "SoftwareFilters": "n/a", "HardwareFilters": "n/a",
            "iEEGReference": "Outer cannula of the microelectrode (article: Recordings were referenced to the outer cannula)",
            "iEEGPlacementScheme": f"Intraoperative microelectrode recording in the {hemi} subthalamic nucleus during DBS surgery (behaviour file name STN{side}); task performed with the {'left' if hand == 'L' else 'right'} hand",
            "Manufacturer": "Alpha Omega", "ManufacturersModelName": "n/a (Alpha Omega MATLAB export; Neuroprobe Sonus shielded tungsten microelectrodes, 1 MOhm, per the article)",
            "ElectricalStimulation": False,
        }, ch_extra=[{"description": ("microelectrode tip, wideband RAW stream" if n.startswith("CRAW") else "macro contact of the same lead, wideband RAW stream")
                      + f"; int16 x BitResolution {float(m[n + '_BitResolution'])} / Gain {float(m[n + '_Gain'])}", "group": f"STN{side}"} for n in names],
            comment=f"b3w10_convert_stnint.py: Dryad 2jm63xsq2 Sub{subn:02d}/Block{blk:02d}/{ao[0]} {','.join(names)}", scale=scale)
        wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "value", "sample", "source_channel"], ev)
        r.update(span_ok=span_ok, sub=sub, block=blk, stem=stem, channels=names, ao_file=ao[0], behaviour=beh[0], hemisphere=hemi, fs_from_times=fs_t, n_events=len(ev), time_begin=t0)
        report["runs"].append(r); scans.setdefault(sub, []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a"])
        print(stem, json.dumps({k: r[k] for k in ("n_samples", "sfreq", "fs_from_times", "channels", "roundtrip_ok", "n_events")}), flush=True)
        del X
    del m, full
    # sourcedata copy (renamed, MAT headers scrubbed)
    for f in files:
        if f.startswith(".") or f == ".DS_Store":
            continue
        src = os.path.join(bd, f); nf = clean_name(f)
        rel = os.path.join(f"Sub{subn:02d}", f"Block{blk:02d}", nf)
        dst = os.path.join(SD, rel); os.makedirs(os.path.dirname(dst), exist_ok=True)
        if nf != f:
            report["renamed"][os.path.relpath(src, X_)] = rel
        if f.endswith(".mat"):
            b = open(src, "rb").read(); b2, ch = scrub_mat_header(b); open(dst, "wb").write(b2)
            if f in beh:
                report["privacy_flags"][rel] = flag_strings(mat_strings(dst))
        else:
            shutil.copy(src, dst)
        report["sourcedata"].append(rel)
for sub, rows in scans.items():
    wtsv(os.path.join(OUT, sub, f"{sub}_scans.tsv"), ["filename", "acq_time"], sorted(rows))
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
for fn in (__file__, "b3w10_common.py", "b3w3_common.py"):
    shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.basename(fn)), os.path.join(OUT, "code", os.path.basename(fn)))
json.dump(report, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
R = report["runs"]
print(json.dumps({"runs": len(R), "subjects": len({r['sub'] for r in R}), "hours": sum(r["duration_s"] for r in R) / 3600, "rt_all": all(r["roundtrip_ok"] for r in R),
                  "renamed": report["renamed"], "privacy_flags": {k: v for k, v in report["privacy_flags"].items() if v}}, indent=1))
print("CONVERT_DONE")
