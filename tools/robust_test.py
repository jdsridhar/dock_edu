"""Robustness + timing: run each case across several seeds; report RMSD + speed.

Reliability is independent of frame count (the search budget is fixed), so we use
a small frame count for the sweep, then time one dock at the app default (1000)."""
import time
import numpy as np
from docking import structures, knowledge, engine

SEEDS = [7, 8, 9, 10, 11]
print(f"{'case':26s} {'rmsd@seeds':28s} {'best':>5s} {'mean':>5s} {'succ':>5s} {'t/run':>6s}")
cases = structures.case_list()
for slug, disp, pdb, fam in cases:
    rec = structures.load_case_study(slug)
    beats = knowledge.get_beats(rec)
    rmsds, times = [], []
    for s in SEEDS:
        t0 = time.time()
        tr = engine.dock(rec, n_frames=400, seed=s, beats=beats)
        times.append(time.time() - t0)
        rmsds.append(tr["rmsd"])
    rmsds = np.array(rmsds, float)
    succ = int((rmsds < 2.0).sum())
    arr = " ".join(f"{r:4.1f}" for r in rmsds)
    print(f"{disp[:26]:26s} [{arr}] {rmsds.min():5.2f} {rmsds.mean():5.2f} "
          f"{succ}/{len(SEEDS):<2d} {np.mean(times):5.2f}s")

# one full-size timing check at the app default
rec = structures.load_case_study("cox2")
t0 = time.time(); tr = engine.dock(rec, n_frames=1000, seed=7, beats=knowledge.get_beats(rec))
print(f"\nCOX-2 @ 1000 frames: {time.time()-t0:.2f}s  rmsd={tr['rmsd']}  "
      f"events={len(tr['events'])}")
print("event timeline:", [(e['phase'], e['frame']) for e in tr['events']])
