"""
Docking trajectory engine.

Generates a frame-by-frame "docking movie" for a structure record: the ligand
approaches, explores the surface, the pocket is detected, a wrong pose is tried
and rejected, the ligand reorients, key interactions switch on, the pose is
refined, and it converges on the real crystallographic pose at the final frame.

All numbers are synthetic but physically motivated; the endpoint is the true
experimental pose, so the movie always lands on a correct, real binding mode.
"""
from __future__ import annotations

import math
import numpy as np

# ordered phases -> (start fraction, display label, colour)
PHASES = [
    ("approach",              0.00, "Approach",            "#3b82f6"),
    ("surface_scan",          0.10, "Surface Scan",        "#22d3ee"),
    ("pocket_detection",      0.26, "Pocket Detection",    "#a78bfa"),
    ("first_candidate",       0.36, "First Candidate",     "#f59e0b"),
    ("rejection",             0.46, "Rejection",           "#ef4444"),
    ("reorientation",         0.56, "Reorientation",       "#eab308"),
    ("interaction_discovery", 0.66, "Interaction Discovery", "#10b981"),
    ("refinement",            0.80, "Refinement",          "#34d399"),
    ("final_pose",            0.93, "Final Pose",          "#4ade80"),
]

# knot waypoints: (frac, (a,b,c) offset in u/v/w axes [A], (pitch,roll,yaw)[deg], quality, clash)
_KNOTS = [
    (0.00, (22.0, 3.0, 0.0),  (160, 120, 40),  0.02, 0.00),
    (0.12, (12.0, 2.0, 4.0),  (90, 60, 120),   0.08, 0.08),
    (0.22, (9.0, -5.0, 3.0),  (40, 150, 220),  0.10, 0.12),
    (0.30, (7.5, 1.0, 0.0),   (110, 60, 150),  0.22, 0.03),
    (0.38, (1.7, 1.3, 0.5),   (150, -40, 165), 0.25, 0.85),
    (0.50, (6.0, 1.0, 0.2),   (150, -40, 165), 0.18, 0.55),
    (0.60, (3.6, -1.0, 0.4),  (45, 20, -30),   0.50, 0.15),
    (0.72, (1.4, 0.2, 0.1),   (12, 8, -10),    0.72, 0.05),
    (0.86, (0.5, 0.1, 0.0),   (4, 3, -3),      0.90, 0.00),
    (1.00, (0.0, 0.0, 0.0),   (0, 0, 0),       1.00, 0.00),
]

POCKET_REACH = 9.0


def _smoothstep(t):
    t = min(1.0, max(0.0, t))
    return t * t * (3 - 2 * t)


def _interp_knots(f, idx):
    """Smoothstep-interpolate knot component `idx` (0=offset vec,1=euler,2=q,3=clash)."""
    for k in range(len(_KNOTS) - 1):
        f0 = _KNOTS[k][0]; f1 = _KNOTS[k + 1][0]
        if f0 <= f <= f1 or k == len(_KNOTS) - 2:
            t = _smoothstep((f - f0) / (f1 - f0 + 1e-9))
            a = _KNOTS[k][idx + 1]; b = _KNOTS[k + 1][idx + 1]
            if isinstance(a, tuple):
                return tuple(a[j] + (b[j] - a[j]) * t for j in range(len(a)))
            return a + (b - a) * t
    return _KNOTS[-1][idx + 1]


def _euler_matrix(pitch, roll, yaw):
    p, r, y = math.radians(pitch), math.radians(roll), math.radians(yaw)
    cx, sx = math.cos(p), math.sin(p)
    cy, sy = math.cos(r), math.sin(r)
    cz, sz = math.cos(y), math.sin(y)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def _basis(u):
    u = u / (np.linalg.norm(u) + 1e-9)
    ref = np.array([0.0, 0.0, 1.0]) if abs(u[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    v = np.cross(u, ref); v /= (np.linalg.norm(v) + 1e-9)
    w = np.cross(u, v)
    return u, v, w


def _phase_at(f):
    name = PHASES[0][0]
    for ph, start, *_ in PHASES:
        if f >= start:
            name = ph
    return name


def generate_trajectory(record, n_frames=1000, seed=7, beats=None):
    rng = np.random.default_rng(seed)
    N = int(n_frames)

    L = np.array(record["ligand_coords"], dtype=float)        # (A,3) crystal pose
    elems = list(record["ligand_elems"])
    A = len(L)
    c0 = L.mean(axis=0)
    local = L - c0

    pc = np.array(record["protein_extent"]["center"], dtype=float)
    P = np.array(record["pocket_center"], dtype=float)
    u_dir = P - pc
    if np.linalg.norm(u_dir) < 2.0:
        u_dir = np.array([1.0, 0.4, 0.2])
    u, v, w = _basis(u_dir)

    atom_phase = rng.uniform(0, 2 * math.pi, size=(A, 3))
    atom_freq = rng.uniform(0.6, 1.4, size=(A, 3))

    centroid = np.zeros((N, 3))
    translation = np.zeros((N, 3))
    euler = np.zeros((N, 3))
    torsions = np.zeros((N, 3))
    occupancy = np.zeros(N)
    shape_fit = np.zeros(N)
    clash = np.zeros(N)
    quality = np.zeros(N)
    frames_xyz = np.zeros((N, A, 3))

    for i in range(N):
        f = i / (N - 1)
        off = _interp_knots(f, 0)
        eul = list(_interp_knots(f, 1))
        q = _interp_knots(f, 2)
        c = _interp_knots(f, 3)

        # world-space centroid offset
        O = off[0] * u + off[1] * v + off[2] * w

        # exploration wander during surface scan
        if 0.10 <= f < 0.30:
            amp = 3.0 * (1 - q)
            O = O + amp * (math.sin(i * 0.20 + 1.0) * v + math.cos(i * 0.16) * w)

        # tumbling spin while searching (fades out by pocket detection)
        if f < 0.30:
            spin = max(0.0, (0.30 - f) / 0.30)
            eul[0] += (i * 5.0) % 360 * spin
            eul[2] += (i * 3.5) % 360 * spin

        # residual jitter, decaying as the pose locks in
        jit = (1 - q)
        O = O + jit * 0.35 * np.array([math.sin(i * 0.5), math.cos(i * 0.43), math.sin(i * 0.37)])
        eul = [eul[0] + jit * 6 * math.sin(i * 0.6),
               eul[1] + jit * 6 * math.cos(i * 0.52),
               eul[2] + jit * 6 * math.sin(i * 0.47)]

        cen = c0 + O
        R = _euler_matrix(*eul)
        coords = (R @ local.T).T + cen

        # intramolecular breathing (keeps bonds, conveys flexibility)
        if jit > 1e-3:
            breath = jit * 0.22 * np.sin(atom_freq * (i * 0.15) + atom_phase)
            coords = coords + breath

        frames_xyz[i] = coords
        centroid[i] = cen
        translation[i] = O
        euler[i] = eul
        quality[i] = q
        clash[i] = c

        pen = min(1.0, max(0.0, 1 - np.linalg.norm(O) / POCKET_REACH))
        occupancy[i] = 100 * pen * (0.75 + 0.25 * q)
        shape_fit[i] = 100 * min(1.0, max(0.0, pen * (0.15 + 0.85 * q) * (1 - 0.85 * c)))
        torsions[i] = [180 * (1 - q) * math.sin(i * 0.05) + 60 * q,
                       120 * (1 - q) * math.cos(i * 0.04) - 30 * q,
                       150 * (1 - q) * math.sin(i * 0.06 + 1) + 15 * q]

    # ---- interactions: assign onset frames ----
    raw = record["interactions"]
    order = {"hydrogen_bond": 0, "salt_bridge": 1, "pi_stacking": 2, "hydrophobic": 3,
             "metal_coordination": 0, "vdw": 4}
    raw_sorted = sorted(raw, key=lambda it: (order.get(it["type"], 5), it.get("dist", 9)))
    n_int = len(raw_sorted)
    onsets = np.linspace(0.66, 0.91, max(1, n_int))
    inter = []
    for j, it in enumerate(raw_sorted):
        ligp = np.array(it["lig"], dtype=float)
        lig_idx = int(np.argmin(np.linalg.norm(L - ligp, axis=1)))
        onset = int(onsets[j] * (N - 1))
        inter.append({
            "type": it["type"], "resn": it["resn"], "resi": it["resi"],
            "dist": it.get("dist"), "prot": [round(x, 2) for x in it["prot"]],
            "lig_xyz": [round(float(x), 2) for x in L[lig_idx]],
            "lig_atom": lig_idx, "onset": onset,
            "label": f"{it['resn']}{it['resi']}",
        })

    # per-frame interaction counts (only count when seated, not during rejection)
    n_hbond = np.zeros(N, dtype=int); n_pi = np.zeros(N, dtype=int)
    n_hydro = np.zeros(N, dtype=int); n_salt = np.zeros(N, dtype=int); n_active = np.zeros(N, dtype=int)
    for it in inter:
        on = it["onset"]
        n_active[on:] += 1
        if it["type"] == "hydrogen_bond":
            n_hbond[on:] += 1
        elif it["type"] == "pi_stacking":
            n_pi[on:] += 1
        elif it["type"] == "salt_bridge":
            n_salt[on:] += 1
        else:
            n_hydro[on:] += 1

    # ---- scoring ----
    occ = occupancy / 100.0
    s_clash = clash * 6.5
    s_vdw = -occ * (0.30 + 0.70 * quality) * 4.6
    s_hbond = -1.35 * n_hbond
    s_hydro = -0.42 * n_hydro - 0.55 * n_pi
    s_elec = -0.95 * n_salt - 0.30 * occ * quality
    s_desolv = -occ * quality * 1.10 + 0.25 * occ
    s_entropy = 0.20 + 1.90 * occ * quality
    total = s_clash + s_vdw + s_hbond + s_hydro + s_elec + s_desolv + s_entropy
    total = _smooth1d(total, 5)

    # ---- clash spheres: pick up to 2 bulky pocket residues ----
    bulk_rank = {"TRP": 5, "PHE": 4, "TYR": 4, "ARG": 4, "HIS": 3, "LEU": 3,
                 "ILE": 3, "MET": 3, "LYS": 3, "GLN": 2, "ASN": 2}
    pocket = record["pocket_residues"]
    clash_res = sorted(pocket, key=lambda r: bulk_rank.get(r["resn"], 1), reverse=True)[:2]
    clash_residues = [{"resn": r["resn"], "resi": r["resi"], "ca": r["ca"]} for r in clash_res]

    # ---- decoy + real pockets for the detection animation ----
    ext = record["protein_extent"]
    radius = 0.5 * float(np.linalg.norm(np.array(ext["max"]) - np.array(ext["min"])))
    hydro_frac = sum(1 for r in pocket if r["resn"] in
                     {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}) / max(1, len(pocket))
    real_pocket = {
        "pos": [round(x, 2) for x in record["pocket_center"]],
        "score": round(0.82 + 0.12 * (1 - hydro_frac) + 0.04 * min(1, len(pocket) / 14), 2),
        "volume": int(28 * len(pocket) + 6 * A),
        "depth": round(9 + 5 * min(1, len(pocket) / 16), 1),
        "drug": round(0.74 + 0.2 * min(1, len(pocket) / 15), 2),
        "hydro": round(hydro_frac, 2),
        "real": True,
    }
    decoys = []
    for d in range(5):
        dirv = rng.normal(size=3); dirv /= (np.linalg.norm(dirv) + 1e-9)
        pos = pc + dirv * radius * rng.uniform(0.45, 0.85)
        decoys.append({
            "pos": [round(float(x), 2) for x in pos],
            "score": round(float(rng.uniform(0.18, 0.6)), 2),
            "volume": int(rng.uniform(80, 320)),
            "depth": round(float(rng.uniform(3, 8)), 1),
            "drug": round(float(rng.uniform(0.15, 0.55)), 2),
            "hydro": round(float(rng.uniform(0.2, 0.7)), 2),
            "real": False,
        })

    # ---- events / annotations ----
    events = _build_events(beats, inter, clash_residues, N)

    bands = [{"phase": ph, "start": int(start * (N - 1)),
              "end": int((PHASES[k + 1][1] if k + 1 < len(PHASES) else 1.0) * (N - 1)),
              "label": label, "color": color}
             for k, (ph, start, label, color) in enumerate(PHASES)]

    def rl(arr, nd=2):
        return [round(float(x), nd) for x in arr]

    return {
        "n_frames": N,
        "ligand_pdb": _multimodel_pdb(frames_xyz, elems),
        "ligand_elems": elems,
        "ligand_natoms": A,
        "centroid": [[round(float(x), 2) for x in c] for c in centroid],
        "translation": [[round(float(x), 2) for x in t] for t in translation],
        "euler": [[round(float(x), 1) for x in e] for e in euler],
        "torsions": [[round(float(x), 1) for x in t] for t in torsions],
        "occupancy": rl(occupancy, 1),
        "shape_fit": rl(shape_fit, 1),
        "clash": rl(clash, 3),
        "scores": {
            "total": rl(total), "hbond": rl(s_hbond), "vdw": rl(s_vdw),
            "hydrophobic": rl(s_hydro), "electrostatic": rl(s_elec),
            "desolvation": rl(s_desolv), "entropy": rl(s_entropy), "clash": rl(s_clash),
        },
        "n_hbond": [int(x) for x in n_hbond],
        "n_active": [int(x) for x in n_active],
        "n_clash": [1 if c > 0.25 else 0 for c in clash],
        "interactions": inter,
        "clash_residues": clash_residues,
        "real_pocket": real_pocket,
        "decoy_pockets": decoys,
        "events": events,
        "phase_bands": bands,
        "pocket_center": [round(float(x), 2) for x in record["pocket_center"]],
        "pocket_radius": round(float(np.std(L) * 2.2 + 3), 2),
        "pocket_residues": pocket,
        "final_total": round(float(total[-1]), 2),
        "best_total": round(float(np.min(total)), 2),
    }


def _smooth1d(arr, w):
    if w < 2:
        return arr
    pad = w // 2
    padded = np.pad(arr, pad, mode="edge")  # replicate edges -> no endpoint attenuation
    k = np.ones(w) / w
    return np.convolve(padded, k, mode="valid")[:len(arr)]


_PHASE_KIND = {
    "approach": "info", "surface_scan": "info", "pocket_detection": "info",
    "first_candidate": "warn", "rejection": "bad", "reorientation": "info",
    "interaction_discovery": "good", "refinement": "good", "final_pose": "good",
}

_FALLBACK_BEATS = {
    "approach": ("Ligand enters from solution", "The ligand drifts toward the protein driven by diffusion, far from any binding site."),
    "surface_scan": ("Scanning the protein surface", "The ligand tumbles across the surface, sampling translations and rotations."),
    "pocket_detection": ("Searching for a druggable pocket", "Candidate cavities are scored; the deepest, most enclosed pocket wins."),
    "first_candidate": ("First candidate pose", "The ligand dives into the pocket in a trial orientation."),
    "rejection": ("Pose rejected — steric clash", "Atoms overlap protein residues; the penalty spikes and the pose is rejected."),
    "reorientation": ("Trying a new orientation", "The ligand flips and re-enters with a more complementary orientation."),
    "interaction_discovery": ("Key interaction discovered", "A hydrogen bond snaps into place and the score drops sharply."),
    "refinement": ("Refining the pose", "Local optimisation tightens contacts and improves shape complementarity."),
    "final_pose": ("Converged docked pose", "The ligand settles into its lowest-energy binding mode."),
}


def _build_events(beats, inter, clash_residues, N):
    bmap = {}
    if beats:
        for b in beats:
            bmap[b.get("phase")] = (b.get("title", ""), b.get("caption", ""))
    events = []
    for ph, start, *_ in PHASES:
        title, caption = bmap.get(ph) or _FALLBACK_BEATS[ph]
        events.append({"frame": int(start * (N - 1)), "phase": ph, "title": title,
                       "caption": caption, "kind": _PHASE_KIND[ph]})
    # micro-events at the hero interaction onsets
    hero = next((it for it in inter if it["type"] == "hydrogen_bond"), None)
    if hero:
        d = f" — {hero['dist']} Å" if hero.get("dist") else ""
        events.append({"frame": hero["onset"], "phase": "interaction_discovery",
                       "title": f"Hydrogen bond formed with {hero['label']}{d}",
                       "caption": "This polar contact anchors the ligand and stabilises the pose.",
                       "kind": "good"})
    pi = next((it for it in inter if it["type"] == "pi_stacking"), None)
    if pi:
        events.append({"frame": pi["onset"], "phase": "refinement",
                       "title": f"π–π stacking with {pi['label']} established",
                       "caption": "Aromatic stacking adds binding energy; the score improves.",
                       "kind": "good"})
    if clash_residues:
        cr = clash_residues[0]
        events.append({"frame": int(0.41 * (N - 1)), "phase": "first_candidate",
                       "title": f"Steric clash with {cr['resn']}{cr['resi']}",
                       "caption": "The trial orientation overlaps the residue — the pose is rejected.",
                       "kind": "bad"})
    events.sort(key=lambda e: e["frame"])
    return events


def _multimodel_pdb(frames_xyz, elems):
    """Build a multi-MODEL PDB string (one MODEL per frame) for 3Dmol trajectory."""
    A = len(elems)
    parts = []
    for i in range(frames_xyz.shape[0]):
        parts.append(f"MODEL{i + 1:>9}")
        fr = frames_xyz[i]
        for j in range(A):
            el = elems[j][:2]
            nm = (el + str(j + 1))[:4].ljust(4)
            x, y, z = fr[j]
            parts.append(f"HETATM{j + 1:>5} {nm} LIG X   1    "
                         f"{x:>8.2f}{y:>8.2f}{z:>8.2f}  1.00  0.00          {el:>2}")
        parts.append("ENDMDL")
    return "\n".join(parts)
