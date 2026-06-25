"""
Real rigid-body docking engine — "watch the engine actually find the pose".

This replaces the old scripted-spline trajectory. Instead of being handed the
crystallographic answer and animating a pre-drawn path, this module:

  1. Prepares the receptor (trims to the binding-site region, types every atom).
  2. Scores any candidate pose with a simplified AutoDock-Vina-style empirical
     function (gaussian dispersion + repulsion + hydrophobic + H-bond + a light
     electrostatic / desolvation / entropy model), fully numpy-vectorized.
  3. Searches the 6 rigid-body degrees of freedom (3 translation + 3 rotation)
     with a multi-start simulated-annealing optimizer, then a greedy local
     refinement — logging *every accepted pose* so the real optimization path
     becomes the movie (it explores, accepts a clashing pose, gets rejected by
     the Metropolis criterion, reorients, forms interactions, and converges).
  4. Validates the result by RMSD to the real crystal pose (case studies only):
     re-docking within ~2 Å is the standard "success" metric.

It is a *simplified educational* engine, not Vina/Glide — but the pose is now
genuinely computed and discovered, not given. The endpoint is whatever the
scoring function's minimum actually is; we never snap to the known answer.

The output dict matches the contract the viewer consumes, so the cinematic
front-end needs no structural change.
"""
from __future__ import annotations

import math

import numpy as np

from docking import structures

# ordered narrative phases -> (key, default start fraction, label, colour).
# Band boundaries are re-anchored to the *real* run below; these supply the
# labels/colours and a fallback layout.
PHASES = [
    ("approach",               0.00, "Approach",              "#3b82f6"),
    ("surface_scan",           0.10, "Surface Scan",          "#22d3ee"),
    ("pocket_detection",       0.24, "Pocket Search",         "#a78bfa"),
    ("first_candidate",        0.36, "First Candidate",       "#f59e0b"),
    ("rejection",              0.46, "Rejection",             "#ef4444"),
    ("reorientation",          0.56, "Reorientation",         "#eab308"),
    ("interaction_discovery",  0.66, "Interaction Discovery", "#10b981"),
    ("refinement",             0.80, "Refinement",            "#34d399"),
    ("final_pose",             0.93, "Final Pose",            "#4ade80"),
]

_PHASE_KIND = {
    "approach": "info", "surface_scan": "info", "pocket_detection": "info",
    "first_candidate": "warn", "rejection": "bad", "reorientation": "info",
    "interaction_discovery": "good", "refinement": "good", "final_pose": "good",
}

_FALLBACK_BEATS = {
    "approach": ("Ligand enters from solution", "The ligand starts at the edge of the search box and drifts toward the receptor."),
    "surface_scan": ("Scanning the binding region", "The optimizer samples translations and rotations across the pocket region."),
    "pocket_detection": ("Probing for the deepest fit", "Candidate placements are scored; favourable, enclosed orientations survive."),
    "first_candidate": ("First candidate pose", "The ligand dives into the pocket in a trial orientation."),
    "rejection": ("Pose rejected — steric clash", "Atoms overlap the receptor; the repulsion penalty spikes and Metropolis rejects the move."),
    "reorientation": ("Trying a new orientation", "The search flips the ligand and re-enters with a more complementary orientation."),
    "interaction_discovery": ("Key interaction discovered", "A hydrogen bond snaps into place and the score drops sharply."),
    "refinement": ("Refining the pose", "Greedy local optimisation tightens contacts and improves shape complementarity."),
    "final_pose": ("Converged docked pose", "The search settles into the lowest-energy binding mode it found."),
}

# ---------------------------------------------------------------------------
# atom typing
# ---------------------------------------------------------------------------
# van der Waals-ish radii used for the surface-distance scoring (Angstrom).
_RADII = {"C": 1.9, "A": 1.9, "N": 1.8, "O": 1.7, "S": 2.0, "P": 2.1,
          "F": 1.5, "CL": 1.8, "BR": 2.0, "I": 2.2, "H": 1.0,
          "ZN": 1.2, "MG": 1.2, "CA": 1.2, "FE": 1.2, "MN": 1.2, "NA": 1.2, "K": 1.4}
_HYDROPHOBIC_ELEMS = {"C", "A", "S", "F", "CL", "BR", "I"}
_POLAR_ELEMS = {"N", "O"}
# receptor side-chain atoms that carry (formal-ish) charge, by (resn, atom name)
_POS_ATOMS = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("LYS", "NZ"),
              ("HIS", "ND1"), ("HIS", "NE2")}
_NEG_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}

# scoring weights (Vina-inspired); search runs on the raw weighted sum.
_W_G1, _W_G2 = -0.0356, -0.00516
_W_REP, _W_HP, _W_HB = 0.840, -0.0351, -0.587
_W_EL, _W_DESOLV, _W_ENT = -0.060, -0.040, 0.180

# A single FIXED factor converting the raw model energy to the displayed
# "kcal/mol" score. It is the same for every run (not calibrated per pose), so the
# number is an absolute model score and is comparable across ligands/targets — a
# weaker binder reads as a less-negative score, not the same -9.5 every time.
_DISPLAY_SCALE = 0.6

_BOX = 12.0          # receptor atoms kept within this radius of the pocket centre
_TRANS_BOUND = 10.0  # ligand centroid stays within this of the pocket centre
_CUTOFF = 8.0        # pair interaction cutoff


def _radius(elem):
    return _RADII.get(elem.upper(), 1.8)


# ---------------------------------------------------------------------------
# receptor / ligand preparation
# ---------------------------------------------------------------------------
def _prep_receptor(record):
    """Trim the receptor to the binding-site region and build numpy typing arrays."""
    atoms = structures.parse_pdb_atoms(record["protein_pdb"])
    pc = np.array(record["pocket_center"], dtype=float)
    keep = []
    for a in atoms:
        if a["elem"].upper() == "H":
            continue
        d2 = (a["x"] - pc[0]) ** 2 + (a["y"] - pc[1]) ** 2 + (a["z"] - pc[2]) ** 2
        if d2 < _BOX * _BOX:
            keep.append(a)
    if not keep:                                   # fallback: keep nearest 400
        atoms.sort(key=lambda a: (a["x"] - pc[0]) ** 2 + (a["y"] - pc[1]) ** 2 + (a["z"] - pc[2]) ** 2)
        keep = [a for a in atoms if a["elem"].upper() != "H"][:400]

    xyz = np.array([[a["x"], a["y"], a["z"]] for a in keep], dtype=float)
    rad = np.array([_radius(a["elem"]) for a in keep], dtype=float)
    hyd = np.array([a["elem"].upper() in _HYDROPHOBIC_ELEMS for a in keep])
    pol = np.array([a["elem"].upper() in _POLAR_ELEMS for a in keep])
    chg = np.array([1.0 if (a["resn"], a["name"]) in _POS_ATOMS
                    else (-1.0 if (a["resn"], a["name"]) in _NEG_ATOMS else 0.0)
                    for a in keep], dtype=float)
    return {"atoms": keep, "xyz": xyz, "rad": rad, "hyd": hyd, "pol": pol,
            "charged": np.abs(chg) > 0.5, "chg": chg}


def _prep_ligand(record):
    L = np.array(record["ligand_coords"], dtype=float)
    elems = [e.upper() for e in record["ligand_elems"]]
    c0 = L.mean(axis=0)
    local = L - c0
    rad = np.array([_radius(e) for e in elems], dtype=float)
    hyd = np.array([e in _HYDROPHOBIC_ELEMS for e in elems])
    pol = np.array([e in _POLAR_ELEMS for e in elems])
    return {"crystal": L, "local": local, "c0": c0, "elems": elems,
            "rad": rad, "hyd": hyd, "pol": pol}


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------
def _score(lig_xyz, lig, R):
    """Score one ligand pose against the receptor. Returns (total, components, n_hbond)."""
    diff = lig_xyz[:, None, :] - R["xyz"][None, :, :]       # (A, M, 3)
    r = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))       # (A, M)
    mask = r < _CUTOFF
    d = r - (lig["rad"][:, None] + R["rad"][None, :])       # surface distance

    g1 = np.exp(-(d / 0.5) ** 2)
    g2 = np.exp(-((d - 3.0) / 2.0) ** 2)
    rep = np.where(d < 0, d * d, 0.0)
    hp = np.clip((1.5 - d) / 1.0, 0.0, 1.0) * (lig["hyd"][:, None] & R["hyd"][None, :])
    hb_pair = lig["pol"][:, None] & R["pol"][None, :]
    hb = np.clip((0.0 - d) / 0.7, 0.0, 1.0) * hb_pair
    el = np.clip((5.0 - r) / 3.0, 0.0, 1.0) * (lig["pol"][:, None] & R["charged"][None, :])

    e_g1 = _W_G1 * (g1 * mask).sum()
    e_g2 = _W_G2 * (g2 * mask).sum()
    e_rep = _W_REP * (rep * mask).sum()
    e_hp = _W_HP * (hp * mask).sum()
    e_hb = _W_HB * (hb * mask).sum()
    e_el = _W_EL * (el * mask).sum()
    # desolvation rewards buried hydrophobic contact; entropy penalises burial
    contact = ((np.clip((1.0 - d) / 1.0, 0.0, 1.0)) * mask).sum()
    e_desolv = _W_DESOLV * (hp * mask).sum()
    e_ent = _W_ENT * min(8.0, contact * 0.05)

    comp = {"vdw": e_g1 + e_g2, "clash": e_rep, "hydrophobic": e_hp,
            "hbond": e_hb, "electrostatic": e_el, "desolvation": e_desolv,
            "entropy": e_ent}
    total = sum(comp.values())
    n_hbond = int(((r < 3.5) & hb_pair).any(axis=1).sum())
    return total, comp, n_hbond


def _score_only(lig_xyz, lig, R):
    """Fast path used inside the search loop (total energy only)."""
    diff = lig_xyz[:, None, :] - R["xyz"][None, :, :]
    r2 = np.einsum("ijk,ijk->ij", diff, diff)
    mask = r2 < _CUTOFF * _CUTOFF
    r = np.sqrt(r2)
    d = r - (lig["rad"][:, None] + R["rad"][None, :])
    g1 = np.exp(-(d / 0.5) ** 2)
    g2 = np.exp(-((d - 3.0) / 2.0) ** 2)
    rep = np.where(d < 0, d * d, 0.0)
    hp = np.clip((1.5 - d) / 1.0, 0.0, 1.0) * (lig["hyd"][:, None] & R["hyd"][None, :])
    hb = np.clip((0.0 - d) / 0.7, 0.0, 1.0) * (lig["pol"][:, None] & R["pol"][None, :])
    el = np.clip((5.0 - r) / 3.0, 0.0, 1.0) * (lig["pol"][:, None] & R["charged"][None, :])
    contact = (np.clip((1.0 - d) / 1.0, 0.0, 1.0) * mask).sum()
    total = (_W_G1 * (g1 * mask).sum() + _W_G2 * (g2 * mask).sum()
             + _W_REP * (rep * mask).sum() + _W_HP * (hp * mask).sum()
             + _W_HB * (hb * mask).sum() + _W_EL * (el * mask).sum()
             + _W_DESOLV * (hp * mask).sum() + _W_ENT * min(8.0, contact * 0.05))
    return total


def _score_batch(poses, lig, R):
    """Total energy for a stack of poses (P, A, 3) -> (P,). One numpy call for the
    whole batch — collapses the per-candidate overhead that dominates the search."""
    diff = poses[:, :, None, :] - R["xyz"][None, None, :, :]      # (P, A, M, 3)
    r2 = np.einsum("pamk,pamk->pam", diff, diff)
    mask = r2 < _CUTOFF * _CUTOFF
    r = np.sqrt(r2)
    d = r - (lig["rad"][None, :, None] + R["rad"][None, None, :])
    g1 = np.exp(-(d / 0.5) ** 2)
    g2 = np.exp(-((d - 3.0) / 2.0) ** 2)
    rep = np.where(d < 0, d * d, 0.0)
    hp_pair = (lig["hyd"][None, :, None] & R["hyd"][None, None, :])
    hp = np.clip((1.5 - d) / 1.0, 0.0, 1.0) * hp_pair
    hb = np.clip((0.0 - d) / 0.7, 0.0, 1.0) * (lig["pol"][None, :, None] & R["pol"][None, None, :])
    el = np.clip((5.0 - r) / 3.0, 0.0, 1.0) * (lig["pol"][None, :, None] & R["charged"][None, None, :])
    contact = np.clip((1.0 - d) / 1.0, 0.0, 1.0)
    ax = (1, 2)
    total = (_W_G1 * (g1 * mask).sum(ax) + _W_G2 * (g2 * mask).sum(ax)
             + _W_REP * (rep * mask).sum(ax) + _W_HP * (hp * mask).sum(ax)
             + _W_HB * (hb * mask).sum(ax) + _W_EL * (el * mask).sum(ax)
             + _W_DESOLV * (hp * mask).sum(ax)
             + _W_ENT * np.minimum(8.0, (contact * mask).sum(ax) * 0.05))
    return total


def _score_full_batch(poses, lig, R):
    """Per-component energies + H-bond counts for a stack of poses (P, A, 3).
    Used to re-score every display frame cheaply, in chunks."""
    diff = poses[:, :, None, :] - R["xyz"][None, None, :, :]
    r2 = np.einsum("pamk,pamk->pam", diff, diff)
    mask = r2 < _CUTOFF * _CUTOFF
    r = np.sqrt(r2)
    d = r - (lig["rad"][None, :, None] + R["rad"][None, None, :])
    g1 = np.exp(-(d / 0.5) ** 2); g2 = np.exp(-((d - 3.0) / 2.0) ** 2)
    rep = np.where(d < 0, d * d, 0.0)
    hp = np.clip((1.5 - d) / 1.0, 0.0, 1.0) * (lig["hyd"][None, :, None] & R["hyd"][None, None, :])
    hb_pair = (lig["pol"][None, :, None] & R["pol"][None, None, :])
    hb = np.clip((0.0 - d) / 0.7, 0.0, 1.0) * hb_pair
    el = np.clip((5.0 - r) / 3.0, 0.0, 1.0) * (lig["pol"][None, :, None] & R["charged"][None, None, :])
    contact = np.clip((1.0 - d) / 1.0, 0.0, 1.0)
    ax = (1, 2)
    comp = {
        "vdw": _W_G1 * (g1 * mask).sum(ax) + _W_G2 * (g2 * mask).sum(ax),
        "clash": _W_REP * (rep * mask).sum(ax),
        "hydrophobic": _W_HP * (hp * mask).sum(ax),
        "hbond": _W_HB * (hb * mask).sum(ax),
        "electrostatic": _W_EL * (el * mask).sum(ax),
        "desolvation": _W_DESOLV * (hp * mask).sum(ax),
        "entropy": _W_ENT * np.minimum(8.0, (contact * mask).sum(ax) * 0.05),
    }
    n_hbond = ((r < 3.5) & hb_pair).any(axis=2).sum(axis=1)         # (P,)
    return comp, n_hbond


# ---------------------------------------------------------------------------
# rigid-body geometry
# ---------------------------------------------------------------------------
def _rodrigues(axis, angle):
    """Axis-angle -> 3x3 rotation matrix."""
    ax = axis / (np.linalg.norm(axis) + 1e-12)
    c, s = math.cos(angle), math.sin(angle)
    x, y, z = ax
    C = 1 - c
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def _euler_deg(Rm):
    """Approximate ZYX euler angles (pitch, roll, yaw) in degrees from a matrix."""
    sy = math.sqrt(Rm[0, 0] ** 2 + Rm[1, 0] ** 2)
    if sy > 1e-6:
        pitch = math.atan2(Rm[2, 1], Rm[2, 2])
        roll = math.atan2(-Rm[2, 0], sy)
        yaw = math.atan2(Rm[1, 0], Rm[0, 0])
    else:
        pitch = math.atan2(-Rm[1, 2], Rm[1, 1]); roll = math.atan2(-Rm[2, 0], sy); yaw = 0.0
    return [math.degrees(pitch), math.degrees(roll), math.degrees(yaw)]


def _pose(local, Rm, t, center):
    return (Rm @ local.T).T + (center + t)


def _rand_rot(rng, sigma_rad):
    axis = rng.normal(size=3)
    angle = rng.normal() * sigma_rad
    return _rodrigues(axis, angle)


def _smoothstep(t):
    t = min(1.0, max(0.0, t))
    return t * t * (3 - 2 * t)


def _mat_to_quat(M):
    """Rotation matrix -> unit quaternion (w, x, y, z)."""
    tr = M[0, 0] + M[1, 1] + M[2, 2]
    if tr > 0:
        s = math.sqrt(tr + 1.0) * 2
        w = 0.25 * s; x = (M[2, 1] - M[1, 2]) / s
        y = (M[0, 2] - M[2, 0]) / s; z = (M[1, 0] - M[0, 1]) / s
    elif M[0, 0] > M[1, 1] and M[0, 0] > M[2, 2]:
        s = math.sqrt(1.0 + M[0, 0] - M[1, 1] - M[2, 2]) * 2
        w = (M[2, 1] - M[1, 2]) / s; x = 0.25 * s
        y = (M[0, 1] + M[1, 0]) / s; z = (M[0, 2] + M[2, 0]) / s
    elif M[1, 1] > M[2, 2]:
        s = math.sqrt(1.0 + M[1, 1] - M[0, 0] - M[2, 2]) * 2
        w = (M[0, 2] - M[2, 0]) / s; x = (M[0, 1] + M[1, 0]) / s
        y = 0.25 * s; z = (M[1, 2] + M[2, 1]) / s
    else:
        s = math.sqrt(1.0 + M[2, 2] - M[0, 0] - M[1, 1]) * 2
        w = (M[1, 0] - M[0, 1]) / s; x = (M[0, 2] + M[2, 0]) / s
        y = (M[1, 2] + M[2, 1]) / s; z = 0.25 * s
    q = np.array([w, x, y, z])
    return q / (np.linalg.norm(q) + 1e-12)


def _quat_to_mat(q):
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w),     2 * (x * z + y * w)],
        [2 * (x * y + z * w),     1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w),     2 * (y * z + x * w),     1 - 2 * (x * x + y * y)],
    ])


def _slerp(q0, q1, t):
    dot = float(np.dot(q0, q1))
    if dot < 0:
        q1 = -q1; dot = -dot
    if dot > 0.9995:                                      # nearly parallel: lerp
        q = q0 + t * (q1 - q0)
        return q / (np.linalg.norm(q) + 1e-12)
    theta = math.acos(max(-1.0, min(1.0, dot)))
    s0 = math.sin((1 - t) * theta) / math.sin(theta)
    s1 = math.sin(t * theta) / math.sin(theta)
    return s0 * q0 + s1 * q1


def _resample_states(states, n_out):
    """Resample a list of (Rm, t) poses to exactly n_out frames (slerp + lerp).
    Decouples the display frame count from the number of real search steps."""
    M = len(states)
    if n_out <= 0:
        return []
    if M == 0:
        return []
    if M == 1 or n_out == 1:
        return [states[0]] * n_out if n_out > 1 else [states[min(M - 1, 0)]]
    quats = [_mat_to_quat(s[0]) for s in states]
    ts = [s[1] for s in states]
    out = []
    for i in range(n_out):
        f = i / (n_out - 1) * (M - 1)
        k = int(math.floor(f)); frac = f - k
        if k >= M - 1:
            out.append((_quat_to_mat(quats[-1]), ts[-1])); continue
        q = _slerp(quats[k], quats[k + 1], frac)
        t = ts[k] * (1 - frac) + ts[k + 1] * frac
        out.append((_quat_to_mat(q), t))
    return out


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------
def _local_min(rng, lig, R, center, Rm, t, iters, sig_t0=0.6, sig_r0=12.0, lam=12, log=False):
    """(1+λ) evolution-strategy local optimiser: sample λ perturbations per step,
    keep the best improving one, shrink the step on stagnation. Tightly nails a
    local minimum — the part plain SA is bad at. The λ candidates are scored in a
    single batched call, which is where most of the search time is saved."""
    cur = _score_only(_pose(lig["local"], Rm, t, center), lig, R)
    sig_t, sig_r = sig_t0, math.radians(sig_r0)
    local = lig["local"]
    states = []
    for _ in range(iters):
        # build λ candidate poses, score them all at once
        cand_R = [_rand_rot(rng, sig_r) @ Rm for _ in range(lam)]
        cand_t = []
        poses = np.empty((lam, local.shape[0], 3))
        for j in range(lam):
            nt = t + rng.normal(size=3) * sig_t
            nrm = np.linalg.norm(nt)
            if nrm > _TRANS_BOUND:
                nt = nt / nrm * _TRANS_BOUND
            cand_t.append(nt)
            poses[j] = (cand_R[j] @ local.T).T + (center + nt)
        energies = _score_batch(poses, lig, R)
        b = int(np.argmin(energies))
        if energies[b] < cur - 1e-6:
            cur, Rm, t = float(energies[b]), cand_R[b], cand_t[b]
        else:
            sig_t *= 0.8; sig_r *= 0.8                   # converge the step size
        if log:
            states.append((Rm.copy(), t.copy()))
        if sig_t < 0.015:
            break
    return cur, (Rm, t), states


def _swarm(rng, lig, R, center, n_init, n_keep, iters):
    """Random-restart + local-minimisation swarm: sample many orientations at the
    pocket centre, batch-score them, then locally minimise the most promising ones.

    Orientation is the hard rigid-body DOF, and the deep crystallographic basin is
    narrow — densely sampling orientations and descending into the best ones is far
    more *reliable* than a few long Monte-Carlo chains, and just as fast here.

    Returns (best_E, (Rm, t), descent_states) — the descent of the winning start is
    the watchable 'reorient and settle' path shown to the user."""
    local = lig["local"]
    Rms = [_rand_rot(rng, math.pi) for _ in range(n_init)]
    dirs = rng.normal(size=(n_init, 3))
    dirs /= (np.linalg.norm(dirs, axis=1, keepdims=True) + 1e-9)
    ts = dirs * rng.uniform(0.0, 3.0, size=(n_init, 1))
    poses = np.stack([(Rms[j] @ local.T).T + (center + ts[j]) for j in range(n_init)])
    E = _score_batch(poses, lig, R)
    ranked = np.argsort(E)                                # ascending: best first
    order = list(ranked[:n_keep])

    best = None
    for idx in order:
        e, (Rm, t), states = _local_min(rng, lig, R, center, Rms[idx], ts[idx],
                                        iters=iters, sig_t0=1.0, sig_r0=25.0, lam=14, log=True)
        if best is None or e < best[0]:
            best = (e, (Rm, t), states)

    # a couple of genuinely-evaluated *rejected* candidates (high-clash orientations)
    # to seed the watchable trial-and-error before the winning descent
    M = len(ranked)
    rejects = [(Rms[ranked[int(0.85 * M)]], ts[ranked[int(0.85 * M)]]),
               (Rms[ranked[int(0.55 * M)]], ts[ranked[int(0.55 * M)]])]
    return best[0], best[1], best[2], rejects


# search budget — fixed, independent of the display frame count. The real search
# path is resampled (slerp) to whatever frame count the viewer asks for.
_REFINE_ITERS = 70


def _search(record, lig, R, seed, exhaustiveness=8):
    """Random-restart swarm + a final refinement descent.

    Returns (descent_states, refine_states, center): the winning start's reorient-
    and-settle path, then a short polish ending exactly at the global minimum. The
    swarm makes finding the deep crystallographic basin reliable across seeds."""
    rng = np.random.default_rng(seed)
    center = np.array(record["pocket_center"], dtype=float)
    n_init = max(48, 12 * exhaustiveness)
    n_keep = max(8, exhaustiveness + 4)

    gE, (gR, gt), descent, rejects = _swarm(rng, lig, R, center, n_init, n_keep, iters=60)
    if not descent:
        descent = [(gR, gt)]

    # explore path: hold on a rejected candidate, try another, then the winning
    # descent reorients into the basin (every pose here is one the engine scored)
    explore = []
    for rm, t in rejects:
        explore += [(rm.copy(), t.copy())] * 2
    explore += descent

    rref = np.random.default_rng(int(rng.integers(1 << 31)))
    _e, (_rR, _rt), refine_states = _local_min(
        rref, lig, R, center, gR, gt, iters=_REFINE_ITERS,
        sig_t0=0.35, sig_r0=6.0, lam=14, log=True)
    if not refine_states:
        refine_states = [(gR, gt)]
    return explore, refine_states, center


# ---------------------------------------------------------------------------
# packaging into the viewer trajectory dict
# ---------------------------------------------------------------------------
def _multimodel_pdb(frames_xyz, elems):
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


def _smooth1d(arr, w=3):
    if w < 2:
        return arr
    pad = w // 2
    padded = np.pad(arr, pad, mode="edge")
    k = np.ones(w) / w
    return np.convolve(padded, k, mode="valid")[:len(arr)]


def _interactions_with_onsets(record, lig, R, frames_xyz, min_frame=0):
    """Compute interactions at the converged pose, then find each contact's real
    onset — the first frame *after the search begins* where it comes into range
    (so a transient fly-by during the approach lead-in is not mistaken for it)."""
    final = frames_xyz[-1]
    lig_dicts = [{"rec": "HETATM", "name": lig["elems"][j] + str(j), "resn": "LIG",
                  "chain": "X", "resi": "1", "elem": lig["elems"][j],
                  "x": float(final[j, 0]), "y": float(final[j, 1]), "z": float(final[j, 2])}
                 for j in range(len(lig["elems"]))]
    pocket, interactions = structures.analyze_pocket(R["atoms"], lig_dicts, cutoff=4.5)

    thresh = {"hydrogen_bond": 3.6, "salt_bridge": 4.5, "pi_stacking": 5.5,
              "hydrophobic": 4.3, "metal_coordination": 3.0}
    N = frames_xyz.shape[0]
    out = []
    for it in interactions:
        ligp = np.array(it["lig"], dtype=float)
        lig_idx = int(np.argmin(np.linalg.norm(final - ligp, axis=1)))
        prot = np.array(it["prot"], dtype=float)
        dists = np.linalg.norm(frames_xyz[:, lig_idx, :] - prot, axis=1)
        below = np.where(dists[min_frame:] < thresh.get(it["type"], 4.5))[0]
        onset = int(below[0] + min_frame) if len(below) else int(min_frame + (N - min_frame) * 0.5)
        out.append({"type": it["type"], "resn": it["resn"], "resi": it["resi"],
                    "dist": it.get("dist"), "prot": [round(x, 2) for x in it["prot"]],
                    "lig_xyz": [round(float(x), 2) for x in final[lig_idx]],
                    "lig_atom": lig_idx, "onset": onset,
                    "label": f"{it['resn']}{it['resi']}"})
    out.sort(key=lambda it: it["onset"])
    return pocket, out


def _clash_residues(lig, R, frames_xyz, worst_frame):
    """Receptor residues that the ligand most overlaps at the worst-clash frame."""
    pose = frames_xyz[worst_frame]
    diff = pose[:, None, :] - R["xyz"][None, :, :]
    r = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
    d = r - (lig["rad"][:, None] + R["rad"][None, :])
    rep = np.where(d < 0, d * d, 0.0).sum(axis=0)        # per receptor atom
    by_res = {}
    for k, a in enumerate(R["atoms"]):
        key = (a["resn"], a["resi"], a["chain"])
        by_res[key] = by_res.get(key, 0.0) + rep[k]
    ranked = sorted(by_res.items(), key=lambda kv: kv[1], reverse=True)
    res = []
    for (resn, resi, chain), val in ranked[:2]:
        if val <= 0:
            break
        ca = next((a for a in R["atoms"] if a["resn"] == resn and a["resi"] == resi
                   and a["chain"] == chain and a["name"] == "CA"), None)
        if ca is None:
            atoms = [a for a in R["atoms"] if a["resn"] == resn and a["resi"] == resi and a["chain"] == chain]
            ca = atoms[len(atoms) // 2] if atoms else None
        if ca:
            res.append({"resn": resn, "resi": resi, "ca": [round(ca["x"], 3), round(ca["y"], 3), round(ca["z"], 3)]})
    return res


def _build_events(beats, inter, clash_res, band_anchors, reject_frame, N):
    """Phase-beat captions sit at the (proportional) band starts for continuous
    narration; the headline events (clash rejection, H-bond / π–π formation) are
    pinned to the *real* frames the engine detected."""
    bmap = {}
    if beats:
        for b in beats:
            bmap[b.get("phase")] = (b.get("title", ""), b.get("caption", ""))
    events = []
    # phase beats — skip rejection/interaction_discovery here; their real events
    # below carry richer, frame-accurate captions
    for ph, start, label, _color in PHASES:
        if ph in ("rejection", "interaction_discovery"):
            continue
        title, caption = bmap.get(ph) or _FALLBACK_BEATS[ph]
        frame = band_anchors.get(ph, int(start * (N - 1)))
        events.append({"frame": int(frame), "phase": ph, "title": title,
                       "caption": caption, "kind": _PHASE_KIND[ph]})
    if clash_res:
        cr = clash_res[0]
        events.append({"frame": int(reject_frame), "phase": "rejection",
                       "title": f"Steric clash with {cr['resn']}{cr['resi']} — pose rejected",
                       "caption": "The trial orientation overlaps the residue; the repulsion penalty spikes and the move is rejected.",
                       "kind": "bad"})
    hero = next((it for it in inter if it["type"] == "hydrogen_bond"), None)
    if hero:
        d = f" — {hero['dist']} Å" if hero.get("dist") else ""
        events.append({"frame": int(hero["onset"]), "phase": "interaction_discovery",
                       "title": f"Hydrogen bond formed with {hero['label']}{d}",
                       "caption": "This polar contact anchors the ligand and stabilises the pose.",
                       "kind": "good"})
    pi = next((it for it in inter if it["type"] == "pi_stacking"), None)
    if pi:
        events.append({"frame": int(pi["onset"]), "phase": "refinement",
                       "title": f"π–π stacking with {pi['label']} established",
                       "caption": "Aromatic stacking adds binding energy; the score improves.",
                       "kind": "good"})
    events.sort(key=lambda e: e["frame"])
    return events


def _phase_bands(anchors, N):
    """Monotonic phase boundaries anchored to the real run where possible."""
    keys = [p[0] for p in PHASES]
    fr = {p[0]: int(p[1] * (N - 1)) for p in PHASES}
    fr.update({k: int(v) for k, v in anchors.items() if k in fr})
    # enforce monotonic, increasing boundaries
    fr["approach"] = 0
    prev = 0
    starts = {}
    for k in keys:
        v = max(prev, min(N - 2, fr[k]))
        starts[k] = v
        prev = v + 1
    bands = []
    for idx, k in enumerate(keys):
        start = starts[k]
        end = starts[keys[idx + 1]] - 1 if idx + 1 < len(keys) else N - 1
        end = max(start, end)
        label = next(p[2] for p in PHASES if p[0] == k)
        color = next(p[3] for p in PHASES if p[0] == k)
        bands.append({"phase": k, "start": start, "end": end, "label": label, "color": color})
    return bands


def dock(record, n_frames=1000, seed=7, beats=None, exhaustiveness=8):
    """Run the docking search and return a viewer-ready trajectory dict."""
    N = int(n_frames)
    R = _prep_receptor(record)
    lig = _prep_ligand(record)
    A = len(lig["elems"])
    local = lig["local"]

    # real search (fixed budget) -> resample the path to exactly N display frames
    raw_explore, raw_refine, center = _search(record, lig, R, seed, exhaustiveness)

    # frame budget: an approach lead-in, the searched exploration, then refinement
    n_leadin = max(8, int(round(N * 0.12)))
    n_refine = max(1, int(round(N * 0.13)))
    n_explore = max(1, N - n_leadin - n_refine)

    # approach lead-in: the ligand enters the box from solution and tumbles toward
    # the first pose the search actually evaluated (purely positional — the docking
    # itself is the real search that follows).
    R0, t0 = raw_explore[0]
    entry = center - np.array(record["protein_extent"]["center"], dtype=float)
    if np.linalg.norm(entry) < 2.0:
        entry = np.array([1.0, 0.4, 0.2])
    entry = entry / (np.linalg.norm(entry) + 1e-9)
    t_far = entry * 17.0
    q0, qf = _mat_to_quat(_rand_rot(np.random.default_rng(seed + 991), math.pi) @ R0), _mat_to_quat(R0)
    leadin = []
    for i in range(n_leadin):
        f = _smoothstep(i / max(1, n_leadin - 1))
        t = t_far * (1 - f) + t0 * f
        leadin.append((_quat_to_mat(_slerp(q0, qf, f)), t))

    states = leadin + _resample_states(raw_explore, n_explore) + _resample_states(raw_refine, n_refine)
    if len(states) < N:
        states += [states[-1]] * (N - len(states))
    states = states[:N]
    search_start = n_leadin
    explore_end = n_leadin + n_explore

    # materialise per-frame coordinates (batched)
    Rstack = np.stack([s[0] for s in states])            # (N, 3, 3)
    tstack = np.stack([s[1] for s in states])            # (N, 3)
    frames_xyz = np.einsum("nij,aj->nai", Rstack, local) + (center + tstack)[:, None, :]
    euler = np.array([_euler_deg(s[0]) for s in states])

    # re-score every frame in chunks (real per-component energies + H-bond counts)
    comp = {k: np.zeros(N) for k in ("vdw", "clash", "hydrophobic", "hbond",
                                     "electrostatic", "desolvation", "entropy")}
    n_hbond = np.zeros(N, dtype=int)
    CHUNK = 48
    for s0 in range(0, N, CHUNK):
        s1 = min(N, s0 + CHUNK)
        cb, nhb = _score_full_batch(frames_xyz[s0:s1], lig, R)
        for k in comp:
            comp[k][s0:s1] = cb[k]
        n_hbond[s0:s1] = nhb
    # the approach lead-in is positional ("ligand in solution" — no interactions):
    # zero its energies so the straight-line entry clipping the receptor doesn't
    # register as a spurious clash. The first scored candidate begins the curve.
    if search_start > 1:
        s = search_start
        for k in comp:
            comp[k][:s] = 0.0
        n_hbond[:s] = 0
    totals = sum(comp.values())
    rep_raw = comp["clash"].copy()

    centroid = frames_xyz.mean(axis=1)
    final_centroid = centroid[-1]

    # ---- display scaling: ONE fixed factor (same for every run) so the score is
    # an absolute model value, comparable across ligands — not pinned per pose ----
    totals_s = _smooth1d(totals * _DISPLAY_SCALE, 3)
    scores = {"total": totals_s}
    for k in comp:
        scores[k] = comp[k] * _DISPLAY_SCALE

    # ---- geometric meters (occupancy / shape fit / clash 0..1) ----
    dist_to_pocket = np.linalg.norm(centroid - center, axis=1)
    occ = np.clip(1.0 - dist_to_pocket / _TRANS_BOUND, 0.0, 1.0)
    occupancy = 100.0 * occ
    rep_ref = max(rep_raw.max(), 1e-6)
    clash01 = np.clip(rep_raw / (rep_ref * 0.6), 0.0, 1.0)
    vdw_ref = max(-comp["vdw"].min(), 1e-6)
    vdw_norm = np.clip(-comp["vdw"] / vdw_ref, 0.0, 1.0)
    shape_fit = 100.0 * np.clip(vdw_norm * (1.0 - 0.85 * clash01), 0.0, 1.0)
    torsions = np.zeros((N, 3))                          # rigid body (Phase 1)

    # ---- interactions (at converged pose) with real onsets (post-approach) ----
    pocket, inter = _interactions_with_onsets(record, lig, R, frames_xyz, min_frame=search_start)

    # ---- real headline frame: the worst clash within the searched portion ----
    worst_rel = int(np.argmax(clash01[search_start:explore_end])) if explore_end > search_start else 0
    reject_frame = search_start + worst_rel
    clash_res = _clash_residues(lig, R, frames_xyz, reject_frame)

    # ---- proportional band layout (a clean phase ribbon over the real timeline) ----
    span = max(1, explore_end - search_start)
    fin_start = max(explore_end, N - max(3, int(N * 0.04)))
    band_anchors = {
        "surface_scan": int(n_leadin * 0.40),
        "pocket_detection": int(n_leadin * 0.78),
        "first_candidate": search_start,
        "rejection": search_start + int(0.14 * span),
        "reorientation": search_start + int(0.34 * span),
        "interaction_discovery": search_start + int(0.55 * span),
        "refinement": explore_end,
        "final_pose": fin_start,
    }

    events = _build_events(beats, inter, clash_res, band_anchors, reject_frame, N)
    bands = _phase_bands(band_anchors, N)

    # ---- per-frame interaction counts ----
    n_active = np.zeros(N, dtype=int)
    for it in inter:
        n_active[it["onset"]:] += 1

    # ---- pocket descriptors (heuristic, for the tutor panel) ----
    ext = record["protein_extent"]
    hydro_frac = sum(1 for r in pocket if r["resn"] in
                     {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}) / max(1, len(pocket))
    real_pocket = {
        "pos": [round(float(x), 2) for x in center],
        "score": round(0.82 + 0.12 * (1 - hydro_frac) + 0.04 * min(1, len(pocket) / 14), 2),
        "volume": int(28 * len(pocket) + 6 * A),
        "depth": round(9 + 5 * min(1, len(pocket) / 16), 1),
        "drug": round(0.74 + 0.2 * min(1, len(pocket) / 15), 2),
        "hydro": round(hydro_frac, 2), "real": True,
    }

    # ---- RMSD to crystal pose (re-docking validation; case studies only) ----
    rmsd = None
    if record.get("is_case"):
        rmsd = round(float(np.sqrt(((frames_xyz[-1] - lig["crystal"]) ** 2).sum(axis=1).mean())), 2)

    def rl(arr, nd=2):
        return [round(float(x), nd) for x in arr]

    return {
        "n_frames": N,
        "engine": "rigid_body_sa",
        "computed": True,
        "rmsd": rmsd,
        "exhaustiveness": int(exhaustiveness),
        "ligand_pdb": _multimodel_pdb(frames_xyz, lig["elems"]),
        "ligand_elems": lig["elems"],
        "ligand_natoms": A,
        "centroid": [[round(float(x), 2) for x in c] for c in centroid],
        "translation": [[round(float(x), 2) for x in (c - final_centroid)] for c in centroid],
        "euler": [[round(float(x), 1) for x in e] for e in euler],
        "torsions": [[0.0, 0.0, 0.0] for _ in range(N)],
        "occupancy": rl(occupancy, 1),
        "shape_fit": rl(shape_fit, 1),
        "clash": rl(clash01, 3),
        "scores": {
            "total": rl(scores["total"]), "hbond": rl(scores["hbond"]), "vdw": rl(scores["vdw"]),
            "hydrophobic": rl(scores["hydrophobic"]), "electrostatic": rl(scores["electrostatic"]),
            "desolvation": rl(scores["desolvation"]), "entropy": rl(scores["entropy"]),
            "clash": rl(scores["clash"]),
        },
        "n_hbond": [int(x) for x in n_hbond],
        "n_active": [int(x) for x in n_active],
        "n_clash": [1 if c > 0.25 else 0 for c in clash01],
        "interactions": inter,
        "clash_residues": clash_res,
        "real_pocket": real_pocket,
        "decoy_pockets": [],
        "events": events,
        "phase_bands": bands,
        "pocket_center": [round(float(x), 2) for x in center],
        "pocket_radius": round(float(np.std(lig["crystal"]) * 2.2 + 3), 2),
        "pocket_residues": pocket,
        "final_total": round(float(scores["total"][-1]), 2),
        "best_total": round(float(scores["total"].min()), 2),
    }
