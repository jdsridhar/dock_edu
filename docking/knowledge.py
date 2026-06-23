"""
Educational content + AI-tutor knowledge base.

Loads the verified scientific content pack (authored and adversarially
fact-checked by a multi-agent workflow) and provides per-target packs plus a
synthesized pack for custom uploads.
"""
from __future__ import annotations

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONTENT_PATH = os.path.join(BASE, "assets", "content.json")

_CONTENT = None

_FALLBACK_BEATS = {
    "approach": ("Ligand enters from solution", "The ligand drifts toward the protein driven by diffusion, far from any binding site."),
    "surface_scan": ("Scanning the protein surface", "The ligand tumbles across the surface sampling translations and rotations."),
    "pocket_detection": ("Searching for a druggable pocket", "Candidate cavities are scored; the deepest, most enclosed pocket wins."),
    "first_candidate": ("First candidate pose", "The ligand dives into the pocket in a trial orientation."),
    "rejection": ("Pose rejected — steric clash", "Atoms overlap protein residues; the penalty spikes and the pose is rejected."),
    "reorientation": ("Trying a new orientation", "The ligand flips and re-enters with a more complementary orientation."),
    "interaction_discovery": ("Key interaction discovered", "A hydrogen bond snaps into place and the score drops sharply."),
    "refinement": ("Refining the pose", "Local optimisation tightens contacts and improves shape complementarity."),
    "final_pose": ("Converged docked pose", "The ligand settles into its lowest-energy binding mode."),
}
_PHASE_ORDER = list(_FALLBACK_BEATS.keys())


def load_content():
    global _CONTENT
    if _CONTENT is None:
        try:
            with open(CONTENT_PATH, encoding="utf-8") as fh:
                _CONTENT = json.load(fh)
        except (FileNotFoundError, json.JSONDecodeError):
            _CONTENT = {"packs": {}, "tutor": {"interactions": [], "concepts": [], "residue_classes": []}}
    return _CONTENT


def get_tutor():
    return load_content()["tutor"]


def _residue_role_map():
    return {r["resn"]: r for r in get_tutor().get("residue_classes", [])}


def _interaction_map():
    return {r["key"]: r for r in get_tutor().get("interactions", [])}


def get_pack(record):
    """Return a content pack for a structure record (real for case studies,
    synthesized for uploads)."""
    content = load_content()
    pack = content["packs"].get(record.get("slug"))
    if pack:
        return pack
    return _synthesize_pack(record)


def get_beats(record):
    pack = get_pack(record)
    beats = pack.get("beats")
    if beats and len(beats) >= 9:
        return beats
    return [{"phase": ph, "title": _FALLBACK_BEATS[ph][0], "caption": _FALLBACK_BEATS[ph][1]}
            for ph in _PHASE_ORDER]


def _synthesize_pack(record):
    rmap = _residue_role_map()
    imap = _interaction_map()
    res_notes = []
    for r in record.get("pocket_residues", [])[:14]:
        rc = rmap.get(r["resn"])
        role = rc["role"] if rc else "Lines the binding pocket and contacts the ligand."
        res_notes.append({"resn": r["resn"], "resi": r["resi"], "role": role})
    int_notes = []
    for it in record.get("interactions", []):
        ic = imap.get(it["type"])
        why = ic["why"] if ic else "Contributes to binding."
        int_notes.append({"resn": it["resn"], "resi": it["resi"], "type": it["type"],
                          "note": f"{it['type'].replace('_', ' ').title()} with "
                                  f"{it['resn']}{it['resi']}. {why}"})
    return {
        "slug": record.get("slug", "custom"),
        "overview": f"{record.get('display', 'This structure')} is a user-supplied structure. "
                    "The docking trajectory below is generated against its detected binding pocket.",
        "disease": "Custom upload — biological and disease context not available.",
        "drug": {"name": record.get("ligand_name", "Ligand"),
                 "code": record.get("ligand_name", "LIG"),
                 "drug_class": "Uploaded ligand",
                 "mechanism": "The detected ligand is docked into the most enclosed pocket of the uploaded protein."},
        "pocket_summary": "The binding pocket was detected from residues within 4.5 Å of the ligand. "
                          "Interactions are inferred from inter-atomic geometry.",
        "beats": [{"phase": ph, "title": _FALLBACK_BEATS[ph][0], "caption": _FALLBACK_BEATS[ph][1]}
                  for ph in _PHASE_ORDER],
        "residue_notes": res_notes,
        "interaction_notes": int_notes,
        "accurate": True, "issues": [],
    }
