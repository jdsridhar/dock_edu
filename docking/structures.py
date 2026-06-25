"""
Structure loading + binding-site analysis.

Provides a single canonical record consumed by the trajectory engine and viewer,
built either from a bundled case study (real RCSB crystal structure) or from a
user upload (PDB / PDBQT protein, optional MOL / SDF / MOL2 / PDB ligand).

No hard dependency on Biopython / OpenBabel — PDB/PDBQT are parsed directly and
RDKit is used opportunistically (with a manual fallback) for MOL/SDF/MOL2.
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT_DIR = os.path.join(BASE, "assets", "structures")

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
STD_AA = set(THREE_TO_ONE)

WATERS = {"HOH", "WAT", "DOD", "H2O"}
IONS = {"NA", "CL", "K", "MG", "ZN", "CA", "MN", "FE", "CU", "NI", "CD", "HG",
        "CO", "SO4", "PO4", "NO3", "BR", "IOD", "F", "CS", "RB", "SR", "BA"}
BUFFERS = {"GOL", "EDO", "PEG", "PG4", "1PE", "2PE", "MPD", "DMS", "TRS", "EPE",
           "MES", "ACT", "ACY", "FMT", "BME", "IMD", "CIT", "TLA", "MLI", "BOG",
           "LDA", "OLC", "OLA", "PLM", "MYR", "SOG", "DIO", "BU3", "P6G", "PE4",
           "12P", "15P", "PGE", "PEU", "MRD", "BEN", "FLC", "CAC", "POP", "SIN"}
SUGARS = {"NAG", "BMA", "MAN", "GAL", "FUC", "GLC", "NDG", "BGC", "XYS", "SIA",
          "FUL", "GLA", "A2G", "GCS"}
COFACTORS = {"HEM", "HEC", "HEA", "HEB", "FAD", "FMN", "NAD", "NAP", "NDP",
             "NAI", "ADP", "ATP", "GDP", "GTP", "AMP", "COA", "PLP", "SAM",
             "SAH", "TPP", "BTN", "B12", "PQQ", "MTE", "MGD", "F3S", "SF4",
             "FES", "PHO", "CLA", "BCL"}
EXCLUDE_HET = WATERS | IONS | BUFFERS | SUGARS | COFACTORS

AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}
POS_RES = {"ARG", "LYS", "HIS"}
NEG_RES = {"ASP", "GLU"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY", "CYS"}

# AutoDock (PDBQT) atom-type -> element
AUTODOCK_ELEM = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O",
                 "OS": "O", "O": "O", "SA": "S", "S": "S", "HD": "H", "HS": "H",
                 "H": "H", "P": "P", "F": "F", "CL": "CL", "BR": "BR", "I": "I",
                 "MG": "MG", "ZN": "ZN", "CA": "CA", "FE": "FE", "MN": "MN"}


# ---------------------------------------------------------------------------
# low-level parsing
# ---------------------------------------------------------------------------
def _elem_from_name(name: str) -> str:
    s = "".join(c for c in name if c.isalpha())
    if not s:
        return "C"
    two = s[:2].upper()
    if two in ("CL", "BR", "NA", "MG", "ZN", "FE", "MN", "CA", "CU", "NI", "CO", "SE"):
        return two
    return s[0].upper()


def parse_pdb_atoms(text: str, pdbqt: bool = False):
    """Parse ATOM/HETATM records from PDB or PDBQT text (first model only)."""
    atoms = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        name = line[12:16].strip()
        altloc = line[16:17].strip()
        if altloc not in ("", "A"):
            continue
        elem = line[76:78].strip().upper()
        if pdbqt:
            elem = AUTODOCK_ELEM.get(elem, _elem_from_name(name))
        elif not elem:
            elem = _elem_from_name(name)
        atoms.append({
            "rec": line[0:6].strip(), "name": name, "resn": line[17:20].strip(),
            "chain": line[21:22].strip(), "resi": line[22:26].strip(),
            "x": x, "y": y, "z": z, "elem": elem,
        })
    return atoms


def parse_ligand(text: str, fmt: str):
    """Parse a small-molecule ligand from various formats -> list of atoms."""
    fmt = fmt.lower().lstrip(".")
    if fmt in ("pdb", "ent"):
        return [a for a in parse_pdb_atoms(text) if a["elem"] != "H"]
    if fmt == "pdbqt":
        return [a for a in parse_pdb_atoms(text, pdbqt=True) if a["elem"] != "H"]
    if fmt in ("mol", "sdf", "mol2"):
        try:
            from rdkit import Chem
            if fmt == "mol2":
                mol = Chem.MolFromMol2Block(text, sanitize=False, removeHs=True)
            elif fmt == "sdf":
                supplier = Chem.SDMolSupplier()
                supplier.SetData(text, sanitize=False, removeHs=True)
                mol = next((m for m in supplier if m is not None), None)
            else:
                mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=True)
            if mol is not None and mol.GetNumConformers():
                conf = mol.GetConformer()
                out = []
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "H":
                        continue
                    p = conf.GetAtomPosition(atom.GetIdx())
                    out.append({"rec": "HETATM", "name": atom.GetSymbol() + str(atom.GetIdx()),
                                "resn": "LIG", "chain": "X", "resi": "1",
                                "x": p.x, "y": p.y, "z": p.z, "elem": atom.GetSymbol().upper()})
                if out:
                    return out
        except Exception:
            pass
        return _parse_molblock_basic(text) if fmt in ("mol", "sdf") else []
    raise ValueError(f"Unsupported ligand format: {fmt}")


def _parse_molblock_basic(text: str):
    lines = text.splitlines()
    if len(lines) < 4:
        return []
    try:
        counts = lines[3]
        natoms = int(counts[0:3])
    except (ValueError, IndexError):
        return []
    out = []
    for i in range(natoms):
        ln = lines[4 + i]
        try:
            x = float(ln[0:10]); y = float(ln[10:20]); z = float(ln[20:30])
            elem = ln[31:34].strip().upper()
        except (ValueError, IndexError):
            continue
        if elem == "H":
            continue
        out.append({"rec": "HETATM", "name": elem + str(i), "resn": "LIG",
                    "chain": "X", "resi": "1", "x": x, "y": y, "z": z, "elem": elem})
    return out


# ---------------------------------------------------------------------------
# geometry / detection
# ---------------------------------------------------------------------------
def _dist(a, b):
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def _centroid(atoms):
    n = len(atoms)
    return [sum(a["x"] for a in atoms) / n, sum(a["y"] for a in atoms) / n,
            sum(a["z"] for a in atoms) / n]


def _extent(atoms):
    xs = [a["x"] for a in atoms]; ys = [a["y"] for a in atoms]; zs = [a["z"] for a in atoms]
    return {"min": [min(xs), min(ys), min(zs)], "max": [max(xs), max(ys), max(zs)],
            "center": [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]}


def detect_ligand(het_atoms, prot_atoms):
    """Pick the most drug-like HETATM residue in contact with the protein."""
    groups = defaultdict(list)
    for a in het_atoms:
        if a["resn"] in EXCLUDE_HET or a["elem"] == "H":
            continue
        groups[(a["chain"], a["resn"], a["resi"])].append(a)
    best, best_key, best_score = None, None, -1
    for key, atoms in groups.items():
        n = len(atoms)
        if n < 6 or n > 120:
            continue
        contacts = sum(1 for la in atoms
                       if any(_dist(la, pa) < 4.5 for pa in prot_atoms
                              if abs(la["x"] - pa["x"]) < 4.5 and abs(la["y"] - pa["y"]) < 4.5))
        score = contacts * 1000 + n
        if contacts and score > best_score:
            best, best_key, best_score = atoms, key, score
    return (best, best_key[1]) if best else (None, None)


def analyze_pocket(prot_atoms, lig_atoms, cutoff=4.5):
    """Return pocket residues + interaction geometry between protein and ligand."""
    res_atoms = defaultdict(list)
    for a in prot_atoms:
        if a["resn"] in STD_AA:
            res_atoms[(a["resn"], a["resi"], a["chain"])].append(a)

    pocket = []
    for key, ratoms in res_atoms.items():
        mind = min((_dist(la, ra) for la in lig_atoms for ra in ratoms), default=999)
        if mind < cutoff:
            ca = next((a for a in ratoms if a["name"] == "CA"), ratoms[0])
            pocket.append({"resn": key[0], "resi": key[1], "chain": key[2], "mind": round(mind, 2),
                           "one": THREE_TO_ONE.get(key[0], "X"),
                           "ca": [round(ca["x"], 3), round(ca["y"], 3), round(ca["z"], 3)]})
    pocket.sort(key=lambda r: r["mind"])

    interactions = []
    for r in pocket[:20]:
        resn, resi = r["resn"], r["resi"]
        ratoms = res_atoms[(resn, resi, r["chain"])]  # match exact chain (multimer-safe)
        hb = None
        for la in lig_atoms:
            if la["elem"] not in ("N", "O"):
                continue
            for ra in ratoms:
                if ra["elem"] not in ("N", "O") or ra["name"] in ("C", "CA"):
                    continue
                d = _dist(la, ra)
                if d < 3.6 and (hb is None or d < hb["dist"]):
                    hb = {"dist": round(d, 2), "lig": [la["x"], la["y"], la["z"]],
                          "prot": [ra["x"], ra["y"], ra["z"]]}
        if hb:
            interactions.append({"type": "hydrogen_bond", "resn": resn, "resi": resi, **hb})
            continue
        if resn in AROMATIC_RES:
            ring = [a for a in ratoms if a["name"] in
                    ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "NE1", "CH2", "CZ2", "CZ3", "CE3", "ND1", "NE2")]
            if ring:
                rc = {"x": _centroid(ring)[0], "y": _centroid(ring)[1], "z": _centroid(ring)[2]}
                mind = min(_dist(la, rc) for la in lig_atoms)
                if mind < 5.5:
                    nearest = min(lig_atoms, key=lambda la: _dist(la, rc))
                    interactions.append({"type": "pi_stacking", "resn": resn, "resi": resi,
                                         "dist": round(mind, 2),
                                         "lig": [nearest["x"], nearest["y"], nearest["z"]],
                                         "prot": [rc["x"], rc["y"], rc["z"]]})
                    continue
        if resn in POS_RES or resn in NEG_RES:
            charged = [a for a in ratoms if a["elem"] in ("N", "O") and a["name"] not in ("N", "O", "CA", "C")]
            best = None
            for la in lig_atoms:
                if la["elem"] not in ("N", "O"):
                    continue
                for ca in charged:
                    d = _dist(la, ca)
                    if d < 4.5 and (best is None or d < best["dist"]):
                        best = {"dist": round(d, 2), "lig": [la["x"], la["y"], la["z"]],
                                "prot": [ca["x"], ca["y"], ca["z"]]}
            if best:
                interactions.append({"type": "salt_bridge", "resn": resn, "resi": resi, **best})
                continue
        if resn in HYDROPHOBIC_RES and r["mind"] < 4.3:
            ra = min(ratoms, key=lambda ra: min(_dist(la, ra) for la in lig_atoms))
            la = min(lig_atoms, key=lambda la: _dist(la, ra))
            interactions.append({"type": "hydrophobic", "resn": resn, "resi": resi,
                                 "dist": round(_dist(la, ra), 2),
                                 "lig": [la["x"], la["y"], la["z"]],
                                 "prot": [ra["x"], ra["y"], ra["z"]]})
    return pocket, interactions


def _protein_pdb_text(prot_atoms):
    """Serialize protein atoms to clean PDB ATOM records for 3Dmol."""
    out = []
    for i, a in enumerate(prot_atoms, 1):
        name = a["name"]
        nm = (" " + name) if (len(name) < 4 and not name[:1].isdigit()) else name
        nm = nm[:4].ljust(4)
        chain = (a["chain"] or "A")[:1]
        out.append(f"ATOM  {i:>5} {nm}{'':1}{a['resn']:>3} {chain}{a['resi']:>4}"
                   f"{'':1}   {a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}"
                   f"{1.0:>6.2f}{0.0:>6.2f}          {a['elem']:>2}")
    out.append("END")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# case studies
# ---------------------------------------------------------------------------
_META_CACHE = None


def case_meta():
    global _META_CACHE
    if _META_CACHE is None:
        with open(os.path.join(STRUCT_DIR, "meta.json"), encoding="utf-8") as fh:
            _META_CACHE = json.load(fh)
    return _META_CACHE


def case_list():
    return [(slug, m["display"], m["source_pdb"], m["family"]) for slug, m in case_meta().items()]


def load_case_study(slug):
    meta = case_meta()[slug]
    with open(os.path.join(STRUCT_DIR, slug + ".pdb"), encoding="utf-8") as fh:
        raw = fh.read()
    protein_pdb = "\n".join(ln for ln in raw.splitlines()
                            if ln.startswith(("ATOM", "TER"))) + "\nEND"
    lig = meta["ligand_atoms"]
    return {
        "slug": slug, "display": meta["display"], "family": meta["family"],
        "source_pdb": meta["source_pdb"], "is_case": True,
        "protein_pdb": protein_pdb,
        "ligand_names": [a["name"] for a in lig],
        "ligand_elems": [a["elem"] for a in lig],
        "ligand_coords": [a["xyz"] for a in lig],
        "ligand_name": meta["ligand_name"],
        "pocket_center": meta["pocket_center"],
        "pocket_residues": meta["pocket_residues"],
        "interactions": meta["interactions"],
        "protein_extent": meta["protein_extent"],
    }


# ---------------------------------------------------------------------------
# uploads
# ---------------------------------------------------------------------------
def build_from_upload(protein_text, protein_fmt, ligand_text=None, ligand_fmt=None, name="Custom"):
    pdbqt = protein_fmt.lower().lstrip(".") == "pdbqt"
    all_atoms = parse_pdb_atoms(protein_text, pdbqt=pdbqt)
    prot_atoms = [a for a in all_atoms if a["resn"] in STD_AA and a["elem"] != "H"]
    if not prot_atoms:
        raise ValueError("No protein atoms (standard amino acids) found in the uploaded structure.")

    # Resolve a ligand to dock: an uploaded one first, otherwise a real ligand bound
    # in the structure. We do NOT fabricate one — docking needs a real molecule.
    lig_atoms, lig_source = None, None
    if ligand_text:
        lig_atoms = parse_ligand(ligand_text, ligand_fmt or "pdb")
        if lig_atoms:
            lig_source = "uploaded"
        else:
            return {"needs_ligand": True, "reason": "ligand_parse_failed", "is_case": False,
                    "display": name, "family": "uploaded structure", "source_pdb": "upload"}
    if not lig_atoms:
        het = [a for a in all_atoms if a["rec"] == "HETATM"]
        lig_atoms, _ = detect_ligand(het, prot_atoms)
        if lig_atoms:
            lig_source = "detected"
    if not lig_atoms:
        # No ligand uploaded and none bound in the structure — ask for one rather
        # than inventing a fake probe molecule and "docking" it.
        het_groups = {(a["resn"]) for a in all_atoms if a["rec"] == "HETATM"
                      and a["resn"] not in WATERS}
        return {"needs_ligand": True, "reason": "no_ligand", "is_case": False,
                "display": name, "family": "uploaded structure", "source_pdb": "upload",
                "het_seen": sorted(het_groups)}

    # If the ligand coordinates do not overlap the protein, dock it into the pocket.
    lc = _centroid(lig_atoms)
    near = any(_dist({"x": lc[0], "y": lc[1], "z": lc[2]}, pa) < 12 for pa in prot_atoms[::10])
    if not near:
        target = _estimate_pocket_center(prot_atoms)
        dx, dy, dz = target[0] - lc[0], target[1] - lc[1], target[2] - lc[2]
        for a in lig_atoms:
            a["x"] += dx; a["y"] += dy; a["z"] += dz

    # Trim to the ligand-bearing chain (uploads can be large multimers).
    warning = None
    chain_contacts = defaultdict(int)
    for la in lig_atoms:
        for pa in prot_atoms:
            if _dist(la, pa) < 6.0:
                chain_contacts[pa["chain"]] += 1
    if chain_contacts:
        keep = max(chain_contacts, key=chain_contacts.get)
        trimmed = [a for a in prot_atoms if a["chain"] == keep]
        if trimmed:
            if len(trimmed) < len(prot_atoms):
                warning = f"Trimmed to chain {keep or '?'} ({len(trimmed)} atoms) around the ligand."
            prot_atoms = trimmed

    # Hard size cap: keep only the binding-site region for very large structures.
    lig_cen = _centroid(lig_atoms)
    if len(prot_atoms) > 7000:
        c = {"x": lig_cen[0], "y": lig_cen[1], "z": lig_cen[2]}
        prot_atoms = [a for a in prot_atoms if _dist(a, c) < 32]
        warning = f"Large structure trimmed to the binding-site region ({len(prot_atoms)} atoms) for performance."

    pocket, interactions = analyze_pocket(prot_atoms, lig_atoms)
    return {
        "slug": "custom", "display": name, "family": "uploaded structure",
        "source_pdb": "upload", "is_case": False, "warning": warning,
        "ligand_source": lig_source,
        "protein_pdb": _protein_pdb_text(prot_atoms),
        "ligand_names": [a["name"] for a in lig_atoms],
        "ligand_elems": [a["elem"] for a in lig_atoms],
        "ligand_coords": [[round(a["x"], 3), round(a["y"], 3), round(a["z"], 3)] for a in lig_atoms],
        "ligand_name": lig_atoms[0]["resn"] if lig_atoms else "LIG",
        "pocket_center": [round(c, 3) for c in lig_cen],
        "pocket_residues": pocket,
        "interactions": interactions,
        "protein_extent": _extent(prot_atoms),
    }


def _estimate_pocket_center(prot_atoms):
    """Lightweight cavity estimate: the buried CA neighbourhood centroid."""
    cas = [a for a in prot_atoms if a["name"] == "CA"] or prot_atoms
    center = _centroid(prot_atoms)
    cc = {"x": center[0], "y": center[1], "z": center[2]}
    scored = []
    sample = cas[:: max(1, len(cas) // 220)]
    for a in sample:
        nb = sum(1 for b in sample if _dist(a, b) < 10)
        scored.append((nb / (1 + _dist(a, cc)), a))
    scored.sort(reverse=True, key=lambda t: t[0])
    top = [a for _, a in scored[:8]] or cas[:1]
    return _centroid(top)
