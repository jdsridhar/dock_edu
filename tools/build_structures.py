"""
Build-time utility: parse raw RCSB PDB files, auto-detect the drug-like ligand,
extract a single protein chain + the ligand, and compute the real binding-pocket
residues and interaction geometry.  Outputs trimmed PDB files plus a metadata
JSON consumed by the application at runtime.

Run once:  python tools/build_structures.py
"""
import json
import math
import os
from collections import defaultdict, OrderedDict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, "assets", "structures_raw")
OUT = os.path.join(BASE, "assets", "structures")

# (raw pdb id, output slug, display name, family)
TARGETS = [
    ("1ACJ", "acetylcholinesterase", "Acetylcholinesterase", "hydrolase"),
    ("1M17", "egfr", "EGFR Kinase", "kinase"),
    ("1FKN", "bace1", "BACE1 (β-Secretase)", "aspartic protease"),
    ("1CX2", "cox2", "COX-2 (Cyclooxygenase-2)", "oxidoreductase"),
]

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
EXCLUDE = WATERS | IONS | BUFFERS | SUGARS | COFACTORS

ELEMENT_MASS = {"H": 1, "C": 12, "N": 14, "O": 16, "P": 31, "S": 32,
                "F": 19, "CL": 35, "BR": 80, "I": 127}
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}
HBOND_DON_ACC = {"SER", "THR", "TYR", "ASN", "GLN", "HIS", "TRP", "ARG", "LYS",
                 "ASP", "GLU", "CYS"}
POS_RES = {"ARG", "LYS", "HIS"}
NEG_RES = {"ASP", "GLU"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO",
                   "GLY", "CYS"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
STD_AA = set(THREE_TO_ONE.keys())


def parse_atom_line(line):
    rec = line[0:6].strip()
    try:
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
    except ValueError:
        return None
    elem = line[76:78].strip().upper()
    name = line[12:16].strip()
    if not elem:  # derive element from atom name
        elem = "".join(c for c in name if c.isalpha())[:2].upper()
        if len(elem) == 2 and elem[1].islower():
            elem = elem[0]
    altloc = line[16:17].strip()
    return {
        "rec": rec,
        "serial": line[6:11].strip(),
        "name": name,
        "altloc": altloc,
        "resn": line[17:20].strip(),
        "chain": line[21:22].strip(),
        "resi": line[22:26].strip(),
        "icode": line[26:27].strip(),
        "x": x, "y": y, "z": z,
        "elem": elem,
        "line": line.rstrip("\n"),
    }


def dist(a, b):
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def centroid(atoms):
    n = len(atoms)
    return [sum(a["x"] for a in atoms) / n,
            sum(a["y"] for a in atoms) / n,
            sum(a["z"] for a in atoms) / n]


def load_first_model(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")):
                a = parse_atom_line(line)
                if a and a["altloc"] in ("", "A"):
                    atoms.append(a)
    return atoms


def build_one(pdb_id, slug, display, family):
    path = os.path.join(RAW, pdb_id + ".pdb")
    atoms = load_first_model(path)
    prot = [a for a in atoms if a["rec"] == "ATOM" and a["resn"] in STD_AA and a["elem"] != "H"]
    het = [a for a in atoms if a["rec"] == "HETATM" and a["elem"] != "H"]

    # group het by (chain, resn, resi)
    groups = defaultdict(list)
    for a in het:
        if a["resn"] in EXCLUDE:
            continue
        groups[(a["chain"], a["resn"], a["resi"])].append(a)

    # score candidate ligands: heavy-atom count in range, in contact with protein
    candidates = []
    for key, gatoms in groups.items():
        n = len(gatoms)
        if n < 8 or n > 90:
            continue
        contacts = 0
        for la in gatoms:
            for pa in prot:
                if abs(la["x"] - pa["x"]) < 4.5 and abs(la["y"] - pa["y"]) < 4.5:
                    if dist(la, pa) < 4.5:
                        contacts += 1
                        break
        candidates.append((contacts, n, key, gatoms))
    candidates.sort(key=lambda t: (t[0], t[1]), reverse=True)

    if not candidates:
        raise RuntimeError(f"No ligand candidate found for {pdb_id}")
    print(f"\n=== {pdb_id} ({display}) candidates ===")
    for contacts, n, key, _ in candidates[:6]:
        print(f"   {key[1]:>4} chain {key[0]} resi {key[2]:>5}  heavy={n:>3}  contacts={contacts}")

    _, lign, ligkey, ligatoms = candidates[0]
    lig_chain = ligkey[0]
    ligname = ligkey[1]

    # choose protein chain that contacts the ligand most
    chain_contacts = defaultdict(int)
    for la in ligatoms:
        for pa in prot:
            if dist(la, pa) < 6.0:
                chain_contacts[pa["chain"]] += 1
    prot_chain = max(chain_contacts, key=chain_contacts.get)
    prot_chain_atoms = [a for a in prot if a["chain"] == prot_chain]

    # pocket residues: protein residues within 4.5 A of any ligand atom
    res_atoms = defaultdict(list)
    for a in prot_chain_atoms:
        res_atoms[(a["resn"], a["resi"])].append(a)

    pocket = []
    for (resn, resi), ratoms in res_atoms.items():
        mind = min(dist(la, ra) for la in ligatoms for ra in ratoms)
        if mind < 4.5:
            pocket.append({"resn": resn, "resi": resi, "mind": round(mind, 2),
                           "ca": next((a for a in ratoms if a["name"] == "CA"), ratoms[0])})
    pocket.sort(key=lambda r: r["mind"])

    lig_cen = centroid(ligatoms)

    # interaction geometry between ligand and pocket
    interactions = []
    for r in pocket:
        resn, resi = r["resn"], r["resi"]
        ratoms = res_atoms[(resn, resi)]
        # hydrogen bonds: polar protein atom (N/O) near polar ligand atom (N/O)
        best_hb = None
        for la in ligatoms:
            if la["elem"] not in ("N", "O"):
                continue
            for ra in ratoms:
                if ra["elem"] not in ("N", "O"):
                    continue
                if ra["name"] in ("C", "CA"):
                    continue
                d = dist(la, ra)
                if d < 3.5 and (best_hb is None or d < best_hb["dist"]):
                    best_hb = {"dist": round(d, 2), "lig": [la["x"], la["y"], la["z"]],
                               "prot": [ra["x"], ra["y"], ra["z"]],
                               "prot_atom": ra["name"], "lig_atom": la["name"]}
        if best_hb:
            interactions.append({"type": "hydrogen_bond", "resn": resn, "resi": resi, **best_hb})
            continue
        # pi-pi / pi-cation for aromatic residues
        if resn in AROMATIC_RES:
            ring = [a for a in ratoms if a["name"] in
                    ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "NE1", "CH2", "CZ2", "CZ3", "CE3", "ND1", "NE2")]
            if ring:
                rc = centroid(ring)
                rc_atom = {"x": rc[0], "y": rc[1], "z": rc[2]}
                mind = min(dist(la, rc_atom) for la in ligatoms)
                if mind < 5.5:
                    nearest = min(ligatoms, key=lambda la: dist(la, rc_atom))
                    interactions.append({"type": "pi_stacking", "resn": resn, "resi": resi,
                                         "dist": round(mind, 2),
                                         "lig": [nearest["x"], nearest["y"], nearest["z"]],
                                         "prot": [round(c, 3) for c in rc]})
                    continue
        # salt bridges
        if resn in POS_RES or resn in NEG_RES:
            charged = [a for a in ratoms if a["elem"] in ("N", "O") and a["name"] not in ("N", "O", "CA", "C")]
            best = None
            for la in ligatoms:
                if la["elem"] not in ("N", "O"):
                    continue
                for ca in charged:
                    d = dist(la, ca)
                    if d < 4.5 and (best is None or d < best["dist"]):
                        best = {"dist": round(d, 2), "lig": [la["x"], la["y"], la["z"]],
                                "prot": [ca["x"], ca["y"], ca["z"]]}
            if best:
                interactions.append({"type": "salt_bridge", "resn": resn, "resi": resi, **best})
                continue
        # hydrophobic contact
        if resn in HYDROPHOBIC_RES and r["mind"] < 4.2:
            ra = min(ratoms, key=lambda ra: min(dist(la, ra) for la in ligatoms))
            la = min(ligatoms, key=lambda la: dist(la, ra))
            interactions.append({"type": "hydrophobic", "resn": resn, "resi": resi,
                                 "dist": round(dist(la, ra), 2),
                                 "lig": [la["x"], la["y"], la["z"]],
                                 "prot": [ra["x"], ra["y"], ra["z"]]})

    # write trimmed PDB: one protein chain + ligand
    out_lines = ["REMARK  Trimmed by Visual Molecular Docking Simulator",
                 f"REMARK  source={pdb_id} protein_chain={prot_chain} ligand={ligname}"]
    serial = 1
    for a in prot_chain_atoms:
        out_lines.append(_fmt_atom(a, serial, "ATOM"))
        serial += 1
    out_lines.append("TER")
    for a in ligatoms:
        out_lines.append(_fmt_atom(a, serial, "HETATM"))
        serial += 1
    out_lines.append("END")
    with open(os.path.join(OUT, slug + ".pdb"), "w") as fh:
        fh.write("\n".join(out_lines) + "\n")

    # ligand element/coord list (heavy atoms) for the trajectory engine
    lig_out = [{"name": a["name"], "elem": a["elem"],
                "xyz": [round(a["x"], 3), round(a["y"], 3), round(a["z"], 3)]}
               for a in ligatoms]

    # protein extent for camera framing
    xs = [a["x"] for a in prot_chain_atoms]
    ys = [a["y"] for a in prot_chain_atoms]
    zs = [a["z"] for a in prot_chain_atoms]
    extent = {"min": [min(xs), min(ys), min(zs)], "max": [max(xs), max(ys), max(zs)],
              "center": [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]}

    meta = {
        "slug": slug, "display": display, "family": family,
        "source_pdb": pdb_id, "protein_chain": prot_chain,
        "ligand_name": ligname, "ligand_natoms": lign,
        "ligand_centroid": [round(c, 3) for c in lig_cen],
        "ligand_atoms": lig_out,
        "pocket_center": [round(c, 3) for c in lig_cen],
        "pocket_residues": [{"resn": r["resn"], "resi": r["resi"], "mind": r["mind"],
                             "ca": [round(r["ca"]["x"], 3), round(r["ca"]["y"], 3), round(r["ca"]["z"], 3)],
                             "one": THREE_TO_ONE.get(r["resn"], "X")} for r in pocket],
        "interactions": interactions,
        "protein_extent": extent,
        "n_protein_atoms": len(prot_chain_atoms),
    }
    print(f"   -> chain {prot_chain}, {len(prot_chain_atoms)} atoms, ligand {ligname} "
          f"({lign} heavy), pocket {len(pocket)} res, {len(interactions)} interactions")
    return meta


def _fmt_atom(a, serial, rec):
    name = a["name"]
    if len(name) < 4 and not name[:1].isdigit():
        name = " " + name
    name = name[:4].ljust(4)
    return (f"{rec:<6}{serial:>5} {name}{'':1}{a['resn']:>3} {a['chain']:1}"
            f"{a['resi']:>4}{'':1}   {a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}          {a['elem']:>2}")


def main():
    meta_all = OrderedDict()
    for pdb_id, slug, display, family in TARGETS:
        meta_all[slug] = build_one(pdb_id, slug, display, family)
    with open(os.path.join(OUT, "meta.json"), "w") as fh:
        json.dump(meta_all, fh, indent=2)
    print("\nWrote", os.path.join(OUT, "meta.json"))
    print("Targets:", ", ".join(meta_all.keys()))


if __name__ == "__main__":
    main()
