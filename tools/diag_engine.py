"""Diagnostic: is poor RMSD a scoring problem or a search problem?"""
import numpy as np
from docking import structures, engine


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


for slug, disp, pdb, fam in structures.case_list():
    rec = structures.load_case_study(slug)
    R = engine._prep_receptor(rec)
    lig = engine._prep_ligand(rec)
    center = np.array(rec["pocket_center"], float)

    # score of the TRUE crystal pose
    cryst = lig["crystal"]
    e_cryst, comp_cryst, nhb_cryst = engine._score(cryst, lig, R)

    # run the search and score the found pose
    states, c2, n_explore = engine._search(rec, lig, R, 600, 7, exhaustiveness=4)
    Rm, t = states[-1]
    found = engine._pose(lig["local"], Rm, t, center)
    e_found, comp_found, nhb_found = engine._score(found, lig, R)

    print(f"== {disp} ({pdb})  lig_atoms={len(lig['elems'])} recept_atoms={len(R['atoms'])}")
    print(f"   crystal raw E = {e_cryst:8.3f}   nHB={nhb_cryst}")
    print(f"   found   raw E = {e_found:8.3f}   nHB={nhb_found}   RMSD={rmsd(found, cryst):.2f}")
    verdict = "SEARCH problem (crystal scores better but not found)" if e_cryst < e_found - 0.2 \
        else ("SCORING problem (found scores >= crystal but wrong geometry)" if rmsd(found, cryst) > 2.5
              else "OK")
    print(f"   verdict: {verdict}")
    # component breakdown
    keys = ["vdw", "clash", "hbond", "hydrophobic", "electrostatic", "desolvation", "entropy"]
    print("   crystal comps:", {k: round(comp_cryst[k], 2) for k in keys})
    print("   found   comps:", {k: round(comp_found[k], 2) for k in keys})
