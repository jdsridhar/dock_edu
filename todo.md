# TODO — Visual Molecular Docking Simulator: Review & Fixes

Status of the review items. **Most are now implemented and browser-verified**
(2026-06-12). Remaining items are lower-priority polish.

> ✅ done & verified · 🟡 partially addressed · ⬜ not started

---

## 🔴 High priority

- ✅ **1. Playback state reset on every rerun.** Fixed two ways: (a) all pure-render
  controls (representation, colour/property map, opacity, surface, surface opacity)
  were **moved into the viewer component**, so changing them no longer triggers a
  Streamlit rerun; (b) frame / play-state / speed are persisted to `sessionStorage`
  and restored on load. *Verified: jumped to frame 999, hit Regenerate (a real
  Streamlit rerun) → restored to 999 instead of resetting to 0.*
- ✅ **2. Unescaped `</script>` injection.** `build_workspace_html` now does
  `json.dumps(...).replace("</", "<\\/")`. *Verified: no raw `</` remains in the data
  payload.*
- ✅ **3. Surface freeze + uploads with no size cap.** Surface is now restricted to
  residues within ~23 Å of the pocket and **built asynchronously** with a "Building
  molecular surface…" note. Uploads are **trimmed to the ligand-bearing chain** and
  hard-capped to the binding-site region (>7000 atoms → within 32 Å), with a user
  warning. *Verified: COX-2 dimer upload trimmed to chain A; surface renders without
  freezing.*
- ✅ **4. No error guard around trajectory generation.** Generation + HTML build are
  wrapped in `try/except → st.error(...)`; the viewer guards empty pockets
  (`zoomTo({})` fallback).

## 🟠 Medium priority

- ✅ **5. "Surface" representation drew cartoon underneath.** `baseProteinStyle()` now
  uses an empty protein style for the `surface` representation (surface + pocket
  sticks only).
- 🟡 **6. Per-frame re-render cost.** Mitigated: renders are gated to integer-frame
  changes (so 0.1× speed renders ~3×/s, not 60×/s), surface is pocket-local, and the
  whole component is paused-idle when not playing. A deeper optimization (diff the
  active-shape set instead of full `removeAllShapes` each frame) is still possible if
  jank shows on low-end GPUs at 1200 frames.
- ✅ **7. Interaction lookup ignored chain.** `analyze_pocket` now carries the chain id
  and matches `(resn, resi, chain)`. *Verified on the COX-2 dimer upload.*
- ✅ **8. Synthetic values presented like measured ones.** Added a visible
  **`SIMULATED`** badge in the viewer header (tooltip explains the endpoint is a real
  crystallographic pose) and a Frame-Inspector note: "translation & rotation are pose
  deltas; torsions are illustrative."
- ✅ **9. "Electrostatic / hydrophobic surface" relabelled.** Colour-map options are now
  **Spectrum / Charge map / Hydropathy (Kyte–Doolittle) / Polarity map** for honesty
  (residue-property colouring, not a computed PB potential).

## 🟢 Low priority / housekeeping

- ✅ **10. `requirements.txt` cleanup.** Removed unused `pandas`; documented that
  `py3Dmol` is not needed (3Dmol.js is bundled) and `rdkit` is upload-only.
- ✅ **11. Repo hygiene.** Added `.gitignore` (`__pycache__/`, `*.log`,
  `assets/structures_raw/`, `tools/content_args.json`); removed the leftover
  `content_args.json` artifact.
- ✅ **12. `assets/structures_raw/` documented** as a build-only input (git-ignored).
- ✅ **13. Keyboard controls.** Space = play/pause, ←/→ step (Shift = jump),
  Home/End. *Verified: Space toggles play.*
- 🟡 **14. Responsive layout.** Added a `max-width:880px` breakpoint that stacks the
  panels under the stage. Functional, but not deeply tested on real mobile; the
  iframe height is still fixed by Streamlit.
- ✅ **15. Always-autoplay.** Now auto-plays only on the **first** visit to a target;
  on later renders it restores the saved play state instead of force-restarting.
- ✅ **16. Broader verification.** Browser-verified: Display controls (representation /
  colour / opacity / surface) apply live with no rerun, surface + hydropathy map
  render, AI tutor, keyboard, and rerun-persistence all work with **zero console
  errors**. All 4 case studies + a multi-chain upload build successfully.

---

## Remaining (optional)
- 🟡 #6 — shape-diff optimization for very high frame counts on weak GPUs.
- 🟡 #14 — true mobile testing / dynamic component height.
- ⬜ Real computed electrostatic potential surface (would need a PB/Coulomb pass) if
  "charge map" should become a genuine potential map.
