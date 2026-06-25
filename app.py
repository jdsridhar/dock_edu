"""
Visual Molecular Docking Simulator — "Watch Molecular Recognition Happen".

A molecular-dynamics-style trajectory player for docking: press PLAY and watch a
ligand search, get rejected, reorient, form interactions, and converge on a real
crystallographic binding pose — with synchronized scoring, shape-fit, occupancy,
a score-evolution timeline, a frame inspector, cinematic event overlays, and an
AI tutor, all in one workspace.
"""
from __future__ import annotations

import hashlib

import streamlit as st
import streamlit.components.v1 as components

from docking import structures, trajectory, knowledge, viewer

st.set_page_config(page_title="Visual Molecular Docking Simulator",
                   page_icon="🧬", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
  .block-container{padding-top:1.1rem;padding-bottom:1rem;max-width:1500px;}
  #MainMenu,footer{visibility:hidden;}
  .vmds-title{font-size:1.7rem;font-weight:800;letter-spacing:-.01em;margin:0;}
  .vmds-sub{color:#7c8db0;font-size:.93rem;margin:.1rem 0 .4rem;}
  .stApp{background:#05070d;}
  section[data-testid="stSidebar"]{background:#0a0f1c;}
  div[data-testid="stExpander"]{border:1px solid #1b2740;border-radius:10px;background:#0a1120;}
</style>
""", unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def _case_record(slug):
    return structures.load_case_study(slug)


@st.cache_data(show_spinner=False)
def _upload_record(protein_text, protein_fmt, ligand_text, ligand_fmt, name, _h):
    return structures.build_from_upload(protein_text, protein_fmt, ligand_text, ligand_fmt, name)


@st.cache_data(show_spinner=False)
def _trajectory(record, n_frames, seed):
    beats = knowledge.get_beats(record)
    return trajectory.generate_trajectory(record, n_frames=n_frames, seed=seed, beats=beats)


# ---------------------------------------------------------------- sidebar
sb = st.sidebar
sb.markdown("### 🧬 Docking Simulator")
sb.caption("Watch molecular recognition happen, frame by frame.")

mode = sb.radio("Structure source", ["Case study", "Upload your own"], horizontal=True)

record, err = None, None
if mode == "Case study":
    cases = structures.case_list()
    labels = [f"{disp}  ·  {pdb}" for slug, disp, pdb, fam in cases]
    idx = sb.selectbox("Target protein", range(len(cases)), format_func=lambda i: labels[i])
    slug = cases[idx][0]
    try:
        record = _case_record(slug)
    except Exception as e:  # noqa: BLE001
        err = f"Could not load case study: {e}"
else:
    sb.caption("Protein: PDB / PDBQT.  Ligand (optional): MOL / SDF / MOL2 / PDB / PDBQT.")
    pfile = sb.file_uploader("Protein structure", type=["pdb", "ent", "pdbqt"])
    lfile = sb.file_uploader("Ligand (optional)", type=["mol", "sdf", "mol2", "pdb", "pdbqt"])
    if pfile is not None:
        try:
            ptext = pfile.getvalue().decode("utf-8", "ignore")
            pfmt = pfile.name.rsplit(".", 1)[-1]
            ltext = lfile.getvalue().decode("utf-8", "ignore") if lfile else None
            lfmt = lfile.name.rsplit(".", 1)[-1] if lfile else None
            h = hashlib.md5((ptext + (ltext or "")).encode()).hexdigest()
            record = _upload_record(ptext, pfmt, ltext, lfmt, pfile.name.rsplit(".", 1)[0], h)
        except Exception as e:  # noqa: BLE001
            err = f"Could not parse upload: {e}"
    else:
        sb.info("Upload a protein structure to generate a docking trajectory.")

sb.markdown("---")
sb.markdown("##### Rendering")
sb.caption("🎛️ Display controls — **representation, colour/property map, opacity and "
           "molecular surface** — live inside the viewer (top of the right-hand panel). "
           "Changing them there will **not** restart the movie.")

sb.markdown("---")
sb.markdown("##### Trajectory")
n_frames = sb.select_slider("Frames", options=[400, 600, 800, 1000, 1200], value=1000)
if "seed" not in st.session_state:
    st.session_state.seed = 7
if sb.button("🔄 Regenerate search path"):
    st.session_state.seed += 1
sb.caption(f"Random search seed: {st.session_state.seed}")

# ---------------------------------------------------------------- header
st.markdown('<p class="vmds-title">Visual Molecular Docking Simulator</p>', unsafe_allow_html=True)
st.markdown('<p class="vmds-sub">Watch Molecular Recognition Happen — a docking trajectory player '
            'that shows <i>how</i> poses are searched, rejected, and accepted.</p>',
            unsafe_allow_html=True)

if err:
    st.error(err)
    st.stop()
if record is None:
    st.info("👈 Pick a built-in case study or upload a structure to begin. "
            "Then press **▶ Play** and slow the speed to **0.1x** to study the search.")
    st.stop()

# no ligand to dock — ask for one instead of inventing a fake molecule
if record.get("needs_ligand"):
    if record.get("reason") == "ligand_parse_failed":
        st.error("Couldn't parse the uploaded ligand file. Supported formats: "
                 "**MOL / SDF / MOL2 / PDB / PDBQT**.")
    else:
        het = record.get("het_seen") or []
        extra = (" The only non-protein groups in this structure are "
                 + ", ".join(het[:8])
                 + " — treated as solvent / ions / cofactors, not drug-like ligands.") if het else ""
        st.warning("**No ligand to dock.** You uploaded a protein but no ligand, and no "
                   "drug-like ligand is bound in this structure." + extra)
        st.info("Docking needs a ligand. Upload one (**MOL / SDF / MOL2 / PDB / PDBQT**) in the "
                "sidebar, or choose a structure that contains a bound inhibitor — e.g. one of the "
                "built-in case studies, or a protein–ligand complex from the RCSB PDB.")
    st.stop()

# ---------------------------------------------------------------- build + render
if record.get("warning"):
    st.warning("⚠️ " + record["warning"])
if record.get("ligand_source") == "detected":
    st.info(f"No ligand uploaded — docking **{record['ligand_name']}** "
            f"({len(record['ligand_coords'])} atoms), the ligand found bound in the uploaded structure.")

try:
    with st.spinner("Generating docking trajectory…"):
        traj = _trajectory(record, n_frames, st.session_state.seed)
        pack = knowledge.get_pack(record)
        tutor = knowledge.get_tutor()
    # initial display defaults; the viewer owns these controls from here on
    options = {"representation": "cartoon", "color_scheme": "spectrum",
               "opacity": 0.85, "show_surface": False, "surface_opacity": 0.55}
    html = viewer.build_workspace_html(record, traj, pack, tutor, options, height=946)
    components.html(html, height=956, scrolling=False)
except Exception as exc:  # noqa: BLE001
    st.error(f"Could not generate the docking trajectory for this structure: {exc}")
    st.caption("If this is a custom upload, check that the protein and ligand parsed correctly "
               "(a clearly-defined binding pocket is needed).")
    st.stop()

# ---------------------------------------------------------------- science + analysis
verified = pack.get("accurate", None)
badge = ""
if record.get("is_case"):
    badge = ("✅ Verified clean by adversarial fact-check"
             if verified else "🛠️ Auto-corrected by adversarial fact-check")

c1, c2 = st.columns([3, 2])
with c1:
    with st.expander(f"📖 The science — {record['display']}", expanded=True):
        if badge:
            st.caption(badge)
        st.markdown(f"**Overview.** {pack['overview']}")
        st.markdown(f"**Disease relevance.** {pack['disease']}")
        d = pack["drug"]
        st.markdown(f"**Inhibitor — {d['name']} ({d['code']}).** *{d['drug_class']}.* {d['mechanism']}")
        st.markdown(f"**Binding pocket.** {pack['pocket_summary']}")
with c2:
    with st.expander("🎓 How to read the movie", expanded=True):
        st.markdown(
            "- **▶ Play**, then drop **Speed → 0.1x** to study each step.\n"
            "- **Blue arrows** pull the ligand toward the pocket; **red clouds/arrows** are steric clashes.\n"
            "- Watch the **Live Scoring** bars and **Score Evolution** spike when a pose is rejected and "
            "drop when interactions form.\n"
            "- **Click any residue or the ligand** in 3D — the **AI Tutor** explains it for the current frame.\n"
            "- The **COMPUTED** badge + **RMSD** chip (top-left of the viewer) compare the engine's "
            "pose to the real crystal structure — under **2 Å** means the search rediscovered it.\n"
            "- The **Frame Inspector** shows translation, rotation, occupancy and clashes."
        )

with st.expander("📈 Trajectory analysis (Plotly)", expanded=False):
    import plotly.graph_objects as go

    frames = list(range(traj["n_frames"]))
    sc = traj["scores"]

    fig = go.Figure()
    for band in traj["phase_bands"]:
        fig.add_vrect(x0=band["start"], x1=band["end"], fillcolor=band["color"], opacity=0.10,
                      line_width=0, annotation_text=band["label"], annotation_position="top",
                      annotation=dict(font_size=9, font_color="#9fb0cc"))
    fig.add_trace(go.Scatter(x=frames, y=sc["total"], name="Total score",
                             line=dict(color="#38bdf8", width=2.4)))
    for ev in traj["events"]:
        col = {"bad": "#ef4444", "good": "#34d399", "warn": "#f59e0b"}.get(ev["kind"], "#93c5fd")
        fig.add_trace(go.Scatter(x=[ev["frame"]], y=[sc["total"][ev["frame"]]], mode="markers",
                                 marker=dict(size=8, color=col, line=dict(color="#0a0f1c", width=1)),
                                 name=ev["title"], hovertext=ev["title"], showlegend=False))
    fig.update_layout(template="plotly_dark", height=320, margin=dict(l=10, r=10, t=24, b=10),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      title="Docking score vs frame (lower = better)",
                      xaxis_title="Frame", yaxis_title="kcal/mol", showlegend=False)
    st.plotly_chart(fig, use_container_width=True)

    comp_names = [("hbond", "H-bond"), ("vdw", "Van der Waals"), ("hydrophobic", "Hydrophobic"),
                  ("electrostatic", "Electrostatic"), ("desolvation", "Desolvation"),
                  ("entropy", "Entropy penalty"), ("clash", "Steric clash")]
    palette = ["#fbbf24", "#34d399", "#94a3b8", "#60a5fa", "#a78bfa", "#f472b6", "#ef4444"]
    fig2 = go.Figure()
    for (k, nm), col in zip(comp_names, palette):
        fig2.add_trace(go.Scatter(x=frames, y=sc[k], name=nm, line=dict(color=col, width=1.6)))
    fig2.update_layout(template="plotly_dark", height=300, margin=dict(l=10, r=10, t=24, b=10),
                       paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                       title="Score components vs frame", xaxis_title="Frame", yaxis_title="kcal/mol",
                       legend=dict(orientation="h", y=-0.25, font=dict(size=10)))
    st.plotly_chart(fig2, use_container_width=True)

    cc1, cc2, cc3, cc4 = st.columns(4)
    cc1.metric("Final docking score", f"{traj['final_total']:.2f} kcal/mol")
    cc2.metric("Best score reached", f"{traj['best_total']:.2f} kcal/mol")
    cc3.metric("Key interactions", f"{len(traj['interactions'])}")
    rmsd = traj.get("rmsd")
    if rmsd is not None:
        verdict = "✓ success" if rmsd < 2 else ("≈ near" if rmsd < 3 else "✗ off-target")
        cc4.metric("RMSD vs crystal", f"{rmsd:.2f} Å", verdict,
                   delta_color="off", help="Deviation of the engine's computed pose from the "
                   "experimental crystal structure. < 2 Å is the standard re-docking success bar.")
    else:
        cc4.metric("RMSD vs crystal", "—", help="No crystal reference for custom uploads.")

st.caption("Poses are **computed** by a simplified rigid-body docking engine (Vina-style scoring, "
           "basin-hopping search) — the movie is the real search and the endpoint is the engine's own "
           "lowest-energy pose, not a pre-set answer. Case studies are validated by **RMSD** against the "
           "RCSB crystal structure. An educational engine, not a production docker. "
           "Built with Streamlit · 3Dmol.js · RDKit · Plotly · NumPy.")
