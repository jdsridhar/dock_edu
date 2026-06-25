"""
Cinematic docking workspace — a single self-contained HTML/JS component.

3Dmol.js is inlined from assets so the viewer works fully offline. All per-frame
data is precomputed in Python; a client-side animation engine drives the 3D
viewer and every synchronized panel (live scoring, shape-fit meter, occupancy,
score-evolution sparkline, frame inspector, interaction list, cinematic event
overlay, and the AI tutor) from one clock, giving smooth variable-speed playback.
"""
from __future__ import annotations

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_JS_PATH = os.path.join(BASE, "assets", "js", "3Dmol-min.js")
_3DMOL_JS = None


def _load_3dmol():
    global _3DMOL_JS
    if _3DMOL_JS is None:
        with open(_JS_PATH, encoding="utf-8") as fh:
            _3DMOL_JS = fh.read()
    return _3DMOL_JS


def build_workspace_html(record, traj, pack, tutor, options, height=940):
    data = {
        "protein_pdb": record["protein_pdb"],
        "traj": traj,
        "pack": pack,
        "tutor": tutor,
        "options": options,
        "meta": {
            "display": record["display"], "source_pdb": record["source_pdb"],
            "family": record["family"], "ligand_name": record["ligand_name"],
            "is_case": record.get("is_case", False),
            "computed": traj.get("computed", False), "rmsd": traj.get("rmsd"),
            "engine": traj.get("engine"), "exhaustiveness": traj.get("exhaustiveness"),
        },
    }
    # escape "</" so structure text containing "</script>" cannot break out of
    # the inline <script> tag (also a minor XSS guard for uploaded files)
    data_json = json.dumps(data, ensure_ascii=True).replace("</", "<\\/")
    html = _TEMPLATE
    html = html.replace("__HEIGHT__", str(height))
    html = html.replace("__3DMOL_JS__", _load_3dmol())
    html = html.replace("__DATA_JSON__", data_json)
    return html


_TEMPLATE = r"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  :root{
    --bg:#070a12; --bg2:#0d1322; --panel:#111a2e; --panel2:#0e1626;
    --line:#1f2b45; --txt:#e6edf7; --mut:#8da2c0; --accent:#38bdf8;
    --good:#34d399; --warn:#f59e0b; --bad:#ef4444; --pur:#a78bfa;
  }
  *{box-sizing:border-box;}
  html,body{margin:0;height:100%;background:var(--bg);color:var(--txt);
    font-family:'Segoe UI',system-ui,-apple-system,sans-serif;overflow:hidden;}
  #app{display:grid;grid-template-columns:1fr 372px;grid-template-rows:1fr auto;
    grid-template-areas:"stage panels" "transport transport";
    height:__HEIGHT__px;gap:8px;padding:8px;}
  #stage{grid-area:stage;position:relative;border-radius:12px;overflow:hidden;
    background:radial-gradient(120% 120% at 50% 0%,#0d1424 0%,#05080f 100%);
    border:1px solid var(--line);}
  #viewer{position:absolute;inset:0;}
  #panels{grid-area:panels;overflow-y:auto;overflow-x:hidden;display:flex;
    flex-direction:column;gap:8px;padding-right:2px;}
  #transport{grid-area:transport;background:var(--panel);border:1px solid var(--line);
    border-radius:12px;padding:10px 14px;}
  .card{background:linear-gradient(180deg,var(--panel) 0%,var(--panel2) 100%);
    border:1px solid var(--line);border-radius:11px;padding:11px 12px;}
  .card h3{margin:0 0 9px;font-size:11px;letter-spacing:.10em;text-transform:uppercase;
    color:var(--mut);font-weight:700;display:flex;justify-content:space-between;align-items:center;}
  .pill{font-size:9.5px;padding:2px 7px;border-radius:999px;background:#1b2740;color:var(--mut);
    letter-spacing:.04em;}
  /* HUD over the stage */
  .hud{position:absolute;pointer-events:none;font-size:12px;}
  #hudTitle{top:12px;left:14px;max-width:60%;}
  #hudTitle .t{font-size:15px;font-weight:700;}
  #hudTitle .s{font-size:11px;color:var(--mut);}
  #hudScore{top:12px;right:14px;text-align:right;}
  #hudScore .v{font-size:30px;font-weight:800;line-height:1;font-variant-numeric:tabular-nums;}
  #hudScore .l{font-size:10px;color:var(--mut);letter-spacing:.08em;text-transform:uppercase;}
  #phaseRibbon{position:absolute;left:14px;bottom:12px;display:flex;align-items:center;gap:8px;
    background:rgba(8,12,22,.7);border:1px solid var(--line);border-radius:999px;padding:5px 12px;
    font-size:12px;backdrop-filter:blur(6px);}
  #phaseDot{width:10px;height:10px;border-radius:50%;}
  #overlay{position:absolute;left:14px;bottom:54px;max-width:62%;background:rgba(9,13,24,.86);
    border-left:3px solid var(--accent);border-radius:8px;padding:10px 13px;backdrop-filter:blur(8px);
    box-shadow:0 8px 30px rgba(0,0,0,.5);opacity:0;transform:translateY(8px);
    transition:opacity .35s,transform .35s;}
  #overlay.show{opacity:1;transform:translateY(0);}
  #overlay .et{font-size:13.5px;font-weight:700;margin-bottom:3px;}
  #overlay .ec{font-size:11.5px;color:#c4d2e8;line-height:1.45;}
  #legend{position:absolute;right:12px;bottom:12px;display:flex;flex-direction:column;gap:3px;
    background:rgba(8,12,22,.62);border:1px solid var(--line);border-radius:8px;padding:7px 9px;
    font-size:10px;color:var(--mut);}
  #legend div{display:flex;align-items:center;gap:6px;}
  #legend i{width:14px;height:3px;border-radius:2px;display:inline-block;}
  /* score bars */
  .srow{display:grid;grid-template-columns:78px 1fr 48px;align-items:center;gap:8px;margin:5px 0;font-size:11px;}
  .srow .nm{color:var(--mut);}
  .track{position:relative;height:9px;background:#0a1120;border-radius:6px;overflow:hidden;}
  .track .center{position:absolute;left:50%;top:0;bottom:0;width:1px;background:#33415c;}
  .track .fill{position:absolute;top:0;bottom:0;border-radius:6px;transition:all .08s linear;}
  .srow .vv{text-align:right;font-variant-numeric:tabular-nums;color:var(--txt);}
  /* meters */
  .meter{height:13px;border-radius:7px;background:#0a1120;position:relative;overflow:hidden;}
  .meter .mf{position:absolute;left:0;top:0;bottom:0;border-radius:7px;transition:width .1s linear;}
  .mrow{display:flex;justify-content:space-between;font-size:11px;color:var(--mut);margin:2px 0 6px;}
  .mrow b{color:var(--txt);font-variant-numeric:tabular-nums;}
  /* inspector */
  .insp{display:grid;grid-template-columns:1fr 1fr;gap:5px 12px;font-size:11px;}
  .insp div{display:flex;justify-content:space-between;border-bottom:1px dashed #1a2540;padding:2px 0;}
  .insp span{color:var(--mut);} .insp b{font-variant-numeric:tabular-nums;font-weight:600;}
  /* interactions */
  #ilist{display:flex;flex-direction:column;gap:4px;max-height:150px;overflow-y:auto;}
  .iitem{display:grid;grid-template-columns:12px 1fr auto;gap:7px;align-items:center;font-size:11px;
    padding:3px 6px;background:#0c1426;border-radius:6px;}
  .iitem .dot{width:9px;height:9px;border-radius:50%;}
  .iitem .ty{color:var(--mut);font-size:9.5px;}
  /* tutor */
  #tutorBtns{display:flex;gap:5px;margin-bottom:8px;flex-wrap:wrap;}
  .tb{font-size:10.5px;padding:4px 9px;border-radius:7px;background:#15203a;color:var(--txt);
    border:1px solid var(--line);cursor:pointer;}
  .tb:hover{background:#1c2b4a;}
  .tb.active{background:var(--accent);color:#04121f;border-color:var(--accent);font-weight:700;}
  #tutorBody{font-size:12px;line-height:1.5;color:#d6e0f0;}
  #tutorBody .th{font-size:14px;font-weight:700;margin-bottom:4px;color:#fff;}
  #tutorBody .tg{font-size:10px;color:var(--mut);text-transform:uppercase;letter-spacing:.06em;margin:7px 0 2px;}
  #tutorHint{font-size:10px;color:var(--mut);margin-top:8px;font-style:italic;}
  /* transport */
  #tlwrap{position:relative;width:100%;height:34px;margin-bottom:8px;cursor:pointer;}
  #tlcanvas{width:100%;height:34px;display:block;}
  #ctrls{display:flex;align-items:center;gap:9px;flex-wrap:wrap;}
  .btn{background:#15203a;border:1px solid var(--line);color:var(--txt);border-radius:8px;
    width:38px;height:34px;font-size:14px;cursor:pointer;display:flex;align-items:center;justify-content:center;}
  .btn:hover{background:#1d2b4a;}
  .btn.play{width:46px;background:var(--accent);color:#04121f;border-color:var(--accent);font-weight:800;}
  #speed{background:#15203a;border:1px solid var(--line);color:var(--txt);border-radius:8px;height:34px;
    padding:0 8px;font-size:12px;}
  #frameReadout{margin-left:auto;font-size:13px;color:var(--mut);font-variant-numeric:tabular-nums;}
  #frameReadout b{color:var(--txt);font-size:16px;}
  .sep{width:1px;height:24px;background:var(--line);}
  /* simulated badge + display controls */
  .badge{display:inline-block;font-size:9px;font-weight:800;letter-spacing:.08em;
    padding:2px 6px;border-radius:5px;background:#3a2a08;color:#fbbf24;border:1px solid #5b430c;
    vertical-align:middle;margin-left:7px;}
  .badge.ok{background:#08291c;color:#34d399;border-color:#13513a;}
  .rmsdchip{display:inline-block;font-size:10px;font-weight:700;padding:2px 7px;border-radius:999px;
    margin-left:7px;vertical-align:middle;background:#0c1426;border:1px solid var(--line);}
  .drow{display:grid;grid-template-columns:96px 1fr;align-items:center;gap:8px;margin:5px 0;font-size:11px;}
  .drow span{color:var(--mut);}
  .drow select,.drow input[type=range]{width:100%;}
  .drow select{background:#0c1426;border:1px solid var(--line);color:var(--txt);border-radius:6px;
    height:26px;font-size:11px;padding:0 6px;}
  .drow input[type=range]{accent-color:var(--accent);}
  .drow .chk{justify-self:start;width:16px;height:16px;accent-color:var(--accent);}
  #surfNote{font-size:10px;color:var(--warn);margin-top:4px;}
  ::-webkit-scrollbar{width:8px;height:8px;}
  ::-webkit-scrollbar-thumb{background:#22304d;border-radius:4px;}
  ::-webkit-scrollbar-track{background:transparent;}
  .loadmsg{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
    color:var(--mut);font-size:13px;z-index:5;}
  @media (max-width:880px){
    #app{grid-template-columns:1fr;grid-template-areas:"stage" "panels" "transport";height:auto;}
    #stage{min-height:420px;height:60vh;} #panels{max-height:none;}
  }
</style>
</head>
<body>
<div id="app">
  <div id="stage">
    <div id="viewer"></div>
    <div class="loadmsg" id="loadmsg">Rendering molecular scene…</div>
    <div class="hud" id="hudTitle"><div class="t"><span id="hudName">—</span><span class="badge" id="modeBadge" title="">COMPUTED</span><span class="rmsdchip" id="hudRmsd" style="display:none"></span></div><div class="s" id="hudSub"></div></div>
    <div class="hud" id="hudScore"><div class="v" id="hudV">—</div><div class="l">Docking Score (kcal/mol)</div></div>
    <div id="phaseRibbon"><span id="phaseDot"></span><span id="phaseLbl">—</span></div>
    <div id="overlay"><div class="et" id="ovT"></div><div class="ec" id="ovC"></div></div>
    <div id="legend">
      <div><i style="background:#3b82f6"></i>Attractive force</div>
      <div><i style="background:#ef4444"></i>Repulsive / clash</div>
      <div><i style="background:#fbbf24"></i>Hydrogen bond</div>
      <div><i style="background:#a78bfa"></i>π–π stacking</div>
    </div>
  </div>

  <div id="panels">
    <div class="card" id="displayCard">
      <h3>Display</h3>
      <div class="drow"><span>Representation</span>
        <select id="ctlRep">
          <option value="cartoon">Cartoon</option>
          <option value="surface">Surface</option>
          <option value="stick">Stick</option>
          <option value="ball_and_stick">Ball &amp; stick</option>
        </select></div>
      <div class="drow"><span>Colour map</span>
        <select id="ctlCol">
          <option value="spectrum">Spectrum (rainbow)</option>
          <option value="electrostatic">Charge map</option>
          <option value="hydrophobic">Hydropathy (Kyte–Doolittle)</option>
          <option value="hydrophilic">Polarity map</option>
        </select></div>
      <div class="drow"><span>Opacity</span>
        <input type="range" id="ctlOp" min="0.2" max="1" step="0.05"></div>
      <div class="drow"><span>Surface</span>
        <input type="checkbox" id="ctlSurf" class="chk"></div>
      <div class="drow"><span>Surface opacity</span>
        <input type="range" id="ctlSurfOp" min="0.1" max="1" step="0.05"></div>
      <div id="surfNote" style="display:none">Building molecular surface…</div>
    </div>
    <div class="card">
      <h3>Live Scoring Engine <span class="pill" id="trendPill">—</span></h3>
      <div id="sbars"></div>
    </div>
    <div class="card">
      <h3>Shape Complementarity</h3>
      <div class="mrow"><span>Shape fit</span><b id="fitV">0%</b></div>
      <div class="meter"><div class="mf" id="fitBar" style="width:0%"></div></div>
      <div class="mrow" style="margin-top:9px"><span>Pocket occupancy</span><b id="occV">0%</b></div>
      <div class="meter"><div class="mf" id="occBar" style="width:0%;background:#38bdf8"></div></div>
    </div>
    <div class="card">
      <h3>Score Evolution</h3>
      <canvas id="spark" height="96"></canvas>
    </div>
    <div class="card">
      <h3>Frame Inspector</h3>
      <div class="insp" id="insp"></div>
      <div style="font-size:9.5px;color:var(--mut);margin-top:7px;font-style:italic">
        Rigid-body docking: translation &amp; rotation are the searched pose offset; the ligand
        conformer is fixed, so torsions stay at 0.</div>
    </div>
    <div class="card">
      <h3>Active Interactions <span class="pill" id="iCount">0</span></h3>
      <div id="ilist"></div>
    </div>
    <div class="card">
      <h3>AI Tutor</h3>
      <div id="tutorBtns">
        <button class="tb active" data-mode="frame">This frame</button>
        <button class="tb" data-mode="pocket">Pocket</button>
        <button class="tb" data-mode="drug">Drug</button>
      </div>
      <div id="tutorBody"></div>
      <div id="tutorHint">Click any residue or the ligand in the 3D view for an explanation.</div>
    </div>
  </div>

  <div id="transport">
    <div id="tlwrap"><canvas id="tlcanvas"></canvas></div>
    <div id="ctrls">
      <button class="btn" id="bRewind" title="Rewind to start">⏮</button>
      <button class="btn" id="bBack" title="Step back">⏪</button>
      <button class="btn play" id="bPlay" title="Play / Pause">▶</button>
      <button class="btn" id="bFwd" title="Step forward">⏩</button>
      <button class="btn" id="bEnd" title="Jump to final pose">⏭</button>
      <button class="btn" id="bStop" title="Stop">⏹</button>
      <div class="sep"></div>
      <span style="font-size:11px;color:var(--mut)">Speed</span>
      <select id="speed">
        <option value="0.1">0.1x</option><option value="0.25">0.25x</option>
        <option value="0.5">0.5x</option><option value="1" selected>1x</option>
        <option value="2">2x</option><option value="5">5x</option>
      </select>
      <div id="frameReadout">Frame <b id="fNow">0</b> / <span id="fTot">0</span></div>
    </div>
  </div>
</div>

<script>__3DMOL_JS__</script>
<script id="docking-data" type="application/json">__DATA_JSON__</script>
<script>
(function(){
  "use strict";
  const DATA = JSON.parse(document.getElementById('docking-data').textContent);
  const T = DATA.traj, PACK = DATA.pack, TUT = DATA.tutor, OPT = DATA.options, META = DATA.meta;
  const N = T.n_frames;
  const BASE_FPS = 30;

  // ---- colours ----
  const ICOLOR = {hydrogen_bond:'#fbbf24', pi_stacking:'#a78bfa', salt_bridge:'#f97316',
    hydrophobic:'#94a3b8', metal_coordination:'#22d3ee', vdw:'#6b7280',
    electrostatic_attraction:'#60a5fa', electrostatic_repulsion:'#ef4444'};
  const ILABEL = {hydrogen_bond:'H-bond', pi_stacking:'π–π', salt_bridge:'salt bridge',
    hydrophobic:'hydrophobic', metal_coordination:'metal', vdw:'van der Waals'};
  const RES_CHARGE = {ARG:1,LYS:1,HIS:0.5,ASP:-1,GLU:-1};
  const RES_HYDRO = {ILE:4.5,VAL:4.2,LEU:3.8,PHE:2.8,CYS:2.5,MET:1.9,ALA:1.8,GLY:-0.4,
    THR:-0.7,SER:-0.8,TRP:-0.9,TYR:-1.3,PRO:-1.6,HIS:-3.2,GLU:-3.5,GLN:-3.5,
    ASP:-3.5,ASN:-3.5,LYS:-3.9,ARG:-4.5};

  // ---- 3Dmol setup ----
  const viewer = $3Dmol.createViewer(document.getElementById('viewer'),
    {backgroundColor:'#070a12', antialias:true});
  const protModel = viewer.addModel(DATA.protein_pdb, 'pdb');

  const POCKET = T.pocket_center;

  // ---- live display state (owned by the component, persisted across reruns) ----
  function loadDisplay(){ try{ return JSON.parse(sessionStorage.getItem('vmds:display')); }catch(e){ return null; } }
  function saveDisplay(){ try{ sessionStorage.setItem('vmds:display', JSON.stringify(DISPLAY)); }catch(e){} }
  const DISPLAY = Object.assign({
    representation: OPT.representation, color_scheme: OPT.color_scheme,
    opacity: OPT.opacity, surface: OPT.show_surface, surface_opacity: OPT.surface_opacity,
  }, loadDisplay() || {});

  // residue -> colour map (reliable across 3Dmol versions via colorscheme map)
  const ALLRES=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',
    'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'];
  function colorFor(r){
    if(DISPLAY.color_scheme==='electrostatic'){
      const c=RES_CHARGE[r];
      if(c>0) return '#3b82f6'; if(c<0) return '#ef4444'; if(c===0.5) return '#c4b5fd';
      return ({SER:1,THR:1,ASN:1,GLN:1,TYR:1,CYS:1,TRP:1}[r])?'#cfd8e8':'#7c8aa5';
    }
    if(DISPLAY.color_scheme==='hydrophobic'){
      const h=RES_HYDRO[r]; if(h===undefined) return '#9aa7bd';
      return h>1.5?'#fb923c':(h>-0.5?'#fcd34d':'#22d3ee');
    }
    if(DISPLAY.color_scheme==='hydrophilic'){
      const h=RES_HYDRO[r]; if(h===undefined) return '#9aa7bd';
      return h<-1.5?'#34d399':(h<0?'#7dd3fc':'#334155');
    }
    return '#5b7bb0';
  }
  function hexInt(h){ return parseInt(h.slice(1),16); }
  let COLORMAP={};
  function rebuildColorMap(){ COLORMAP={}; ALLRES.forEach(r=>COLORMAP[r]=hexInt(colorFor(r))); }
  function colSpec(isCartoon){
    if(DISPLAY.color_scheme==='spectrum') return isCartoon?{color:'spectrum'}:{};
    return {colorscheme:{prop:'resn', map:COLORMAP}};
  }
  function baseProteinStyle(){
    const op=DISPLAY.opacity; let s={};
    if(DISPLAY.representation==='stick'){ s={stick:Object.assign({radius:0.13}, colSpec(false))}; }
    else if(DISPLAY.representation==='ball_and_stick'){
      s={stick:Object.assign({radius:0.14}, colSpec(false)),
         sphere:Object.assign({scale:0.23}, colSpec(false))};
    } else if(DISPLAY.representation==='surface'){ s={}; }  // surface drawn separately
    else { s={cartoon:Object.assign({opacity:op}, colSpec(true))}; }
    protModel.setStyle({}, s);
  }

  const pocketResis = T.pocket_residues.map(r=>parseInt(r.resi)).filter(n=>!isNaN(n));
  protModel.setClickable({resi: pocketResis}, true, function(atom){ showResidue(atom.resn, ''+atom.resi); });

  // ligand as multi-frame trajectory model
  const ligModel = viewer.addModelsAsFrames(T.ligand_pdb, 'pdb');
  ligModel.setStyle({}, {stick:{radius:0.22, colorscheme:'yellowCarbon'},
    sphere:{scale:0.30, colorscheme:'yellowCarbon'}});
  ligModel.setClickable({}, true, function(){ setMode('drug'); });

  // ---- molecular surface (restricted to the pocket region, built async) ----
  let _surfId=null, _nearby=null;
  function clearSurfaces(){
    try{ if(typeof viewer.removeAllSurfaces==='function') viewer.removeAllSurfaces();
         else if(_surfId!=null && viewer.removeSurface) viewer.removeSurface(_surfId); }catch(e){}
    _surfId=null;
  }
  function nearbyResisForSurface(){
    if(_nearby) return _nearby;
    const atoms=protModel.selectedAtoms({}); const set={};
    for(let k=0;k<atoms.length;k++){ const a=atoms[k];
      const dx=a.x-POCKET[0], dy=a.y-POCKET[1], dz=a.z-POCKET[2];
      if(dx*dx+dy*dy+dz*dz < 23*23){ const ri=parseInt(a.resi); if(!isNaN(ri)) set[ri]=1; } }
    _nearby=Object.keys(set).map(Number); if(!_nearby.length) _nearby=pocketResis; return _nearby;
  }
  function buildSurface(){
    const note=document.getElementById('surfNote'); if(note) note.style.display='block';
    setTimeout(function(){
      try{
        const ss=(DISPLAY.color_scheme==='spectrum')?{color:'#5b7bb0'}
          :{colorscheme:{prop:'resn', map:COLORMAP}};
        const ret=viewer.addSurface($3Dmol.SurfaceType.VDW,
          Object.assign({opacity:DISPLAY.surface_opacity}, ss), {resi:nearbyResisForSurface()});
        if(typeof ret==='number') _surfId=ret;
        viewer.render();
      }catch(e){ console.warn('surface failed', e); }
      if(note) note.style.display='none';
    }, 30);
  }

  function applyProteinOnly(){
    rebuildColorMap(); baseProteinStyle();
    viewer.addStyle({resi: pocketResis}, {stick:{radius:0.18, colorscheme:'cyanCarbon'}});
    viewer.render();
  }
  function applyDisplay(){
    applyProteinOnly();
    clearSurfaces();
    if(DISPLAY.surface || DISPLAY.representation==='surface') buildSurface();
    viewer.render();
  }

  function wireDisplay(){
    const rep=document.getElementById('ctlRep'), col=document.getElementById('ctlCol'),
      op=document.getElementById('ctlOp'), su=document.getElementById('ctlSurf'),
      so=document.getElementById('ctlSurfOp');
    rep.value=DISPLAY.representation; col.value=DISPLAY.color_scheme;
    op.value=DISPLAY.opacity; su.checked=!!DISPLAY.surface; so.value=DISPLAY.surface_opacity;
    rep.onchange=()=>{ DISPLAY.representation=rep.value; saveDisplay(); applyDisplay(); };
    col.onchange=()=>{ DISPLAY.color_scheme=col.value; saveDisplay(); applyDisplay(); };
    op.oninput =()=>{ DISPLAY.opacity=parseFloat(op.value); saveDisplay(); applyProteinOnly(); };
    su.onchange=()=>{ DISPLAY.surface=su.checked; saveDisplay(); applyDisplay(); };
    so.onchange=()=>{ DISPLAY.surface_opacity=parseFloat(so.value); saveDisplay(); applyDisplay(); };
  }

  wireDisplay();
  applyDisplay();
  if(pocketResis.length) viewer.zoomTo({resi: pocketResis}); else viewer.zoomTo({});
  viewer.zoom(0.82);
  viewer.render();
  document.getElementById('loadmsg').style.display='none';

  // ---- dynamic shapes per frame ----
  const CFINAL = T.centroid[N-1];
  function vadd(a,b){return {x:a[0]+b.x,y:a[1]+b.y,z:a[2]+b.z};}
  function frameLig(it, cen){
    return {x: it.lig_xyz[0]+(cen[0]-CFINAL[0]),
            y: it.lig_xyz[1]+(cen[1]-CFINAL[1]),
            z: it.lig_xyz[2]+(cen[2]-CFINAL[2])};
  }
  function bandOf(i){ for(const b of T.phase_bands){ if(i>=b.start && i<=b.end) return b;} return T.phase_bands[T.phase_bands.length-1]; }

  function drawShapes(i){
    const cen = T.centroid[i];
    const cenV = {x:cen[0],y:cen[1],z:cen[2]};
    const pocV = {x:POCKET[0],y:POCKET[1],z:POCKET[2]};
    const total = T.scores.total[i];
    const favor = Math.max(0, Math.min(1, -total/12));
    const band = bandOf(i);

    // pocket cavity glow (brighter during detection)
    const glow = band.phase==='pocket_detection' ? 0.22 : 0.07;
    viewer.addSphere({center:pocV, radius:T.pocket_radius, color:'#a78bfa',
      opacity:glow, wireframe: band.phase!=='pocket_detection'});

    // decoy cavities during detection
    if(band.phase==='pocket_detection'){
      const prog = (i-band.start)/Math.max(1,(band.end-band.start));
      T.decoy_pockets.forEach(d=>{
        const fade = Math.max(0, 0.5*(1-prog));
        if(fade>0.02) viewer.addSphere({center:{x:d.pos[0],y:d.pos[1],z:d.pos[2]},
          radius:3.0, color:'#64748b', opacity:fade, wireframe:true});
      });
    }

    // attractive force toward pocket when outside
    const dx=POCKET[0]-cen[0], dy=POCKET[1]-cen[1], dz=POCKET[2]-cen[2];
    const dlen=Math.sqrt(dx*dx+dy*dy+dz*dz)||1;
    if(dlen>1.3){
      const L = Math.min(dlen-0.5, 2.0+6.0*favor);
      const end={x:cen[0]+dx/dlen*L, y:cen[1]+dy/dlen*L, z:cen[2]+dz/dlen*L};
      viewer.addArrow({start:cenV, end:end, radius:0.16, color:'#3b82f6', mid:0.8});
    }

    // active interactions
    let nactive=0;
    for(const it of T.interactions){
      if(i < it.onset) continue;
      nactive++;
      const lp = frameLig(it, cen);
      const pp = {x:it.prot[0],y:it.prot[1],z:it.prot[2]};
      const col = ICOLOR[it.type]||'#cbd5e1';
      viewer.addLine({start:lp, end:pp, color:col, dashed:true, linewidth:2.5});
      if(it.type==='hydrogen_bond' || it.type==='pi_stacking'){
        viewer.addLabel(it.label+(it.dist?(' '+it.dist+'Å'):''),
          {position:{x:(lp.x+pp.x)/2,y:(lp.y+pp.y)/2,z:(lp.z+pp.z)/2},
           fontSize:9, fontColor:'#06121f', backgroundColor:col, backgroundOpacity:0.85,
           borderThickness:0, inFront:true});
      }
    }

    // clash visualisation
    const clash = T.clash[i];
    if(clash>0.25){
      const pen = T.scores.clash[i];
      T.clash_residues.forEach((cr,idx)=>{
        const mid={x:(cen[0]+cr.ca[0])/2,y:(cen[1]+cr.ca[1])/2,z:(cen[2]+cr.ca[2])/2};
        viewer.addSphere({center:mid, radius:0.7+1.4*clash, color:'#ef4444', opacity:0.30});
        // repulsive arrow pushing ligand away from residue
        const rx=cen[0]-cr.ca[0], ry=cen[1]-cr.ca[1], rz=cen[2]-cr.ca[2];
        const rl=Math.sqrt(rx*rx+ry*ry+rz*rz)||1;
        viewer.addArrow({start:{x:cr.ca[0],y:cr.ca[1],z:cr.ca[2]},
          end:{x:cr.ca[0]+rx/rl*(1.5+2*clash),y:cr.ca[1]+ry/rl*(1.5+2*clash),z:cr.ca[2]+rz/rl*(1.5+2*clash)},
          radius:0.16, color:'#ef4444', mid:0.8});
        if(idx===0) viewer.addLabel('Steric clash  +'+pen.toFixed(1),
          {position:mid, fontSize:10, fontColor:'#fff', backgroundColor:'#b91c1c',
           backgroundOpacity:0.9, borderThickness:0, inFront:true});
      });
    }
    return nactive;
  }

  // ---- panels ----
  const COMPS = [['hbond','H-Bond'],['vdw','Van der Waals'],['hydrophobic','Hydrophobic'],
    ['electrostatic','Electrostatic'],['desolvation','Desolvation'],
    ['entropy','Entropy Penalty'],['clash','Steric Clash']];
  const sbars = document.getElementById('sbars');
  const barEls = {};
  COMPS.forEach(([k,nm])=>{
    const row=document.createElement('div'); row.className='srow';
    row.innerHTML='<div class="nm">'+nm+'</div><div class="track"><div class="center"></div>'
      +'<div class="fill"></div></div><div class="vv">0.0</div>';
    sbars.appendChild(row);
    barEls[k]={fill:row.querySelector('.fill'), val:row.querySelector('.vv')};
  });
  function setBar(k, value){
    const el=barEls[k]; const scale=7.0;
    const frac=Math.max(-1,Math.min(1,value/scale));
    const w=Math.abs(frac)*50;
    el.fill.style.width=w+'%';
    if(value<=0){ el.fill.style.right='50%'; el.fill.style.left=''; el.fill.style.background='#34d399'; }
    else { el.fill.style.left='50%'; el.fill.style.right=''; el.fill.style.background='#ef4444'; }
    el.val.textContent=value.toFixed(1);
  }

  const insp=document.getElementById('insp');
  const inspRows=[['Frame','fr'],['Phase','ph'],['Translation','tr'],['Rotation','ro'],
    ['Torsions','to'],['Occupancy','oc'],['H-bonds','hb'],['Clashes','cl'],['Score','sc']];
  const inspEls={};
  inspRows.forEach(([lbl,k])=>{
    const d=document.createElement('div'); d.innerHTML='<span>'+lbl+'</span><b></b>';
    insp.appendChild(d); inspEls[k]=d.querySelector('b');
  });

  function fitColor(v){ return v<40?'#ef4444':(v<70?'#f59e0b':'#34d399'); }

  // sparkline (score evolution)
  const spark=document.getElementById('spark');
  const sctx=spark.getContext('2d');
  function sizeSpark(){ spark.width=spark.clientWidth*devicePixelRatio; spark.height=96*devicePixelRatio;
    sctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0); }
  const totals=T.scores.total;
  let tmin=Math.min.apply(null,totals), tmax=Math.max.apply(null,totals);
  const pad=(tmax-tmin)*0.1+0.5; tmin-=pad; tmax+=pad;
  function drawSpark(i){
    const W=spark.clientWidth, H=96;
    sctx.clearRect(0,0,W,H);
    // phase band tints
    T.phase_bands.forEach(b=>{
      sctx.fillStyle=b.color+'18';
      const x0=b.start/(N-1)*W, x1=b.end/(N-1)*W;
      sctx.fillRect(x0,0,x1-x0,H);
    });
    // zero line
    const yz=H-(0-tmin)/(tmax-tmin)*H;
    sctx.strokeStyle='#33415c'; sctx.lineWidth=1; sctx.setLineDash([3,3]);
    sctx.beginPath(); sctx.moveTo(0,yz); sctx.lineTo(W,yz); sctx.stroke(); sctx.setLineDash([]);
    // curve
    sctx.beginPath();
    for(let f=0; f<N; f++){
      const x=f/(N-1)*W, y=H-(totals[f]-tmin)/(tmax-tmin)*H;
      f===0?sctx.moveTo(x,y):sctx.lineTo(x,y);
    }
    sctx.strokeStyle='#38bdf8'; sctx.lineWidth=1.8; sctx.stroke();
    // filled area under curve up to i
    const xi=i/(N-1)*W;
    sctx.lineTo(xi,H); sctx.lineTo(0,H); sctx.closePath();
    sctx.fillStyle='rgba(56,189,248,.10)'; sctx.fill();
    // marker
    const ym=H-(totals[i]-tmin)/(tmax-tmin)*H;
    sctx.strokeStyle='#fcd34d'; sctx.lineWidth=1.4;
    sctx.beginPath(); sctx.moveTo(xi,0); sctx.lineTo(xi,H); sctx.stroke();
    sctx.fillStyle='#fcd34d'; sctx.beginPath(); sctx.arc(xi,ym,3.2,0,7); sctx.fill();
  }

  // timeline canvas
  const tlc=document.getElementById('tlcanvas'); const tctx=tlc.getContext('2d');
  function sizeTL(){ tlc.width=tlc.clientWidth*devicePixelRatio; tlc.height=34*devicePixelRatio;
    tctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0); }
  function drawTL(i){
    const W=tlc.clientWidth, H=34;
    tctx.clearRect(0,0,W,H);
    T.phase_bands.forEach(b=>{
      const x0=b.start/(N-1)*W, x1=b.end/(N-1)*W;
      tctx.fillStyle=b.color+'55'; tctx.fillRect(x0,6,x1-x0,16);
      tctx.fillStyle=b.color; tctx.fillRect(x0,6,1.5,16);
    });
    // event ticks
    T.events.forEach(e=>{
      const x=e.frame/(N-1)*W;
      tctx.fillStyle = e.kind==='bad'?'#ef4444':(e.kind==='good'?'#34d399':(e.kind==='warn'?'#f59e0b':'#93c5fd'));
      tctx.beginPath(); tctx.moveTo(x,3); tctx.lineTo(x-3,0); tctx.lineTo(x+3,0); tctx.closePath(); tctx.fill();
    });
    // progress + playhead
    const xi=i/(N-1)*W;
    tctx.fillStyle='rgba(56,189,248,.18)'; tctx.fillRect(0,6,xi,16);
    tctx.fillStyle='#fcd34d'; tctx.fillRect(xi-1,2,2.4,28);
  }

  // ---- overlay events ----
  let lastEventIdx=-1;
  function updateOverlay(i){
    let cur=-1;
    for(let k=0;k<T.events.length;k++){ if(T.events[k].frame<=i) cur=k; else break; }
    const ov=document.getElementById('overlay');
    if(cur<0){ ov.classList.remove('show'); return; }
    const e=T.events[cur];
    const linger=Math.max(18, Math.round(N*0.08));
    if(i-e.frame>linger){ ov.classList.remove('show'); return; }
    if(cur!==lastEventIdx){
      document.getElementById('ovT').textContent=e.title;
      document.getElementById('ovC').textContent=e.caption;
      const col=e.kind==='bad'?'#ef4444':(e.kind==='good'?'#34d399':(e.kind==='warn'?'#f59e0b':'#38bdf8'));
      ov.style.borderLeftColor=col;
      lastEventIdx=cur;
    }
    ov.classList.add('show');
  }

  // ---- interactions list ----
  const ilist=document.getElementById('ilist');
  function updateIList(i){
    const act=T.interactions.filter(it=>i>=it.onset);
    document.getElementById('iCount').textContent=act.length;
    ilist.innerHTML='';
    act.forEach(it=>{
      const d=document.createElement('div'); d.className='iitem';
      d.innerHTML='<span class="dot" style="background:'+(ICOLOR[it.type]||'#ccc')+'"></span>'
        +'<div>'+it.label+'<div class="ty">'+(ILABEL[it.type]||it.type)+'</div></div>'
        +'<div style="color:var(--mut)">'+(it.dist?it.dist+' Å':'')+'</div>';
      d.onclick=()=>showResidue(it.resn, it.resi);
      ilist.appendChild(d);
    });
    if(!act.length) ilist.innerHTML='<div style="color:var(--mut);font-size:11px;padding:4px">No interactions formed yet.</div>';
  }

  // ---- AI tutor ----
  const RESCLASS={}; (TUT.residue_classes||[]).forEach(r=>RESCLASS[r.resn]=r);
  const RESNOTE={}; (PACK.residue_notes||[]).forEach(r=>RESNOTE[r.resn+r.resi]=r.role);
  const INTNOTE={}; (PACK.interaction_notes||[]).forEach(r=>INTNOTE[r.resn+r.resi]=r);
  const ICONCEPT={}; (TUT.interactions||[]).forEach(r=>ICONCEPT[r.key]=r);
  let tutorMode='frame';
  function setMode(m){ tutorMode=m;
    document.querySelectorAll('.tb').forEach(b=>b.classList.toggle('active', b.dataset.mode===m));
    renderTutor(); }
  document.querySelectorAll('.tb').forEach(b=>b.onclick=()=>setMode(b.dataset.mode));

  function showResidue(resn, resi){
    tutorMode='residue:'+resn+':'+resi;
    document.querySelectorAll('.tb').forEach(b=>b.classList.remove('active'));
    renderTutor();
  }
  function esc(s){ return (s||'').replace(/</g,'&lt;'); }
  function renderTutor(){
    const body=document.getElementById('tutorBody'); const i=cur;
    if(tutorMode==='drug'){
      const d=PACK.drug;
      body.innerHTML='<div class="th">'+esc(d.name)+'  <span style="font-size:11px;color:var(--mut)">('+esc(d.code)+')</span></div>'
        +'<div class="tg">Class</div>'+esc(d.drug_class)
        +'<div class="tg">Mechanism</div>'+esc(d.mechanism); return;
    }
    if(tutorMode==='pocket'){
      const p=T.real_pocket;
      body.innerHTML='<div class="th">Binding Pocket</div>'+esc(PACK.pocket_summary)
        +'<div class="tg">Metrics</div>'
        +'<div class="insp" style="margin-top:4px">'
        +row2('Pocket score',p.score)+row2('Volume',p.volume+' Å³')+row2('Depth',p.depth+' Å')
        +row2('Druggability',p.drug)+row2('Hydrophobicity',p.hydro)+'</div>'; return;
    }
    if(tutorMode.startsWith('residue:')){
      const [,resn,resi]=tutorMode.split(':');
      const rc=RESCLASS[resn]||{name:resn,property:'',role:''};
      const note=RESNOTE[resn+resi];
      const inote=INTNOTE[resn+resi];
      let h='<div class="th">'+esc(resn)+esc(resi)+' · '+esc(rc.name)+'</div>'
        +'<div class="tg">Property</div>'+esc(rc.property)
        +'<div class="tg">General role</div>'+esc(rc.role);
      if(note){ h+='<div class="tg">In this binding site</div>'+esc(note); }
      if(inote){ h+='<div class="tg">'+esc((ILABEL[inote.type]||inote.type))+' interaction</div>'+esc(inote.note); }
      body.innerHTML=h; return;
    }
    // default: explain current frame
    const band=bandOf(i);
    const beat=(PACK.beats||[]).find(b=>b.phase===band.phase) || {title:band.label, caption:''};
    const act=T.interactions.filter(it=>i>=it.onset);
    let trend='';
    if(i>3){ const d=T.scores.total[i]-T.scores.total[i-Math.min(i,8)];
      trend = d<-0.2?'The score is <b style="color:#34d399">improving</b> as binding strengthens.'
        : (d>0.2?'The score is <b style="color:#ef4444">worsening</b> — the pose is unfavourable here.'
        :'The score is roughly stable.'); }
    let h='<div class="th">'+esc(beat.title)+'</div>'+esc(beat.caption);
    h+='<div class="tg">What is happening</div>'+trend;
    if(act.length){ h+='<div class="tg">Interactions holding the ligand</div>'
      + act.map(a=>a.label+' ('+(ILABEL[a.type]||a.type)+')').join(', '); }
    if(T.clash[i]>0.25){ h+='<div class="tg" style="color:#ef4444">Warning</div>'
      +'Steric clash detected — atoms overlap, so this pose is penalised and likely rejected.'; }
    body.innerHTML=h;
  }
  function row2(l,v){ return '<div><span>'+l+'</span><b>'+v+'</b></div>'; }

  // ---- master frame render ----
  let cur=0, lastRendered=-1;
  function updatePanels(i){
    const sc=T.scores;
    setBar('hbond',sc.hbond[i]); setBar('vdw',sc.vdw[i]); setBar('hydrophobic',sc.hydrophobic[i]);
    setBar('electrostatic',sc.electrostatic[i]); setBar('desolvation',sc.desolvation[i]);
    setBar('entropy',sc.entropy[i]); setBar('clash',sc.clash[i]);
    const tot=sc.total[i];
    const hv=document.getElementById('hudV'); hv.textContent=tot.toFixed(1);
    hv.style.color = tot<-6?'#34d399':(tot<0?'#a3e635':(tot<3?'#f59e0b':'#ef4444'));
    // trend pill
    const tp=document.getElementById('trendPill');
    if(i>3){ const d=tot-sc.total[i-Math.min(i,8)];
      tp.textContent=d<-0.2?'improving ▼':(d>0.2?'worsening ▲':'stable');
      tp.style.color=d<-0.2?'#34d399':(d>0.2?'#ef4444':'#8da2c0'); }
    // meters
    const fit=T.shape_fit[i], occ=T.occupancy[i];
    const fb=document.getElementById('fitBar'); fb.style.width=fit+'%'; fb.style.background=fitColor(fit);
    document.getElementById('fitV').textContent=Math.round(fit)+'%';
    document.getElementById('occBar').style.width=occ+'%';
    document.getElementById('occV').textContent=Math.round(occ)+'%';
    // inspector
    const tr=T.translation[i], eu=T.euler[i], to=T.torsions[i], band=bandOf(i);
    inspEls.fr.textContent=i; inspEls.ph.textContent=band.label;
    inspEls.tr.textContent=tr[0]+', '+tr[1]+', '+tr[2]+' Å';
    inspEls.ro.textContent=eu[0]+'°, '+eu[1]+'°, '+eu[2]+'°';
    inspEls.to.textContent=to[0]+'°, '+to[1]+'°, '+to[2]+'°';
    inspEls.oc.textContent=Math.round(occ)+'%';
    inspEls.hb.textContent=T.n_hbond[i];
    inspEls.cl.textContent=T.n_clash[i];
    inspEls.sc.textContent=tot.toFixed(2);
    // phase ribbon
    document.getElementById('phaseLbl').textContent=band.label;
    document.getElementById('phaseDot').style.background=band.color;
    // readouts
    document.getElementById('fNow').textContent=i;
    updateIList(i); updateOverlay(i); drawSpark(i); drawTL(i);
    if(tutorMode==='frame') renderTutor();
  }
  function renderFrame(f){
    const i=Math.max(0,Math.min(N-1,Math.round(f)));
    cur=i;
    if(i===lastRendered) return;
    lastRendered=i;
    ligModel.setFrame(i);
    viewer.removeAllShapes(); viewer.removeAllLabels();
    drawShapes(i);
    viewer.render();
    updatePanels(i);
  }

  // ---- playback engine ----
  const STATE_KEY='vmds:state:'+META.display+':'+N;
  function loadPlayState(){ try{ return JSON.parse(sessionStorage.getItem(STATE_KEY)); }catch(e){ return null; } }
  function savePlayState(){ try{ sessionStorage.setItem(STATE_KEY,
    JSON.stringify({frac:cur/(N-1), playing:playing, speed:speed})); }catch(e){} }
  let playing=false, dir=1, speed=1, fpos=0, last=null, _saveTick=0;
  function setPlay(p){ playing=p; document.getElementById('bPlay').textContent=p?'⏸':'▶';
    if(p) last=null; savePlayState(); }
  function togglePlay(){ if(fpos>=N-1){ fpos=0; } dir=1; setPlay(!playing); }
  function tick(ts){
    if(playing){
      if(last==null) last=ts;
      const dt=Math.min(0.1,(ts-last)/1000); last=ts;
      fpos+=dir*speed*BASE_FPS*dt;
      if(fpos>=N-1){ fpos=N-1; setPlay(false); }
      if(fpos<=0){ fpos=0; if(dir<0) setPlay(false); }
      renderFrame(fpos);
      if((++_saveTick % 15)===0) savePlayState();
    }
    requestAnimationFrame(tick);
  }
  function seek(f){ fpos=Math.max(0,Math.min(N-1,f)); renderFrame(fpos); savePlayState(); }

  document.getElementById('bPlay').onclick=togglePlay;
  document.getElementById('bStop').onclick=()=>{ setPlay(false); seek(0); };
  document.getElementById('bRewind').onclick=()=>{ seek(0); };
  document.getElementById('bEnd').onclick=()=>{ setPlay(false); seek(N-1); };
  document.getElementById('bBack').onclick=()=>{ setPlay(false); seek(fpos-Math.max(1,Math.round(N*0.01))); };
  document.getElementById('bFwd').onclick=()=>{ setPlay(false); seek(fpos+Math.max(1,Math.round(N*0.01))); };
  document.getElementById('speed').onchange=(e)=>{ speed=parseFloat(e.target.value); savePlayState(); };

  // keyboard controls: space = play/pause, ←/→ step (shift = jump), Home/End
  window.addEventListener('keydown',function(e){
    const tag=(e.target&&e.target.tagName)||'';
    if(tag==='SELECT'||tag==='INPUT'||tag==='TEXTAREA') return;
    if(e.code==='Space'){ e.preventDefault(); togglePlay(); }
    else if(e.code==='ArrowRight'){ setPlay(false); seek(fpos+(e.shiftKey?Math.round(N*0.05):1)); }
    else if(e.code==='ArrowLeft'){ setPlay(false); seek(fpos-(e.shiftKey?Math.round(N*0.05):1)); }
    else if(e.code==='Home'){ setPlay(false); seek(0); }
    else if(e.code==='End'){ setPlay(false); seek(N-1); }
  });

  // timeline scrubbing
  let scrubbing=false;
  function tlSeek(ev){
    const r=tlc.getBoundingClientRect();
    const x=(ev.touches?ev.touches[0].clientX:ev.clientX)-r.left;
    seek(Math.round(x/r.width*(N-1)));
  }
  const tlwrap=document.getElementById('tlwrap');
  tlwrap.addEventListener('mousedown',e=>{scrubbing=true; setPlay(false); tlSeek(e);});
  window.addEventListener('mousemove',e=>{ if(scrubbing) tlSeek(e); });
  window.addEventListener('mouseup',()=>{scrubbing=false;});
  tlwrap.addEventListener('touchstart',e=>{scrubbing=true; setPlay(false); tlSeek(e); e.preventDefault();},{passive:false});
  tlwrap.addEventListener('touchmove',e=>{ if(scrubbing){tlSeek(e); e.preventDefault();} },{passive:false});
  tlwrap.addEventListener('touchend',()=>{scrubbing=false;});

  // ---- header + init ----
  document.getElementById('hudName').textContent=META.display;
  document.getElementById('hudSub').textContent=
    (META.is_case?('PDB '+META.source_pdb+' · '):'')+'ligand '+META.ligand_name+' · '+META.family;

  // mode badge: COMPUTED by the real engine vs (legacy) SIMULATED
  (function(){
    const badge=document.getElementById('modeBadge'), chip=document.getElementById('hudRmsd');
    if(META.computed){
      badge.textContent='COMPUTED'; badge.classList.add('ok');
      badge.title='Pose computed by the rigid-body docking engine (simplified Vina-style scoring, '
        +'basin-hopping search'+(META.exhaustiveness?(' ×'+META.exhaustiveness):'')+'). '
        +'The endpoint is the engine’s own lowest-energy pose, not a given answer.';
    } else {
      badge.textContent='SIMULATED';
      badge.title='Synthetic, physically-motivated trajectory for teaching.';
    }
    if(META.rmsd!=null){
      const r=META.rmsd, ok=r<2.0, near=r<3.0;
      chip.style.display='inline-block';
      chip.textContent='RMSD '+r.toFixed(2)+' Å '+(ok?'✓':(near?'≈':'✗'))+' vs crystal';
      chip.style.color=ok?'#34d399':(near?'#f59e0b':'#ef4444');
      chip.style.borderColor=ok?'#13513a':(near?'#5b430c':'#5b1414');
      chip.title='Root-mean-square deviation of the computed pose from the experimental '
        +'crystal structure. < 2 Å is the standard re-docking success criterion.';
    }
  })();
  document.getElementById('fTot').textContent=N-1;
  document.getElementById('fNow').textContent='0';

  function resizeAll(){ try{viewer.resize();}catch(e){} sizeSpark(); sizeTL(); drawSpark(cur); drawTL(cur); }
  window.addEventListener('resize',resizeAll);
  sizeSpark(); sizeTL();

  // restore prior position/speed so changing display options never loses your place
  const savedState=loadPlayState();
  if(savedState){
    speed=savedState.speed||1;
    document.getElementById('speed').value=String(speed);
    fpos=Math.max(0,Math.min(N-1,(savedState.frac||0)*(N-1)));
    renderFrame(fpos);
  } else {
    renderFrame(0);
  }
  setMode('frame');
  requestAnimationFrame(tick);

  // auto-play only on the first visit to this target; otherwise restore play state
  if(savedState){ if(savedState.playing) { dir=1; setPlay(true); } }
  else { setTimeout(()=>{ dir=1; setPlay(true); }, 900); }
})();
</script>
</body>
</html>
"""
