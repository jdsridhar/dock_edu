#!/usr/bin/env python3
"""
Precomputes trajectories for the 4 case studies and builds static HTML files
suitable for direct hosting on GitHub Pages, along with a landing page.
"""
import os
import sys

# Ensure visualizer root is on the path
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

from docking import structures, trajectory, knowledge, viewer


def main():
    # Support printing unicode on Windows terminals
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
        
    print("Starting static site build...")
    
    # 1. Output directory
    out_dir = os.path.join(BASE, "static")
    os.makedirs(out_dir, exist_ok=True)
    
    # Cases to export
    cases = structures.case_list()
    print(f"Found {len(cases)} case studies: {[c[0] for c in cases]}")
    
    for slug, disp, pdb, fam in cases:
        print(f"Exporting case study: {disp} ({slug})...")
        record = structures.load_case_study(slug)
        beats = knowledge.get_beats(record)
        
        # Dock ligand (seed=7, frames=1000)
        print("  Docking ligand (this may take a few seconds)...")
        traj = trajectory.generate_trajectory(record, n_frames=1000, seed=7, beats=beats)
        
        pack = knowledge.get_pack(record)
        tutor = knowledge.get_tutor()
        
        # Build viewer workspace HTML
        options = {
            "representation": "cartoon", 
            "color_scheme": "spectrum",
            "opacity": 0.85, 
            "show_surface": False, 
            "surface_opacity": 0.55
        }
        # Force height to occupy full screen/iframe beautifully
        html_content = viewer.build_workspace_html(record, traj, pack, tutor, options, height=940)
        
        # Save HTML file
        out_path = os.path.join(out_dir, f"{slug}.html")
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(html_content)
        print(f"  Successfully wrote {out_path} ({len(html_content)} bytes)")

    print("Static case study pages built successfully!")


if __name__ == "__main__":
    main()
