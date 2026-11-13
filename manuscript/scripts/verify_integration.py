"""Targeted main-result asset contract; not whole-manuscript scientific QA.
Runs two deterministic asset builds; writes only the six declared main assets.
Prints its report to stdout, leaving historical integration/build logs untouched.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import subprocess
import sys

sys.dont_write_bytecode = True
from build_assets import ROOT, MAN, OUTPUTS, SOURCES, RULES, SCENARIOS, POPS, METRICS, GRID, load_authority

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def inventory():
    return {p.relative_to(ROOT).as_posix(): digest(p) for p in ROOT.rglob("*")
            if p.is_file() and p not in OUTPUTS}

def main():
    before = inventory()
    data = load_authority()
    checks = {}
    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            raise AssertionError(name)
    snapshots = []
    for number in (1, 2):
        proc = subprocess.run([sys.executable, "-B", str(MAN/"scripts/build_assets.py")],
                              cwd=ROOT, capture_output=True, text=True)
        check(f"build_{number}_succeeded", proc.returncode == 0)
        snapshots.append({p.relative_to(ROOT).as_posix(): digest(p) for p in OUTPUTS})
    check("two_consecutive_builds_byte_identical", snapshots[0] == snapshots[1])
    check("all_other_repository_files_unchanged", before == inventory())
    t1 = OUTPUTS[-2].read_text(encoding="utf-8")
    t2 = OUTPUTS[-1].read_text(encoding="utf-8")
    rows1 = []
    for line in t1.splitlines():
        if line.startswith(("Primary:", "Secondary:")):
            fields = [x.strip() for x in line.removesuffix(r" \\").split("&")]
            rows1.append(fields)
    check("four_medium_classification_rows", len(rows1) == 4)
    for rule, label in zip(RULES, ("Primary: two-sided", "Secondary: directional")):
        for pop in POPS:
            name = "EU Member States" if pop == POPS[0] else r"EU27\_2020"
            rows = [r for r in rows1 if r[:2] == [label, name]]
            check(f"unique_medium_row_{rule}_{pop}", len(rows) == 1)
            original, counts, ua, ub = data["medium"][rule, pop]
            check(f"medium_numbers_{rule}_{pop}", tuple(map(int, rows[0][2:])) == counts+(ua, ub))
    widths_a = [float(r["A_ci_high"])-float(r["A_ci_low"]) for r in data["medium"][RULES[0],POPS[0]][0]]
    widths_b = [float(r["B1_ci_high"])-float(r["B1_ci_low"]) for r in data["medium"][RULES[0],POPS[0]][0]]
    check("all_27_B1_intervals_wider_independent_recount", sum(b>a for a,b in zip(widths_a,widths_b)) == 27)
    a,b,ratio,wider = data["precision"][RULES[0],POPS[0]]
    check("precision_table_display_matches_authority", all(x in t1 for x in
        (f"Median A interval width: {a:.4f}", f"median B1 width: {b:.4f}",
         f"Median paired B1/A width ratio: {ratio:.3f}", f"B1 wider: {wider}/27")))
    from build_assets import LABELS
    rows2=[]
    for line in t2.splitlines():
        if line.startswith(("Member States &", r"EU27\_2020 &")):
            rows2.append([x.strip() for x in line.removesuffix(r" \\").split("&")])
    check("twelve_short_window_rows", len(rows2) == 12)
    for pop in POPS:
        for metric in METRICS:
            for scenario, label in zip(SCENARIOS, ("U","B")):
                name = "Member States" if pop == POPS[0] else r"EU27\_2020"
                rr = [r for r in rows2 if r[:3] == [name, LABELS[metric], label]]
                check(f"unique_short_row_{pop}_{metric}_{label}", len(rr) == 1)
                frozen = data["short"][scenario,pop,metric,.95]
                check(f"short_numbers_{pop}_{metric}_{label}", tuple(map(int,rr[0][3:])) ==
                    tuple(int(frozen[k]) for k in ("agreement","resolved_class_disagreement","unresolved","total")))
    check("all_figures_nonempty", all(p.is_file() and p.stat().st_size>1000 for p in OUTPUTS[:4]))
    # Inspect the actual plotted x/y arrays independently against serialized N06 rows.
    import build_assets
    captured = []
    save = build_assets.savefig
    def capture(fig, name):
        captured.append((fig,name))
    build_assets.savefig = capture
    try:
        build_assets.draw_fig1(data)
        build_assets.draw_fig2(data)
    finally:
        build_assets.savefig = save
    fig1,fig2 = captured[0][0],captured[1][0]
    check("fig1_two_panels", len(fig1.axes) == 2)
    check("fig1_primary_class_group_sizes", [len(c.get_offsets()) for c in fig1.axes[0].collections] == [12,15,1])
    check("fig1_precision_27_pairs", len(fig1.axes[1].lines) == 27 and len(fig1.axes[1].collections) == 54)
    for line in fig1.axes[1].lines:
        check(f"fig1_pair_{line.get_ydata()[0]}_B1_wider", line.get_xdata()[1] > line.get_xdata()[0])
    check("fig2_three_panels",len(fig2.axes)==3)
    for ax,metric in zip(fig2.axes,METRICS):
        check("fig2_exact_ticks_"+metric, tuple(ax.get_xticks())==GRID)
        check("fig2_points_only_"+metric, len(ax.lines)==1)  # only 50% reference line
        for collection,scenario in zip(ax.collections,SCENARIOS):
            offsets=collection.get_offsets()
            check("fig2_exact_x_"+metric+scenario,tuple(float(x) for x in offsets[:,0])==GRID)
            expected=[100*int(data["short"][scenario,POPS[0],metric,t]["unresolved"])/int(data["short"][scenario,POPS[0],metric,t]["total"]) for t in GRID]
            check("fig2_exact_y_"+metric+scenario,all(abs(float(y)-z)<1e-12 for y,z in zip(offsets[:,1],expected)))
    build_assets.plt.close(fig1);build_assets.plt.close(fig2)
    check("verification_plot_inspection_writes_nothing", before == inventory())
    print(json.dumps(dict(scope="P3-C3A main-result assets only",checks=checks,all_asset_checks_pass=True,
        sources=data["source_sha256"],output_sha256=snapshots[-1],
        scientific_validation_rerun=False,whole_manuscript_QA=False,
        staged_representation_note="Protected Results prose/captions and supplements still await separate alignment."),indent=2))

if __name__ == "__main__":
    main()
