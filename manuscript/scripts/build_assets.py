"""Generate manuscript figures and LaTeX tables from frozen result CSVs.

Requires reportlab for vector PDFs; optionally uses pdftoppm for PNG previews.
No scientific model or classification is recomputed here.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor

ROOT = Path(__file__).resolve().parents[2]
MAN = ROOT / "manuscript"
DATA = ROOT / "data" / "processed"
FIG = MAN / "figures"
TAB = MAN / "tables"
SUP = MAN / "supplementary"
BLUE = HexColor("#1766a3")
ORANGE = HexColor("#d87216")
GRAY = HexColor("#666666")
LIGHT = HexColor("#d8d8d8")


def read(name):
    with (DATA / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def text_file(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8", newline="\n")


def preview(pdf):
    tool = shutil.which("pdftoppm")
    if tool:
        subprocess.run([tool, "-f", "1", "-singlefile", "-r", "220", "-png", str(pdf), str(pdf.with_suffix(""))], check=True)


def draw_fig1(A, B, C):
    b = {r["geo"]: r for r in B}
    cls = {r["geo"]: r["stability_class"] for r in C if r["metric"] == "medium_term_slope"}
    assert set(b) == {r["geo"] for r in A} and len(A) == 28 and len(cls) == 28
    points = [(r["geo"], float(r["slope_pp_per_year"]), float(b[r["geo"]]["slope_pp_per_year"]), cls[r["geo"]]) for r in A]
    values = [v for _, x, y, _ in points for v in (x, y)]
    lo, hi = min(values), max(values)
    pad = (hi-lo)*.09
    lo -= pad; hi += pad
    pdf = FIG / "fig1_medium_term_slopes.pdf"
    w, h = 540, 455
    c = canvas.Canvas(str(pdf), pagesize=(w,h), invariant=1)
    left, bottom, width, height = 68, 64, 440, 320
    xmap = lambda v: left + (v-lo)/(hi-lo)*width
    ymap = lambda v: bottom + (v-lo)/(hi-lo)*height
    c.setFont("Helvetica-Bold", 12)
    c.drawCentredString(w/2, h-23, "Medium-term CMUR slopes, 2010-2022")
    c.setStrokeColor(LIGHT); c.setLineWidth(.5)
    for v in (-1,-.5,0,.5,1):
        if lo<=v<=hi:
            c.line(xmap(v),bottom,xmap(v),bottom+height)
            c.line(left,ymap(v),left+width,ymap(v))
            c.setFont("Helvetica",8); c.setFillColor(GRAY)
            c.drawCentredString(xmap(v),bottom-14,f"{v:g}")
            c.drawRightString(left-7,ymap(v)-3,f"{v:g}")
    c.setStrokeColor(GRAY); c.setLineWidth(.8)
    c.rect(left,bottom,width,height)
    c.setStrokeColor(BLUE); c.setLineWidth(1.2)
    c.line(xmap(lo),ymap(lo),xmap(hi),ymap(hi))
    for geo,x,y,category in points:
        c.setFillColor(BLUE if category=="AGREEMENT" else ORANGE)
        if geo=="EU27_2020":
            c.setStrokeColor(HexColor("#222222")); c.setLineWidth(1.4)
            c.rect(xmap(x)-4,ymap(y)-4,8,8,fill=1,stroke=1)
        else:
            c.circle(xmap(x),ymap(y),3.3,fill=1,stroke=0)
    c.setFillColor(HexColor("#222222")); c.setFont("Helvetica",10)
    c.drawCentredString(left+width/2,33,"Released annual slope (percentage points/year)")
    c.saveState();c.translate(19,bottom+height/2);c.rotate(90)
    c.drawCentredString(0,0,"Anchor-only slope (percentage points/year)");c.restoreState()
    c.setFont("Helvetica",8.5)
    c.setFillColor(BLUE);c.circle(80,h-43,3,fill=1,stroke=0)
    c.setFillColor(HexColor("#222222"));c.drawString(89,h-46,"Agreement")
    c.setFillColor(ORANGE);c.circle(171,h-43,3,fill=1,stroke=0)
    c.setFillColor(HexColor("#222222"));c.drawString(180,h-46,"Unresolved")
    c.save();preview(pdf)


def draw_fig2(C):
    labels = [("annual_change","Annual change"),("turning","Turning"),("acceleration","Acceleration")]
    summary = []
    for code,label in labels:
        rows=[r for r in C if r["metric"]==code and r["geo"]!="EU27_2020" and int(r["year"])<=2022]
        count=Counter(r["stability_class"] for r in rows)
        assert count["CHANGE"]==0 and count["AGREEMENT"]+count["UNRESOLVED"]==len(rows)
        summary.append((label,count["AGREEMENT"],count["UNRESOLVED"],len(rows)))
    pdf=FIG/"fig2_short_window_support.pdf"
    w,h=555,284
    c=canvas.Canvas(str(pdf),pagesize=(w,h),invariant=1)
    c.setFont("Helvetica-Bold",12)
    c.drawCentredString(w/2,h-25,"Short-window support through 2022")
    c.setFont("Helvetica",10);c.setFillColor(BLUE);c.rect(90,h-51,10,10,fill=1,stroke=0)
    c.setFillColor(HexColor("#222222"));c.drawString(105,h-49,"Agreement")
    c.setFillColor(ORANGE);c.rect(194,h-51,10,10,fill=1,stroke=0)
    c.setFillColor(HexColor("#222222"));c.drawString(209,h-49,"Unresolved")
    bar_x,bar_w=135,320
    for i,(label,agree,unresolved,total) in enumerate(summary):
        y=h-99-i*65
        c.setFont("Helvetica",10);c.setFillColor(HexColor("#222222"))
        c.drawString(20,y+8,label)
        c.setFillColor(ORANGE);c.rect(bar_x,y,bar_w,20,fill=1,stroke=0)
        c.setFillColor(BLUE);c.rect(bar_x,y,bar_w*agree/total,20,fill=1,stroke=0)
        c.setStrokeColor(LIGHT);c.rect(bar_x,y,bar_w,20,fill=0,stroke=1)
        c.setFillColor(HexColor("#222222"));c.setFont("Helvetica",10)
        c.drawRightString(w-18,y+6,f"{agree} / {total} agree")
        c.drawString(bar_x,y-13,f"{unresolved} unresolved")
    c.save();preview(pdf)


def latex_tables(A,B,C):
    medium=[r for r in C if r["metric"]=="medium_term_slope"]
    ms=Counter(r["stability_class"] for r in medium if r["geo"]!="EU27_2020")
    eu=Counter(r["stability_class"] for r in medium if r["geo"]=="EU27_2020")
    assert sum(ms.values())==27 and sum(eu.values())==1
    table1=r"""\begin{table}[htbp]
\centering
\caption{Medium-term direction classification for 2010--2022. Member-State counts compare the released annual series with the reporting-year-anchor analysis.}
\label{tab:mediumterm}
\begin{tabular}{lrrrr}
\toprule
Population & Agreement & Resolved-class disagreement & Unresolved & Total \\
\midrule
"""+f"EU Member States & {ms['AGREEMENT']} & {ms['CHANGE']} & {ms['UNRESOLVED']} & 27 "+"\\"*2+"\n"+f"EU27\\_2020 & {eu['AGREEMENT']} & {eu['CHANGE']} & {eu['UNRESOLVED']} & 1 "+"\\"*2+"\n"+r"""\bottomrule
\end{tabular}
\end{table}
"""
    text_file(TAB/"table1_medium_term_summary.tex",table1)
    lines=[]
    for code,label in (("annual_change","Annual change"),("turning","Turning interpretation"),("acceleration","Acceleration")):
        rr=[r for r in C if r["metric"]==code and r["geo"]!="EU27_2020" and int(r["year"])<=2022]
        cc=Counter(r["stability_class"] for r in rr)
        assert cc["AGREEMENT"]+cc["CHANGE"]+cc["UNRESOLVED"]==len(rr)
        lines.append(f"{label} & {cc['AGREEMENT']} & {cc['CHANGE']} & {cc['UNRESOLVED']} & {len(rr)} "+"\\"*2)
    table2=r"""\begin{table}[htbp]
\centering
\caption{Short-window Member-State classifications through 2022 under the specified cross-variant consensus rule. Counts are geography-by-window claim instances.}
\label{tab:shortwindow}
\begin{tabular}{lrrrr}
\toprule
Metric & Agreement & Resolved-class disagreement & Unresolved & Total \\
\midrule
"""+"\n".join(lines)+"\n"+r"""\bottomrule
\end{tabular}
\end{table}
"""
    text_file(TAB/"table2_short_window_summary.tex",table2)
    b={r["geo"]:r for r in B}
    assert len(A)==len(b)==28 and {r["geo"] for r in A}==set(b)
    lines=[]
    for r in sorted(A,key=lambda x:x["geo"]):
        other=b[r["geo"]]
        x,y=float(r["slope_pp_per_year"]),float(other["slope_pp_per_year"])
        geo=r["geo"].replace("_",r"\_")
        lines.append(f"{geo} & {x:.4f} & {y:.4f} & {y-x:.4f} & {other['n']} & {r['direction']} & {other['direction']} "+"\\"*2)
    supp=r"""\begin{longtable}{lrrrrll}
\caption{Geography-level medium-term CMUR slope comparison, 2010--2022. A uses the released annual series; B1 uses admissible reporting-year anchors. Slopes and Difference are in percentage points per year. Difference is B1 minus A, calculated from unrounded estimates; B1 $n$ counts admissible anchors.}\label{tab:S1}\\
\toprule
Geo & A slope & B1 slope & Difference & B1 $n$ & A class & B1 class \\
\midrule
\endfirsthead
\toprule
Geo & A slope & B1 slope & Difference & B1 $n$ & A class & B1 class \\
\midrule
\endhead
"""+"\n".join(lines)+"\n"+r"""\bottomrule
\end{longtable}
"""
    text_file(SUP/"tableS1_country_slopes.tex",supp)
    interval_lines=[]
    for r in sorted(A,key=lambda x:x["geo"]):
        other=b[r["geo"]]
        geo=r["geo"].replace("_",r"\_")
        a_interval=f"[{float(r['ci_low']):.4f}, {float(r['ci_high']):.4f}]"
        b_interval=f"[{float(other['ci_low']):.4f}, {float(other['ci_high']):.4f}]"
        interval_lines.append(f"{geo} & {a_interval} & {b_interval} & {float(other['p_positive']):.6f} "+"\\"*2)
    supp_intervals=r"""\begin{longtable}{lccc}
\caption{Country-level medium-term slope uncertainty and sign support, 2010--2022. A and B1 are defined in Table~\ref{tab:S1}; intervals are standard two-sided 95\% Student-$t$ OLS intervals in percentage points per year. B1 $p_+$ is the classical one-sided Student-$t$ sign-support quantity defined in Section~\ref{sec:mediumterm}.}\label{tab:S2}\\
\toprule
Geo & A 95\% interval & B1 95\% interval & B1 $p_+$ \\
\midrule
\endfirsthead
\toprule
Geo & A 95\% interval & B1 95\% interval & B1 $p_+$ \\
\midrule
\endhead
"""+"\n".join(interval_lines)+"\n"+r"""\bottomrule
\end{longtable}
"""
    text_file(SUP/"tableS2_country_intervals.tex",supp_intervals)


def main():
    FIG.mkdir(parents=True,exist_ok=True)
    A=read("branchA_medium_term.csv")
    B=read("branchB1_anchor_medium_term.csv")
    C=read("temporal_stability_comparison.csv")
    draw_fig1(A,B,C)
    draw_fig2(C)
    latex_tables(A,B,C)
    print("Generated manuscript figures and tables from frozen result CSVs")


if __name__=="__main__":
    main()
