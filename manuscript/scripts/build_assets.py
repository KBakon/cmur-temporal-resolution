"""Main assets from frozen N01/N05/N06/N07; no scientific recomputation."""
from pathlib import Path
from statistics import median
from collections import Counter
import csv, hashlib, json, re, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2]
MAN=ROOT/"manuscript"; FIG=MAN/"figures"; TAB=MAN/"tables"
SOURCES=dict(N01=ROOT/"analysis/wp_p1a_temporal_support_diagnostics/medium_term_harmonized_rules.csv",
 N05=ROOT/"analysis/wp_p1c_calibration_robustness/medium_term_revision_candidate.csv",
 N06=ROOT/"analysis/wp_p1c_calibration_robustness/threshold_survival.csv",
 N07=ROOT/"analysis/wp_p1c_calibration_robustness/threshold_claim_classifications.csv")
RULES=("COMMON_TWO_SIDED_95","COMMON_DIRECTIONAL_95")
SCENARIOS=("UNSCALED_STRESS","BROWNIAN_SD_1_SQRT2")
POPS=("Member States","EU27_2020")
METRICS=("annual_change","turning","acceleration")
GRID=(.80,.85,.90,.925,.95,.975,.99)
LABELS=dict(annual_change="Annual change",turning="Turning",acceleration="Acceleration")
COLORS=("#1766A3","#C76B17")
OUTPUTS=[FIG/f"{n}.{e}" for n in ("fig1_medium_term_slopes","fig2_short_window_support") for e in ("pdf","png")]
OUTPUTS += [TAB/"table1_medium_term_summary.tex",TAB/"table2_short_window_summary.tex"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
 "axes.spines.right":False,"pdf.fonttype":42,"savefig.dpi":240})
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok,msg):
 if not ok: raise ValueError(msg)
def read(p):
 with p.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def counts(rr):
 c=Counter(r["stability_class"] for r in rr)
 return tuple(c[k] for k in ("AGREEMENT","CHANGE","UNRESOLVED"))
def load_authority():
 freeze=(ROOT/"artifacts/SCIENTIFIC_FREEZE_HANDOFF.md").read_text(encoding="utf-8")
 hashes={}
 for name,p in SOURCES.items():
  # Parse the exact indexed table row, including path and SHA-256.
  line=next((l for l in freeze.splitlines() if l.startswith("| "+name+" |")),None)
  require(line is not None,f"Missing freeze index {name}")
  fields=line.split("|"); indexed=fields[2].strip().strip(chr(96)); h=fields[-2].strip().strip(chr(96))
  require(indexed==p.relative_to(ROOT).as_posix(),f"{name} path mismatch")
  require(p.is_file() and digest(p)==h,f"SCIENTIFIC_AUTHORITY_CONFLICT: {name} SHA-256")
  hashes[name]=h
 n01,n05,n06,n07=[read(SOURCES[k]) for k in ("N01","N05","N06","N07")]
 require(len(n01)==56 and len({(r["rule"],r["geo"]) for r in n01})==56,"N01 complete unique grid")
 geos={r["geo"] for r in n01 if r["population"]==POPS[0]}
 require(len(geos)==27 and "EU27_2020" not in geos,"MS population")
 require(len(n05)==4,"N05 complete grid")
 medium={};precision={}
 for rule,lock in zip(RULES,((12,0,15),(15,0,12))):
  for pop in POPS:
   rr=[r for r in n01 if r["rule"]==rule and r["population"]==pop]
   require({r["geo"] for r in rr}==(geos if pop==POPS[0] else {"EU27_2020"}),"N01 missing geography")
   cc=counts(rr);require(cc==(lock if pop==POPS[0] else (1,0,0)),"Medium lock")
   pp=[r for r in n05 if r["rule"]==rule and r["population"]==pop]
   require(len(pp)==1,"N05 missing/duplicate"); pp=pp[0]
   require(cc==tuple(int(pp[k]) for k in ("agreement","resolved_class_disagreement","unresolved")),"N01/N05 classes")
   ua=sum(r["A_class"]=="UNRESOLVED" for r in rr);ub=sum(r["B1_class"]=="UNRESOLVED" for r in rr)
   require((ua,ub)==(int(pp["A_unresolved"]),int(pp["B1_unresolved"])),"Branch unresolved")
   wa=[float(r["A_ci_high"])-float(r["A_ci_low"]) for r in rr]
   wb=[float(r["B1_ci_high"])-float(r["B1_ci_low"]) for r in rr]
   stats=(median(wa),median(wb),median(b/a for a,b in zip(wa,wb)),sum(b>a for a,b in zip(wa,wb)))
   for x,k in zip(stats[:3],("median_A_CI_width","median_B1_CI_width","median_paired_B1_A_CI_width_ratio")):
    require(math.isclose(x,float(pp[k]),rel_tol=0,abs_tol=1e-12),"N05 precision")
   medium[rule,pop]=(rr,cc,ua,ub);precision[rule,pop]=stats
 require(medium[RULES[0],POPS[0]][2:]==(8,15),"Primary branch lock")
 stats=precision[RULES[0],POPS[0]]
 require(tuple(round(x,d) for x,d in zip(stats[:3],(4,4,3)))==(.3386,.6219,1.858) and stats[3]==27,"27/27 precision lock")
 require(len(n06)==84 and len(n07)==13328,"Short-window row grid")
 grouped={}
 for r in n07:grouped.setdefault((r["scenario"],r["population"],r["metric"],float(r["threshold"])),[]).append(r)
 locks={SCENARIOS[0]:((10,0,314),(2,0,295),(1,0,296)),SCENARIOS[1]:((29,0,295),(13,0,284),(7,0,290))}
 short={}
 for r in n06:
  key=(r["scenario"],r["population"],r["metric"],float(r["threshold"]))
  require(key not in short,"N06 duplicate")
  s,p,m,t=key;require(s in SCENARIOS and p in POPS and m in METRICS and t in GRID,"Unexpected N06 key")
  rr=grouped.get(key,[]);n=(324 if m==METRICS[0] else 297) if p==POPS[0] else (12 if m==METRICS[0] else 11)
  require(len(rr)==n==int(r["total"]) and len({(x["geo"],x["year"],x["metric"]) for x in rr})==n,"N07 complete claims")
  cc=tuple(int(r[k]) for k in ("agreement","resolved_class_disagreement","unresolved"))
  require(counts(rr)==cc and sum(cc)==n,"N06/N07 reconciliation")
  require(math.isclose(cc[2]/n,float(r["unresolved_proportion"]),abs_tol=1e-15),"N06 share")
  if p==POPS[0]:require(cc[2]>n/2,"Majority unresolved lock")
  if t==.95:require(cc==(locks[s][METRICS.index(m)] if p==POPS[0] else (0,0,n)),"Reference count lock")
  short[key]=r
 require(set(short)=={(s,p,m,t) for s in SCENARIOS for p in POPS for m in METRICS for t in GRID},"Missing threshold group")
 adverse={(r["scenario"],r["geo"],r["year"],r["metric"],float(r["threshold"])) for r in n07 if r["stability_class"]=="CHANGE"}
 require(adverse=={(SCENARIOS[1],"BE","2019",METRICS[0],t) for t in (.8,.85)},"Adverse-case lock")
 return dict(medium=medium,precision=precision,short=short,source_sha256=hashes)
def savefig(fig,name):
 fig.savefig(FIG/f"{name}.pdf",bbox_inches="tight",metadata={"CreationDate":None,"ModDate":None,"Creator":"CMUR deterministic asset builder"})
 fig.savefig(FIG/f"{name}.png",bbox_inches="tight",metadata={"Software":"CMUR deterministic asset builder"})
 plt.close(fig)
def draw_fig1(d):
 rr=d["medium"][RULES[0],POPS[0]][0];eu=d["medium"][RULES[0],POPS[1]][0][0]
 fig,(ax,bx)=plt.subplots(1,2,figsize=(11.3,6.6),gridspec_kw={"width_ratios":[1.08,1]},layout="constrained")
 for cls,color in zip(("AGREEMENT","UNRESOLVED"),COLORS):
  rows=[r for r in rr if r["stability_class"]==cls]
  ax.scatter([float(r["A_slope_pp_per_year"]) for r in rows],[float(r["B1_slope_pp_per_year"]) for r in rows],
   s=36,c=color,edgecolors="white",linewidths=.5,label=f"MS {cls.lower()} (n={len(rows)})",zorder=3)
 ax.scatter(float(eu["A_slope_pp_per_year"]),float(eu["B1_slope_pp_per_year"]),s=65,marker="s",c=COLORS[0],edgecolors="#222222",label="EU27_2020 (agreement)",zorder=4)
 vals=[float(r[k]) for r in rr+[eu] for k in ("A_slope_pp_per_year","B1_slope_pp_per_year")]
 lo,hi=min(vals),max(vals);pad=(hi-lo)*.1
 ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],c="#999999",ls="--",lw=1)
 ax.axhline(0,c="#dddddd",lw=.8);ax.axvline(0,c="#dddddd",lw=.8)
 ax.set(xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad),xlabel="A: released annual slope (pp/year)",ylabel="B1: source-anchor slope (pp/year)",title="A   Direction: primary common two-sided 95%")
 ax.set_aspect("equal",adjustable="box");ax.legend(frameon=False,fontsize=9,loc="upper left")
 ordered=sorted(rr,key=lambda r:(float(r["B1_ci_high"])-float(r["B1_ci_low"]),r["geo"]))
 for i,r in enumerate(ordered):
  a=float(r["A_ci_high"])-float(r["A_ci_low"]);b=float(r["B1_ci_high"])-float(r["B1_ci_low"])
  bx.plot([a,b],[i,i],c="#bbbbbb",lw=1.2)
  bx.scatter(a,i,c=COLORS[0],s=24,marker="o",zorder=3);bx.scatter(b,i,c=COLORS[1],s=24,marker="D",zorder=3)
 bx.set_yticks(range(27),[r["geo"] for r in ordered],fontsize=8)
 bx.set(xlabel="Classical 95% OLS interval width (pp/year)",title=f"B   Model-based precision: B1 wider in {d['precision'][RULES[0],POPS[0]][3]}/27",ylim=(-1,27),xlim=(0,None))
 bx.grid(axis="x",alpha=.2)
 bx.legend(handles=[Line2D([],[],marker="o",ls="",c=COLORS[0],label="A: released annual"),Line2D([],[],marker="D",ls="",c=COLORS[1],label="B1: source anchors")],frameon=False,fontsize=9,loc="lower right")
 fig.suptitle("Medium-term CMUR direction and classical OLS precision, 2010–2022",fontsize=13,fontweight="bold")
 fig.supxlabel("Unresolved denotes an unmet support criterion; interval widths are conditional on the classical OLS model.",fontsize=9)
 savefig(fig,"fig1_medium_term_slopes")
def draw_fig2(d):
 fig,axes=plt.subplots(1,3,figsize=(11.8,4.5),sharey=True,layout="constrained")
 for ax,m,panel in zip(axes,METRICS,"ABC"):
  for s,color,marker in zip(SCENARIOS,COLORS,("o","s")):
   ys=[100*int(d["short"][s,POPS[0],m,t]["unresolved"])/int(d["short"][s,POPS[0],m,t]["total"]) for t in GRID]
   ax.scatter(GRID,ys,s=50 if marker=="o" else 27,marker=marker,facecolors="none" if marker=="o" else color,edgecolors=color,linewidths=1.5,label="U: unscaled stress" if s==SCENARIOS[0] else "B: Brownian amplitude",zorder=3)
  ax.set_xticks(GRID,[str(t) for t in GRID],rotation=55,ha="right",fontsize=9)
  ax.set(xlim=(.792,1.0),ylim=(0,104),xlabel="Sign-support threshold",title=f"{panel}   {LABELS[m]}")
  ax.axhline(50,c="#777777",lw=.9,ls="--");ax.grid(axis="y",alpha=.2)
 axes[0].set_ylabel("Member-State unresolved share (%)")
 axes[1].legend(frameon=False,loc="lower center",fontsize=9)
 fig.suptitle("Short-window CMUR support through 2022 across the tested threshold grid",fontsize=13,fontweight="bold")
 fig.supxlabel("Points only: no interpolation. Dashed line marks 50% unresolved; shares describe classifications, not probabilities of truth.",fontsize=9)
 savefig(fig,"fig2_short_window_support")
def table_texts(d):
 end=r" \\"
 lines=[r"\begin{table}[htbp]",r"\centering",r"\small",
 r"\caption{Medium-term CMUR comparison, 2010--2022: primary common two-sided 95\% and secondary common directional 95\% rules. A is the released annual series and B1 the admissible source-anchor series. Both rules use the same fixed fits. A/D/U denotes agreement/resolved-class disagreement/unresolved. Precision summaries concern Member States and classical model-based 95\% OLS intervals.}",
 r"\label{tab:mediumterm}",r"\begin{tabular}{llrrrrr}",r"\toprule",
 r"Rule & Population & A & D & U & \shortstack{A\\unresolved} & \shortstack{B1\\unresolved}"+end,r"\midrule",r"\multicolumn{7}{l}{\textit{Part A: direction classification}}"+end]
 for rule,label in zip(RULES,("Primary: two-sided","Secondary: directional")):
  for pop in POPS:
   rr,cc,ua,ub=d["medium"][rule,pop];p="EU Member States" if pop==POPS[0] else r"EU27\_2020"
   lines.append(f"{label} & {p} & {cc[0]} & {cc[1]} & {cc[2]} & {ua} & {ub}"+end)
 a,b,ratio,wider=d["precision"][RULES[0],POPS[0]]
 lines += [r"\midrule",r"\multicolumn{7}{l}{\textit{Part B: Member-State classical OLS precision}}"+end,
 rf"\multicolumn{{7}}{{l}}{{Median A interval width: {a:.4f} pp/year; median B1 width: {b:.4f} pp/year.}}"+end,
 rf"\multicolumn{{7}}{{l}}{{Median paired B1/A width ratio: {ratio:.3f}; B1 wider: {wider}/27 Member States.}}"+end,
 r"\bottomrule",r"\end{tabular}",r"\end{table}",""]
 one="\n".join(lines)
 lines=[r"\begin{table}[htbp]",r"\centering",r"\small",r"\setlength{\tabcolsep}{4pt}",
 r"\caption{Source-aware short-window comparisons through 2022 at $\tau=0.95$, including the deterministic vote and minimum-three unanimity rule. U is the primary unscaled empirical stress scenario; B is the secondary Brownian-motivated $1/\sqrt{2}$ residual-amplitude sensitivity.}",
 r"\label{tab:shortwindow}",r"\begin{tabular}{lllrrrr}",r"\toprule",
 r"Population & Metric & Scenario & \shortstack{Agree-\\ment} & \shortstack{Resolved-class\\disagreement} & Unresolved & Total"+end,r"\midrule"]
 for pop in POPS:
  for m in METRICS:
   for s,label in zip(SCENARIOS,("U","B")):
    r=d["short"][s,pop,m,.95];p="Member States" if pop==POPS[0] else r"EU27\_2020"
    lines.append(f"{p} & {LABELS[m]} & {label} & {r['agreement']} & {r['resolved_class_disagreement']} & {r['unresolved']} & {r['total']}"+end)
  if pop==POPS[0]:lines.append(r"\midrule")
 lines += [r"\bottomrule",r"\end{tabular}",r"\end{table}",""]
 return one,"\n".join(lines)
def main():
 d=load_authority();FIG.mkdir(exist_ok=True);TAB.mkdir(exist_ok=True)
 draw_fig1(d);draw_fig2(d)
 for p,t in zip(OUTPUTS[-2:],table_texts(d)):p.write_text(t,encoding="utf-8",newline="\n")
 require(all(digest(SOURCES[k])==h for k,h in d["source_sha256"].items()),"Scientific inputs changed")
 print(json.dumps(dict(source_sha256=d["source_sha256"],output_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in OUTPUTS},primary_MS=d["medium"][RULES[0],POPS[0]][1],secondary_MS=d["medium"][RULES[1],POPS[0]][1],precision_MS=d["precision"][RULES[0],POPS[0]],short_groups_checked=len(d["short"]),supplementary_files_written=False),indent=2))
if __name__=="__main__":main()
