import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import numpy as np
from sklearn.metrics import roc_curve, auc

sns.set_style("whitegrid")
sns.set_palette("dark")
COLORS  = {"nt100m": "#2E75B6", "nt500m": "#C00000"}
LABELS  = {"nt100m": "NT-100M", "nt500m": "NT-500M"}
os_import = __import__('os')
os_import.makedirs("figures/outputs", exist_ok=True)

# ── Figure 1: Position Sensitivity Heatmap ────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 4))
for ax, model in zip(axes, ["nt100m", "nt500m"]):
    df  = pd.read_csv(f"results/perturbation_{model}.csv")
    # Average flip rate per position across all sequences
    pss = df.groupby("mut_position")["flipped"].mean()
    ax.bar(pss.index, pss.values, color=COLORS[model], alpha=0.7, width=1.0)
    ax.axvline(x=200, color="red", linestyle="--", linewidth=1.5, label="GT/AG position")
    ax.axvline(x=201, color="red", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Position in sequence")
    ax.set_ylabel("Flip rate")
    ax.set_title(f"{LABELS[model]} — Position Sensitivity")
    ax.legend()
plt.tight_layout()
plt.savefig("figures/outputs/fig1_pss_heatmap.png", dpi=300, bbox_inches="tight")
print("Saved fig1_pss_heatmap.png")

# ── Figure 2: FFR vs SFR bar chart ───────────────────────────
with open("results/metrics_summary.json") as f:
    metrics = json.load(f)

models = [m for m in ["nt100m", "nt500m"] if m in metrics]
x      = np.arange(len(models))
width  = 0.25

fig, ax = plt.subplots(figsize=(9, 5))
ffrs = [metrics[m]["FFR"] for m in models]
esrs = [metrics[m]["ESR"] for m in models]
sfrs = [metrics[m]["SFR"] for m in models]
bsrs = [metrics[m]["BSR"] for m in models]

ax.bar(x - width, ffrs, width, label="FFR — Consensus (Type1)", color="#2E75B6")
ax.bar(x,         esrs, width, label="ESE region (Type2)",      color="#70AD47")
ax.bar(x + width, sfrs, width, label="SFR — Silent (Type3)",    color="#C00000")

for i, bsr in enumerate(bsrs):
    ax.text(i, max(ffrs[i], sfrs[i]) + 0.03,
            f"BSR={bsr:.1f}", ha="center", fontsize=11, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels([LABELS[m] for m in models])
ax.set_ylabel("Flip Rate")
ax.set_title("Functional vs Silent Flip Rate\n(higher FFR and lower SFR = better biological specificity)")
ax.legend()
ax.set_ylim(0, 1.15)
plt.tight_layout()
plt.savefig("figures/outputs/fig2_ffr_sfr.png", dpi=300, bbox_inches="tight")
print("Saved fig2_ffr_sfr.png")

# ── Figure 3: ROC curves ──────────────────────────────────────
cv  = pd.read_csv("data/processed/clinvar_gencode_merged.csv")
fig, ax = plt.subplots(figsize=(7, 6))

for model in models:
    df   = pd.read_csv(f"results/perturbation_{model}.csv")
    t1   = df[df["perturbation_type"] == "Type1"]
    merged = pd.merge(cv, t1, on=["chrom","pos"], how="inner")
    if len(merged) > 10:
        fpr, tpr, _ = roc_curve(merged["label"], merged["conf_change"])
        auroc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=COLORS[model],
                label=f"{LABELS[model]} (AUROC={auroc:.2f})", linewidth=2)

ax.plot([0,1],[0,1], "k--", label="Random (0.50)", linewidth=1)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ClinVar Pathogenicity Prediction\n(using model confidence change as score)")
ax.legend()
plt.tight_layout()
plt.savefig("figures/outputs/fig3_roc.png", dpi=300, bbox_inches="tight")
print("Saved fig3_roc.png")

# ── Figure 4: Confidence change by pathogenicity ──────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, model in zip(axes, models):
    df     = pd.read_csv(f"results/perturbation_{model}.csv")
    t1     = df[df["perturbation_type"] == "Type1"]
    merged = pd.merge(cv, t1, on=["chrom","pos"], how="inner")
    if len(merged) > 10:
        merged["Significance"] = merged["label"].map(
            {1: "Pathogenic", 0: "Benign"}
        )
        sns.violinplot(
            data=merged, x="Significance", y="conf_change",
            palette={"Pathogenic": "#C00000", "Benign": "#2E75B6"},
            ax=ax
        )
        ax.set_title(f"{LABELS[model]}\nConfidence Change by Pathogenicity")
        ax.set_ylabel("Confidence Change (baseline - perturbed)")
        ax.set_xlabel("")
        ax.axhline(y=0, color="black", linestyle="--", linewidth=1)

plt.tight_layout()
plt.savefig("figures/outputs/fig4_clinvar_violin.png", dpi=300, bbox_inches="tight")
print("Saved fig4_clinvar_violin.png")
print("\nAll figures saved to figures/outputs/")
