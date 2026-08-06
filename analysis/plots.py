# Import required libraries
import pandas as pd
import matplotlib.pyplot as plt

# Read in summmary results
df = pd.read_csv("results/scoring/summary_metrics.csv")

pivot = df.pivot(index="task",
                 columns="prompt",
                 values="Recall")

# Plot recall
pivot.plot(kind="bar")

plt.ylabel("Recall")
plt.tight_layout()
plt.show()             


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------
# Load data
# -----------------------
df = pd.read_csv("results/scoring/summary_metrics.csv")

# Keep only code-review tasks
code_tasks = ["bulkRNA", "seurat_qc", "ml_disease"]
code = df[df["task"].isin(code_tasks)]

# Order prompts
prompt_order = ["baseline", "evidence_based", "blind"]
code["prompt"] = pd.Categorical(
    code["prompt"],
    categories=prompt_order,
    ordered=True
)

pivot = (
    code.pivot(index="task", columns="prompt", values="Recall")
    .loc[code_tasks]
)

# -----------------------
# Plot styling
# -----------------------
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 18,
    "axes.labelsize": 14,
    "figure.facecolor": "white"
})

colors = ["#4E79A7", "#F28E2B", "#59A14F"]

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(pivot.index))
width = 0.24

for i, prompt in enumerate(prompt_order):
    bars = ax.bar(
        x + (i-1)*width,
        pivot[prompt],
        width,
        label=prompt.replace("_", " ").title(),
        color=colors[i]
    )

    # Value labels
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,
            h + 0.02,
            f"{h:.2f}",
            ha="center",
            fontsize=10
        )

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(["Bulk RNA", "Seurat QC", "ML Disease"])
ax.set_ylim(0, 1.1)

ax.set_ylabel("Recall")
ax.set_xlabel("Benchmark Task")
ax.set_title(
    "Prompt Engineering Significantly Affects Error Detection Recall",
    pad=20,
    weight="bold"
)

ax.grid(axis="y", linestyle="--", alpha=0.35)
ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(title="Prompt Strategy", frameon=False)

plt.tight_layout()

plt.savefig(
    "results/figures/figure1_recall_by_prompt.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()