# Plots for the coupled environmental-modal runs.

import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data = {}
with open("history.csv") as f:
    for row in csv.DictReader(f):
        c = row["case"]
        if c not in data:
            data[c] = {"t": [], "T": [], "s": [], "f": []}
        data[c]["t"].append(float(row["t_hours"]))
        data[c]["T"].append(float(row["T"]))
        data[c]["s"].append(float(row["sigma"]))
        data[c]["f"].append(float(row["f1"]))

fig, ax = plt.subplots(3, 1, figsize=(7, 9), sharex=True)

for c in sorted(data):
    ax[0].plot(data[c]["t"], data[c]["T"], label=c)
    ax[1].plot(data[c]["t"], [s / 1e3 for s in data[c]["s"]], label=c)
    ax[2].plot(data[c]["t"], data[c]["f"], label=c)

ax[0].set_ylabel("temperature (K)")
ax[1].set_ylabel("axial stress (kPa)")
ax[2].set_ylabel("f1 (Hz)")
ax[2].set_xlabel("time (hours)")

for a in ax:
    a.grid(alpha=0.3)
ax[0].legend()

plt.tight_layout()
plt.savefig("docs/figures/path_dependence.png", dpi=150)
print("wrote docs/figures/path_dependence.png")
