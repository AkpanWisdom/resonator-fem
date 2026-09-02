import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

x, m1, m2, m3 = [], [], [], []
with open("modes.csv") as f:
    for row in csv.DictReader(f):
        x.append(float(row["x"]))
        m1.append(float(row["mode1"]))
        m2.append(float(row["mode2"]))
        m3.append(float(row["mode3"]))

plt.figure(figsize=(7, 4))
plt.plot(x, m1, label="mode 1")
plt.plot(x, m2, label="mode 2")
plt.plot(x, m3, label="mode 3")
plt.axhline(0, color="k", lw=0.5)
plt.xlabel("x / L")
plt.ylabel("normalized amplitude")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("docs/figures/mode_shapes.png", dpi=150)
print("wrote docs/figures/mode_shapes.png")
