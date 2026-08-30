import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("rose-pine-matplotlib/themes/rose-pine-moon.mplstyle")
from pyobis import occurrences

query = occurrences.search(taxonid=125612, size=50_000)
df = query.execute()

col_map = {}
for col in df.columns:
    low = col.lower()
    if low == "decimallongitude":
        col_map[col] = "decimalLongitude"
    elif low == "decimallatitude":
        col_map[col] = "decimalLatitude"
    elif low == "species":
        col_map[col] = "species"
df = df.rename(columns=col_map)

if "species" not in df.columns:
    for candidate in ["scientificName", "scientificname", "taxonName", "taxonname"]:
        if candidate in df.columns:
            df = df.rename(columns={candidate: "species"})
            break

df = df.dropna(subset=["decimalLatitude", "species"])
df = df[df["decimalLatitude"].between(-90, 90)]
df["species"] = df["species"].str.strip()

BIN_WIDTH = 5
bins = np.arange(-90, 91, BIN_WIDTH)
mids = (bins[:-1] + bins[1:]) / 2

df["lat_bin"] = pd.cut(df["decimalLatitude"], bins=bins, labels=mids, include_lowest=True)

richness = (
    df.groupby("lat_bin", observed=True)["species"]
    .nunique()
    .reset_index()
    .rename(columns={"lat_bin": "lat", "species": "S"})
)
richness["lat"] = richness["lat"].astype(float)
richness = richness[richness["S"] > 0].sort_values("lat").reset_index(drop=True)

x = richness["lat"].values
y = richness["S"].values

X = np.column_stack([x**2, x, np.ones_like(x)])
beta = np.linalg.solve(X.T @ X, X.T @ y)
a, b, c = beta

y_pred = X @ beta
ss_res = np.sum((y - y_pred) ** 2)
ss_tot = np.sum((y - y.mean()) ** 2)
R2 = 1 - ss_res / ss_tot
rmse = np.sqrt(ss_res / len(y))

sigma2 = ss_res / (len(y) - 3)
cov = sigma2 * np.linalg.inv(X.T @ X)
se = np.sqrt(np.diag(cov))

x_smooth = np.linspace(-90, 90, 500)
X_smooth = np.column_stack([x_smooth**2, x_smooth, np.ones_like(x_smooth)])
y_smooth = X_smooth @ beta

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle("Tetraodontidae Latitudinal Gradient – Quadratic Fit (Normal Equations)",
             fontsize=13, fontweight="bold")

ax = axes[0]
ax.scatter(x, y, color="black", s=50, zorder=5, label="Observed", alpha=0.85)
ax.plot(x_smooth, y_smooth, color="steelblue", linewidth=2.5,
        label=f"$S = {a:.4f}x^2 + {b:.4f}x + {c:.2f}$\n$R^2 = {R2:.3f}$")
ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
ax.axvline(-b / (2*a), color="steelblue", linestyle=":", linewidth=1.5,
           alpha=0.7, label=f"Peak at {-b/(2*a):.1f}°")
ax.set_xlabel("Latitude (°)", fontsize=12)
ax.set_ylabel("Species Richness", fontsize=12)
ax.set_title("Quadratic Fit")
ax.legend(fontsize=9)
ax.set_xlim(-90, 90)

residuals = y - y_pred
axes[1].bar(x, residuals, width=BIN_WIDTH * 0.85,
            color=["#eb6f92" if r < 0 else "#9ccfd8" for r in residuals],
            edgecolor="white", alpha=0.85)
axes[1].axhline(0, color="black", linewidth=1)
axes[1].set_xlabel("Latitude (°)", fontsize=12)
axes[1].set_ylabel("Residual (obs − pred)", fontsize=12)
axes[1].set_title("Residuals")
axes[1].set_xlim(-90, 90)

plt.tight_layout()
plt.savefig("Tetraodontidae_quadratic_fit.png", dpi=150, bbox_inches="tight")
plt.show()

richness["S_pred"] = y_pred
richness["residual"] = residuals
richness.to_csv("Tetraodontidae_quadratic_fit.csv", index=False)