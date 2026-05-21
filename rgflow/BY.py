import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D 

# --- Limits ---
points = 1500 
upp = 70
upNf = 30

# --- PARAMETERS for gauge system WITH NO SCALARS ---
def b0(Nf, p): return -(2/3) * (-6 + 9*Nf - 2*p)

# --- CAF conditions ---
def F1(Nf, p): return b0(Nf, p)

# --- Setting up the Grid for Plotting ---
Nf_vals = np.linspace(0.1, upNf, points)
p_vals = np.linspace(0, upp, points)
Nf, p = np.meshgrid(Nf_vals, p_vals)

# Evaluate conditions
CAF_region_BY = (F1(Nf, p) < 0)

# --- Plotting ---
plt.rcParams.update({
    "text.usetex": True,           
    "font.family": "serif",
    "axes.linewidth": 1.5          
})

fig, ax = plt.subplots(figsize=(10, 8))

# Colors
color_free = '#1f78b4'  # Light blue 

# 1. Plot the single region
ax.contourf(Nf, p, CAF_region_BY, levels=[0.5, 1.5], colors=[color_free], antialiased=True)

# 2. Plot the b0 = 0 boundary line directly
ax.contour(Nf, p, b0(Nf, p), levels=[0], colors='black', linewidths=2.5, linestyles='dashed')

# Formatting axes
ax.set_xlabel('$N$', fontsize=32)
ax.set_ylabel('$p$', fontsize=32)
ax.set_xlim(3, upNf)
ax.set_ylim(0, upp)
ax.set_xticks([5, 10, 15, 20, 25, 30])
ax.set_yticks([0,10, 20, 30, 40, 50, 60, 70])

ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)

# Legend formatting
legend_elements = [
    Patch(facecolor=color_free, label='CAF region'),
    Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b_0=0$')
]
# Placed in upper left so it doesn't overlap the colored area (which fills the bottom right)
ax.legend(handles=legend_elements, loc='lower right', fontsize=22, frameon=True, 
          edgecolor='black', fancybox=False, title=r'\textbf{BY model}', title_fontsize=24)

plt.tight_layout()
plt.savefig("CAFBY_noscalar.pdf", format='pdf', bbox_inches='tight')
plt.show()