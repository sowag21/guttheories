import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fractions import Fraction

# ==========================================
# 1. GLOBAL STYLING
# ==========================================
plt.rcParams.update({
    "text.usetex": True,           
    "font.family": "serif",
    "axes.linewidth": 1.5          
})

COLOR_FLOW = "#1380c9"  # Solid dark blue for stream lines

# ==========================================
# 2. THEORETICAL PARAMETERS
# ==========================================
BYN = 1.0  

# Float versions for vector field calculations
def b0(Nf, p): return -(-13 + 18*Nf - 4*p) / 3
def c1(Nf): return (3*(5 - 2*Nf - 3*Nf**2)) / Nf
def c2(Nf, p): return (3*Nf + 5) / 2 * BYN
def d1(Nf): return 4*(4 + Nf)
def d2(Nf): return 6/Nf - 6*Nf
def d4(Nf): return (3*(2 - 4*Nf + Nf**2 + Nf**3)) / (4*Nf**2)
def k(Nf, p): return (b0(Nf, p) - d2(Nf))**2 - 4*d1(Nf)*d4(Nf)

# Exact fractions for the legends
def b0_frac(Nf, p): return Fraction(-(-13 + 18*Nf - 4*p), 3)
def d1_frac(Nf): return Fraction(4*(4 + Nf), 1)
def d2_frac(Nf): return Fraction(6, Nf) - Fraction(6*Nf, 1)
def d4_frac(Nf): return Fraction(3*(2 - 4*Nf + Nf**2 + Nf**3), 4*Nf**2)
def k_frac(Nf, p): return (b0_frac(Nf, p) - d2_frac(Nf))**2 - 4 * d1_frac(Nf) * d4_frac(Nf)

# ==========================================
# 3. BETA FUNCTIONS
# ==========================================
def betaG(ag, Nf, p):
    return b0(Nf, p) * ag**2

def betaH(ag, ah, Nf, p):
    return ah * (c1(Nf) * ag + c2(Nf, p) * ah)

def betaLambda(ag, al, Nf, p):
    return al * (d2(Nf) * ag + d1(Nf) * al) + d4(Nf) * ag**2

# Helper to format fractions beautifully in LaTeX
def format_latex_frac(f):
    if f.denominator == 1:
        return f"{f.numerator}"
    elif f.numerator < 0:
        return f"-\\frac{{{-f.numerator}}}{{{f.denominator}}}"
    else:
        return f"\\frac{{{f.numerator}}}{{{f.denominator}}}"

# ==========================================
# 4. PLOTTING FUNCTIONS
# ==========================================

def plot_gauge_yukawa_flow(Nf, p, ag_max=1.0, ah_max=1.0):
    """Plots the Gauge vs. Yukawa (alpha_g vs alpha_y) flow diagram."""
    ag_vals = np.linspace(0.001, ag_max, 200)
    ah_vals = np.linspace(0.001, ah_max, 200)
    AG, AH = np.meshgrid(ag_vals, ah_vals)

    U = betaG(AG, Nf, p)
    V = betaH(AG, AH, Nf, p)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.streamplot(AG, AH, U, V, color=COLOR_FLOW, linewidth=1.5, density=1.5, arrowsize=1.5)

    ax.set_xlabel(r'$\alpha_g$', fontsize=32)
    ax.set_ylabel(r'$\alpha_y$', fontsize=32)
    ax.set_xlim(0, ag_max)
    ax.set_ylim(0, ah_max)
    
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which='minor', width=1.0, length=4, direction='in', top=True, right=True)

    b0_exact = b0_frac(Nf, p)
    legend_elements = [
        Line2D([0], [0], color='none', label=fr'$b_0 = {format_latex_frac(b0_exact)}$'),
        Line2D([0], [0], color='none', label=fr'$c_1 = {int(c1(Nf))}$'),
        Line2D([0], [0], color='none', label=fr'$c_2 = {int(c2(Nf, p))}$')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=22, 
              title=r'\textbf{Parameters}', title_fontsize=24,
              frameon=True, edgecolor='black', fancybox=False, handlelength=0, handletextpad=0)

    plt.tight_layout()
    plt.savefig(f"Flowplot_Yukawa_N{Nf}p{p}.pdf", format='pdf', bbox_inches='tight')
    plt.show()


def plot_gauge_scalar_flow(Nf, p, ag_max=1.0, al_max=1.0):
    """Plots the Gauge vs. Scalar (alpha_g vs alpha_lambda) flow diagram."""
    ag_vals = np.linspace(0.001, ag_max, 200)
    al_vals = np.linspace(-al_max, al_max, 200)
    AG, AL = np.meshgrid(ag_vals, al_vals)

    U = betaG(AG, Nf, p)
    V = betaLambda(AG, AL, Nf, p)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.streamplot(AG, AL, U, V, color=COLOR_FLOW, linewidth=1.5, density=1.5, arrowsize=1.5)

    # Fixed flow lines
    k_val = k(Nf, p)
    if k_val >= 0:
        r1 = (b0(Nf, p) - d2(Nf) + np.sqrt(k_val)) / (2 * d1(Nf))
        r2 = (b0(Nf, p) - d2(Nf) - np.sqrt(k_val)) / (2 * d1(Nf))
        if abs(r1) > 0.001:
            ax.plot([0, ag_max], [0, r1 * ag_max], color='green', linewidth=2.5, linestyle='dashed')
        if abs(r2) > 0.001:
            ax.plot([0, ag_max], [0, r2 * ag_max], color='green', linewidth=2.5, linestyle='dashed')

    ax.set_xlabel(r'$\alpha_g$', fontsize=32)
    ax.set_ylabel(r'$\alpha_\lambda$', fontsize=32)
    ax.set_xlim(0, ag_max)
    ax.set_ylim(-al_max, al_max)
    
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which='minor', width=1.0, length=4, direction='in', top=True, right=True)

    b0_exact = b0_frac(Nf, p)
    k_exact = k_frac(Nf, p)
    legend_elements = [
        Line2D([0], [0], color='none', label=fr'$b_0 = {format_latex_frac(b0_exact)}$'),
        Line2D([0], [0], color='none', label=fr'$k = {format_latex_frac(k_exact)}$')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=22, 
              title=r'\textbf{Parameters}', title_fontsize=24,
              frameon=True, edgecolor='black', fancybox=False, handlelength=0, handletextpad=0)

    plt.tight_layout()
    plt.savefig(f"Flowplot_Scalar_N{Nf}p{p}.pdf", format='pdf', bbox_inches='tight')
    plt.show()

# ==========================================
# 5. EXECUTION
# ==========================================

# Generate the Gauge-Yukawa Flow Plot
plot_gauge_yukawa_flow(Nf=5, p=30)

# Generate the Gauge-Scalar Flow Plot
plot_gauge_scalar_flow(Nf=10, p=12)