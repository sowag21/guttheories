import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D 
import math

# --- Limits ---
points = 1500 
BYN = 1
upp = 70
upNf = 30

# --- PARAMETERS for GG gauge-scalar system WITH FIXED FLOW ---
def b0(Nf, p): return -(1/3) * (11 + 18*Nf - 4*p)
def c1(Nf): return (3 * (5 + 2*Nf - 3*Nf**2))/Nf
def c2(Nf, p): return (3*Nf - 1)/2 
def d1(Nf): return 4*(4 + Nf)
def d2(Nf): return 6/Nf - 6*Nf
def d3(Nf): return 2*(-1 + Nf) 
def d4(Nf): return (3*(2 - 4*Nf + Nf**2 + Nf**3))/(4*Nf**2)
def d5(Nf): return -((Nf - 1)/2) 

def ahFixed(Nf, p): return (b0(Nf, p) - c1(Nf))/c2(Nf, p)

def d2prime(Nf, p): return d2(Nf) + d3(Nf)*(b0(Nf, p) - c1(Nf))/c2(Nf, p)
def d4prime(Nf, p): return d4(Nf) + d5(Nf)*((b0(Nf, p) - c1(Nf))/c2(Nf, p))**2

def k(Nf, p): return (b0(Nf, p) - d2(Nf))**2 - 4*d1(Nf)*d4(Nf)
def kprime(Nf, p): return (b0(Nf, p) - d2prime(Nf, p))**2 - 4*d1(Nf)*d4prime(Nf, p)

# --- CAF conditions ---
def F1(Nf, p): return b0(Nf, p)
def F2(Nf, p): return b0(Nf, p) - c1(Nf)
def F3(Nf, p): return k(Nf, p)
def F4(Nf, p): 
    k_val = k(Nf, p)
    safe_sqrt = np.where(k_val >= 0, np.sqrt(np.maximum(k_val, 0)), -1) 
    return b0(Nf, p) - d2(Nf) + safe_sqrt

def Fl1(Nf, p): return b0(Nf, p)
def Fl2(Nf, p): return b0(Nf, p) - c1(Nf)
def Fl3(Nf, p): return kprime(Nf, p)
def Fl4(Nf, p): 
    k_val = kprime(Nf, p)
    safe_sqrt = np.where(k_val >= 0, np.sqrt(np.maximum(k_val, 0)), -1)
    return b0(Nf, p) - d2prime(Nf, p) + safe_sqrt

# --- Find the smallest N and p ---
def find_absolute_minimums():
    min_N_flow = None
    min_p_flow = None
    
    # Start scanning from N=1 upwards
    for N_int in range(5, int(upNf) + 1):
        for p_int in range(0, int(upp) + 1):
            
            try:
                # Check fixed flow conditions (dark blue region)
                if (Fl1(N_int, p_int) < 0 and 
                    Fl2(N_int, p_int) > 0 and 
                    Fl3(N_int, p_int) >= 0 and 
                    Fl4(N_int, p_int) > 0):
                    
                    min_N_flow = N_int
                    min_p_flow = p_int
                    return min_N_flow, min_p_flow
            except ZeroDivisionError:
                pass
                
    return None, None

min_N, min_p = find_absolute_minimums()

print("\n--- MINIMUM COORDINATES ---")
if min_N is not None:
    print(f"The absolute smallest valid integer N is: {min_N}")
    print(f"The smallest valid integer p for N={min_N} is: {min_p}")
else:
    print("No valid integer coordinates found in the specified range.")
print("---------------------------\n")

# --- Setting up the Grid for Plotting ---
Nf_vals = np.linspace(0.1, upNf, points)
p_vals = np.linspace(0, upp, points)
Nf, p = np.meshgrid(Nf_vals, p_vals)

# Evaluate conditions
CAF_region_BY = (F1(Nf, p) < 0) & (F2(Nf, p) > 0) & (F3(Nf, p) >= 0) & (F4(Nf, p) > 0)
CAF_region_BY_Flow = (Fl1(Nf, p) < 0) & (Fl2(Nf, p) > 0) & (Fl3(Nf, p) >= 0) & (Fl4(Nf, p) > 0)

# Apply bounds to the smaller region
BY_bounds = (Nf >= 2.5) & (p >= 0)
CAF_region_BY = CAF_region_BY & BY_bounds

# --- Plotting ---
plt.rcParams.update({
    "text.usetex": True,           
    "font.family": "serif",
    "axes.linewidth": 1.5,
    "hatch.linewidth": 10.0   
})

fig, ax = plt.subplots(figsize=(10, 8))

# Colors
color_flow = "#4e97c7"  # Darker blue 
color_free = '#a6cee3'  # Light blue 

# 1. Plot the regions (Order and Style Swapped)
# Base layer: Solid dark blue (On fixed flow)
ax.contourf(Nf, p, CAF_region_BY_Flow, levels=[0.5, 1.5], colors=[color_flow], antialiased=True)

# MATHEMATICAL STRIPES 
stripe_mask = ((p + 3.0 * Nf) % 30) < 16 

# Keep only the stripes that fall inside the Off Fixed Flow region
CAF_region_BY_Striped = CAF_region_BY & stripe_mask
ax.contourf(Nf, p, CAF_region_BY_Striped, levels=[0.5, 1.5], colors=[color_free], antialiased=True)

# 2. Plot the b0 = 0 and k = 0 boundary lines directly
ax.contour(Nf, p, b0(Nf, p), levels=[0], colors='black', linewidths=2.5, linestyles='dashed')
ax.contour(Nf, p, k(Nf, p), levels=[0], colors='black', linewidths=2.5, linestyles='dashdot')
ax.contour(Nf, p, Fl4(Nf,p), levels=[0], colors='black', linewidths=2.5, linestyles='dotted')


# =========================================================
# NEW: MATHEMATICA p FUNCTIONS (Updated formulas)
# =========================================================
def p_branches(Nf):
    with np.errstate(divide='ignore', invalid='ignore'):
        # Denominator: 16 * (Nf - 1) * Nf * (Nf + 4)
        denom = 16.0 * (Nf - 1.0) * Nf * (Nf + 4.0)
        
        # Polynomial base: -720 + 4*Nf*(19 + Nf*(168 + Nf*(2 - 9*Nf)))
        base = -720.0 + 4.0 * Nf * (19.0 + Nf * (168.0 + Nf * (2.0 - 9.0 * Nf)))
        
        # Inner sqrt: (1 - 3*Nf)^2 * (Nf - 1)^2 * (Nf + 4)^2 * (Nf^2 + 2*Nf - 2)
        inner_sqrt = (1.0 - 3.0 * Nf)**2 * (Nf - 1.0)**2 * (Nf + 4.0)**2 * (Nf**2 + 2.0 * Nf - 2.0)
        
        # Safe sqrt calculation
        safe_inner = np.where(inner_sqrt >= 0, inner_sqrt, np.nan)
        sqrt_term = 3.0 * np.sqrt(6.0) * np.sqrt(safe_inner)
        
        # The two branches
        p_minus = (base - sqrt_term) / denom
        p_plus = (base + sqrt_term) / denom
        
        # Filter extremes just in case plotting strays
        p_minus = np.where(np.abs(p_minus) > 500, np.nan, p_minus)
        p_plus = np.where(np.abs(p_plus) > 500, np.nan, p_plus)
        
    return p_minus, p_plus

# Generate 1D array specifically for plotting the lines smoothly
Nf_line = np.linspace(3, upNf, points)
p_minus_vals, p_plus_vals = p_branches(Nf_line)

# Plot the two branches
#ax.plot(Nf_line, p_minus_vals, color='#e41a1c', linewidth=2.5, zorder=5)  # Red line
#ax.plot(Nf_line, p_plus_vals, color='#984ea3', linewidth=2.5, zorder=5)   # Purple line
# =========================================================


# Formatting axes
ax.set_xlabel('$N$', fontsize=32)
ax.set_ylabel('$p$', fontsize=32)
ax.set_xlim(5, upNf)
ax.set_xticks([5, 10, 15, 20, 25, 30])
ax.set_ylim(0, upp)

ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True)

# Legend formatting
legend_elements = [
    Patch(facecolor=color_free, label='Off fixed flow'),
    Patch(facecolor=color_flow, label='On fixed flow'),
    Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b=0$'),
    Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashdot', label=r'$k = 0$'),
    Line2D([0], [0], color='black', linewidth=2.5, linestyle='dotted', label=r'$b - d_2^\prime + \sqrt{k^\prime} = 0$'),
    #Line2D([0], [0], color='#e41a1c', linewidth=2.5, linestyle='solid', label=r'$p_1$'),
    #Line2D([0], [0], color='#984ea3', linewidth=2.5, linestyle='solid', label=r'$p_2$')
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=22, frameon=True, 
          edgecolor='black', fancybox=False, title=r'\textbf{GG model}', title_fontsize=24)

plt.tight_layout()
plt.savefig("CAFGG_fund_Ng1.pdf", format='pdf', bbox_inches='tight')
plt.show()