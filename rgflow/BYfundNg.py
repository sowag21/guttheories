import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D 
import math

# --- Global Plot Settings ---
plt.rcParams.update({
    "text.usetex": True,           
    "font.family": "serif",
    "axes.linewidth": 1.5          
})

color_flow = '#4e97c7'  # Darker blue
color_free = '#a6cee3'  # Light blue

# --- Limits ---
points = 1500 
upp = 70
upNf = 30

def check_both_flows_BY(Nf, p, Ng):
    """
    Evaluates the BY model math for a specific (Nf, p, Ng) coordinate.
    Returns two booleans: (is_on_fixed_flow, is_off_fixed_flow).
    """
    M_val = Ng * Nf + 4 * Ng + p
    Kmin_val = min(Ng, M_val) 
    
    # Base equations
    b0 = -((-1 + 22*Nf - 4*Ng*Nf - 12*Ng - 4*p)/3)
    c1 = (3*(5 - 2*Nf - 3*Nf**2))/Nf
    c2 = Kmin_val + Kmin_val*Nf + (3 + Nf)/2
    d1 = 4*(4 + Nf)
    d2 = 6/Nf - 6*Nf
    d3 = 2*Kmin_val*(1 + Nf)
    d4 = (3*(2 - 4*Nf + Nf**2 + Nf**3))/(4*Nf**2)
    d5 = -(Kmin_val * (Nf + 3)/2)

    # ---------------------------------------------------------
    # 1. OFF FIXED FLOW (CAF_region_BY)
    # ---------------------------------------------------------
    k_val = (b0 - d2)**2 - 4*d1*d4
    
    F1 = b0
    F2 = b0 - c1
    F3 = k_val
    F4 = b0 - d2 + math.sqrt(k_val) if k_val >= 0 else -1
    
    is_off_flow = (F1 < 0) and (F2 > 0) and (F3 >= 0) and (F4 > 0)

    # ---------------------------------------------------------
    # 2. ON FIXED FLOW (CAF_region_BY_Flow)
    # ---------------------------------------------------------
    # Prevent division by zero just in case
    if c2 == 0:
        return False, is_off_flow
        
    ahFixed = (b0 - c1) / c2
    d2prime = d2 + d3*ahFixed
    d4prime = d4 + d5*(ahFixed**2)
    
    kprime_val = (b0 - d2prime)**2 - 4*d1*d4prime
    
    Fl1 = b0
    Fl2 = b0 - c1
    Fl3 = kprime_val
    Fl4 = b0 - d2prime + math.sqrt(kprime_val) if kprime_val >= 0 else -1
    
    is_on_flow = (Fl1 < 0) and (Fl2 > 0) and (Fl3 >= 0) and (Fl4 > 0)

    return is_on_flow, is_off_flow

def find_absolute_minimums_BY():
    max_Ng = 12 # Scans Ng=1 through 12
    results = {}
    
    for Ng in range(1, max_Ng + 1):
        found_fixed = False
        found_free = False
        
        results[Ng] = {
            'fixed': {'N': None, 'p': None},
            'free': {'N': None, 'p': None}
        }
        
        # Start scanning from N=5 upwards
        for N_int in range(3, upNf + 1):
            if found_fixed and found_free:
                break # Both found, move to next Ng
                
            for p_int in range(0, upp + 1):
                try:
                    is_on, is_off = check_both_flows_BY(N_int, p_int, Ng)
                    
                    if is_on and not found_fixed:
                        results[Ng]['fixed'] = {'N': N_int, 'p': p_int}
                        found_fixed = True
                        
                    if is_off and not found_free:
                        results[Ng]['free'] = {'N': N_int, 'p': p_int}
                        found_free = True
                        
                except ZeroDivisionError:
                    pass
                    
    return results

# --- Setting up the Grid for Plotting ---
Nf_vals = np.linspace(0.1, upNf, points)
p_vals = np.linspace(0, upp, points)
Nf, p = np.meshgrid(Nf_vals, p_vals)

# ==========================================================
# LOOP OVER Ng VALUES (1 through 5)
# ==========================================================
for Ng in range(1, 6):
    print(f"Calculating and plotting for Ng = {Ng}...")
    
    # --- PARAMETERS (Defined as functions just like before) ---
    def M(Nf, p): return Ng * Nf + 4 * Ng + p
    def Kmin(Nf, p): return np.minimum(Ng, M(Nf, p))
    
    def b0(Nf, p): return -((-1 + 22*Nf - 4*Ng*Nf - 12*Ng - 4*p)/3)
    def c1(Nf): return (3*(5 - 2*Nf - 3*Nf**2))/Nf
    def c2(Nf, p): return (Kmin(Nf, p) + Kmin(Nf, p)*Nf + (3 + Nf)/2)
    def d1(Nf): return 4*(4 + Nf)
    def d2(Nf): return 6/Nf - 6*Nf
    def d3(Nf, p): return 2*Kmin(Nf, p)*(1 + Nf)
    def d4(Nf): return (3*(2 - 4*Nf + Nf**2 + Nf**3))/(4*Nf**2)
    def d5(Nf, p): return -(Kmin(Nf, p) * (Nf + 3)/2)

    def ahFixed(Nf, p): return (b0(Nf, p) - c1(Nf))/c2(Nf, p)

    def d2prime(Nf, p): return d2(Nf) + d3(Nf, p)*(b0(Nf, p) - c1(Nf))/c2(Nf, p)
    def d4prime(Nf, p): return d4(Nf) + d5(Nf, p)*((b0(Nf, p) - c1(Nf))/c2(Nf, p))**2
    
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

    # Evaluate conditions on the grid
    CAF_region_BY = (F1(Nf, p) < 0) & (F2(Nf, p) > 0) & (F3(Nf, p) >= 0) & (F4(Nf, p) > 0)
    CAF_region_BY_Flow = (Fl1(Nf, p) < 0) & (Fl2(Nf, p) > 0) & (Fl3(Nf, p) >= 0) & (Fl4(Nf, p) > 0)

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(10, 8))

    # Base layer: Solid dark blue (On fixed flow)
    ax.contourf(Nf, p, CAF_region_BY_Flow, levels=[0.5, 1.5], colors=[color_flow], antialiased=True)

    # ---------------------------------------------------------
    # MATHEMATICAL STRIPES 
    # '30' = Spacing, '16' = Line thickness
    stripe_mask = ((p + 3.0 * Nf) % 30) < 16 

    # Mask the stripes so they only appear in the 'Off fixed flow' region
    CAF_region_BY_Striped = CAF_region_BY & stripe_mask

    # Plot the stripes as solid light blue shapes
    ax.contourf(Nf, p, CAF_region_BY_Striped, levels=[0.5, 1.5], colors=[color_free], antialiased=True)
    # ---------------------------------------------------------

    # Plot the boundary lines
    ax.contour(Nf, p, b0(Nf, p), levels=[0], colors='black', linewidths=2.5, linestyles='dashed')
    ax.contour(Nf, p, k(Nf, p), levels=[0], colors='black', linewidths=2.5, linestyles='dashdot')
    
    # Red lines specific to the BY plot
    ax.contour(Nf, p, b0(Nf, p)-c1(Nf), levels=[0], colors='red', linewidths=2.5, linestyles='solid')
    ax.contour(Nf, p, F4(Nf,p), levels=[0], colors='red', linewidths=2.5, linestyles='solid')

    # Formatting axes
    ax.set_xlabel('$N$', fontsize=32)
    ax.set_ylabel('$p$', fontsize=32)
    ax.set_xlim(3, upNf)
    ax.set_ylim(0, upp)
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    ax.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True, pad=12)

    # Legend formatting (Restored to original solid colors)
    legend_elements = [
        Patch(facecolor=color_free, label='Off fixed flow'),
        Patch(facecolor=color_flow, label='On fixed flow'),
        Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b_0=0$'),
        Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashdot', label=r'$k=0$')
    ]
    
    # Restored the title to the original top position, dynamically updating Ng
    ax.legend(handles=legend_elements, loc='upper left', fontsize=22, frameon=True, 
              edgecolor='black', fancybox=False, title=fr'\textbf{{BY model ($N_g={Ng}$)}}', title_fontsize=24)

    plt.tight_layout()
    
    # Save dynamically named file
    filename = f"CAFplotBY_Ng{Ng}.pdf"
    plt.savefig(filename, format='pdf', bbox_inches='tight')
    
    # Close the figure to clear memory before the next loop iteration
    plt.close(fig)

print("All plots generated and saved successfully!")

# ==========================================================
# EXECUTION CODE
# ==========================================================
print("\nCalculating absolute minimums for 'On fixed flow' and 'off fixed flow' (BY Model)...")
minimums = find_absolute_minimums_BY()

print("--- Absolute Minimums (N >= 5) ---")
for Ng, vals in minimums.items():
    print(f"\nModel Ng = {Ng}:")
    
    f_res = vals['fixed']
    if f_res['N'] is not None:
        print(f"  [On fixed flow]        Minimum N = {f_res['N']:<2}, p = {f_res['p']}")
    else:
        print("  [On fixed flow]        No valid region found.")
        
    n_res = vals['free']
    if n_res['N'] is not None:
        print(f"  [alpha_y off fixed flow] Minimum N = {n_res['N']:<2}, p = {n_res['p']}")
    else:
        print("  [alpha_y off fixed flow] No valid region found.")