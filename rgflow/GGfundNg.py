import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D 
import math

# --- Global Plotting Settings ---
plt.rcParams.update({
    "text.usetex": True,           
    "font.family": "serif",
    "axes.linewidth": 1.5,
    "hatch.linewidth": 3.0   
})

# Colors
color_flow = "#4e97c7"  # Darker blue 
color_free = '#a6cee3'  # Light blue 

# --- Limits & Grid Setup ---
points = 1500 
upp = 70
upNf = 30

Nf_vals = np.linspace(0.1, upNf, points)
p_vals = np.linspace(0, upp, points)
Nf, p = np.meshgrid(Nf_vals, p_vals)

def check_both_flows_GG_fund(Nf, p, Ng):
    """
    Evaluates the GG model math for a specific (Nf, p, Ng) coordinate.
    Returns two booleans: (is_on_fixed_flow, is_off_fixed_flow).
    """
    M_val = Ng * Nf - 4 * Ng + p
    Kmin_val = min(Ng, M_val) 
    
    # Base equations
    b0 = -((-1 + 22*Nf - 4*Ng*Nf + 12*Ng - 4*p)/3)
    c1 = (3*(5 + 2*Nf - 3*Nf**2))/Nf
    c2 = -Kmin_val + Kmin_val*Nf + (1 + Nf)/2
    d1 = 4*(4 + Nf)
    d2 = 6/Nf - 6*Nf
    d3 = 2*Kmin_val*(-1 + Nf)
    d4 = (3*(2 - 4*Nf + Nf**2 + Nf**3))/(4*Nf**2)
    d5 = -(Kmin_val * (Nf - 1)/2)

    # ---------------------------------------------------------
    # 1. OFF FIXED FLOW (CAF_region_BY / free flow)
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


def find_absolute_minimums_GG_fund():
    max_Ng = 14 # Scans up to Ng=14
    results = {}
    
    for Ng in range(1, max_Ng + 1):
        found_fixed = False
        found_free = False
        
        results[Ng] = {
            'fixed': {'N': None, 'p': None},
            'free': {'N': None, 'p': None}
        }
        
        # Start scanning from N=5 upwards
        for N_int in range(5, upNf + 1):
            if found_fixed and found_free:
                break # Move to next Ng if both are found
                
            for p_int in range(0, upp + 1):
                try:
                    is_on, is_off = check_both_flows_GG_fund(N_int, p_int, Ng)
                    
                    # Log the FIRST valid coordinate found for fixed flow
                    if is_on and not found_fixed:
                        results[Ng]['fixed'] = {'N': N_int, 'p': p_int}
                        found_fixed = True
                        
                    # Log the FIRST valid coordinate found for off fixed flow
                    if is_off and not found_free:
                        results[Ng]['free'] = {'N': N_int, 'p': p_int}
                        found_free = True
                        
                except ZeroDivisionError:
                    pass
                    
    return results

def find_max_ng():
    last_valid_Ng = 0
    last_valid_points = []
    
    for test_Ng in range(1, 100): 
        current_points = []
        
        # Loop through integers to find all valid points for this test_Ng
        for N_int in range(5, int(upNf) + 1): 
            for p_int in range(0, int(upp) + 1):
                
                # Inline definitions using test_Ng
                def M(Nf, p): return test_Ng * Nf - 4 * test_Ng + p
                def Kmin(Nf, p): return min(test_Ng, M(Nf, p))
                def b0(Nf, p): return -((-1 + 22*Nf - 4*test_Ng*Nf + 12*test_Ng - 4*p)/3)
                def c1(Nf): return (3*(5 + 2*Nf - 3*Nf**2))/Nf
                def c2(Nf, p): return (-Kmin(Nf, p) + Kmin(Nf, p)*Nf + (1 + Nf)/2) 
                def d1(Nf): return 4*(4 + Nf)
                def d2(Nf): return 6/Nf - 6*Nf
                def d3(Nf, p): return 2*Kmin(Nf, p)*(-1 + Nf) 
                def d4(Nf): return (3*(2 - 4*Nf + Nf**2 + Nf**3))/(4*Nf**2)
                def d5(Nf, p): return -(Kmin(Nf, p) * (Nf - 1)/2) 
                def d2prime(Nf, p): return d2(Nf) + d3(Nf, p)*(b0(Nf, p) - c1(Nf))/c2(Nf, p)
                def d4prime(Nf, p): return d4(Nf) + d5(Nf, p)*((b0(Nf, p) - c1(Nf))/c2(Nf, p))**2
                def kprime(Nf, p): return (b0(Nf, p) - d2prime(Nf, p))**2 - 4*d1(Nf)*d4prime(Nf, p)

                try:
                    # Check fixed flow conditions
                    if b0(N_int, p_int) < 0 and (b0(N_int, p_int) - c1(N_int)) > 0:
                        k_val = kprime(N_int, p_int)
                        if k_val >= 0:
                            fl4_val = b0(N_int, p_int) - d2prime(N_int, p_int) + np.sqrt(k_val)
                            if fl4_val > 0:
                                current_points.append((N_int, p_int))
                except ZeroDivisionError:
                    pass
        
        # If we found points, save them as the latest valid ones
        if len(current_points) > 0:
            last_valid_Ng = test_Ng
            last_valid_points = current_points
        else:
            # The moment we find NO points, return the data from the previous loop
            return last_valid_Ng, last_valid_points

# --- Run the analysis and print cleanly ---
max_ng_val, surviving_coords = find_max_ng()

print(f"\n--- MATHEMATICAL ANALYSIS ---")
print(f"Maximum Ng with valid integer solutions (for N >= 5): Ng = {max_ng_val}")
print(f"Surviving coordinates (N, p) at this limit: {surviving_coords}")
print(f"-----------------------------\n")

# ==========================================================
# LOOP OVER Ng VALUES (1 through 13)
# ==========================================================
for Ng in range(1, 14):
    print(f"Calculating and plotting for Ng = {Ng}...")
    
    # --- PARAMETERS (Defined as functions) ---
    def M(Nf, p): return Ng * Nf - 4 * Ng + p
    def Kmin(Nf, p): return np.minimum(Ng, M(Nf, p))
    
    def b0(Nf, p): return -((-1 + 22*Nf - 4*Ng*Nf + 12*Ng - 4*p)/3)
    def c1(Nf): return (3*(5 + 2*Nf - 3*Nf**2))/Nf
    def c2(Nf, p): return (-Kmin(Nf, p) + Kmin(Nf, p)*Nf + (1 + Nf)/2) 
    def d1(Nf): return 4*(4 + Nf)
    def d2(Nf): return 6/Nf - 6*Nf
    def d3(Nf, p): return 2*Kmin(Nf, p)*(-1 + Nf) 
    def d4(Nf): return (3*(2 - 4*Nf + Nf**2 + Nf**3))/(4*Nf**2)
    def d5(Nf, p): return -(Kmin(Nf, p) * (Nf - 1)/2) 

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

    # ---------------------------------------------------------
    # PAINTER's ALGORITHM (Stripes ONLY in the overlap)
    # ---------------------------------------------------------
    
    # 1. Base layer: Solid dark blue (On fixed flow)
    ax.contourf(Nf, p, CAF_region_BY_Flow, levels=[0.5, 1.5], colors=[color_flow], antialiased=True)
    
    # 2. Base layer: Solid light blue (Off fixed flow)
    # Must be plotted before intersection layers
    ax.contourf(Nf, p, CAF_region_BY, levels=[0.5, 1.5], colors=[color_free], antialiased=True)

    # 3. Intersection Base: Re-plot the overlap area as solid dark blue
    # This covers up the light blue from step 2 where they intersect
    intersection = CAF_region_BY_Flow & CAF_region_BY
    ax.contourf(Nf, p, intersection, levels=[0.5, 1.5], colors=[color_flow], antialiased=True)
    
    # 4. Intersection Stripes: Overlay tight mathematical stripes in the intersection
    # (If the intersection in GG is very large and you prefer the wide stripes from before,
    # simply change this line to: stripe_mask_tight = ((p + 3.0 * Nf) % 30) < 16)
    stripe_mask_tight = ((p + 3.0 * Nf) % 30) < 16
    striped_intersection = intersection & stripe_mask_tight
    ax.contourf(Nf, p, striped_intersection, levels=[0.5, 1.5], colors=[color_free], antialiased=True)
    # ---------------------------------------------------------
    # =========================================================
    # NEW: MATHEMATICA p FUNCTIONS (Depends on Ng)
    # =========================================================
    def p_branches(Nf, Ng):
        with np.errstate(divide='ignore', invalid='ignore'):
            # Terms outside the fraction
            base_terms = -0.25 + Nf + 3.0 * Ng - Nf * Ng
            
            # Polynomial inside the square root: 8 - 14*Nf + 5*Nf^3 + Nf^4
            inner_sqrt = 8.0 - 14.0 * Nf + 5.0 * Nf**3 + Nf**4
            
            # Safe sqrt calculation
            safe_inner = np.where(inner_sqrt >= 0, inner_sqrt, np.nan)
            sqrt_term = 3.0 * np.sqrt(3.0) * np.sqrt(safe_inner)
            
            # Calculate the minus and plus branches
            p_minus = base_terms + (9.0 - sqrt_term) / (2.0 * Nf)
            p_plus  = base_terms + (9.0 + sqrt_term) / (2.0 * Nf)
            
            # Filter extreme values just to keep the plot clean
            p_minus = np.where(np.abs(p_minus) > 500, np.nan, p_minus)
            p_plus = np.where(np.abs(p_plus) > 500, np.nan, p_plus)
            
        return p_minus, p_plus

    # Generate 1D array specifically for plotting the lines smoothly
    Nf_line = np.linspace(3, upNf, points)
    p_minus_vals, p_plus_vals = p_branches(Nf_line, Ng)

    # Plot the two branches
    ax.plot(Nf_line, p_minus_vals, color='#e41a1c', linewidth=2.5, zorder=5)  # Red line
    ax.plot(Nf_line, p_plus_vals, color='#984ea3', linewidth=2.5, zorder=5)   # Purple line
    # =========================================================

    # Plot the boundary lines
    ax.contour(Nf, p, b0(Nf, p), levels=[0], colors='black', linewidths=2.5, linestyles='dashed')
    ax.contour(Nf, p, k(Nf, p), levels=[0], colors='black', linewidths=2.5, linestyles='dashdot')
    ax.contour(Nf, p, Fl4(Nf,p), levels=[0], colors='black', linewidths=2.5, linestyles='dotted')

    # Formatting axes
    ax.set_xlabel('$N$', fontsize=32)
    ax.set_ylabel('$p$', fontsize=32)
    ax.set_xlim(5, upNf)
    ax.set_ylim(0, upp)
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    ax.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True, pad=12)

    # Legend formatting (Solid colors only)
    legend_elements = [
        Patch(facecolor=color_free, label='Off fixed flow'),
        Patch(facecolor=color_flow, label='On fixed flow'),
        Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b_0=0$'),
        Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashdot', label=r'$k=0$'),
        Line2D([0], [0], color='black', linewidth=2.5, linestyle='dotted', label=r'$b - d_2^\prime + \sqrt{k^\prime} = 0$')
    ]
    
    # Restored the title to the original top position, dynamically updating Ng
    ax.legend(handles=legend_elements, loc='upper left', fontsize=22, frameon=True, 
              edgecolor='black', fancybox=False, title=fr'\textbf{{GG model ($N_g={Ng}$)}}', title_fontsize=24)
    plt.tight_layout()
    
    # Save dynamically named file
    filename = f"CAFGG_Ng{Ng}.pdf"
    plt.savefig(filename, format='pdf', bbox_inches='tight')
    
    # Close the figure to clear memory before the next loop iteration
    plt.close(fig)

print("All plots generated and saved successfully!")

# ==========================================================
# EXECUTION CODE
# ==========================================================
print("\nCalculating absolute minimums for 'On fixed flow' and 'off fixed flow' (GG Model - Fundamental)...")
minimums = find_absolute_minimums_GG_fund()

print("--- Absolute Minimums (N >= 5) ---")
for Ng, vals in minimums.items():
    print(f"\nModel Ng = {Ng}:")
    
    f_res = vals['fixed']
    if f_res['N'] is not None:
        print(f"  [On fixed flow]          Minimum N = {f_res['N']:<2}, p = {f_res['p']}")
    else:
        print("  [On fixed flow]          No valid region found.")
        
    n_res = vals['free']
    if n_res['N'] is not None:
        print(f"  [alpha_y off fixed flow] Minimum N = {n_res['N']:<2}, p = {n_res['p']}")
    else:
        print("  [alpha_y off fixed flow] No valid region found.")