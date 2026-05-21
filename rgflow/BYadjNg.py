import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import math

# --- Vectorized Descartes' Rule of Signs ---
def count_sign_changes(coeffs):
    """
    Properly skips zeros to count sign changes, perfectly mimicking 
    Mathematica's DeleteCases[list, 0] behavior.
    """
    signs = np.sign(coeffs)
    changes = np.zeros(coeffs.shape[:-1], dtype=int)
    last_nonzero_sign = np.zeros(coeffs.shape[:-1], dtype=int)
    
    for i in range(coeffs.shape[-1]):
        current_sign = signs[..., i]
        # Valid means it is neither zero nor a NaN 
        valid = (current_sign != 0) & ~np.isnan(current_sign)
        
        # A change occurs if we have a valid sign, a previous valid sign, and they differ
        is_change = valid & (last_nonzero_sign != 0) & (current_sign != last_nonzero_sign)
        changes[is_change] += 1
        
        # Update the last seen non-zero sign
        last_nonzero_sign = np.where(valid, current_sign, last_nonzero_sign)
        
    return changes

def check_both_flows(Nf, p, Ng):
    """
    Evaluates the math for a specific (Nf, p, Ng) coordinate 
    and returns two booleans: (is_fixed_flow, is_noY_flow).
    """
    # Lower sign evaluated: K = min(p, (N + 4)*Ng + p)
    M_val = Ng * Nf + 4 * Ng + p
    KminVal = min(p, M_val) 
    
    b0 = -((-4*Ng*Nf - 12*Ng + 21*Nf - 4*p) / 3)
    c1 = 6/Nf - 6*Nf
    c2 = 2*KminVal + Nf - 3/Nf
    
    # 1. Gauge Condition
    if not (b0 < 0 and (b0 - c1) > 0):
        return False, False
        
    alphag = -1 / b0
    a = (4 * (Nf**2 - 9)) / Nf
    c_val = 12                                     
    e = Nf**2 + 7
    f = (12 * (3 + Nf**2)) / Nf**2
    g = (4 * (2*Nf**2 - 3)) / Nf
    i = 18 * alphag**2

    # --- Fixed Flow Calculations ---
    alphay_fixed = (c1 - b0) / (b0 * c2)
    b_coeff_fixed = 1 + 4*KminVal*alphay_fixed - 12*Nf*alphag
    d_coeff_fixed = 3*Nf*alphag**2 - 2*KminVal*alphay_fixed**2
    
    alpha_fixed = e*a**2 + f*c_val**2 - g*c_val*a
    beta_fixed = 2*e*a*b_coeff_fixed - g*c_val*b_coeff_fixed - b_coeff_fixed*c_val*a
    gamma_fixed = e*(b_coeff_fixed**2 + 2*a*d_coeff_fixed) - g*c_val*d_coeff_fixed - b_coeff_fixed*c_val*b_coeff_fixed + i*c_val**2
    delta_fixed = 2*e*b_coeff_fixed*d_coeff_fixed - b_coeff_fixed*c_val*d_coeff_fixed
    epsilon_fixed = e * d_coeff_fixed**2
    
    Delta_fixed = (256 * alpha_fixed**3 * epsilon_fixed**3 
             - 192 * alpha_fixed**2 * beta_fixed * delta_fixed * epsilon_fixed**2 
             - 128 * alpha_fixed**2 * gamma_fixed**2 * epsilon_fixed**2 
             + 144 * alpha_fixed**2 * gamma_fixed * delta_fixed**2 * epsilon_fixed 
             - 27 * alpha_fixed**2 * delta_fixed**4 
             + 144 * alpha_fixed * beta_fixed**2 * gamma_fixed * epsilon_fixed**2 
             - 6 * alpha_fixed * beta_fixed**2 * delta_fixed**2 * epsilon_fixed 
             - 80 * alpha_fixed * beta_fixed * gamma_fixed**2 * delta_fixed * epsilon_fixed 
             + 18 * alpha_fixed * beta_fixed * gamma_fixed * delta_fixed**3 
             + 16 * alpha_fixed * gamma_fixed**4 * epsilon_fixed 
             - 4 * alpha_fixed * gamma_fixed**3 * delta_fixed**2 
             - 27 * beta_fixed**4 * epsilon_fixed**2 
             + 18 * beta_fixed**3 * gamma_fixed * delta_fixed * epsilon_fixed 
             - 4 * beta_fixed**3 * delta_fixed**3 
             - 4 * beta_fixed**2 * gamma_fixed**3 * epsilon_fixed 
             + beta_fixed**2 * gamma_fixed**2 * delta_fixed**2)
             
    P_fixed = 8*alpha_fixed*gamma_fixed - 3*beta_fixed**2
    D0_fixed = 64*alpha_fixed**3*epsilon_fixed - 16*alpha_fixed**2*gamma_fixed**2 + 16*alpha_fixed*beta_fixed**2*gamma_fixed - 16*alpha_fixed**2*beta_fixed*delta_fixed - 3*beta_fixed**4

    roots_4_fixed = (Delta_fixed > 0) and (D0_fixed < 0) and (P_fixed < 0)
    roots_2_fixed = (Delta_fixed < 0)
    has_real_fixed = roots_4_fixed or roots_2_fixed

    # --- No Yukawa (alpha_y = 0) Flow Calculations ---
    b_coeff_noY = 1 - 12*Nf*alphag
    d_coeff_noY = 3*Nf*alphag**2
    
    alpha_noY = e*a**2 + f*c_val**2 - g*c_val*a
    beta_noY = 2*e*a*b_coeff_noY - g*c_val*b_coeff_noY - b_coeff_noY*c_val*a
    gamma_noY = e*(b_coeff_noY**2 + 2*a*d_coeff_noY) - g*c_val*d_coeff_noY - b_coeff_noY*c_val*b_coeff_noY + i*c_val**2
    delta_noY = 2*e*b_coeff_noY*d_coeff_noY - b_coeff_noY*c_val*d_coeff_noY
    epsilon_noY = e * d_coeff_noY**2
    
    Delta_noY = (256 * alpha_noY**3 * epsilon_noY**3 
             - 192 * alpha_noY**2 * beta_noY * delta_noY * epsilon_noY**2 
             - 128 * alpha_noY**2 * gamma_noY**2 * epsilon_noY**2 
             + 144 * alpha_noY**2 * gamma_noY * delta_noY**2 * epsilon_noY 
             - 27 * alpha_noY**2 * delta_noY**4 
             + 144 * alpha_noY * beta_noY**2 * gamma_noY * epsilon_noY**2 
             - 6 * alpha_noY * beta_noY**2 * delta_noY**2 * epsilon_noY 
             - 80 * alpha_noY * beta_noY * gamma_noY**2 * delta_noY * epsilon_noY 
             + 18 * alpha_noY * beta_noY * gamma_noY * delta_noY**3 
             + 16 * alpha_noY * gamma_noY**4 * epsilon_noY 
             - 4 * alpha_noY * gamma_noY**3 * delta_noY**2 
             - 27 * beta_noY**4 * epsilon_noY**2 
             + 18 * beta_noY**3 * gamma_noY * delta_noY * epsilon_noY 
             - 4 * beta_noY**3 * delta_noY**3 
             - 4 * beta_noY**2 * gamma_noY**3 * epsilon_noY 
             + beta_noY**2 * gamma_noY**2 * delta_noY**2)

    P_noY = 8*alpha_noY*gamma_noY - 3*beta_noY**2
    D0_noY = 64*alpha_noY**3*epsilon_noY - 16*alpha_noY**2*gamma_noY**2 + 16*alpha_noY*beta_noY**2*gamma_noY - 16*alpha_noY**2*beta_noY*delta_noY - 3*beta_noY**4

    roots_4_noY = (Delta_noY > 0) and (D0_noY < 0) and (P_noY < 0)
    roots_2_noY = (Delta_noY < 0)
    has_real_noY = roots_4_noY or roots_2_noY

    return has_real_fixed, has_real_noY

def find_absolute_minimums():
    upNf = 30
    upp = 70
    
    results = {}
    
    for Ng in range(1, 6):
        found_fixed = False
        found_noY = False
        
        results[Ng] = {
            'fixed': {'N': None, 'p': None},
            'noY': {'N': None, 'p': None}
        }
        
        # Start scanning from N=5 upwards
        for N_int in range(2, int(upNf) + 1):
            if found_fixed and found_noY:
                break 
                
            for p_int in range(0, int(upp) + 1):
                try:
                    is_fixed, is_noY = check_both_flows(N_int, p_int, Ng)
                    
                    if is_fixed and not found_fixed:
                        results[Ng]['fixed'] = {'N': N_int, 'p': p_int}
                        found_fixed = True
                        
                    if is_noY and not found_noY:
                        results[Ng]['noY'] = {'N': N_int, 'p': p_int}
                        found_noY = True
                        
                except ZeroDivisionError:
                    pass
                    
    return results

# --- Limits ---
points = 1500 
upp = 70
upNf = 30

# --- Setting up the Grid for Plotting ---
Nf_vals = np.linspace(0.1, upNf, points)
p_vals = np.linspace(0, upp, points)
Nf, p = np.meshgrid(Nf_vals, p_vals)

# --- Global Plot Settings ---
plt.rcParams.update({
    "text.usetex": True,           
    "font.family": "serif",
    "axes.linewidth": 1.5,
    "hatch.linewidth": 3.0  
})

# --- Color Palettes ---
color_caf_flow = '#1f78b4'  
color_caf_free = '#a6cee3'  
color_roots_4  = '#33a02c'  
color_roots_3  = '#4daf4a'  
color_roots_2  = '#b2df8a'  
color_roots_1  = '#fdbf6f'  
color_roots_0  = '#e41a1c'  

# ==========================================================
# LOOP OVER Ng VALUES (1 through 5)
# ==========================================================
for Ng in range(1, 6):
    print(f"Calculating and plotting for Ng = {Ng}...")
    
    with np.errstate(divide='ignore', invalid='ignore'):
        # ---------------------------------------------------------
        # 1. SHARED BASE PARAMETERS (Lower Sign Evaluated)
        # ---------------------------------------------------------
        M_val = Ng * Nf + 4 * Ng + p
        KminVal = np.minimum(p, M_val) 
        
        b0 = -((-4*Ng*Nf - 12*Ng + 21*Nf - 4*p) / 3)
        c1 = 6/Nf - 6*Nf
        c2 = 2*KminVal + Nf - 3/Nf
        
        alphag = -1 / b0
        
        a = (4 * (Nf**2 - 9)) / Nf
        c_val = 12                                     
        e = Nf**2 + 7
        f = (12 * (3 + Nf**2)) / Nf**2
        g = (4 * (2*Nf**2 - 3)) / Nf
        i = 18 * alphag**2

        # ---------------------------------------------------------
        # 2. STANDARD FIXED FLOW (alphay != 0)
        # ---------------------------------------------------------
        alphay_fixed = (c1 - b0) / (b0 * c2)
        
        b_coeff_fixed = 1 + 4*KminVal*alphay_fixed - 12*Nf*alphag
        d_coeff_fixed = 3*Nf*alphag**2 - 2*KminVal*alphay_fixed**2
        
        alpha_fixed = e*a**2 + f*c_val**2 - g*c_val*a
        beta_fixed = 2*e*a*b_coeff_fixed - g*c_val*b_coeff_fixed - b_coeff_fixed*c_val*a
        gamma_fixed = e*(b_coeff_fixed**2 + 2*a*d_coeff_fixed) - g*c_val*d_coeff_fixed - b_coeff_fixed*c_val*b_coeff_fixed + i*c_val**2
        delta_fixed = 2*e*b_coeff_fixed*d_coeff_fixed - b_coeff_fixed*c_val*d_coeff_fixed
        epsilon_fixed = e * d_coeff_fixed**2
        
        Delta_fixed = (256 * alpha_fixed**3 * epsilon_fixed**3 
                 - 192 * alpha_fixed**2 * beta_fixed * delta_fixed * epsilon_fixed**2 
                 - 128 * alpha_fixed**2 * gamma_fixed**2 * epsilon_fixed**2 
                 + 144 * alpha_fixed**2 * gamma_fixed * delta_fixed**2 * epsilon_fixed 
                 - 27 * alpha_fixed**2 * delta_fixed**4 
                 + 144 * alpha_fixed * beta_fixed**2 * gamma_fixed * epsilon_fixed**2 
                 - 6 * alpha_fixed * beta_fixed**2 * delta_fixed**2 * epsilon_fixed 
                 - 80 * alpha_fixed * beta_fixed * gamma_fixed**2 * delta_fixed * epsilon_fixed 
                 + 18 * alpha_fixed * beta_fixed * gamma_fixed * delta_fixed**3 
                 + 16 * alpha_fixed * gamma_fixed**4 * epsilon_fixed 
                 - 4 * alpha_fixed * gamma_fixed**3 * delta_fixed**2 
                 - 27 * beta_fixed**4 * epsilon_fixed**2 
                 + 18 * beta_fixed**3 * gamma_fixed * delta_fixed * epsilon_fixed 
                 - 4 * beta_fixed**3 * delta_fixed**3 
                 - 4 * beta_fixed**2 * gamma_fixed**3 * epsilon_fixed 
                 + beta_fixed**2 * gamma_fixed**2 * delta_fixed**2)
                 
        P_fixed = 8*alpha_fixed*gamma_fixed - 3*beta_fixed**2
        D0_fixed = 64*alpha_fixed**3*epsilon_fixed - 16*alpha_fixed**2*gamma_fixed**2 + 16*alpha_fixed*beta_fixed**2*gamma_fixed - 16*alpha_fixed**2*beta_fixed*delta_fixed - 3*beta_fixed**4

        roots_4_fixed = (Delta_fixed > 0) & (D0_fixed < 0) & (P_fixed < 0)
        roots_2_fixed = (Delta_fixed < 0)
        has_real_roots_fixed = roots_4_fixed | roots_2_fixed

        # ---------------------------------------------------------
        # 3. MODIFIED FLOW (alphay = 0)
        # ---------------------------------------------------------
        b_coeff_noY = 1 - 12*Nf*alphag
        d_coeff_noY = 3*Nf*alphag**2
        
        alpha_noY = e*a**2 + f*c_val**2 - g*c_val*a
        beta_noY = 2*e*a*b_coeff_noY - g*c_val*b_coeff_noY - b_coeff_noY*c_val*a
        gamma_noY = e*(b_coeff_noY**2 + 2*a*d_coeff_noY) - g*c_val*d_coeff_noY - b_coeff_noY*c_val*b_coeff_noY + i*c_val**2
        delta_noY = 2*e*b_coeff_noY*d_coeff_noY - b_coeff_noY*c_val*d_coeff_noY
        epsilon_noY = e * d_coeff_noY**2
        
        Delta_noY = (256 * alpha_noY**3 * epsilon_noY**3 
                 - 192 * alpha_noY**2 * beta_noY * delta_noY * epsilon_noY**2 
                 - 128 * alpha_noY**2 * gamma_noY**2 * epsilon_noY**2 
                 + 144 * alpha_noY**2 * gamma_noY * delta_noY**2 * epsilon_noY 
                 - 27 * alpha_noY**2 * delta_noY**4 
                 + 144 * alpha_noY * beta_noY**2 * gamma_noY * epsilon_noY**2 
                 - 6 * alpha_noY * beta_noY**2 * delta_noY**2 * epsilon_noY 
                 - 80 * alpha_noY * beta_noY * gamma_noY**2 * delta_noY * epsilon_noY 
                 + 18 * alpha_noY * beta_noY * gamma_noY * delta_noY**3 
                 + 16 * alpha_noY * gamma_noY**4 * epsilon_noY 
                 - 4 * alpha_noY * gamma_noY**3 * delta_noY**2 
                 - 27 * beta_noY**4 * epsilon_noY**2 
                 + 18 * beta_noY**3 * gamma_noY * delta_noY * epsilon_noY 
                 - 4 * beta_noY**3 * delta_noY**3 
                 - 4 * beta_noY**2 * gamma_noY**3 * epsilon_noY 
                 + beta_noY**2 * gamma_noY**2 * delta_noY**2)

        P_noY = 8*alpha_noY*gamma_noY - 3*beta_noY**2
        D0_noY = 64*alpha_noY**3*epsilon_noY - 16*alpha_noY**2*gamma_noY**2 + 16*alpha_noY*beta_noY**2*gamma_noY - 16*alpha_noY**2*beta_noY*delta_noY - 3*beta_noY**4

        roots_4_noY = (Delta_noY > 0) & (D0_noY < 0) & (P_noY < 0)
        roots_2_noY = (Delta_noY < 0)
        has_real_roots_noY = roots_4_noY | roots_2_noY

        # ---------------------------------------------------------
        # 4. DESCARTES' RULE OF SIGNS 
        # ---------------------------------------------------------
        coeffs_noY = np.stack([epsilon_noY, delta_noY, gamma_noY, beta_noY, alpha_noY], axis=-1)
        sign_changes_noY = count_sign_changes(coeffs_noY)

        # ---------------------------------------------------------
        # 5. NEW ANALYTICAL ROOTS (Lower Sign Hardcoded)
        # ---------------------------------------------------------
        # Lower sign: \mp -> + , \pm -> -
        P1 = -18 + Nf**2 * (-3 + 4*Ng) + 4*Nf * (3*Ng + p)
        P2 = -3 + Nf**2 + 2*Nf * KminVal
        P3 = 9 - 20*Nf**2 + Nf**4
        P4 = 567 + 315*Nf**2 - 35*Nf**4 + Nf**6

        B_term = (-3 + Nf**2) * (Nf*(15 + 4*Ng) + 4*(3*Ng + p)) + \
                 (72 + Nf**2*(42 - 8*Ng) - 8*Nf*(3*Ng + p)) * KminVal #mistake on 3Ng should be pos
        
        C_term = 2 * P1**2 * KminVal - 27 * Nf * P2**2
        
        D_term = 23328 * Nf * P2**2 - 12 * Nf * B_term**2 - 48 * (3 - 2*Nf**2) * C_term + \
                 (7 + Nf**2) * (Nf * B_term**2 + 8 * (9 - Nf**2) * C_term)
        Denom = Nf**6*(Nf*(-21+4*Ng)+4*(3*Ng+p))**12*P2**12

        # Note: Assumed a minus sign before 432 based on polynomial alternation
        Delta0_new = 1/Denom*(
            - B_term**4 * D_term**2 * Nf**2 * (1 + Nf**2)**2 * P3**2 
            + B_term**2 * D_term**3 * Nf * (7 + Nf**2) * P3**2 
            - 64 * C_term * B_term**6 * Nf**3 * (1 + Nf**2)**3 * P3**3 
            + 72 * C_term * B_term**4 * D_term * Nf**2 * (1 + Nf**2) * (7 + Nf**2) * P3**3 
            + 432 * C_term**2 * B_term**4 * Nf**2 * (7 + Nf**2)**2 * P3**4 
            + B_term**2 * D_term**3 * Nf * (1 + Nf**2)**2 * P4 
            - D_term**4 * (7 + Nf**2) * P4 
            + 72 * C_term * B_term**4 * D_term * Nf**2 * (1 + Nf**2)**3 * P3 * P4 
            - 80 * C_term * B_term**2 * D_term**2 * Nf * (1 + Nf**2) * (7 + Nf**2) * P3 * P4 
            + 96 * C_term**2 * B_term**4 * Nf**2 * (1 + Nf**2)**2 * (7 + Nf**2) * P3**2 * P4 
            - 576 * C_term**2 * B_term**2 * D_term * Nf * (7 + Nf**2)**2 * P3**2 * P4 
            + 432 * C_term**2 * B_term**4 * Nf**2 * (1 + Nf**2)**4 * P4**2 
            - 576 * C_term**2 * B_term**2 * D_term * Nf * (1 + Nf**2)**2 * (7 + Nf**2) * P4**2 
            + 128 * C_term**2 * D_term**2 * (7 + Nf**2)**2 * P4**2 
            - 3072 * C_term**3 * B_term**2 * Nf * (1 + Nf**2) * (7 + Nf**2)**2 * P3 * P4**2 
            - 4096 * C_term**4 * (7 + Nf**2)**3 * P4**3
        )
        # ---------------------------------------------------------
        # NEW ANALYTICAL ROOTS: NO YUKAWA (Lower Sign Hardcoded)
        # ---------------------------------------------------------
        # Pre-compute powers of Nf to avoid redundant array calculations
        Nf2 = Nf**2
        Nf3 = Nf**3
        Nf4 = Nf2**2
        Nf6 = Nf4 * Nf2
        Nf8 = Nf4**2
        Nf10 = Nf8 * Nf2
        Nf12 = Nf6**2
        
        # Lower sign evaluated for \mp: becomes +
        A_term = Nf * (-21 + 4 * Ng) + 4 * (3 * Ng + p)
        
        P0_noY = -217326564 - 56077596*Nf2 + 23223267*Nf4 - 2630538*Nf6 + 166727*Nf8 - 9548*Nf10 + 188*Nf12
        P1_noY = -33434856 + 10958328*Nf2 + 7518420*Nf4 - 1773243*Nf6 + 171254*Nf8 - 9917*Nf10 + 350*Nf12
        P2_noY = -64297800 + 83041848*Nf2 + 15587964*Nf4 - 10068291*Nf6 + 1385624*Nf8 - 104921*Nf10 + 4472*Nf12
        P3_noY = -571536 + 1910304*Nf2 - 245412*Nf4 - 262236*Nf6 + 52675*Nf8 - 5338*Nf10 + 247*Nf12
        P4_noY = -1714608 + 14362272*Nf2 - 8214156*Nf4 - 2139588*Nf6 + 725257*Nf8 - 98878*Nf10 + 4837*Nf12
        P5_noY = 35964 - 46278*Nf2 - 2949*Nf4 + 3476*Nf6 - 629*Nf8 + 32*Nf10
        P6_noY = 107892 - 365634*Nf2 + 31473*Nf4 + 27116*Nf6 - 6143*Nf8 + 320*Nf10
        P7_noY = (-5 + Nf2)**2 * (-81 - 18*Nf2 + 2*Nf4)

        Delta0_noY_analytical = (
            - 2176782336 * Nf4 * P0_noY 
            - 725594112 * A_term * Nf3 * P1_noY 
            - 10077696 * A_term**2 * Nf2 * P2_noY 
            - 15116544 * A_term**3 * Nf * P3_noY 
            - 34992 * A_term**4 * P4_noY 
            - 139968 * A_term**5 * Nf * P5_noY 
            - 216 * A_term**6 * P6_noY 
            - 288 * A_term**7 * Nf * P7_noY 
            - A_term**8 * P7_noY
        )

        # ---------------------------------------------------------
        # 6. BOOLEAN LOGIC FOR FINAL COMBINED PLOT
        # ---------------------------------------------------------
        gauge_cond = (b0 < 0) & ((b0 - c1) > 0)
        
        CAF_fixed_flow = gauge_cond & has_real_roots_fixed
        CAF_noY_flow = gauge_cond & has_real_roots_noY

    # ==========================================================
    # PLOT 1: ROOTS ON FIXED FLOW
    # ==========================================================
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    ax1.contourf(Nf, p, roots_2_fixed, levels=[0.5, 1.5], colors=[color_roots_2], antialiased=True)
    ax1.contourf(Nf, p, roots_4_fixed, levels=[0.5, 1.5], colors=[color_roots_4], antialiased=True)
    ax1.contour(Nf, p, b0, levels=[0], colors='black', linewidths=2.5, linestyles='dashed')
    #ax1.contour(Nf, p, Delta_fixed, levels=[0], colors='#e41a1c', linewidths=2.5, linestyles='solid') 

    ax1.set_xlabel('$N$', fontsize=32); ax1.set_ylabel('$p$', fontsize=32)
    ax1.set_xlim(3, upNf); ax1.set_ylim(0, upp)
    ax1.set_xticks([5, 10, 15, 20, 25, 30]); ax1.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax1.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True, pad=12)

    leg1 = [Patch(facecolor=color_roots_2, label='2 real, 2 complex roots'),
            Patch(facecolor='white', edgecolor='black', label='4 complex roots'),
            Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b_0=0$'),
            #Line2D([0], [0], color='#e41a1c', linewidth=2.5, linestyle='solid', label=r'$\Delta=0$')
    ]
    ax1.legend(handles=leg1, loc='lower right', fontsize=22, frameon=True, 
               edgecolor='black', fancybox=False, title=fr'\textbf{{BY model ($N_g={Ng}$)}}', title_fontsize=24)
    
    plt.tight_layout()
    plt.savefig(f"CAFplot_1_RootsFixed_Ng{Ng}.pdf", format='pdf', bbox_inches='tight')
    plt.close(fig1)

    # ==========================================================
    # PLOT 2: ROOTS WITH NO YUKAWA
    # ==========================================================
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    ax2.contourf(Nf, p, roots_2_noY, levels=[0.5, 1.5], colors=[color_roots_2], antialiased=True)
    ax2.contourf(Nf, p, roots_4_noY, levels=[0.5, 1.5], colors=[color_roots_4], antialiased=True)
    ax2.contour(Nf, p, b0, levels=[0], colors='black', linewidths=2.5, linestyles='dashed')
    #ax2.contour(Nf, p, Delta_noY, levels=[0], colors='#984ea3', linewidths=2.5, linestyles='solid') 

    ax2.set_xlabel('$N$', fontsize=32); ax2.set_ylabel('$p$', fontsize=32)
    ax2.set_xlim(3, upNf); ax2.set_ylim(0, upp)
    ax2.set_xticks([5, 10, 15, 20, 25, 30]); ax2.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax2.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True, pad=12)

    leg2 = [Patch(facecolor=color_roots_4, label='4 real roots'),
            Patch(facecolor=color_roots_2, label='2 real, 2 complex roots'),
            Patch(facecolor='white', edgecolor='black', label='4 complex roots'),
            Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b_0=0$'),
            #Line2D([0], [0], color='#984ea3', linewidth=2.5, linestyle='solid', label=r'$\Delta=0$')
    ]
    ax2.legend(handles=leg2, loc='lower right', fontsize=22, frameon=True, 
               edgecolor='black', fancybox=False, title=fr'\textbf{{BY model ($N_g={Ng}$)}}', title_fontsize=24)
    
    plt.tight_layout()
    plt.savefig(f"CAFplot_2_RootsNoY_Ng{Ng}.pdf", format='pdf', bbox_inches='tight')
    plt.close(fig2)

    # ==========================================================
    # PLOT 3: FINAL COMBINED CAF (With New Analytical Line)
    # ==========================================================
    fig3, ax3 = plt.subplots(figsize=(10, 8))
    
    ax3.contourf(Nf, p, CAF_fixed_flow, levels=[0.5, 1.5], colors=[color_caf_flow], antialiased=True)
    ax3.contourf(Nf, p, CAF_noY_flow, levels=[0.5, 1.5], colors=[color_caf_free], antialiased=True)

    intersection = CAF_fixed_flow & CAF_noY_flow
    ax3.contourf(Nf, p, CAF_fixed_flow, levels=[0.5, 1.5], colors=[color_caf_flow], antialiased=True)
    
    stripe_mask_tight = ((p + 3.0 * Nf) % 0.20) < 0.8  
    striped_intersection = intersection & stripe_mask_tight

    ax3.contourf(Nf, p, striped_intersection, levels=[0.5, 1.5], colors=[color_caf_free], antialiased=True)
    
    ax3.contour(Nf, p, b0, levels=[0], colors='black', linewidths=2.5, linestyles='dashed')
    #ax3.contour(Nf, p, Delta_fixed, levels=[0], colors='#e41a1c', linewidths=2.5, linestyles='solid') # Red
    ax3.contour(Nf, p, Delta_noY, levels=[0], colors="black", linewidths=2.5, linestyles='solid')   # Purple
    
    # NEW: Analytical Delta0 line
    ax3.contour(Nf, p, Delta0_new, levels=[0], colors='black', linewidths=3.0, linestyles='dashdot') # Orange/Gold Dotted
    #ax3.contour(Nf, p, Delta0_noY_analytical, levels=[0], colors="#187526", linewidths=3.0, linestyles='solid')

    ax3.set_xlabel('$N$', fontsize=32); ax3.set_ylabel('$p$', fontsize=32)
    ax3.set_xlim(3, upNf); ax3.set_ylim(0, upp)
    ax3.set_xticks([5, 10, 15, 20, 25, 30]); ax3.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax3.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True, pad=12)

    leg3 = [
        Patch(facecolor=color_caf_flow, label='On fixed flow'),
        Patch(facecolor=color_caf_free, label=r'$\alpha_y$ off fixed flow'),
        Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b_0=0$'),
        #Line2D([0], [0], color='#e41a1c', linewidth=2.5, linestyle='solid', label=r'$\Delta = 0$ (On fixed flow)'),
        Line2D([0], [0], color='black', linewidth=2.5, linestyle='solid', label=r'$\Delta = 0$ $\alpha_y$ off fixed flow '),
        Line2D([0], [0], color='black', linewidth=3.0, linestyle='dashdot', label=r'$\Delta = 0$ on fixed flow')
    ]
    ax3.legend(handles=leg3, loc='upper left', fontsize=22, frameon=True, 
               edgecolor='black', fancybox=False, title=fr'\textbf{{BY model ($N_g={Ng}$)}}', title_fontsize=24)

    plt.tight_layout()
    plt.savefig(f"CAFplot_3_FinalCombined_Ng{Ng}.pdf", format='pdf', bbox_inches='tight')
    plt.close(fig3)

    

    # ==========================================================
    # PLOT 4: DESCARTES RULE OF SIGNS (No Yukawa)
    # ==========================================================
    fig4, ax4 = plt.subplots(figsize=(10, 8))
    
    ax4.contourf(Nf, p, sign_changes_noY == 4, levels=[0.5, 1.5], colors=[color_roots_4], antialiased=True)
    ax4.contourf(Nf, p, sign_changes_noY == 3, levels=[0.5, 1.5], colors=[color_roots_3], antialiased=True)
    ax4.contourf(Nf, p, sign_changes_noY == 2, levels=[0.5, 1.5], colors=[color_roots_2], antialiased=True)
    ax4.contourf(Nf, p, sign_changes_noY == 1, levels=[0.5, 1.5], colors=[color_roots_1], antialiased=True)
    ax4.contourf(Nf, p, sign_changes_noY == 0, levels=[0.5, 1.5], colors=[color_roots_0], antialiased=True)
    
    ax4.contour(Nf, p, b0, levels=[0], colors='black', linewidths=2.5, linestyles='dashed')
    ax4.contour(Nf, p, Delta_noY, levels=[0], colors='#984ea3', linewidths=2.5, linestyles='solid')

    ax4.set_xlabel('$N$', fontsize=32); ax4.set_ylabel('$p$', fontsize=32)
    ax4.set_xlim(3, upNf); ax4.set_ylim(0, upp)
    ax4.set_xticks([5, 10, 15, 20, 25, 30]); ax4.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
    ax4.tick_params(axis='both', which='major', labelsize=26, width=1.5, length=8, direction='in', top=True, right=True, pad=12)

    leg4 = [Patch(facecolor='none', edgecolor='none', label=fr'\textbf{{BY Model}} $(N_g = {Ng})$'),
            Patch(facecolor=color_roots_4, label='4 sign changes'),
            Patch(facecolor=color_roots_3, label='3 sign changes'),
            Patch(facecolor=color_roots_2, label='2 sign changes'),
            Patch(facecolor=color_roots_1, label='1 sign change'),
            Patch(facecolor=color_roots_0, label='0 sign changes'),
            Line2D([0], [0], color='black', linewidth=2.5, linestyle='dashed', label=r'$b_0=0$'),
            Line2D([0], [0], color='#984ea3', linewidth=2.5, linestyle='solid', label=r'$\Delta=0$')]
    
    ax4.legend(handles=leg4, loc='lower right', fontsize=18, frameon=True, edgecolor='black')

    plt.tight_layout()
    plt.savefig(f"CAFplot_4_DescartesNoY_Ng{Ng}.pdf", format='pdf', bbox_inches='tight')
    plt.close(fig4)

print("All 20 plots generated successfully with final formatting!")

# ==========================================================
# EXECUTION CODE
# ==========================================================
print("\nCalculating absolute minimums for 'On fixed flow' and 'alpha_y off fixed flow'...")
minimums = find_absolute_minimums()

print("--- Absolute Minimums (N >= 5) ---")
for Ng, vals in minimums.items():
    print(f"\nModel Ng = {Ng}:")
    
    f_res = vals['fixed']  
    if f_res['N'] is not None:
        print(f"  [On fixed flow]        Minimum N = {f_res['N']:<2}, p = {f_res['p']}")
    else:
        print("  [On fixed flow]        No valid region found.")
        
    n_res = vals['noY']
    if n_res['N'] is not None:
        print(f"  [alpha_y off fixed flow] Minimum N = {n_res['N']:<2}, p = {n_res['p']}")
    else:
        print("  [alpha_y off fixed flow] No valid region found.")