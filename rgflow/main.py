""" Paper: Conformal Phase Diagram of Complete Asymptotically Free Theories"

To try and plot something similar to figure one to better understand the graph and how it is produced. 

"""
import numpy as np
import matplotlib.pyplot as plt


#Function which takes initial conditions and makes the fixed plot

def plot_fixed_flow(b0, c1, c2):
    #Require condition of equation 11
    if b0-c1>0:
        c_fixed = (b0 - c1) / 2

        # Range of gauge coupling (x - axis)
        alpha_g = np.linspace(0.001,2,1000)

        # Compute alpha_h values
        alpha_H_fixed = c_fixed * alpha_g

        # Create figure
        plt.figure(figsize=(8,6))
        plt.plot(alpha_g, alpha_H_fixed, "g-", linewidth = 2, label = "Fixed flow")
        plt.grid()
        plt.show()

    else: print("The constants didn't live up to the condition of eq. 11")



def plot_rg_flows_h_g(b0=-1, c1=-2, c2=1):
    """
    Plot 1-loop RG flow for gauge (α_g) and Yukawa (α_H) couplings.
    dα_g/dlnμ = b0 * α_g²
    dα_H/dlnμ = (c1 + c2 * α_H/α_g) * α_H * α_g
    """
    # Fixed ratio (α_H / α_g)
    c_fixed = (b0 - c1) / c2
    if c_fixed <= 0:
        print("The constants didn't live up to the condition of eq. 11")
        return
    
    # Grid for flow field - every point in the meshgrid is possible values of the coupling constants (α_g, α_H)
    αg = np.linspace(0.001, 1, 500)
    αH = np.linspace(0.001, 1, 500)
    g, H = np.meshgrid(αg, αH)

    # RG flow equations using the points of the meshgrid
    dαg = - b0 * g**2
    dαH = - (c1 + c2 * (H / g)) * H * g

    # Normalize arrows for clarity
    norm = np.sqrt(dαg**2 + dαH**2)
    dαg /= norm
    dαH /= norm

    # --- Plot ---
    plt.figure(figsize=(8, 6))
    plt.streamplot(g, H, dαg, dαH, color='lightsteelblue',
                   density=1.4, linewidth=1.2, arrowsize=1)

    # Fixed flow line
    αg_line = np.linspace(0.001, 1, 500)
    αH_fixed = c_fixed * αg_line
    plt.plot(αg_line, αH_fixed, 'g-', lw=2, label=fr'Fixed flow')
    plt.plot(0,0, 'ro', markersize = 10, label='n fixed point')

    # Axis labels and layout
    plt.xlabel(r'Gauge coupling $\alpha_g$', fontsize=13)
    plt.ylabel(r'Yukawa coupling $\alpha_H$', fontsize=13)
    plt.title(fr'1-loop RG Flow  ($b_0={b0}$, $c_1={c1}$, $c_2={c2}$)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('H_g_rg_flow_plot.png', dpi=300)
    plt.show()



def plot_rg_flows_l_g(b0=-1, d1=2, d2=-4, d4=1/2):
    """
    Plot 1-loop RG flow for gauge (α_g) and scalar (α_lambda) couplings.
    dα_g/dlnμ = b0 * α_g²
    dα_lambda/dlnμ = (1 / b0) * (α_lambda / α_g) * ( d1 * (α_lambda / α_g) + d2 ) + (d4 / b0)
    """
    # Fixed ratio (α_lambda / α_g)
    k = (b0 - d2)**2 - 4 * d1 * d4
    if k < 0:
        print("The constants didn't live up to the condition k>0")
        return
    c_fixed1 = (b0 - d2 + np.sqrt(k)) / (2 * d1)
    c_fixed2 = (b0 - d2 - np.sqrt(k)) / (2 * d1)
    
    # Grid for flow field - every point in the meshgrid is possible values of the coupling constants (α_g, α_H)
    αg = np.linspace(0.001, 1, 500)
    α_lambda = np.linspace(-1, 1, 500)
    g, l = np.meshgrid(αg, α_lambda)

    # RG flow equations using the points of the meshgrid
    dαg = - b0 * g**2
    dα_lambda =  (1 / b0) * (l / g) * ( d1 * (l / g) + d2 ) + (d4 / b0)

    # Normalize arrows for clarity
    norm = np.sqrt(dαg**2 + dα_lambda**2)
    dαg /= norm
    dα_lambda /= norm

    # --- Plot ---
    plt.figure(figsize=(8, 6))
    plt.streamplot(g, l, dαg, dα_lambda, color='lightsteelblue',
                   density=1.4, linewidth=1.2, arrowsize=1)

    # Fixed flow line
    αg_line = np.linspace(0, 1, 500)
    α_lambda_fixed1 = c_fixed1 * αg_line
    α_lambda_fixed2 = c_fixed2 * αg_line
    plt.plot(αg_line, α_lambda_fixed1, 'g-', lw=2, label=fr'Fixed flow')
    plt.plot(αg_line, α_lambda_fixed2, 'g--', lw=2, label=fr'Fixed flow')
    plt.plot(0,0, 'ro', markersize = 10, label='n fixed point')

    # Axis labels and layout
    plt.xlabel(r'Gauge coupling $\alpha_g$', fontsize=13)
    plt.ylabel(r'Scalar coupling $\alpha_\lambda$', fontsize=13)
    plt.ylim(-1, 1)
    plt.title(fr'1-loop RG Flow  ($b_0={b0}$, $d_1={d1}$, $d_2={d2}$, $d_4={d4}$)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('l_g_rg_flow_plot.png', dpi=300)
    plt.show()


def plot_rg_flows_g_chiral():
    """ I have considered the constant b0 = 0 and then isolated for p the critical number of flavors for asymptotic freedom. That means that for p < p_critical the theory is AF and for p > p_critical the theory is IR free (I think)."""
    # Linspace of N
    nn = np.linspace(0, 20, 400)
    p_v = 11 * nn * 3 / 4
    p_BY = (6 * nn + 4) * 3 / 4 # From BY model
    p_GG = (6 * nn - 4) * 3 / 4  # From GG model

    # --- Plot ---
    plt.figure(figsize=(8, 6))
    plt.plot(nn, p_v, 'b-', lw=2, label='Vector-like fermions')
    plt.plot(nn, p_BY, 'r--', lw=2, label='Chiral fermions (BY model)')
    plt.plot(nn, p_GG, 'g-.', lw=2, label='Chiral fermions (GG model)')
    # Axis labels and layout
    plt.xlabel(r'Number of colors $N$', fontsize=13)
    plt.ylabel(r'Critical number of flavors $p$', fontsize=13)
    plt.xlim(5, 20)
    plt.ylim(0, 120)
    plt.title(r'Critical number of flavors for asymptotic freedom chiral vs vector like SU(N)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('critical_flow_plot_c_v.png', dpi=300)
    plt.show()

    
def plot_rg_flows_h_gchiral(N = 5, p = 25):
    """
    Plot 1-loop RG flow for gauge (α_g - chiral) and yukawa (α_H) couplings.
    dα_g/dlnμ = b0 * α_g² = - (6N - 4/3p + 4) * α_g² (GG)
    dα_g/dlnμ = b0 * α_g² = - (6N + 4/3p - 4) * α_g² (BY)
    dα_H/dlnμ = α_H * ( (c_1 * α_g) + (c_2 * α_H) ) = α_g * α_H * (-9 * N - 6 + 15/N ) + α_H² * (3*N / 2 + 1 + 3/2) (GG)
    dα_H/dlnμ = α_H * ( (c_1 * α_g) + (c_2 * α_H) ) = α_g * α_H * (-9 * N + 6 + 15/N ) + α_H² * (3*N / 2 + 1 - 3/2) (BY)
    """
    #Define constants for the models
    # GG model
    b0_GG =  (6*N - 4/3*p + 4)
    c1_GG = -9 * N - 6 + 15/N
    c2_GG = 3*N / 2 + 1 + 3/2
    # BY model
    b0_BY =  (6*N + 4/3*p - 4)
    c1_BY = -9 * N + 6 + 15/N
    c2_BY = 3*N / 2 + 1 - 3/2
    
    # Fixed ratio (α_H / α_g) 
    c_fixed_GG = (-b0_GG - c1_GG) / c2_GG
    c_fixed_BY = (-b0_BY - c1_BY) / c2_BY

    # If function ensuring that the model is asymptotically free (p<(6N -+4)3/4) so when p> (6N -+4)3/4 we break the loop - less than or larger than?
    #if p > b0_BY: # From BY model
     #   return print("The constants didn't live up to the condition of AF - BY")
    #if p > b0_GG:  # From GG model 
     #   return print("The constants didn't live up to the condition of AF - GG")

    # Grid for flow field - every point in the meshgrid is possible values of the coupling constants (α_g, α_H)
    αg = np.linspace(0.001, 1, 500)
    α_H = np.linspace(0.001, 1, 500)
    g, H = np.meshgrid(αg, α_H)

    # RG flow equations using the points of the meshgrid for GG model
    dαg_GG =  b0_GG * g**2 
    dα_H_GG = - g * H * c1_GG - H**2 * c2_GG

    # RG flow equations using the points of the meshgrid for BY model
    dαg_BY =  b0_BY * g**2 
    dα_H_BY = - g * H * c1_BY - H**2 * c2_BY

    # Normalize arrows for clarity GG
    norm_GG = np.sqrt(dαg_GG**2 + dα_H_GG**2)
    dαg_GG /= norm_GG
    dα_H_GG /= norm_GG

    # Normalize arrows for clarity BY
    norm_BY = np.sqrt(dαg_BY**2 + dα_H_BY**2)
    dαg_BY /= norm_BY
    dα_H_BY /= norm_BY

    # --- Plot ---
    plt.figure(figsize=(8, 6))
    plt.streamplot(g, H, dαg_GG, dα_H_GG, color='lightsteelblue',
                   density=1.4, linewidth=1.2, arrowsize=1)


    # Fixed flow line
    αg_line = np.linspace(0.001, 1, 500)
    αH_fixed_GG = c_fixed_GG * αg_line
    plt.plot(αg_line, αH_fixed_GG, 'g--', lw=2, label=fr'Fixed flow')
    plt.plot(0,0, 'ro', markersize = 10, label='n fixed point')

    # Axis labels and layout
    plt.xlabel(r'Gauge coupling $\alpha_g$', fontsize=13)
    plt.ylabel(r'Yukawa coupling $\alpha_H$', fontsize=13)
    plt.ylim(-0.1, 1)
    plt.title(fr'1-loop RG Flow GG chiral theory ($N={N}$, $p={p}$)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'H_gCC_rg_flow_plotN{N}p{str(p).replace(".", "_")}.png', dpi=300)
    plt.show()

    # --- Plot ---
    plt.figure(figsize=(8, 6))
    plt.streamplot(g, H, dαg_BY, dα_H_BY, color='lightsteelblue',
                   density=1.4, linewidth=1.2, arrowsize=1)
    plt.plot(0,0, 'ro', markersize = 10, label='n fixed point')
    
    # Fixed flow line
    αg_line = np.linspace(0.001, 1, 500)
    αH_fixed_BY = c_fixed_BY * αg_line
    plt.plot(αg_line, αH_fixed_BY, 'g--', lw=2, label=fr'Fixed flow')

    # Axis labels and layout
    plt.xlabel(r'Gauge coupling $\alpha_g$', fontsize=13)
    plt.ylabel(r'Scalar coupling $\alpha_H$', fontsize=13)
    plt.ylim(-0.1, 1)
    plt.title(fr'1-loop RG Flow BY chiral theory ($N={N}$, $p={p}$)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('H_gBY_rg_flow_plot.png', dpi=300)
    plt.show()  



# Run the plots
if __name__ == "__main__":
    #plot_fixed_flow(b0 = -1, c1 = -2, c2 = 1)
    #plot_rg_flows_h_g(b0=-1, c1=-2, c2=1)
    #plot_rg_flows_l_g(b0=-1, d1=2, d2=-4, d4=1/2)
    #plot_rg_flows_g_chiral()
    plot_rg_flows_h_gchiral(N = 7, p = 35) # Change N and p to see different flows , N=7, p=35 makes p>p_critical for GG
    plot_rg_flows_h_gchiral(N = 7, p = 34)
    plot_rg_flows_h_gchiral(N = 7, p = 20)
