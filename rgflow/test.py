import numpy as np
import matplotlib.pyplot as plt
def plot_rg_flows_g_chiral():
    """ Plot fix point for 1 loop bet function for non-abelian gauge (α_g) couplings. Comparison of chiral and vector-like fermions."""
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
    # Fixed ratio (α_H / α_g) ?

    # If function ensuring that the model is asymptotically free - less than or larger than?
    if p < (6 * N + 4 ) * 3 / 4: # From BY model
        return print("The constants didn't live up to the condition of AF")
    if p < (6 * N - 4) * 3 / 4:  # From GG model 
        return print("The constants didn't live up to the condition of AF")

    # Grid for flow field - every point in the meshgrid is possible values of the coupling constants (α_g, α_H)
    αg = np.linspace(0.001, 1, 500)
    α_H = np.linspace(0.001, 1, 500)
    g, H = np.meshgrid(αg, α_H)

    # RG flow equations using the points of the meshgrid for GG model
    dαg_GG = - (6*N - 4/3*p + 4) * g**2 
    dα_H_GG =  g * H * (-9 * N - 6 + 15/N ) + H**2 * (3*N / 2 + 1 + 3/2)

    # RG flow equations using the points of the meshgrid for BY model
    dαg_BY = - (6*N + 4/3*p - 4) * g**2 
    dα_H_BY =  g * H * (-9 * N + 6 + 15/N ) + H**2 * (3*N / 2 + 1 - 3/2)

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

    # Axis labels and layout
    plt.xlabel(r'Gauge coupling $\alpha_g$', fontsize=13)
    plt.ylabel(r'Scalar coupling $\alpha_H$', fontsize=13)
    plt.ylim(-1, 1)
    plt.title(fr'1-loop RG Flow GG chiral theory ($N={N}$, $p={p}$)', fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('H_gCC_rg_flow_plot.png', dpi=300)
    plt.show()    
# Run the plots
if __name__ == "__beta_chiral_vs_vector__":
    plot_rg_flows_g_chiral()
    plot_rg_flows_h_gchiral(N = 5, p = 25)
    