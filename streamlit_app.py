import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

def titration_pH_strong_acid_strong_base(V_titrant, V_initial, C_initial, C_titrant):
    """
    Calculate the pH for a strong acid (e.g., HCl) titrated with a strong base (e.g., NaOH)
    using a simple 1:1 neutralization reaction.
    
    Parameters:
      V_titrant : numpy array
          Array of titrant volumes (in liters)
      V_initial : float
          Initial acid volume (in liters)
      C_initial : float
          Acid concentration (M)
      C_titrant : float
          Titrant (base) concentration (M)
          
    Returns:
      pH_vals : numpy array
          Calculated pH values at each titrant volume.
    """
    moles_acid = C_initial * V_initial
    pH_vals = []
    
    for Vb in V_titrant:
        moles_base = C_titrant * Vb
        V_total = V_initial + Vb
        
        if moles_base < moles_acid:
            # Before equivalence: excess acid
            H_conc = (moles_acid - moles_base) / V_total
            pH = -np.log10(H_conc)
        elif np.isclose(moles_base, moles_acid, rtol=1e-6):
            # At equivalence: pH ~7 for strong acid-strong base titration
            pH = 7.0
        else:
            # After equivalence: excess base
            OH_conc = (moles_base - moles_acid) / V_total
            pOH = -np.log10(OH_conc)
            pH = 14 - pOH
        pH_vals.append(pH)
    
    return np.array(pH_vals)

def main():
    st.title("Titration Curve: Strong Acid with Strong Base")
    
    # Set parameters for the titration:
    # For example, 50 mL of 0.1 M acid titrated with 0.1 M base.
    V_initial = 0.05       # Acid volume in liters (50 mL)
    C_initial = 0.1        # Acid concentration in M
    C_titrant = 0.1        # Base concentration in M
    
    # Generate an array of titrant volumes from 0 to 0.1 L (0 to 100 mL)
    V_titrant = np.linspace(0, 0.1, 500)
    pH_vals = titration_pH_strong_acid_strong_base(V_titrant, V_initial, C_initial, C_titrant)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(V_titrant * 1000, pH_vals, label="Titration Curve")
    ax.axvline(x=V_initial * 1000, color='grey', linestyle='--', label="Equivalence Point")
    ax.set_xlabel("Volume of Base added (mL)")
    ax.set_ylabel("pH")
    ax.set_title("Titration Curve: Strong Acid with Strong Base")
    ax.legend()
    ax.grid(True)
    
    st.pyplot(fig)

if __name__ == "__main__":
    main()
