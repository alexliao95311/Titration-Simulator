import numpy as np
import matplotlib.pyplot as plt

# Acid & Base Dissociation Constants (Ka/Kb)
acids = {
    "HCl": 1e6,   # Strong acid
    "CH3COOH": 1.75e-5,  # Weak acid (Acetic Acid)
    "H2SO4": 1.2e3   # Strong Acid (Sulfuric)
}

bases = {
    "NaOH": 1e6,   # Strong base
    "NH3": 1.8e-5   # Weak base (Ammonia)
}

# User selects acids/bases
print("Available acids:", list(acids.keys()))
print("Available bases:", list(bases.keys()))

titrant = input("Enter the titrant (acid or base): ")
analyte = input("Enter the analyte (acid or base): ")

if titrant in acids and analyte in bases:
    titrant_type = "acid"
    analyte_type = "base"
    Ka_or_Kb = acids[titrant]
elif titrant in bases and analyte in acids:
    titrant_type = "base"
    analyte_type = "acid"
    Ka_or_Kb = bases[titrant]
else:
    print("Invalid choices! Please restart.")
    exit()

# User inputs known values
titrant_conc = float(input(f"Enter {titrant} concentration (M): "))
titrant_vol = float(input(f"Enter {titrant} volume (mL): "))
analyte_vol = float(input(f"Enter {analyte} volume (mL): "))

# Calculate the analyte concentration (since we have known titrant data)
analyte_conc = (titrant_conc * titrant_vol) / analyte_vol
print(f"Calculated {analyte} concentration: {analyte_conc:.4f} M")

# Set up titration calculations
volume_titrant_added = np.linspace(0, titrant_vol * 2, 100)
pH_values = []

for Vt in volume_titrant_added:
    # Moles of analyte & titrant
    moles_analyte = analyte_conc * (analyte_vol / 1000)
    moles_titrant = titrant_conc * (Vt / 1000)

    if titrant_type == "acid":
        if moles_titrant < moles_analyte:
            remaining_base = moles_analyte - moles_titrant
            OH_conc = remaining_base / (analyte_vol + Vt) * 1000
            pH_values.append(14 + np.log10(OH_conc))
        else:
            excess_acid = moles_titrant - moles_analyte
            H_conc = excess_acid / (analyte_vol + Vt) * 1000
            pH_values.append(-np.log10(H_conc))
    else:
        if moles_titrant < moles_analyte:
            remaining_acid = moles_analyte - moles_titrant
            H_conc = np.sqrt(Ka_or_Kb * (remaining_acid / (analyte_vol + Vt) * 1000)) if Ka_or_Kb < 1e5 else (remaining_acid / (analyte_vol + Vt) * 1000)
            pH_values.append(-np.log10(H_conc))
        else:
            excess_base = moles_titrant - moles_analyte
            OH_conc = excess_base / (analyte_vol + Vt) * 1000
            pH_values.append(14 + np.log10(OH_conc))

# Plot titration curve
plt.plot(volume_titrant_added, pH_values, label=f"Titration: {titrant} → {analyte}")
plt.axhline(y=7, color='gray', linestyle='--', label="Neutral pH")
plt.xlabel("Volume of Titrant Added (mL)")
plt.ylabel("pH")
plt.title("Titration Curve")
plt.legend()
plt.grid()
plt.show()
