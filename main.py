import numpy as np
import matplotlib.pyplot as plt
import csv

# Acid & Base Dissociation Constants (Ka/Kb)
acids = {
    "HCl": 1e6,         # Strong acid
    "H2SO4": 1e3,       # Strong acid (first dissociation)
    "H3PO4": 7.5e-3,    # Phosphoric acid (first dissociation)
    "CH3COOH": 1.75e-5, # Acetic acid
    "H2CO3": 4.3e-7,    # Carbonic acid
    "NH4+": 5.6e-10     # Ammonium ion
}

bases = {
    "NaOH": 1e6,        # Strong base
    "KOH": 1e6,         # Strong base
    "Ca(OH)2": 1e6,     # Strong base
    "NH3": 1.8e-5,      # Ammonia
    "HCO3-": 5.6e-11,   # Bicarbonate
    "CH3COO-": 5.7e-10  # Acetate
}

# Corresponding conjugate pairs
conjugate_pairs = {
    "CH3COOH": "CH3COO-",
    "CH3COO-": "CH3COOH",
    "NH3": "NH4+",
    "NH4+": "NH3",
    "H2CO3": "HCO3-",
    "HCO3-": "H2CO3"
}

# Constants
Kw = 1.0e-14  # Water dissociation constant

def get_user_input():
    """Get user input for titration parameters"""
    print("\n--- TITRATION SIMULATOR ---")
    print("Available acids:", list(acids.keys()))
    print("Available bases:", list(bases.keys()))
    
    titrant = input("Enter the titrant: ")
    analyte = input("Enter the analyte: ")
    
    # Determine if this is an acid-base or base-acid titration
    if titrant in acids and analyte in bases:
        titrant_type = "acid"
        analyte_type = "base"
        Ka_titrant = acids[titrant]
        Kb_analyte = bases[analyte]
    elif titrant in bases and analyte in acids:
        titrant_type = "base"
        analyte_type = "acid"
        Kb_titrant = bases[titrant]
        Ka_analyte = acids[analyte]
    else:
        print("Invalid choices! Please restart.")
        exit()
    
    # Get concentrations and volumes
    titrant_conc = float(input(f"Enter {titrant} concentration (M): "))
    analyte_conc = float(input(f"Enter {analyte} concentration (M): "))
    analyte_vol = float(input(f"Enter {analyte} volume (mL): "))
    max_titrant_vol = float(input(f"Enter maximum {titrant} volume to simulate (mL): "))
    
    return {
        'titrant': titrant,
        'analyte': analyte,
        'titrant_type': titrant_type,
        'analyte_type': analyte_type,
        'titrant_conc': titrant_conc,
        'analyte_conc': analyte_conc,
        'analyte_vol': analyte_vol,
        'max_titrant_vol': max_titrant_vol
    }

def calculate_pH_strong_acid(concentration):
    """Calculate pH of a strong acid solution"""
    if concentration <= 0:
        return 7.0
    elif concentration < 1e-6:  # Very dilute solutions need water dissociation
        h_conc = (concentration + np.sqrt(concentration**2 + 4*Kw))/2
        return -np.log10(h_conc)
    else:
        return -np.log10(concentration)

def calculate_pH_strong_base(concentration):
    """Calculate pH of a strong base solution"""
    if concentration <= 0:
        return 7.0
    elif concentration < 1e-6:  # Very dilute solutions need water dissociation
        oh_conc = (concentration + np.sqrt(concentration**2 + 4*Kw))/2
        h_conc = Kw / oh_conc
        return -np.log10(h_conc)
    else:
        h_conc = Kw / concentration
        return -np.log10(h_conc)

def calculate_pH_weak_acid(concentration, Ka):
    """Calculate pH of a weak acid solution"""
    if concentration <= 0:
        return 7.0
    # Solve quadratic equation for [H+]: [H+]^2 + Ka[H+] - Ka*concentration = 0
    h_conc = (-Ka + np.sqrt(Ka**2 + 4*Ka*concentration)) / 2
    return -np.log10(h_conc)

def calculate_pH_weak_base(concentration, Kb):
    """Calculate pH of a weak base solution"""
    if concentration <= 0:
        return 7.0
    # Calculate [OH-] first
    oh_conc = (-Kb + np.sqrt(Kb**2 + 4*Kb*concentration)) / 2
    h_conc = Kw / oh_conc
    return -np.log10(h_conc)

def calculate_buffer_pH(acid_conc, base_conc, Ka):
    """Calculate pH of a buffer solution using Henderson-Hasselbalch"""
    if acid_conc <= 0 or base_conc <= 0:
        if acid_conc > 0:
            return calculate_pH_weak_acid(acid_conc, Ka)
        elif base_conc > 0:
            Kb = Kw / Ka
            return calculate_pH_weak_base(base_conc, Kb)
        else:
            return 7.0
            
    pKa = -np.log10(Ka)
    return pKa + np.log10(base_conc / acid_conc)

def simulate_titration(params):
    """Simulate titration and calculate pH values"""
    titrant = params['titrant']
    analyte = params['analyte']
    titrant_type = params['titrant_type']
    analyte_type = params['analyte_type']
    titrant_conc = params['titrant_conc']
    analyte_conc = params['analyte_conc']
    analyte_vol = params['analyte_vol']
    max_titrant_vol = params['max_titrant_vol']
    
    # Set up titration calculations
    volume_titrant_added = np.linspace(0, max_titrant_vol, 200)  # More points for smoother curve
    pH_values = []
    
    # Calculate equivalence point volume
    equiv_vol = (analyte_conc * analyte_vol) / titrant_conc
    
    # Get Ka/Kb values
    if titrant_type == "acid":
        Ka_titrant = acids[titrant]
        is_titrant_strong = Ka_titrant > 1e0
        
        Kb_analyte = bases[analyte]
        is_analyte_strong = Kb_analyte > 1e0
        
        if not is_analyte_strong and analyte in conjugate_pairs:
            conjugate_acid = conjugate_pairs[analyte]
            Ka_conjugate = acids.get(conjugate_acid, Kw / Kb_analyte)
        else:
            Ka_conjugate = Kw / Kb_analyte
    else:  # titrant is base
        Kb_titrant = bases[titrant]
        is_titrant_strong = Kb_titrant > 1e0
        
        Ka_analyte = acids[analyte]
        is_analyte_strong = Ka_analyte > 1e0
        
        if not is_analyte_strong and analyte in conjugate_pairs:
            conjugate_base = conjugate_pairs[analyte]
            Kb_conjugate = bases.get(conjugate_base, Kw / Ka_analyte)
        else:
            Kb_conjugate = Kw / Ka_analyte
    
    # Initial pH of analyte solution before any titrant is added
    if analyte_type == "acid":
        if is_analyte_strong:
            initial_pH = calculate_pH_strong_acid(analyte_conc)
        else:
            initial_pH = calculate_pH_weak_acid(analyte_conc, Ka_analyte)
    else:  # analyte is base
        if is_analyte_strong:
            initial_pH = calculate_pH_strong_base(analyte_conc)
        else:
            initial_pH = calculate_pH_weak_base(analyte_conc, Kb_analyte)
    
    # For each volume of titrant added, calculate pH
    for Vt in volume_titrant_added:
        # Calculate moles
        moles_analyte = analyte_conc * (analyte_vol / 1000)  # Convert to L
        moles_titrant = titrant_conc * (Vt / 1000)           # Convert to L
        total_volume = (analyte_vol + Vt) / 1000             # Total volume in L
        
        if Vt == 0:
            pH_values.append(initial_pH)
            continue
            
        # ACID-BASE TITRATION
        if titrant_type == "acid" and analyte_type == "base":
            if moles_titrant < moles_analyte:  # Before equivalence
                if is_analyte_strong and is_titrant_strong:
                    # Strong acid + Strong base: excess base
                    excess_base = (moles_analyte - moles_titrant) / total_volume
                    pH = calculate_pH_strong_base(excess_base)
                elif not is_analyte_strong and is_titrant_strong:
                    # Strong acid + Weak base: buffer region
                    remaining_base = moles_analyte - moles_titrant
                    formed_conjugate = moles_titrant
                    pH = calculate_buffer_pH(formed_conjugate / total_volume, 
                                            remaining_base / total_volume, 
                                            Ka_conjugate)
                elif is_analyte_strong and not is_titrant_strong:
                    # Weak acid + Strong base: excess base with weak acid partially neutralized
                    excess_base = (moles_analyte - moles_titrant) / total_volume
                    pH = calculate_pH_strong_base(excess_base)
                else:
                    # Weak acid + Weak base: complex buffer
                    remaining_base = moles_analyte - moles_titrant
                    formed_conjugate = moles_titrant
                    pH = calculate_buffer_pH(formed_conjugate / total_volume, 
                                            remaining_base / total_volume, 
                                            Ka_conjugate)
            
            elif moles_titrant > moles_analyte:  # After equivalence
                if is_titrant_strong:
                    # Excess strong acid
                    excess_acid = (moles_titrant - moles_analyte) / total_volume
                    pH = calculate_pH_strong_acid(excess_acid)
                else:
                    # Excess weak acid
                    excess_acid = (moles_titrant - moles_analyte) / total_volume
                    pH = calculate_pH_weak_acid(excess_acid, Ka_titrant)
            
            else:  # At equivalence point
                if is_titrant_strong and is_analyte_strong:
                    # Strong acid + Strong base = neutral salt
                    pH = 7.0
                elif not is_titrant_strong and is_analyte_strong:
                    # Weak acid + Strong base = basic salt
                    salt_conc = moles_analyte / total_volume
                    Kb_salt = Kw / Ka_titrant
                    pH = calculate_pH_weak_base(salt_conc, Kb_salt)
                elif is_titrant_strong and not is_analyte_strong:
                    # Strong acid + Weak base = acidic salt
                    salt_conc = moles_analyte / total_volume
                    pH = calculate_pH_weak_acid(salt_conc, Ka_conjugate)
                else:
                    # Weak acid + Weak base = depends on relative Ka/Kb
                    salt_conc = moles_analyte / total_volume
                    if Ka_titrant > Kb_analyte:
                        pH = calculate_pH_weak_acid(salt_conc, Ka_conjugate)
                    else:
                        Kb_salt = Kw / Ka_titrant
                        pH = calculate_pH_weak_base(salt_conc, Kb_salt)
        
        # BASE-ACID TITRATION
        else:  # titrant is base, analyte is acid
            if moles_titrant < moles_analyte:  # Before equivalence
                if is_analyte_strong and is_titrant_strong:
                    # Strong base + Strong acid: excess acid
                    excess_acid = (moles_analyte - moles_titrant) / total_volume
                    pH = calculate_pH_strong_acid(excess_acid)
                elif not is_analyte_strong and is_titrant_strong:
                    # Strong base + Weak acid: buffer region
                    remaining_acid = moles_analyte - moles_titrant
                    formed_conjugate = moles_titrant
                    pH = calculate_buffer_pH(remaining_acid / total_volume, 
                                           formed_conjugate / total_volume, 
                                           Ka_analyte)
                elif is_analyte_strong and not is_titrant_strong:
                    # Weak base + Strong acid: excess acid
                    excess_acid = (moles_analyte - moles_titrant) / total_volume
                    pH = calculate_pH_strong_acid(excess_acid)
                else:
                    # Weak base + Weak acid: complex buffer
                    remaining_acid = moles_analyte - moles_titrant
                    formed_conjugate = moles_titrant
                    pH = calculate_buffer_pH(remaining_acid / total_volume, 
                                           formed_conjugate / total_volume, 
                                           Ka_analyte)
            
            elif moles_titrant > moles_analyte:  # After equivalence
                if is_titrant_strong:
                    # Excess strong base
                    excess_base = (moles_titrant - moles_analyte) / total_volume
                    pH = calculate_pH_strong_base(excess_base)
                else:
                    # Excess weak base
                    excess_base = (moles_titrant - moles_analyte) / total_volume
                    pH = calculate_pH_weak_base(excess_base, Kb_titrant)
            
            else:  # At equivalence point
                if is_titrant_strong and is_analyte_strong:
                    # Strong base + Strong acid = neutral salt
                    pH = 7.0
                elif not is_titrant_strong and is_analyte_strong:
                    # Weak base + Strong acid = acidic salt
                    salt_conc = moles_analyte / total_volume
                    pH = calculate_pH_weak_acid(salt_conc, Kw / Kb_titrant)
                elif is_titrant_strong and not is_analyte_strong:
                    # Strong base + Weak acid = basic salt
                    salt_conc = moles_analyte / total_volume
                    Kb_salt = Kw / Ka_analyte
                    pH = calculate_pH_weak_base(salt_conc, Kb_salt)
                else:
                    # Weak base + Weak acid = depends on relative Ka/Kb
                    salt_conc = moles_analyte / total_volume
                    if Ka_analyte > Kb_titrant:
                        pH = calculate_pH_weak_acid(salt_conc, Kw / Kb_titrant)
                    else:
                        Kb_salt = Kw / Ka_analyte
                        pH = calculate_pH_weak_base(salt_conc, Kb_salt)
        
        pH_values.append(pH)
    
    # Calculate half-equivalence point pH
    half_equiv_vol = equiv_vol / 2
    if titrant_type == "acid" and not is_analyte_strong:
        half_equiv_pH = -np.log10(Ka_conjugate)
    elif titrant_type == "base" and not is_analyte_strong:
        half_equiv_pH = -np.log10(Ka_analyte)
    else:
        # Interpolate from the data
        half_equiv_idx = np.abs(volume_titrant_added - half_equiv_vol).argmin()
        half_equiv_pH = pH_values[half_equiv_idx]
    
    # Calculate equivalence point pH
    equiv_idx = np.abs(volume_titrant_added - equiv_vol).argmin()
    equiv_pH = pH_values[equiv_idx]
    
    return {
        'volume_titrant_added': volume_titrant_added,
        'pH_values': pH_values,
        'half_equiv_vol': half_equiv_vol,
        'half_equiv_pH': half_equiv_pH,
        'equiv_vol': equiv_vol,
        'equiv_pH': equiv_pH
    }

def save_data_to_csv(results, params):
    """Save titration data to CSV file"""
    csv_filename = f"titration_{params['titrant']}_{params['analyte']}.csv"
    with open(csv_filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Volume of Titrant (mL)", "pH Value"])
        for i in range(len(results['volume_titrant_added'])):
            writer.writerow([results['volume_titrant_added'][i], results['pH_values'][i]])
    
    print(f"Titration data exported to {csv_filename}")
    return csv_filename

def plot_titration_curve(results, params):
    """Plot titration curve with key points"""
    plt.figure(figsize=(10, 6))
    
    # Plot the titration curve
    plt.plot(results['volume_titrant_added'], results['pH_values'], 'b-', linewidth=2)
    
    # Add horizontal line at pH 7
    plt.axhline(y=7, color='gray', linestyle='--', label="Neutral pH")
    
    # Add vertical lines at half-equivalence and equivalence points
    plt.axvline(x=results['half_equiv_vol'], color='green', linestyle='--', 
                label=f"Half-Equivalence Point (pH={results['half_equiv_pH']:.2f})")
    plt.axvline(x=results['equiv_vol'], color='red', linestyle='--', 
                label=f"Equivalence Point (pH={results['equiv_pH']:.2f})")
    
    # Mark buffer region if applicable
    if ((params['titrant_type'] == 'acid' and acids[params['titrant']] < 1e0) or
        (params['titrant_type'] == 'base' and acids[params['analyte']] < 1e0)):
        buffer_start = max(0, results['equiv_vol'] * 0.2)
        buffer_end = min(results['equiv_vol'] * 0.8, results['max_titrant_vol'])
        plt.axvspan(buffer_start, buffer_end, alpha=0.2, color='yellow', label="Buffer Region")
    
    # Add labels and title
    plt.xlabel("Volume of Titrant Added (mL)")
    plt.ylabel("pH")
    title = f"Titration Curve: {params['titrant']} ({params['titrant_conc']} M) → "
    title += f"{params['analyte']} ({params['analyte_conc']} M, {params['analyte_vol']} mL)"
    plt.title(title)
    
    # Add grid and legend
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best')
    
    return plt

def main():
    """Main function to run the titration simulator"""
    # Get user input
    params = get_user_input()
    
    # Simulate titration
    results = simulate_titration(params)
    
    # Save data to CSV
    csv_filename = save_data_to_csv(results, params)
    
    # Print key results
    print("\n--- TITRATION RESULTS ---")
    print(f"Initial pH: {results['pH_values'][0]:.2f}")
    print(f"Half-Equivalence Point: {results['half_equiv_vol']:.2f} mL (pH = {results['half_equiv_pH']:.2f})")
    print(f"Equivalence Point: {results['equiv_vol']:.2f} mL (pH = {results['equiv_pH']:.2f})")
    print(f"Final pH: {results['pH_values'][-1]:.2f}")
    
    # Plot titration curve
    plt = plot_titration_curve(results, params)
    plt.show()

if __name__ == "__main__":
    main()