import numpy as np

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

# Stoichiometric factors
acid_protons = {
    "HCl": 1,
    "H2SO4": 1,  # Considering only first dissociation in typical conditions
    "H3PO4": 1,  # Considering only first dissociation in typical conditions
    "CH3COOH": 1,
    "H2CO3": 1,  # Considering only first dissociation in typical conditions
    "NH4+": 1
}

base_oh = {
    "NaOH": 1,
    "KOH": 1,
    "Ca(OH)2": 2,
    "NH3": 1,
    "HCO3-": 1,
    "CH3COO-": 1
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
    if concentration <= 0 or Ka <= 0:
        return 7.0
    # Solve quadratic equation for [H+]: [H+]^2 + Ka[H+] - Ka*concentration = 0
    # Using the quadratic formula: x = (-b ± √(b² - 4ac)) / 2a
    # Where a=1, b=Ka, c=-Ka*concentration
    h_conc = (-Ka + np.sqrt(Ka**2 + 4*Ka*concentration)) / 2
    return -np.log10(h_conc)

def calculate_pH_weak_base(concentration, Kb):
    """Calculate pH of a weak base solution"""
    if concentration <= 0 or Kb <= 0:
        return 7.0
    # Calculate [OH-] first using the quadratic formula
    oh_conc = (-Kb + np.sqrt(Kb**2 + 4*Kb*concentration)) / 2
    h_conc = Kw / oh_conc
    return -np.log10(h_conc)

def calculate_buffer_pH(acid_conc, conj_base_conc, Ka):
    """Calculate pH of a buffer solution using Henderson-Hasselbalch"""
    if acid_conc <= 0 or conj_base_conc <= 0 or Ka <= 0:
        if acid_conc > 0:
            return calculate_pH_weak_acid(acid_conc, Ka)
        elif conj_base_conc > 0:
            Kb = Kw / Ka
            return calculate_pH_weak_base(conj_base_conc, Kb)
        else:
            return 7.0
            
    pKa = -np.log10(Ka)
    return pKa + np.log10(conj_base_conc / acid_conc)

def calculate_weak_acid_salt_pH(salt_conc, Ka):
    """Calculate pH of salt solution from weak acid (conjugate base)"""
    if salt_conc <= 0 or Ka <= 0:
        return 7.0
    
    # For the conjugate base of a weak acid
    Kb = Kw / Ka
    
    # Calculate [OH-] using the approximation for salt of weak acid
    # OH- comes from hydrolysis: A- + H2O ⇌ HA + OH-
    # [OH-] = √(Kb × salt_conc)
    oh_conc = np.sqrt(Kb * salt_conc)
    h_conc = Kw / oh_conc
    
    return -np.log10(h_conc)

def calculate_weak_base_salt_pH(salt_conc, Kb):
    """Correctly calculate pH of solution formed by conjugate acid of a weak base"""
    if salt_conc <= 0 or Kb <= 0:
        return 7.0
    
    Ka = Kw / Kb
    return calculate_pH_weak_acid(salt_conc, Ka)

def get_conjugate(species):
    """Get the conjugate acid or base of a species if available"""
    return conjugate_pairs.get(species, None)

def is_strong_acid(species):
    """Check if a species is a strong acid"""
    return species in acids and acids[species] > 1e0

def is_strong_base(species):
    """Check if a species is a strong base"""
    return species in bases and bases[species] > 1e0

def simulate_titration(params):
    print("🔸 simulate_titration loaded from:", __file__)
    """Simulate titration and return results"""
    titrant = params['titrant']
    analyte = params['analyte']
    titrant_type = params['titrant_type']
    analyte_type = params['analyte_type']
    titrant_conc = params['titrant_conc']
    analyte_conc = params['analyte_conc']
    analyte_vol = params['analyte_vol']
    max_titrant_vol = params['max_titrant_vol']
    resolution = params.get('resolution', 300)  # Default resolution if not provided
    
    # Get the number of protons/OH from the species
    n_H = acid_protons.get(titrant, 1) if titrant_type == "acid" else base_oh.get(titrant, 1)
    n_OH = base_oh.get(analyte, 1) if analyte_type == "base" else acid_protons.get(analyte, 1)
    
    # Initial moles of the analyte
    if analyte_type == "base":
        moles_analyte = analyte_conc * (analyte_vol / 1000) * n_OH
    else:  # acid
        moles_analyte = analyte_conc * (analyte_vol / 1000) * acid_protons.get(analyte, 1)
    
    # Equivalence volume in mL
    equiv_vol = (moles_analyte / (titrant_conc * n_H)) * 1000
    
    # Create volume array with more points near equivalence for better resolution
    # Generate more points around equivalence point for smoother curve
    vol_before = np.linspace(0, max(0.95 * equiv_vol, 0), int(resolution * 0.4))
    vol_around = np.linspace(0.95 * equiv_vol, min(1.05 * equiv_vol, max_titrant_vol), int(resolution * 0.2))
    vol_after = np.linspace(min(1.05 * equiv_vol, max_titrant_vol), max_titrant_vol, int(resolution * 0.4))
    volume_titrant_added = np.unique(np.concatenate([vol_before, vol_around, vol_after]))
    
    # Get dissociation constants
    Ka_titrant = acids.get(titrant, None)
    Kb_titrant = bases.get(titrant, None)
    Ka_analyte = acids.get(analyte, None)
    Kb_analyte = bases.get(analyte, None)
    
    # Check strength (strong if Ka or Kb > 1)
    is_titrant_strong_acid = is_strong_acid(titrant)
    is_titrant_strong_base = is_strong_base(titrant)
    is_analyte_strong_acid = is_strong_acid(analyte)
    is_analyte_strong_base = is_strong_base(analyte)
    
    # Calculate initial pH
    pH_values = []
    
    if analyte_type == "acid":
        if is_analyte_strong_acid:
            initial_pH = calculate_pH_strong_acid(analyte_conc)
        else:
            initial_pH = calculate_pH_weak_acid(analyte_conc, Ka_analyte)
    else:  # analyte is base
        if is_analyte_strong_base:
            initial_pH = calculate_pH_strong_base(analyte_conc)
        else:
            initial_pH = calculate_pH_weak_base(analyte_conc, Kb_analyte)
    
    # Calculate pH for each volume of titrant
    for Vt in volume_titrant_added:
        total_vol = (analyte_vol + Vt) / 1000  # Total volume in liters
        
        # moles of titrant added
        moles_titrant = titrant_conc * (Vt / 1000) * n_H
        
        # Calculate pH based on titration type
        if titrant_type == "acid" and analyte_type == "base":
            # Acid titrating base
            moles_base_remaining = moles_analyte - moles_titrant
            
            if moles_base_remaining > 0:
                # Before equivalence point
                if is_analyte_strong_base:
                    # Strong base remaining
                    conc_base = moles_base_remaining / total_vol
                    pH = calculate_pH_strong_base(conc_base)
                else:
                    # Weak base remaining and forming a buffer with its conjugate acid
                    moles_conj_acid = moles_titrant  # conjugate acid formed
                    if moles_conj_acid > 0:  # Buffer region
                        conc_base = moles_base_remaining / total_vol
                        conc_conj_acid = moles_conj_acid / total_vol
                        Ka_conj_acid = Kw / Kb_analyte  # Ka of conjugate acid
                        pH = calculate_buffer_pH(conc_conj_acid, conc_base, Ka_conj_acid)
                    else:  # Just weak base
                        conc_base = moles_base_remaining / total_vol
                        pH = calculate_pH_weak_base(conc_base, Kb_analyte)
            
            elif moles_base_remaining < 0:
                # After equivalence point: excess acid
                moles_excess_acid = abs(moles_base_remaining)
                conc_excess_acid = moles_excess_acid / total_vol
                
                if is_titrant_strong_acid:
                    pH = calculate_pH_strong_acid(conc_excess_acid)
                else:
                    pH = calculate_pH_weak_acid(conc_excess_acid, Ka_titrant)
            
            else:
                # Exactly at equivalence point: salt solution
                salt_conc = moles_analyte / total_vol
                
                if is_titrant_strong_acid and is_analyte_strong_base:
                    # Salt of strong acid and strong base: neutral
                    pH = 7.0
                elif is_titrant_strong_acid and not is_analyte_strong_base:
                    # Salt of strong acid and weak base: acidic
                    # Correct calculation for salt of weak base with strong acid
                    pH = calculate_weak_base_salt_pH(salt_conc, Kb_analyte)
                elif not is_titrant_strong_acid and is_analyte_strong_base:
                    # Salt of weak acid and strong base: basic
                    # Correct calculation for salt of weak acid with strong base
                    pH = calculate_weak_acid_salt_pH(salt_conc, Ka_titrant)
                else:
                    # Exactly at equivalence: salt of conjugate acid and conjugate base
                    conj_acid = get_conjugate(analyte)    # e.g., CH3COOH
                    conj_base = get_conjugate(titrant)    # e.g., NH3
                    Ka_eq = acids[conj_acid]
                    Kb_eq = bases[conj_base]
                    pKa = -np.log10(Ka_eq)
                    pKb = -np.log10(Kb_eq)
                    pH = 7 + 0.5 * (pKa - pKb)
        
        elif titrant_type == "base" and analyte_type == "acid":
            # Base titrating acid
            moles_acid_remaining = moles_analyte - moles_titrant
            
            if moles_acid_remaining > 0:
                # Before equivalence point
                if is_analyte_strong_acid:
                    # Strong acid remaining
                    conc_acid = moles_acid_remaining / total_vol
                    pH = calculate_pH_strong_acid(conc_acid)
                else:
                    # Weak acid remaining and forming a buffer with its conjugate base
                    moles_conj_base = moles_titrant  # conjugate base formed
                    if moles_conj_base > 0:  # Buffer region
                        conc_acid = moles_acid_remaining / total_vol
                        conc_conj_base = moles_conj_base / total_vol
                        pH = calculate_buffer_pH(conc_acid, conc_conj_base, Ka_analyte)
                    else:  # Just weak acid
                        conc_acid = moles_acid_remaining / total_vol
                        pH = calculate_pH_weak_acid(conc_acid, Ka_analyte)
            
            elif moles_acid_remaining < 0:
                # After equivalence point: excess base
                moles_excess_base = abs(moles_acid_remaining)
                conc_excess_base = moles_excess_base / total_vol
                
                if is_titrant_strong_base:
                    pH = calculate_pH_strong_base(conc_excess_base)
                else:
                    pH = calculate_pH_weak_base(conc_excess_base, Kb_titrant)
            
            else:
                # Exactly at equivalence point: salt solution
                salt_conc = moles_analyte / total_vol
                
                if is_titrant_strong_base and is_analyte_strong_acid:
                    # Salt of strong base and strong acid: neutral
                    pH = 7.0
                elif is_titrant_strong_base and not is_analyte_strong_acid:
                    # Salt of strong base and weak acid: basic
                    # Correct calculation for salt of weak acid with strong base
                    pH = calculate_weak_acid_salt_pH(salt_conc, Ka_analyte)
                elif not is_titrant_strong_base and is_analyte_strong_acid:
                    # Salt of weak base and strong acid: acidic
                    # Correct calculation for salt of weak base with strong acid
                    pH = calculate_weak_base_salt_pH(salt_conc, Kb_titrant)
                else:
                    # Exactly at equivalence: salt of conjugate acid and conjugate base
                    conj_base = get_conjugate(analyte)    # e.g., CH3COO-
                    conj_acid = get_conjugate(titrant)    # e.g., NH4+
                    Ka_eq = acids[conj_acid]
                    Kb_eq = bases[conj_base]
                    pKa = -np.log10(Ka_eq)
                    pKb = -np.log10(Kb_eq)
                    pH = 7 + 0.5 * (pKa - pKb)
        
        else:
            # Default to neutral (this shouldn't happen with proper validation)
            pH = 7.0
        
        # Guard against unrealistic pH values
        pH = max(0, min(14, pH))  # Keep pH between 0 and 14
        pH_values.append(pH)
    
    # Calculate half-equivalence point (where pH = pKa for weak acid titrations)
    half_equiv_vol = equiv_vol / 2

    # For different scenarios, half-equivalence pH has different meanings
    if titrant_type == "base" and analyte_type == "acid" and not is_analyte_strong_acid:
        # For weak acid titrated with base, pH at half-equivalence = pKa of weak acid
        half_equiv_pH = -np.log10(Ka_analyte)
    elif titrant_type == "acid" and analyte_type == "base" and not is_analyte_strong_base:
        # For weak base titrated with acid, pH at half-equivalence = 14 - pKb of weak base = pKa of conjugate acid
        half_equiv_pH = 14 + np.log10(Kb_analyte)
    else:
        # For other cases, just find the pH at the half equivalence volume
        idx_half = np.abs(volume_titrant_added - half_equiv_vol).argmin()
        half_equiv_pH = pH_values[idx_half]
    
    # Find pH at true equivalence (should already be calculated but we'll ensure accuracy)
    idx_eq = np.abs(volume_titrant_added - equiv_vol).argmin()
    equiv_pH = pH_values[idx_eq]
    
    return {
        'volume_titrant_added': volume_titrant_added,
        'pH_values': pH_values,
        'half_equiv_vol': half_equiv_vol,
        'half_equiv_pH': half_equiv_pH,
        'equiv_vol': equiv_vol,
        'equiv_pH': equiv_pH,
        'initial_pH': pH_values[0],
        'final_pH': pH_values[-1]
    }