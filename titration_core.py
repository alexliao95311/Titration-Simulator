import numpy as np

# Constants
Kw = 1.0e-14  # Water dissociation constant

# Acid & Base Dissociation Constants (Ka/Kb)
acids = {
    "HCl": 1e6,         # Strong acid
    "H2SO4": 1e3,       # Strong acid (first dissociation)
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

# Helper functions for identifying species
def is_strong_acid(species):
    """Check if a species is a strong acid"""
    return species in acids and acids[species] > 1e0

def is_strong_base(species):
    """Check if a species is a strong base"""
    return species in bases and bases[species] > 1e0

def get_conjugate(species):
    """Get the conjugate acid or base of a species if available"""
    return conjugate_pairs.get(species, None)

# Basic pH calculations
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

def simulate_titration(params):
    """Simulate titration and return results"""
    # Extract parameters
    titrant = params['titrant']
    analyte = params['analyte']
    # auto-detect acid/base type
    if titrant in acids and titrant not in bases:
        titrant_type = "acid"
    elif titrant in bases and titrant not in acids:
        titrant_type = "base"
    else:
        titrant_type = params.get('titrant_type')

    if analyte in acids and analyte not in bases:
        analyte_type = "acid"
    elif analyte in bases and analyte not in acids:
        analyte_type = "base"
    else:
        analyte_type = params.get('analyte_type')
    print(f"[DEBUG] simulate_titration: titrant={titrant} ({titrant_type}), analyte={analyte} ({analyte_type})")
    titrant_conc = params['titrant_conc']
    analyte_conc = params['analyte_conc']
    analyte_vol = params['analyte_vol']
    max_titrant_vol = params['max_titrant_vol']
    resolution = params.get('resolution', 300)  # Default resolution if not provided
    
    # Get the number of protons/OH from the species
    if titrant_type == "acid":
        n_H = acid_protons.get(titrant, 1)
    else:
        n_H = base_oh.get(titrant, 1)
        
    if analyte_type == "base":
        n_OH = base_oh.get(analyte, 1)
    else:
        n_OH = acid_protons.get(analyte, 1)
    
    # Initial moles of the analyte - Fixed to use correct stoichiometric factor
    if analyte_type == "base":
        moles_analyte = analyte_conc * (analyte_vol / 1000) * n_OH
    else:  # acid
        moles_analyte = analyte_conc * (analyte_vol / 1000) * n_OH  # FIXED: should be acid_protons
    
    # Equivalence volume in mL
    equiv_vol = (moles_analyte / (titrant_conc * n_H)) * 1000
    
    # Create volume array with more points near equivalence for better resolution
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
    
    # Handle weak–weak titrations via analytic piecewise formulas
    weak_weak = False
    if titrant_type=="acid" and analyte_type=="base" and not (is_titrant_strong_acid or is_analyte_strong_base):
        weak_weak = True
        Ka1 = acids[get_conjugate(analyte)]        # HA from analyte
        Ka2 = acids[titrant]                      # HA from titrant
    elif titrant_type=="base" and analyte_type=="acid" and not (is_titrant_strong_base or is_analyte_strong_acid):
        weak_weak = True
        # For weak-weak: get Ka of titrant's conjugate acid and Kb of analyte's conjugate base
        Ka_titrant_conj = acids[get_conjugate(titrant)]
        Kb_analyte_conj = bases[get_conjugate(analyte)]
    if weak_weak:
        print(f"[DEBUG] weak-weak logic triggered for titrant={titrant}/{titrant_type}, analyte={analyte}/{analyte_type}")
        pH_values = []
        half_equiv_pH = None
        # precompute eq pH
        if titrant_type == "acid":
            pKb = -np.log10(bases[get_conjugate(titrant)])
            pKa = -np.log10(acids[get_conjugate(analyte)])
            eq_pH = 7 + 0.5*(pKb - pKa)
            print(f"[DEBUG] weak-weak eq_pH computed: {eq_pH}")
        else:
            # weak-weak: base titrant and acid analyte
            Ka_titrant_conj = acids[get_conjugate(titrant)]
            Kb_analyte_conj = bases[get_conjugate(analyte)]
            eq_pH = 0.5 * (
                -np.log10(Ka_titrant_conj)
                + (14 + np.log10(Kb_analyte_conj))
            )
            print(f"[DEBUG] weak-weak eq_pH computed: {eq_pH}")
        for Vt_ml in volume_titrant_added:
            Vt_L = Vt_ml/1000
            Vtot = (analyte_vol + Vt_ml)/1000
            m_t = titrant_conc * Vt_L * (acid_protons.get(titrant,1) if titrant_type=="acid" else base_oh.get(titrant,1))
            if Vt_ml == 0:
                # initial pH
                if analyte_type == "base":
                    pH = calculate_pH_weak_base(analyte_conc, bases[analyte])
                else:
                    pH = calculate_pH_weak_acid(analyte_conc, acids[analyte])
            elif Vt_ml < equiv_vol:
                # buffer region
                if titrant_type == "acid":
                    m_B = moles_analyte - m_t
                    m_BH = m_t
                    pH = calculate_buffer_pH(m_BH/Vtot, m_B/Vtot, Ka1)
                else:
                    m_A = moles_analyte - m_t
                    m_A_conj = m_t
                    pH = calculate_buffer_pH(m_A/Vtot, m_A_conj/Vtot, acids[analyte])
            elif abs(Vt_ml - equiv_vol) < 1e-6:
                pH = eq_pH
            else:
                if titrant_type == "acid":
                    # post-equivalence: excess titrant (weak acid)
                    excess_moles = m_t - moles_analyte
                    conc_excess = excess_moles / Vtot
                    pH = calculate_pH_weak_acid(conc_excess, Ka2)
                else:
                    # post-equivalence: excess titrant (weak base)
                    excess_moles = m_t - moles_analyte
                    conc_excess = excess_moles / Vtot
                    pH = calculate_pH_weak_base(conc_excess, bases[titrant])
            pH_values.append(max(0.0, min(14.0, pH)))
        # recompute half-equivalence pH
        half_equiv_pH = -np.log10(Ka1) if titrant_type=="acid" else -np.log10(acids[analyte])
        # Override pH at the equivalence point to our manual calculation for weak-weak titrations
        idx_eq = np.abs(volume_titrant_added - equiv_vol).argmin()
        pH_values[idx_eq] = eq_pH
        return {
            'volume_titrant_added': volume_titrant_added,
            'pH_values': pH_values,
            'half_equiv_vol': equiv_vol/2,
            'half_equiv_pH': half_equiv_pH,
            'equiv_vol': equiv_vol,
            'equiv_pH': eq_pH,
            'initial_pH': pH_values[0],
            'final_pH': pH_values[-1]
        }
    # --- End weak-weak special case ---
    # Fallback: original code for all other titrations
    # Calculate initial pH
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
    # Make sure initial pH is within a realistic range
    initial_pH = max(0, min(14, initial_pH))
    pH_values = [initial_pH]
    # Calculate pH for each volume of titrant
    for Vt in volume_titrant_added[1:]:  # Skip first point (already calculated)
        tol = 1e-12  # numerical tolerance for equivalence checks
        total_vol = (analyte_vol + Vt) / 1000  # Total volume in liters
        # moles of titrant added
        moles_titrant = titrant_conc * (Vt / 1000) * n_H
        # Calculate pH based on titration type
        if titrant_type == "acid" and analyte_type == "base":
            # Acid titrating base
            moles_base_remaining = moles_analyte - moles_titrant
            # Determine titration region
            if abs(moles_base_remaining) <= tol:
                # Equivalence point
                if (is_titrant_strong_acid and is_analyte_strong_base):
                    # Strong acid - strong base makes a neutral salt
                    pH = 7.0
                else:
                    # Weak acid or weak base creates a salt solution
                    # Get conjugate species
                    conj_acid_analyte = get_conjugate(analyte)
                    conj_base_titrant = get_conjugate(titrant)
                    if conj_acid_analyte in acids and conj_base_titrant in bases:
                        # Calculate salt pH from conjugate species
                        Ka_conj_acid = acids[conj_acid_analyte]
                        Kb_conj_base = bases[conj_base_titrant]
                        pH = 7 + 0.5 * (np.log10(Kb_conj_base) - np.log10(Ka_conj_acid))
                    else:
                        # Default to neutral if conjugates not found
                        pH = 7.0
            elif moles_base_remaining > tol:
                # Before equivalence: excess base or buffer
                if is_analyte_strong_base:
                    # Strong base - direct calculation
                    conc_base = moles_base_remaining / total_vol
                    pH = calculate_pH_strong_base(conc_base)
                else:
                    # Weak base - buffer region or excess base
                    moles_conj_acid = moles_titrant
                    if moles_conj_acid > 0 and moles_base_remaining > 0:
                        # Buffer region - Henderson-Hasselbalch
                        conc_base = moles_base_remaining / total_vol
                        conc_conj_acid = moles_conj_acid / total_vol
                        Ka_conj_acid = Kw / Kb_analyte
                        pH = calculate_buffer_pH(conc_conj_acid, conc_base, Ka_conj_acid)
                    else:
                        # Just excess base
                        conc_base = moles_base_remaining / total_vol
                        pH = calculate_pH_weak_base(conc_base, Kb_analyte)
            else:
                # After equivalence: excess acid
                moles_excess_acid = abs(moles_base_remaining)
                conc_excess_acid = moles_excess_acid / total_vol
                if is_titrant_strong_acid:
                    pH = calculate_pH_strong_acid(conc_excess_acid)
                else:
                    pH = calculate_pH_weak_acid(conc_excess_acid, Ka_titrant)
        elif titrant_type == "base" and analyte_type == "acid":
            # Base titrating acid
            moles_acid_remaining = moles_analyte - moles_titrant
            # Determine titration region
            if abs(moles_acid_remaining) <= tol:
                # Equivalence point
                if (is_titrant_strong_base and is_analyte_strong_acid):
                    # Strong base - strong acid makes a neutral salt
                    pH = 7.0
                else:
                    # Weak acid or weak base creates a salt solution
                    # Get conjugate species
                    conj_acid_titrant = get_conjugate(titrant)
                    conj_base_analyte = get_conjugate(analyte)
                    if conj_acid_titrant in acids and conj_base_analyte in bases:
                        # Calculate salt pH from conjugate species
                        Ka_conj_acid = acids[conj_acid_titrant]
                        Kb_conj_base = bases[conj_base_analyte]
                        pH = 7 + 0.5 * (np.log10(Kb_conj_base) - np.log10(Ka_conj_acid))
                    else:
                        # Default to neutral if conjugates not found
                        pH = 7.0
            elif moles_acid_remaining > tol:
                # Before equivalence: excess acid or buffer
                if is_analyte_strong_acid:
                    # Strong acid - direct calculation
                    conc_acid = moles_acid_remaining / total_vol
                    pH = calculate_pH_strong_acid(conc_acid)
                else:
                    # Weak acid - buffer region or excess acid
                    moles_conj_base = moles_titrant
                    if moles_conj_base > 0 and moles_acid_remaining > 0:
                        # Buffer region - Henderson-Hasselbalch
                        conc_acid = moles_acid_remaining / total_vol
                        conc_conj_base = moles_conj_base / total_vol
                        pH = calculate_buffer_pH(conc_acid, conc_conj_base, Ka_analyte)
                    else:
                        # Just excess acid
                        conc_acid = moles_acid_remaining / total_vol
                        pH = calculate_pH_weak_acid(conc_acid, Ka_analyte)
            else:
                # After equivalence: excess base
                moles_excess_base = abs(moles_acid_remaining)
                conc_excess_base = moles_excess_base / total_vol
                if is_titrant_strong_base:
                    pH = calculate_pH_strong_base(conc_excess_base)
                else:
                    pH = calculate_pH_weak_base(conc_excess_base, Kb_titrant)
        else:
            # Default to neutral (this shouldn't happen with proper validation)
            pH = 7.0
        # Guard against unrealistic pH values
        pH = max(0, min(14, pH))
        pH_values.append(pH)
    # Calculate half-equivalence point (where pH = pKa for weak acid titrations)
    half_equiv_vol = equiv_vol / 2
    # Handle half-equivalence point pH calculations properly
    if titrant_type == "base" and analyte_type == "acid" and not is_analyte_strong_acid:
        # For weak acid titrated with base, pH at half-equivalence = pKa of weak acid
        half_equiv_pH = -np.log10(Ka_analyte)
    elif titrant_type == "acid" and analyte_type == "base" and not is_analyte_strong_base:
        # For weak base titrated with acid, pH at half-equivalence = pKw - pKb of weak base
        half_equiv_pH = 14 - (-np.log10(Kb_analyte))
    else:
        # For other cases, find the pH at half equivalence volume from the data
        idx_half = np.abs(volume_titrant_added - half_equiv_vol).argmin()
        half_equiv_pH = pH_values[idx_half]
    # Determine the index closest to the equivalence volume
    idx_eq = np.abs(volume_titrant_added - equiv_vol).argmin()
    # Manual pH at equivalence for weak–weak titrations
    if titrant_type == "acid" and analyte_type == "base" and not (is_titrant_strong_acid or is_analyte_strong_base):
        conj_acid = get_conjugate(analyte)
        conj_base = get_conjugate(titrant)
        if conj_acid in acids and conj_base in bases:
            Ka_c = acids[conj_acid]
            Kb_c = bases[conj_base]
            equiv_pH = 7 + 0.5 * (np.log10(Kb_c) - np.log10(Ka_c))
        else:
            equiv_pH = pH_values[idx_eq]
    elif titrant_type == "base" and analyte_type == "acid" and not (is_titrant_strong_base or is_analyte_strong_acid):
        conj_acid = get_conjugate(titrant)
        conj_base = get_conjugate(analyte)
        if conj_acid in acids and conj_base in bases:
            Ka_c = acids[conj_acid]
            Kb_c = bases[conj_base]
            equiv_pH = 7 + 0.5 * (np.log10(Kb_c) - np.log10(Ka_c))
        else:
            equiv_pH = pH_values[idx_eq]
    else:
        # Strong–strong or mixed cases
        if (titrant_type == "acid" and analyte_type == "base" and is_titrant_strong_acid and is_analyte_strong_base) \
        or (titrant_type == "base" and analyte_type == "acid" and is_titrant_strong_base and is_analyte_strong_acid):
            equiv_pH = 7.0
        else:
            equiv_pH = pH_values[idx_eq]
    # Override the curve value at equivalence to our manual calculation
    pH_values[idx_eq] = equiv_pH
    return {
        'volume_titrant_added': volume_titrant_added,
        'pH_values': pH_values,
        'half_equiv_vol': half_equiv_vol,
        'half_equiv_pH': half_equiv_pH,
        'equiv_vol': equiv_vol,
        'equiv_pH': equiv_pH,
        'initial_pH': initial_pH,
        'final_pH': pH_values[-1]
    }
    