import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import math

class TitrationSimulator:
    """
    A class to simulate acid-base titrations and generate titration curves.
    """
    
    # Constants
    Kw = 1e-14  # Water dissociation constant
    
    # Define acid/base properties: Ka or Kb values (use None for strong acids/bases)
    ACIDS = {
        "HCl": None,  # Strong acid
        "H2SO4": None,  # Strong acid (first proton)
        "HNO3": None,  # Strong acid
        "CH3COOH": 1.8e-5,  # Acetic acid
        "HCOOH": 1.8e-4,  # Formic acid
        "H3PO4": 7.5e-3,  # Phosphoric acid (first proton)
        "H2CO3": 4.3e-7,  # Carbonic acid
        "HF": 6.8e-4,  # Hydrofluoric acid
        "H3BO3": 5.8e-10,  # Boric acid
    }
    
    BASES = {
        "NaOH": None,  # Strong base
        "KOH": None,  # Strong base
        "Ba(OH)2": None,  # Strong base
        "NH3": 1.8e-5,  # Ammonia
        "C5H5N": 1.7e-9,  # Pyridine
        "CH3NH2": 4.4e-4,  # Methylamine
    }
    
    # Indicators with pH ranges and color changes
    INDICATORS = {
        "Methyl violet": (0.0, 1.6, "yellow", "violet"),
        "Thymol blue (acid)": (1.2, 2.8, "red", "yellow"),
        "Methyl orange": (3.1, 4.4, "red", "orange/yellow"),
        "Methyl red": (4.4, 6.2, "red", "yellow"),
        "Bromothymol blue": (6.0, 7.6, "yellow", "blue"),
        "Phenol red": (6.4, 8.0, "yellow", "red"),
        "Phenolphthalein": (8.3, 10.0, "colorless", "pink"),
        "Thymolphthalein": (9.3, 10.5, "colorless", "blue"),
        "Alizarin yellow R": (10.1, 12.0, "yellow", "red"),
    }
    
    def __init__(self, analyte, titrant, analyte_conc, titrant_conc, analyte_vol):
        """
        Initialize the titration simulation with the analyte and titrant.
        
        Parameters:
        -----------
        analyte : str
            The analyte (acid or base)
        titrant : str
            The titrant (base or acid)
        analyte_conc : float
            Concentration of the analyte in mol/L
        titrant_conc : float
            Concentration of the titrant in mol/L
        analyte_vol : float
            Volume of the analyte in mL
        """
        self.analyte = analyte
        self.titrant = titrant
        self.analyte_conc = analyte_conc
        self.titrant_conc = titrant_conc
        self.analyte_vol = analyte_vol
        
        # Determine if analyte is acid or base
        self.analyte_is_acid = analyte in self.ACIDS
        self.titrant_is_acid = titrant in self.ACIDS
        
        # Get Ka or Kb values
        if self.analyte_is_acid:
            self.analyte_K = self.ACIDS[analyte]
        else:
            self.analyte_K = self.BASES[analyte]
            
        if self.titrant_is_acid:
            self.titrant_K = self.ACIDS[titrant]
        else:
            self.titrant_K = self.BASES[titrant]
            
        # Validate that one is acid and one is base
        if self.analyte_is_acid == self.titrant_is_acid:
            raise ValueError("One reactant must be an acid and the other must be a base.")
    
    def calculate_pH(self, titrant_vol):
        """
        Calculate the pH for a given volume of titrant added.
        
        Parameters:
        -----------
        titrant_vol : float
            Volume of titrant added in mL
            
        Returns:
        --------
        float
            The pH value
        """
        # Calculate moles of analyte and titrant
        moles_analyte = self.analyte_conc * (self.analyte_vol / 1000)  # Convert mL to L
        moles_titrant = self.titrant_conc * (titrant_vol / 1000)  # Convert mL to L
        
        # Total volume in L
        total_vol = (self.analyte_vol + titrant_vol) / 1000
        
        # Different cases based on whether we have acid-base or base-acid titration
        if self.analyte_is_acid and not self.titrant_is_acid:  # Acid analyte, Base titrant
            return self._calculate_pH_acid_analyte(moles_analyte, moles_titrant, total_vol)
        elif not self.analyte_is_acid and self.titrant_is_acid:  # Base analyte, Acid titrant
            return self._calculate_pH_base_analyte(moles_analyte, moles_titrant, total_vol)
    
    def _calculate_pH_acid_analyte(self, moles_acid, moles_base, total_vol):
        """Helper method to calculate pH for acid analyte titrated with base"""
        # Excess of acid or base
        excess = moles_acid - moles_base
        
        # Strong acid - Strong base
        if self.analyte_K is None and self.titrant_K is None:
            if excess > 0:  # Excess acid
                conc_H = excess / total_vol
                return -math.log10(conc_H)
            elif excess < 0:  # Excess base
                conc_OH = abs(excess) / total_vol
                conc_H = self.Kw / conc_OH
                return -math.log10(conc_H)
            else:  # Equivalence point - neutral solution
                return 7.0
                
        # Weak acid - Strong base
        elif self.analyte_K is not None and self.titrant_K is None:
            Ka = self.analyte_K
            
            if excess > 0:  # Before equivalence point - buffer region
                # [A-]/[HA] ratio
                if moles_base > 0:
                    ratio = moles_base / excess
                    return -math.log10(Ka) + math.log10(ratio)
                else:  # No base added yet
                    # Use quadratic formula for weak acid
                    c = moles_acid / total_vol
                    h = (-Ka + math.sqrt(Ka**2 + 4*Ka*c)) / 2
                    return -math.log10(h)
            elif excess < 0:  # After equivalence point - excess strong base
                conc_OH = abs(excess) / total_vol
                conc_H = self.Kw / conc_OH
                return -math.log10(conc_H)
            else:  # Equivalence point - salt of weak acid
                # For the salt of a weak acid, pH is determined by the hydrolysis of the conjugate base
                # pKa + pKb = 14 → pKb = 14 - pKa
                # At equivalence point, pH = 7 + (pKa - pKb)/2 = 7 + (pKa - (14 - pKa))/2 = 7 + (2*pKa - 14)/2 = pKa
                return -math.log10(Ka)
                
        # Strong acid - Weak base
        elif self.analyte_K is None and self.titrant_K is not None:
            Kb = self.titrant_K
            
            if excess > 0:  # Excess strong acid
                conc_H = excess / total_vol
                return -math.log10(conc_H)
            elif excess < 0:  # After equivalence - mixture of weak base and its salt
                # [B]/[BH+] ratio
                ratio = abs(excess) / moles_acid
                pOH = -math.log10(Kb) + math.log10(ratio)
                return 14 - pOH
            else:  # Equivalence point - salt of weak base
                # pH determined by hydrolysis of conjugate acid
                Ka_conj = self.Kw / Kb
                # At equivalence point for salt of weak base, pH = pKa of conjugate acid = 14 - pKb of original base
                return -math.log10(Ka_conj)
        
        return 7.0  # Default fallback
    
    def _calculate_pH_base_analyte(self, moles_base, moles_acid, total_vol):
        """Helper method to calculate pH for base analyte titrated with acid"""
        # Excess of base or acid
        excess = moles_base - moles_acid
        
        # Strong base - Strong acid
        if self.analyte_K is None and self.titrant_K is None:
            if excess > 0:  # Excess base
                conc_OH = excess / total_vol
                conc_H = self.Kw / conc_OH
                return -math.log10(conc_H)
            elif excess < 0:  # Excess acid
                conc_H = abs(excess) / total_vol
                return -math.log10(conc_H)
            else:  # Equivalence point - neutral solution
                return 7.0
                
        # Strong base - Weak acid
        elif self.analyte_K is None and self.titrant_K is not None:
            Ka = self.titrant_K
            
            if excess > 0:  # Excess strong base
                conc_OH = excess / total_vol
                conc_H = self.Kw / conc_OH
                return -math.log10(conc_H)
            elif excess < 0:  # After equivalence - excess weak acid
                # Use Henderson-Hasselbalch for buffer solution
                c_acid = abs(excess) / total_vol
                c_salt = moles_base / total_vol
                pH = -math.log10(Ka) + math.log10(c_salt / c_acid)
                return pH
            else:  # Equivalence point - salt of weak acid
                # At equivalence point, pH = pKa
                return -math.log10(Ka)
                
        # Weak base - Strong acid
        elif self.analyte_K is not None and self.titrant_K is None:
            Kb = self.analyte_K
            Ka_conj = self.Kw / Kb  # Ka of conjugate acid
            
            if excess > 0:  # Before equivalence point - buffer region
                # [B]/[BH+] ratio
                ratio = excess / moles_acid if moles_acid > 0 else float('inf')
                if ratio == float('inf'):  # No acid added yet
                    # Use quadratic formula for weak base
                    c = moles_base / total_vol
                    oh = (-Kb + math.sqrt(Kb**2 + 4*Kb*c)) / 2
                    h = self.Kw / oh
                    return -math.log10(h)
                pOH = -math.log10(Kb) + math.log10(ratio)
                return 14 - pOH
            elif excess < 0:  # After equivalence - excess strong acid
                conc_H = abs(excess) / total_vol
                return -math.log10(conc_H)
            else:  # Equivalence point - salt of weak base
                # At equivalence point, pH = 14 - pKb = pKa of conjugate acid
                return -math.log10(Ka_conj)
        
        # Weak base - Weak acid (if this ever gets implemented)
        # Need to consider relative strengths for mixed weak acid/base reactions
        
        return 7.0  # Default fallback
    
    def generate_titration_curve(self, max_vol=None, points=100):
        """
        Generate data for a titration curve.
        
        Parameters:
        -----------
        max_vol : float, optional
            Maximum volume of titrant to add (mL). If None, calculated automatically.
        points : int, optional
            Number of points to calculate
            
        Returns:
        --------
        tuple
            (volumes, pH values) as numpy arrays
        """
        if max_vol is None:
            # Calculate volume needed for complete reaction plus some excess
            moles_analyte = self.analyte_conc * (self.analyte_vol / 1000)
            vol_equivalence = moles_analyte * 1000 / self.titrant_conc  # in mL
            max_vol = vol_equivalence * 2  # Go to twice the equivalence volume
        
        # Generate volumes with more points around the equivalence point
        vol_equivalence = self.get_equivalence_point()[0]
        
        # Create non-uniform distribution of points with higher density near equivalence point
        # 40% of points before equivalence, 20% around equivalence, 40% after equivalence
        vols_before = np.linspace(0, 0.95*vol_equivalence, int(0.4*points))
        vols_around = np.linspace(0.95*vol_equivalence, 1.05*vol_equivalence, int(0.2*points))
        vols_after = np.linspace(1.05*vol_equivalence, max_vol, int(0.4*points))
        volumes = np.unique(np.concatenate([vols_before, vols_around, vols_after]))
        
        # Calculate pH for each volume
        pH_values = np.array([self.calculate_pH(vol) for vol in volumes])
        
        return volumes, pH_values
    
    def get_equivalence_point(self):
        """
        Calculate the equivalence point of the titration.
        
        Returns:
        --------
        tuple
            (volume at equivalence point (mL), pH at equivalence point)
        """
        # Calculate volume at equivalence point
        moles_analyte = self.analyte_conc * (self.analyte_vol / 1000)
        vol_equivalence = moles_analyte * 1000 / self.titrant_conc  # in mL
        
        # Calculate pH at equivalence point
        pH_equivalence = self.calculate_pH(vol_equivalence)
        
        return vol_equivalence, pH_equivalence
    
    def get_half_equivalence_point(self):
        """
        Calculate the half-equivalence point of the titration.
        
        Returns:
        --------
        tuple
            (volume at half-equivalence point (mL), pH at half-equivalence point)
        """
        vol_equivalence, _ = self.get_equivalence_point()
        vol_half_equivalence = vol_equivalence / 2
        
        # At half-equivalence point for weak acid/base titrations:
        # For weak acid + strong base: pH = pKa
        # For weak base + strong acid: pH = 14 - pKb = pKa of conjugate acid
        pH_half_equivalence = self.calculate_pH(vol_half_equivalence)
        
        return vol_half_equivalence, pH_half_equivalence
    
    def get_buffer_region(self):
        """
        Calculate the buffer region of the titration (relevant for weak acid/base titrations).
        
        Returns:
        --------
        tuple
            (start volume, end volume) of buffer region in mL
        """
        vol_equivalence, _ = self.get_equivalence_point()
        
        # Buffer region is typically considered to be from 10% to 90% neutralization
        buffer_start = vol_equivalence * 0.1
        buffer_end = vol_equivalence * 0.9
        
        return buffer_start, buffer_end
    
    def find_inflection_point(self, volumes, pH_values):
        """
        Find the inflection point in the titration curve, which corresponds to the equivalence point.
        
        Parameters:
        -----------
        volumes : numpy.ndarray
            Array of titrant volumes
        pH_values : numpy.ndarray
            Array of corresponding pH values
            
        Returns:
        --------
        tuple
            (volume at inflection point, pH at inflection point)
        """
        # Apply Savitzky-Golay filter to smooth the pH values
        if len(pH_values) > 10:  # Need enough points for smoothing
            try:
                pH_smooth = savgol_filter(pH_values, min(11, len(pH_values) - (len(pH_values) % 2 - 1)), 3)
            except:
                pH_smooth = pH_values
        else:
            pH_smooth = pH_values
        
        # Calculate first derivative (dpH/dV)
        dpH = np.gradient(pH_smooth)
        dV = np.gradient(volumes)
        first_derivative = dpH / dV
        
        # Find the maximum of the first derivative
        max_idx = np.argmax(first_derivative)
        
        return volumes[max_idx], pH_values[max_idx]
    
    def recommend_indicator(self):
        """
        Recommend an appropriate indicator for the titration based on the pH range at the equivalence point.
        
        Returns:
        --------
        dict
            Information about the recommended indicator
        """
        _, pH_eq = self.get_equivalence_point()
        
        best_indicator = None
        min_distance = float('inf')
        
        for name, (low_pH, high_pH, color1, color2) in self.INDICATORS.items():
            # Calculate the midpoint of the indicator's pH range
            mid_pH = (low_pH + high_pH) / 2
            
            # Calculate how close this midpoint is to our equivalence point pH
            distance = abs(mid_pH - pH_eq)
            
            # If this indicator is better than our current best, update
            if distance < min_distance:
                min_distance = distance
                best_indicator = {
                    "name": name,
                    "pH_range": (low_pH, high_pH),
                    "color_change": f"{color1} to {color2}",
                    "suitability": "Excellent" if low_pH <= pH_eq <= high_pH else "Usable"
                }
        
        # If no suitable indicator found within range, mark the closest as "Usable"
        if best_indicator["suitability"] != "Excellent":
            best_indicator["suitability"] = "Usable, but transition may not be ideal"
            
        return best_indicator
    
    def plot_titration_curve(self, max_vol=None, save_path=None):
        """
        Plot the titration curve.
        
        Parameters:
        -----------
        max_vol : float, optional
            Maximum volume of titrant to add (mL). If None, calculated automatically.
        save_path : str, optional
            Path to save the plot. If None, the plot is displayed instead.
            
        Returns:
        --------
        tuple
            (fig, ax) Matplotlib figure and axis objects
        """
        volumes, pH_values = self.generate_titration_curve(max_vol)
        
        # Calculate key points
        eq_vol, eq_pH = self.get_equivalence_point()
        half_eq_vol, half_eq_pH = self.get_half_equivalence_point()
        
        # Get inflection point (should be close to equivalence point, but more accurate)
        infl_vol, infl_pH = self.find_inflection_point(volumes, pH_values)
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot the titration curve
        ax.plot(volumes, pH_values, 'b-', linewidth=2, label='Titration Curve')
        
        # Mark equivalence point and half-equivalence point
        ax.plot(eq_vol, eq_pH, 'ro', markersize=8, label=f'Equivalence Point (pH = {eq_pH:.2f})')
        ax.plot(half_eq_vol, half_eq_pH, 'go', markersize=8, label=f'Half-Equivalence Point (pH = {half_eq_pH:.2f})')
        
        # Add vertical line at equivalence point
        ax.axvline(x=eq_vol, color='r', linestyle='--', alpha=0.5)
        
        # Add horizontal line at half-equivalence point pH (= pKa for weak acid)
        if (self.analyte_is_acid and self.analyte_K is not None) or \
           (not self.analyte_is_acid and self.titrant_K is not None):
            ax.axhline(y=half_eq_pH, color='g', linestyle='--', alpha=0.5)
        
        # Get buffer region if applicable
        if (self.analyte_is_acid and self.analyte_K is not None) or \
           (not self.analyte_is_acid and self.titrant_K is not None):
            buffer_start, buffer_end = self.get_buffer_region()
            ax.axvspan(buffer_start, buffer_end, alpha=0.2, color='green', label='Buffer Region')
        
        # Get recommended indicator
        indicator = self.recommend_indicator()
        
        # Add indicator transition range
        low_pH, high_pH = indicator["pH_range"]
        ax.axhspan(low_pH, high_pH, alpha=0.2, color='orange', label=f'Indicator: {indicator["name"]}')
        
        # Set labels and title
        ax.set_xlabel('Volume of Titrant (mL)', fontsize=12)
        ax.set_ylabel('pH', fontsize=12)
        
        # Create title based on acid-base combination
        acid_name = self.analyte if self.analyte_is_acid else self.titrant
        base_name = self.titrant if self.analyte_is_acid else self.analyte
        acid_type = "Strong" if (self.analyte_is_acid and self.analyte_K is None) or \
                              (self.titrant_is_acid and self.titrant_K is None) else "Weak"
        base_type = "Strong" if (not self.analyte_is_acid and self.analyte_K is None) or \
                              (not self.titrant_is_acid and self.titrant_K is None) else "Weak"
        
        title = f"{acid_type} Acid ({acid_name}) - {base_type} Base ({base_name}) Titration"
        ax.set_title(title, fontsize=14)
        
        # Add a grid
        ax.grid(True, alpha=0.3)
        
        # Add legend
        ax.legend(loc='best')
        
        # Make y-axis start from 0
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(max(0, ymin), min(14, ymax + 1))
        
        # Set x-axis to start from 0
        ax.set_xlim(0, ax.get_xlim()[1])
        
        # Add text annotations with key information
        info_text = (
            f"Analyte: {self.analyte} ({self.analyte_conc:.3f} M, {self.analyte_vol:.1f} mL)\n"
            f"Titrant: {self.titrant} ({self.titrant_conc:.3f} M)\n"
            f"Equivalence Point: {eq_vol:.2f} mL (pH = {eq_pH:.2f})\n"
            f"Half-Equivalence Point: {half_eq_vol:.2f} mL (pH = {half_eq_pH:.2f})\n"
            f"Recommended Indicator: {indicator['name']} ({indicator['color_change']})"
        )
        
        # Place text box in upper left or bottom right depending on curve shape
        if eq_pH > 7:  # Basic titration, place text in bottom right
            ax.text(0.95, 0.05, info_text, transform=ax.transAxes, fontsize=10,
                    verticalalignment='bottom', horizontalalignment='right',
                    bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.8})
        else:  # Acidic titration, place text in upper left
            ax.text(0.05, 0.95, info_text, transform=ax.transAxes, fontsize=10,
                    verticalalignment='top', horizontalalignment='left',
                    bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.8})
        
        # Adjust figure layout
        plt.tight_layout()
        
        # Save or show the plot
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig, ax


def run_titration_simulation():
    """
    Run an interactive titration simulation with user input.
    """
    print("Acid-Base Titration Simulator")
    print("-----------------------------")
    
    # Get available acids and bases
    acids = list(TitrationSimulator.ACIDS.keys())
    bases = list(TitrationSimulator.BASES.keys())
    
    # Print available acids and bases
    print("\nAvailable acids:")
    for i, acid in enumerate(acids, 1):
        acid_type = "Strong" if TitrationSimulator.ACIDS[acid] is None else "Weak"
        print(f"{i}. {acid} ({acid_type})")
    
    print("\nAvailable bases:")
    for i, base in enumerate(bases, 1):
        base_type = "Strong" if TitrationSimulator.BASES[base] is None else "Weak"
        print(f"{i}. {base} ({base_type})")
    
    # Get user input
    try:
        # Analyte (the solution in the flask)
        print("\nSelect analyte type:")
        print("1. Acid")
        print("2. Base")
        analyte_type = int(input("Enter your choice (1 or 2): "))
        
        if analyte_type == 1:  # Acid analyte
            acid_idx = int(input(f"Select acid (1-{len(acids)}): ")) - 1
            analyte = acids[acid_idx]
            
            base_idx = int(input(f"Select base titrant (1-{len(bases)}): ")) - 1
            titrant = bases[base_idx]
        else:  # Base analyte
            base_idx = int(input(f"Select base (1-{len(bases)}): ")) - 1
            analyte = bases[base_idx]
            
            acid_idx = int(input(f"Select acid titrant (1-{len(acids)}): ")) - 1
            titrant = acids[acid_idx]
        
        # Get concentrations and volume
        analyte_conc = float(input(f"Enter {analyte} concentration (mol/L): "))
        titrant_conc = float(input(f"Enter {titrant} concentration (mol/L): "))
        analyte_vol = float(input(f"Enter {analyte} volume (mL): "))
        
        # Create simulator
        simulator = TitrationSimulator(analyte, titrant, analyte_conc, titrant_conc, analyte_vol)
        
        # Generate titration curve
        fig, ax = simulator.plot_titration_curve()
        
        # Display key information
        eq_vol, eq_pH = simulator.get_equivalence_point()
        half_eq_vol, half_eq_pH = simulator.get_half_equivalence_point()
        indicator = simulator.recommend_indicator()
        
        print("\nResults:")
        print("--------")
        print(f"Equivalence Point: {eq_vol:.2f} mL (pH = {eq_pH:.2f})")
        print(f"Half-Equivalence Point: {half_eq_vol:.2f} mL (pH = {half_eq_pH:.2f})")
        print(f"Recommended Indicator: {indicator['name']} ({indicator['color_change']})")
        print(f"Indicator Suitability: {indicator['suitability']}")
        
        # Show the plot
        plt.show()
        
    except (ValueError, IndexError) as e:
        print(f"Error: {e}")
        print("Please enter valid numerical values.")


if __name__ == "__main__":
    # Example usage
    # You can run this directly for an interactive session or use the sample code below
    
    # Example: Manually create and run a titration simulation
    # Uncomment and modify these examples as needed
    
    # Example 1: Weak acid (acetic acid) titrated with strong base (NaOH)
    print("Example 1: Weak acid (acetic acid) titrated with strong base (NaOH)")
    sim1 = TitrationSimulator("CH3COOH", "NaOH", 0.1, 0.1, 15)
    eq_vol1, eq_pH1 = sim1.get_equivalence_point()
    half_eq_vol1, half_eq_pH1 = sim1.get_half_equivalence_point()
    indicator1 = sim1.recommend_indicator()
    print(f"Equivalence Point: {eq_vol1:.2f} mL (pH = {eq_pH1:.2f})")
    print(f"Half-Equivalence Point: {half_eq_vol1:.2f} mL (pH = {half_eq_pH1:.2f})")
    print(f"Recommended Indicator: {indicator1['name']} ({indicator1['color_change']})\n")
    
    # Example 2: Weak base (ammonia) titrated with strong acid (HCl)
    print("Example 2: Weak base (ammonia) titrated with strong acid (HCl)")
    sim2 = TitrationSimulator("NH3", "HCl", 0.1, 0.1, 15)
    eq_vol2, eq_pH2 = sim2.get_equivalence_point()
    half_eq_vol2, half_eq_pH2 = sim2.get_half_equivalence_point()
    indicator2 = sim2.recommend_indicator()
    print(f"Equivalence Point: {eq_vol2:.2f} mL (pH = {eq_pH2:.2f})")
    print(f"Half-Equivalence Point: {half_eq_vol2:.2f} mL (pH = {half_eq_pH2:.2f})")
    print(f"Recommended Indicator: {indicator2['name']} ({indicator2['color_change']})\n")
    
    # Example 3: Weak acid (phosphoric acid) titrated with weak base (methylamine)
    print("Example 3: Weak acid (phosphoric acid) titrated with weak base (methylamine)")
    sim3 = TitrationSimulator("H3PO4", "CH3NH2", 0.1, 0.1, 15)
    eq_vol3, eq_pH3 = sim3.get_equivalence_point()
    half_eq_vol3, half_eq_pH3 = sim3.get_half_equivalence_point()
    indicator3 = sim3.recommend_indicator()
    print(f"Equivalence Point: {eq_vol3:.2f} mL (pH = {eq_pH3:.2f})")
    print(f"Half-Equivalence Point: {half_eq_vol3:.2f} mL (pH = {half_eq_pH3:.2f})")
    print(f"Recommended Indicator: {indicator3['name']} ({indicator3['color_change']})\n")
    
    # You can generate and display plots for these examples
    # Uncomment to see plots (one at a time to avoid overlapping windows)
    #sim1.plot_titration_curve()
    #plt.show()
    
    sim2.plot_titration_curve()
    plt.show()
    
    #sim3.plot_titration_curve()
    #plt.show()
    
    # For interactive mode with user input:
    # run_titration_simulation()