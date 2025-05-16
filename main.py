import io
import csv
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import titration_core as tc  # Import the corrected core functionality

# Configure the page
st.set_page_config(
    page_title="Titration Simulator", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better appearance
st.markdown("""
<style>
    .stApp {
        max-width: 1200px;
        margin: 0 auto;
    }
    .info-box {
        background-color: #f0f8ff;
        padding: 10px;
        border-radius: 5px;
        border-left: 3px solid #4B8BBE;
    }
    .results-header {
        padding: 5px 10px;
        border-radius: 5px;
        margin-bottom: 10px;
    }
</style>
""", unsafe_allow_html=True)

# Title and description
st.title("Titration Simulator")

# Create sidebar for input parameters
st.sidebar.header("Titration Parameters")

# Set up tabs for different types of input methods
tab_mode = st.sidebar.radio(
    "Select Input Mode",
    ["Standard Setup", "Advanced Options"]
)

if tab_mode == "Standard Setup":
    # Create two columns for titrant and analyte selection
    col_t, col_a = st.sidebar.columns(2)
    
    with col_t:
        st.subheader("Titrant")
        # Group acid and base options with headings
        titrant_options = ([f"🧪 ACIDS"] + 
                          list(tc.acids.keys()) + 
                          [f"🧫 BASES"] + 
                          list(tc.bases.keys()))
        
        # Default to HCl (index 1, first acid)
        titrant = st.selectbox(
            "Select Titrant",
            options=titrant_options,
            index=1,  # HCl
            key="titrant"
        )
        
        # Disable the headers from being selected
        if titrant in ["🧪 ACIDS", "🧫 BASES"]:
            st.error("Please select an actual chemical species.")
            st.stop()
            
        # Display titrant info
        if titrant in tc.acids:
            st.markdown(f"**Type**: Acid")
            st.markdown(f"**Ka**: {tc.acids[titrant]:.2e}")
            is_strong = "Yes" if tc.is_strong_acid(titrant) else "No"
            st.markdown(f"**Strong Acid**: {is_strong}")
            titrant_type = "acid"
        else:
            st.markdown(f"**Type**: Base")
            st.markdown(f"**Kb**: {tc.bases[titrant]:.2e}")
            is_strong = "Yes" if tc.is_strong_base(titrant) else "No"
            st.markdown(f"**Strong Base**: {is_strong}")
            titrant_type = "base"
            
    with col_a:
        st.subheader("Analyte")
        # Group acid and base options with headings
        analyte_options = ([f"🧪 ACIDS"] + 
                          list(tc.acids.keys()) + 
                          [f"🧫 BASES"] + 
                          list(tc.bases.keys()))
        
        # Default to NaOH for acid titrant, HCl for base titrant
        default_analyte_index = titrant_options.index("NaOH") if titrant in tc.acids else 1
        
        analyte = st.selectbox(
            "Select Analyte",
            options=analyte_options,
            index=default_analyte_index,
            key="analyte"
        )
        
        # Disable the headers from being selected
        if analyte in ["🧪 ACIDS", "🧫 BASES"]:
            st.error("Please select an actual chemical species.")
            st.stop()
            
        # Display analyte info
        if analyte in tc.acids:
            st.markdown(f"**Type**: Acid")
            st.markdown(f"**Ka**: {tc.acids[analyte]:.2e}")
            is_strong = "Yes" if tc.is_strong_acid(analyte) else "No"
            st.markdown(f"**Strong Acid**: {is_strong}")
            analyte_type = "acid"
        else:
            st.markdown(f"**Type**: Base")
            st.markdown(f"**Kb**: {tc.bases[analyte]:.2e}")
            is_strong = "Yes" if tc.is_strong_base(analyte) else "No"
            st.markdown(f"**Strong Base**: {is_strong}")
            analyte_type = "base"
    
    # Validate acid-base combination
    if (titrant_type == "acid" and analyte_type == "acid") or (titrant_type == "base" and analyte_type == "base"):
        st.sidebar.warning("⚠️ For typical titrations, titrant and analyte should be of opposite types (acid-base or base-acid).")
    
    # Concentrations and volumes
    st.sidebar.subheader("Concentrations & Volumes")
    
    conc_t = st.sidebar.number_input(
        f"{titrant} concentration (M)",
        min_value=0.001,
        max_value=10.0,
        value=0.1,
        step=0.01,
        format="%.3f",
        key="conc_t"
    )
    
    conc_a = st.sidebar.number_input(
        f"{analyte} concentration (M)",
        min_value=0.001,
        max_value=10.0,
        value=0.1,
        step=0.01,
        format="%.3f",
        key="conc_a"
    )
    
    vol_a = st.sidebar.number_input(
        f"{analyte} volume (mL)",
        min_value=1.0,
        max_value=1000.0,
        value=25.0,
        step=5.0,
        key="vol_a"
    )
    
    # Calculate theoretical equivalence point
    if titrant in tc.acids:
        n_H = tc.acid_protons.get(titrant, 1)
    else:
        n_H = tc.base_oh.get(titrant, 1)
        
    if analyte in tc.bases:
        n_OH = tc.base_oh.get(analyte, 1)
    else:
        n_OH = tc.acid_protons.get(analyte, 1)
        
    theoretical_equiv = vol_a * conc_a * n_OH / (conc_t * n_H)
    
    # Set max volume with sensible default based on equivalence point
    max_vt = st.sidebar.number_input(
        f"Max {titrant} volume (mL)",
        min_value=theoretical_equiv * 0.5,
        max_value=theoretical_equiv * 5,
        value=theoretical_equiv * 2,
        step=5.0,
        help=f"Theoretical equivalence point is around {theoretical_equiv:.1f} mL",
        key="max_vt"
    )

else:  # Advanced Options
    # Direct input of all parameters
    st.sidebar.subheader("Titrant")
    titrant = st.sidebar.selectbox(
        "Select Titrant",
        options=list(tc.acids.keys()) + list(tc.bases.keys()),
        index=0,  # HCl
        key="titrant_adv"
    )
    
    if titrant in tc.acids:
        titrant_type = "acid"
    else:
        titrant_type = "base"
        
    st.sidebar.subheader("Analyte")
    analyte = st.sidebar.selectbox(
        "Select Analyte",
        options=list(tc.acids.keys()) + list(tc.bases.keys()),
        index=len(tc.acids),  # NaOH (first base)
        key="analyte_adv"
    )
    
    if analyte in tc.acids:
        analyte_type = "acid"
    else:
        analyte_type = "base"
        
    # Show warning for same-type titration
    if (titrant_type == "acid" and analyte_type == "acid") or (titrant_type == "base" and analyte_type == "base"):
        st.sidebar.warning("⚠️ For typical titrations, titrant and analyte should be of opposite types (acid-base or base-acid).")

    # Advanced concentration and volume inputs
    st.sidebar.subheader("Concentrations & Volumes")
    
    conc_t = st.sidebar.number_input(
        f"{titrant} concentration (M)",
        min_value=0.0001,
        max_value=100.0,
        value=0.1,
        step=0.01,
        format="%.4f",
        key="conc_t_adv"
    )
    
    conc_a = st.sidebar.number_input(
        f"{analyte} concentration (M)",
        min_value=0.0001,
        max_value=100.0,
        value=0.1,
        step=0.01,
        format="%.4f",
        key="conc_a_adv"
    )
    
    vol_a = st.sidebar.number_input(
        f"{analyte} volume (mL)",
        min_value=0.1,
        max_value=10000.0,
        value=25.0,
        step=1.0,
        key="vol_a_adv"
    )
    
    max_vt = st.sidebar.number_input(
        f"Max {titrant} volume (mL)",
        min_value=0.1,
        max_value=10000.0,
        value=50.0,
        step=1.0,
        key="max_vt_adv"
    )

# Add slider for plot resolution
resolution = st.sidebar.slider(
    "Plot Resolution",
    min_value=100,
    max_value=1000,
    value=300,
    step=100,
    help="Higher values give smoother curves but may slow down the app"
)

# Simulate button
simulate = st.sidebar.button("🔬 Simulate Titration", key="simulate_btn", use_container_width=True)

# When simulate button is clicked
if simulate:
    # Build params dict based on selection
    if titrant in tc.acids and analyte in tc.bases:
        titrant_type = "acid"
        analyte_type = "base"
    elif titrant in tc.bases and analyte in tc.acids:
        titrant_type = "base"
        analyte_type = "acid"
    else:
        if titrant in ["🧪 ACIDS", "🧫 BASES"] or analyte in ["🧪 ACIDS", "🧫 BASES"]:
            st.error("Please select valid chemical species for titrant and analyte.")
        else:
            st.error(f"Invalid titrant/analyte combination: {titrant}/{analyte}")
        st.stop()

    # Display titration info
    st.markdown(f"""
    <div class='results-header'>
    <h3>Titration: {titrant} ({conc_t} M) → {analyte} ({conc_a} M, {vol_a} mL)</h3>
    </div>
    """, unsafe_allow_html=True)

    # Create loading spinner for simulation
    with st.spinner("Calculating titration curve..."):
        params = {
            'titrant': titrant,
            'analyte': analyte,
            'titrant_type': titrant_type,
            'analyte_type': analyte_type,
            'titrant_conc': conc_t,
            'analyte_conc': conc_a,
            'analyte_vol': vol_a,
            'max_titrant_vol': max_vt,
            'resolution': resolution
        }
        
        # Run the simulation
        results = tc.simulate_titration(params)
    
    # Create two columns layout for displaying results
    col1, col2 = st.columns([3, 2])
    
    with col1:
        # Plot the curve with matplotlib
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot the titration curve
        ax.plot(results['volume_titrant_added'], results['pH_values'], 'b-', linewidth=2.5)
        
        # Add horizontal line at pH 7
        ax.axhline(y=7, color='gray', linestyle='--', alpha=0.7, label="Neutral pH")
        
        # Add vertical lines at half-equivalence and equivalence points
        ax.axvline(x=results['half_equiv_vol'], color='green', linestyle='--', alpha=0.7,
                   label=f"Half-Equiv Point (pH={results['half_equiv_pH']:.2f})")
        ax.axvline(x=results['equiv_vol'], color='red', linestyle='--', alpha=0.7,
                   label=f"Equiv Point (pH={results['equiv_pH']:.2f})")
        
        # Mark buffer region if applicable for weak acid/base titrations
        if ((titrant_type == 'acid' and not tc.is_strong_acid(titrant)) or
            (analyte_type == 'acid' and not tc.is_strong_acid(analyte))):
            buffer_start = max(0, results['equiv_vol'] * 0.25)
            # use the sidebar’s max_vt value instead of results[...]
            buffer_end = min(results['equiv_vol'] * 0.75, max_vt)
            ax.axvspan(buffer_start, buffer_end, alpha=0.2, color='yellow', label="Buffer Region")
        
        # Improve the appearance
        ax.set_xlabel("Volume of Titrant Added (mL)", fontsize=12)
        ax.set_ylabel("pH", fontsize=12)
        
        # Set appropriate y-axis limits
        min_pH = max(0, min(results['pH_values']) - 0.5)
        max_pH = min(14, max(results['pH_values']) + 0.5)
        ax.set_ylim(min_pH, max_pH)
        
        # Add grid and legend
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='best', fontsize=10)
        
        # Show inflection points more clearly
        equiv_idx = np.abs(results['volume_titrant_added'] - results['equiv_vol']).argmin()
        half_equiv_idx = np.abs(results['volume_titrant_added'] - results['half_equiv_vol']).argmin()
        
        ax.plot(results['equiv_vol'], results['pH_values'][equiv_idx], 'ro', markersize=8)
        ax.plot(results['half_equiv_vol'], results['pH_values'][half_equiv_idx], 'go', markersize=8)
        
        # Add title with chemical equation
        if titrant_type == "acid" and analyte_type == "base":
            title = f"{titrant} + {analyte} → Salt + H₂O"
        else:
            title = f"{analyte} + {titrant} → Salt + H₂O"
        ax.set_title(title, fontsize=14)
        
        # Display the plot
        st.pyplot(fig)
        
        # Create CSV data for download
        csv_data = io.StringIO()
        writer = csv.writer(csv_data)
        writer.writerow(["Volume (mL)", "pH"])
        for v, p in zip(results['volume_titrant_added'], results['pH_values']):
            writer.writerow([f"{v:.3f}", f"{p:.3f}"])
        
        st.download_button(
            label="📥 Download Data as CSV",
            data=csv_data.getvalue(),
            file_name=f"titration_{titrant}_{analyte}.csv",
            mime="text/csv"
        )
    
    with col2:
        # Get pH values from the actual data points
        half_equiv_idx = np.abs(results['volume_titrant_added'] - results['half_equiv_vol']).argmin()
        equiv_idx = np.abs(results['volume_titrant_added'] - results['equiv_vol']).argmin()
        
        # For strong acid-strong base titrations, the equivalence pH should be 7.0
        # Adjust the displayed value if needed
        if ((titrant_type == "acid" and analyte_type == "base" and 
            tc.is_strong_acid(titrant) and tc.is_strong_base(analyte)) or
            (titrant_type == "base" and analyte_type == "acid" and 
            tc.is_strong_base(titrant) and tc.is_strong_acid(analyte))):
            # Strong acid-strong base titration should have pH 7 at equivalence
            equiv_display_pH = 7.0
        else:
            # Use the pH from the simulation for all other cases
            equiv_display_pH = results['pH_values'][equiv_idx]
        
        # Format results nicely with the correct pH values
        results_data = {
            "Parameter": ["Initial pH", "Half-Equivalence Point", 
                        "Equivalence Point", "Final pH"],
            "Value": [
                f"{results['pH_values'][0]:.2f}",
                f"{results['half_equiv_vol']:.2f} mL (pH {results['pH_values'][half_equiv_idx]:.2f})",
                f"{results['equiv_vol']:.2f} mL (pH {equiv_display_pH:.2f})",
                f"{results['pH_values'][-1]:.2f}"
            ]
        }

        # Display results
        st.table(results_data)

    # Create a new row of columns for Titration Analysis and pH Indicator Guide
    analysis_col, indicator_col = st.columns(2)

    with analysis_col:
        # Add information about the titration
        st.subheader("Titration Analysis")
        
        # Generate appropriate analysis based on the titration type
        if titrant_type == "acid" and analyte_type == "base":
            if tc.is_strong_acid(titrant) and tc.is_strong_base(analyte):
                st.info("**Strong Acid - Strong Base Titration**\n\n"
                    "The equivalence point is at pH 7 (neutral). "
                    "This titration has a sharp endpoint with a steep curve near the equivalence point.")
            elif tc.is_strong_acid(titrant) and not tc.is_strong_base(analyte):
                st.info("**Strong Acid - Weak Base Titration**\n\n"
                    "The equivalence point is acidic (pH < 7). "
                    "The curve is less steep near the equivalence point compared to "
                    "strong acid-strong base titrations.")
            elif not tc.is_strong_acid(titrant) and tc.is_strong_base(analyte):
                st.info("**Weak Acid - Strong Base Titration**\n\n"
                    "The equivalence point is basic (pH > 7). "
                    "There's a buffer region before the equivalence point "
                    "where pH changes slowly with added titrant.")
            else:
                st.info("**Weak Acid - Weak Base Titration**\n\n"
                    "The equivalence point depends on the relative strengths of the acid and base. "
                    "The titration curve is less steep with a less distinct endpoint.")
        else:  # base titrating acid
            if tc.is_strong_base(titrant) and tc.is_strong_acid(analyte):
                st.info("**Strong Base - Strong Acid Titration**\n\n"
                    "The equivalence point is at pH 7 (neutral). "
                    "This titration has a sharp endpoint with a steep curve near the equivalence point.")
            elif tc.is_strong_base(titrant) and not tc.is_strong_acid(analyte):
                st.info("**Strong Base - Weak Acid Titration**\n\n"
                    "The equivalence point is basic (pH > 7). "
                    "The curve is less steep near the equivalence point compared to "
                    "strong base-strong acid titrations.")
            elif not tc.is_strong_base(titrant) and tc.is_strong_acid(analyte):
                st.info("**Weak Base - Strong Acid Titration**\n\n"
                    "The equivalence point is acidic (pH < 7). "
                    "There's a buffer region before the equivalence point "
                    "where pH changes slowly with added titrant.")
            else:
                st.info("**Weak Base - Weak Acid Titration**\n\n"
                    "The equivalence point depends on the relative strengths of the base and acid. "
                    "The titration curve is less steep with a less distinct endpoint.")
        # Suggest appropriate indicator
        equiv_pH = results['equiv_pH']
        if equiv_pH < 3.0:
            st.success("**Suggested indicator**: Methyl violet (pH range: 0-2)")
        elif equiv_pH < 5.0:
            st.success("**Suggested indicator**: Methyl orange (pH range: 3.1-4.4)")
        elif equiv_pH < 7.0:
            st.success("**Suggested indicator**: Bromothymol blue (pH range: 6.0-7.6)")
        elif equiv_pH < 9.0:
            st.success("**Suggested indicator**: Phenolphthalein (pH range: 8.0-10.0)")
        else:
            st.success("**Suggested indicator**: Alizarin yellow (pH range: 10.0-12.0)")
            
    with indicator_col:
        # Add indicator color for pH
        st.subheader("pH Indicator Guide")
        
        # Use expander to save space
        with st.expander("Show pH Indicator Colors", expanded=True):
            # Display a gradient of pH indicators
            indicators = {
                "Red (pH < 3)": "Methyl violet, Congo red",
                "Orange (pH 3-5)": "Methyl orange, Bromocresol green",
                "Yellow (pH 5-6)": "Methyl red, Bromothymol blue",
                "Green (pH 6-8)": "Bromothymol blue, Phenol red",
                "Blue (pH 8-10)": "Phenolphthalein, Thymolphthalein",
                "Purple (pH > 10)": "Alizarin yellow, Indigo carmine"
            }
            
            for color, indicator in indicators.items():
                st.markdown(f"**{color}**: {indicator}")
        
        