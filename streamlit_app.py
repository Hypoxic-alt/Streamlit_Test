import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import random

# -------------------------------
# Titration Calculation and Plotting Functions

def titration_pH(V_titrant, titration_type, p_value, V_initial=0.05, C_initial=0.1, C_titrant=0.1):
    moles_initial = C_initial * V_initial
    pH = np.zeros_like(V_titrant)
    for i, Vt in enumerate(V_titrant):
        moles_titrant = C_titrant * Vt
        V_total = V_initial + Vt
        if titration_type == 'Weak Acid with Strong Base':
            pKa = p_value
            Ka = 10**(-pKa)
            if moles_titrant < moles_initial:
                if moles_titrant == 0:
                    C_eff = moles_initial / V_total
                    H_conc = np.sqrt(Ka * C_eff)
                    pH[i] = -np.log10(H_conc)
                else:
                    ratio = moles_titrant / (moles_initial - moles_titrant)
                    pH[i] = pKa + np.log10(ratio)
            elif np.isclose(moles_titrant, moles_initial, rtol=1e-6):
                C_conjugate = moles_initial / V_total
                Kb = 1e-14 / Ka
                OH_conc = np.sqrt(Kb * C_conjugate)
                pOH = -np.log10(OH_conc)
                pH[i] = 14 - pOH
            else:
                excess_OH = (moles_titrant - moles_initial) / V_total
                pOH = -np.log10(excess_OH)
                pH[i] = 14 - pOH
        elif titration_type == 'Weak Base with Strong Acid':
            pKa = p_value  # pKₐ of the conjugate acid
            Ka = 10**(-pKa)
            Kb = 1e-14 / Ka
            if moles_titrant < moles_initial:
                if moles_titrant == 0:
                    C_eff = moles_initial / V_total
                    OH_conc = np.sqrt(Kb * C_eff)
                    pOH = -np.log10(OH_conc)
                    pH[i] = 14 - pOH
                else:
                    moles_base = moles_initial - moles_titrant
                    moles_conjugate = moles_titrant
                    ratio = moles_base / moles_conjugate
                    pH[i] = pKa + np.log10(ratio)
            elif np.isclose(moles_titrant, moles_initial, rtol=1e-6):
                C_conjugate = moles_initial / V_total
                H_conc = np.sqrt(Ka * C_conjugate)
                pH[i] = -np.log10(H_conc)
            else:
                excess_H = (moles_titrant - moles_initial) / V_total
                pH[i] = -np.log10(excess_H)
        elif titration_type == 'Strong Acid with Strong Base':
            if moles_titrant < moles_initial:
                H_conc = (moles_initial - moles_titrant) / V_total
                pH[i] = -np.log10(H_conc)
            elif np.isclose(moles_titrant, moles_initial, rtol=1e-6):
                pH[i] = 7.0
            else:
                OH_conc = (moles_titrant - moles_initial) / V_total
                pOH = -np.log10(OH_conc)
                pH[i] = 14 - pOH
        elif titration_type == 'Strong Base with Strong Acid':
            if moles_titrant < moles_initial:
                OH_conc = (moles_initial - moles_titrant) / V_total
                pOH = -np.log10(OH_conc)
                pH[i] = 14 - pOH
            elif np.isclose(moles_titrant, moles_initial, rtol=1e-6):
                pH[i] = 7.0
            else:
                H_conc = (moles_titrant - moles_initial) / V_total
                pH[i] = -np.log10(H_conc)
    return pH

def plot_titration(titration_type, p_value, V_initial, show_eq_line, show_eq_ph_line, show_half_eq_line, show_description):
    V_titrant = np.linspace(0, 0.1, 500)
    pH_vals = titration_pH(V_titrant, titration_type, p_value, V_initial, 0.1, 0.1)
    V_eq = V_initial
    pH_at_eq = titration_pH(np.array([V_eq]), titration_type, p_value, V_initial, 0.1, 0.1)[0]
    V_half_eq = V_eq / 2
    pH_half_eq = titration_pH(np.array([V_half_eq]), titration_type, p_value, V_initial, 0.1, 0.1)[0]
    
    if titration_type == 'Weak Acid with Strong Base':
        specific_xlabel = 'Volume of Strong Base added (mL)'
        specific_title = f'Titration of Weak Acid (pKₐ = {p_value}) with Strong Base'
    elif titration_type == 'Weak Base with Strong Acid':
        specific_xlabel = 'Volume of Strong Acid added (mL)'
        specific_title = f'Titration of Weak Base (pKₐ = {p_value}) with Strong Acid'
    elif titration_type == 'Strong Acid with Strong Base':
        specific_xlabel = 'Volume of Strong Base added (mL)'
        specific_title = 'Titration of Strong Acid with Strong Base'
    elif titration_type == 'Strong Base with Strong Acid':
        specific_xlabel = 'Volume of Strong Acid added (mL)'
        specific_title = 'Titration of Strong Base with Strong Acid'
    
    xlabel = specific_xlabel if show_description else "Volume of Titrant added (mL)"
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(V_titrant * 1000, pH_vals, label="Titration Curve")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("pH")
    if show_description:
        ax.set_title(specific_title)
    if show_eq_line:
        ax.axvline(x=V_eq * 1000, color="grey", linestyle="--", label="Equivalence Volume")
    if show_eq_ph_line:
        ax.axhline(y=pH_at_eq, color="purple", linestyle="--", label=f"Equiv. pH: {pH_at_eq:.2f}")
    if show_half_eq_line:
        ax.axhline(y=pH_half_eq, color="green", linestyle=":", label=f"Half-Equiv. pH: {pH_half_eq:.2f}")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 14)
    ax.legend()
    ax.grid(True)
    return fig

# -------------------------------
# Quiz Data and Logic

descriptions = {
    'Weak Acid with Strong Base': (
        "Titration of Weak Acid with Strong Base:\n"
        "• The pH starts higher than expected for a strong acid because the weak acid is only partially dissociated.\n"
        "• A buffer region is present due to significant amounts of both the weak acid and its conjugate base.\n"
        "• The pH at the equivalence point is greater than 7 because the conjugate base hydrolyzes to form OH⁻."
    ),
    'Weak Base with Strong Acid': (
        "Titration of Weak Base with Strong Acid:\n"
        "• The pH starts lower than expected for a strong base because the weak base is only partially ionized.\n"
        "• A buffer region exists due to the presence of both the weak base and its conjugate acid.\n"
        "• The pH at the equivalence point is less than 7 because the conjugate acid donates H⁺."
    ),
    'Strong Acid with Strong Base': (
        "Titration of Strong Acid with Strong Base:\n"
        "• The pH starts very low due to the complete ionization of the strong acid.\n"
        "• No buffer region is observed as the reaction goes to completion.\n"
        "• The equivalence point is neutral (around pH 7) due to complete neutralization."
    ),
    'Strong Base with Strong Acid': (
        "Titration of Strong Base with Strong Acid:\n"
        "• The pH starts very high due to the complete ionization of the strong base.\n"
        "• No buffer region is observed because the reaction goes to completion.\n"
        "• The equivalence point is neutral (around pH 7) due to complete neutralization."
    )
}

def generate_scenario():
    scenario_type = random.choice([
        "Weak Acid with Strong Base",
        "Weak Base with Strong Acid",
        "Strong Acid with Strong Base",
        "Strong Base with Strong Acid"
    ])
    if scenario_type == "Weak Acid with Strong Base":
        possible_values = [4.0, 4.5, 5.0, 5.5]
        correct_pKa = random.choice(possible_values)
    elif scenario_type == "Weak Base with Strong Acid":
        possible_values = [9.0, 9.5, 10.0, 10.5]
        correct_pKa = random.choice(possible_values)
    else:
        possible_values = ["N/A"]
        correct_pKa = "N/A"
    if scenario_type == "Weak Acid with Strong Base":
        eq_pH = "greater than 7"
        buffer_ans = "Yes"
    elif scenario_type == "Weak Base with Strong Acid":
        eq_pH = "less than 7"
        buffer_ans = "Yes"
    else:
        eq_pH = "7"
        buffer_ans = "No"
    V_acid = random.choice([0.03, 0.04, 0.05, 0.06])
    return {
        "type": scenario_type,
        "pKa": correct_pKa,
        "pKa_options": possible_values,
        "eq_pH": eq_pH,
        "buffer": buffer_ans,
        "V_acid": V_acid
    }

# -------------------------------
# Streamlit Session State Initialization

if "scenario" not in st.session_state:
    st.session_state.scenario = generate_scenario()
if "submitted" not in st.session_state:
    st.session_state.submitted = False
if "feedback" not in st.session_state:
    st.session_state.feedback = ""

# -------------------------------
# Streamlit App Layout

st.title("Titration Quiz")

# Get the current scenario
scenario = st.session_state.scenario
p_val = scenario["pKa"] if scenario["pKa"] != "N/A" else 0

# Display the titration plot
if st.session_state.submitted:
    fig = plot_titration(scenario["type"], p_val, scenario["V_acid"], True, True, True, True)
else:
    fig = plot_titration(scenario["type"], p_val, scenario["V_acid"], False, False, False, False)
st.pyplot(fig)

# If the answer has been submitted, show the descriptive text.
if st.session_state.submitted:
    st.markdown(descriptions[scenario["type"]])

st.header("Answer the following:")

with st.form(key="quiz_form"):
    user_type = st.radio("Titration Type:", 
                         options=["Weak Acid with Strong Base", "Weak Base with Strong Acid",
                                  "Strong Acid with Strong Base", "Strong Base with Strong Acid"])
    user_pH = st.radio("Equivalence pH:", options=["7", "less than 7", "greater than 7"])
    user_buffer = st.radio("Buffer Region:", options=["Yes", "No"])
    
    # For pKₐ: if scenario is weak acid/base, generate numeric options; otherwise, use "N/A".
    if scenario["type"] in ["Weak Acid with Strong Base", "Weak Base with Strong Acid"]:
        numeric_options = list(scenario["pKa_options"])
        correct = scenario["pKa"]
        if correct not in numeric_options:
            numeric_options.append(correct)
        others = [opt for opt in numeric_options if opt != correct]
        if len(others) >= 2:
            chosen = random.sample(others, 2)
        else:
            chosen = others
        final_numeric = chosen + [correct]
        final_numeric = [str(x) for x in final_numeric]
        final_options = list(dict.fromkeys(final_numeric + ["N/A"]))  # Unique options
    else:
        final_options = ["N/A", "4.0", "7.0", "10.0"]
    user_pKa = st.selectbox("pKₐ:", options=final_options)
    
    # Equivalence Volume (in mL) based on the acid volume (in liters)
    correct_vol = int(scenario["V_acid"] * 1000)
    candidates = sorted({x for x in [correct_vol - 10, correct_vol, correct_vol + 10, correct_vol + 20] if 30 <= x <= 60})
    final_vol_options = [str(x) for x in candidates]
    user_vol = st.selectbox("Equivalence Volume (mL):", options=final_vol_options)
    
    submitted = st.form_submit_button(label="Submit Answer")

# -------------------------------
# Process Form Submission

if submitted:
    feedback_text = ""
    if user_type == scenario["type"]:
        feedback_text += "Titration Type: Correct!\n"
    else:
        feedback_text += f"Titration Type: Incorrect. Correct answer: {scenario['type']}\n"
    if user_pH == scenario["eq_pH"]:
        feedback_text += "Equivalence pH: Correct!\n"
    else:
        feedback_text += f"Equivalence pH: Incorrect. Correct answer: {scenario['eq_pH']}\n"
    if user_buffer == scenario["buffer"]:
        feedback_text += "Buffer Region: Correct!\n"
    else:
        feedback_text += f"Buffer Region: Incorrect. Correct answer: {scenario['buffer']}\n"
    if scenario["pKa"] != "N/A":
        try:
            if abs(float(user_pKa) - float(scenario["pKa"])) < 0.1:
                feedback_text += "pKₐ: Correct!\n"
            else:
                feedback_text += f"pKₐ: Incorrect. Correct answer: {scenario['pKa']}\n"
        except:
            feedback_text += f"pKₐ: Incorrect. Correct answer: {scenario['pKa']}\n"
    else:
        if user_pKa == "N/A":
            feedback_text += "pKₐ: Correct!\n"
        else:
            feedback_text += "pKₐ: Incorrect. Correct answer: N/A\n"
    if user_vol == str(correct_vol):
        feedback_text += "Equivalence Volume: Correct!\n"
    else:
        feedback_text += f"Equivalence Volume: Incorrect. Correct answer: {correct_vol}\n"
    
    st.session_state.feedback = feedback_text
    st.session_state.submitted = True
    st.experimental_rerun()  # Rerun to update the display with feedback and updated plot

# Display feedback if submitted
if st.session_state.submitted:
    st.subheader("Feedback:")
    st.text(st.session_state.feedback)
    if st.button("Next Question"):
        st.session_state.scenario = generate_scenario()
        st.session_state.submitted = False
        st.session_state.feedback = ""
        st.experimental_rerun()
