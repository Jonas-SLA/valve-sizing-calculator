// ===================================================================================
// Valve Cv Calculator Script
// ===================================================================================
// This script calculates the required flow coefficient (Cv) for a control valve
// based on fluid properties, process conditions, and valve type.
//
// Units used in calculations:
// Q = Flow Rate in US Gallons per Minute (GPM)
// P = Pressure in Pounds per Square Inch (PSI)
// D = Diameter in Inches
// ===================================================================================

/*******
    # preconditions: user inputs read, units known (psig/psia)
    P1 = toAbsPressure(input.inletPressure, input.isGauge)
    P2 = toAbsPressure(input.outletPressure, input.isGauge)
    Q  = input.flowRate  # GPM
    T  = input.temperature
    D  = input.pipeDiameter # in

    # fluid props
    if input.Pv_provided:
    Pv = input.Pv  # must be psia
    else:
    Pv = computeVaporPressure(fluid, T)  # Antoine or table

    rho = input.density || densityFromTemp(fluid, T)
    SG = rho / 1000.0
    mu_cP = input.viscosity || estimateViscosity(fluid, T)  # cP
    nu_cSt = mu_cP / (rho/1000)  # convert dynamic->kinematic

    FL = valveDefaults[ valveType ].Fl
    Pc = input.Pc || lookupCriticalPressure(fluid) # psia if available

    FF = calculateFF(Pv, Pc)
    deltaP = P1 - P2
    deltaP_choked = FL*FL * (P1 - FF*Pv)
    choked = (deltaP >= deltaP_choked)

    # initial Cv
    if (choked) Cv0 = Q / (N1 * FL) * sqrt(SG / (P1 - FF*Pv))
    else Cv0 = Q * sqrt(SG / deltaP)

    # Rev and FR iterations
    Cv = Cv0
    Re = calcRev(Cv, nu_cSt, ...); 
    if Re <= 10000:
    for i in range(maxIter):
        Re = calcRev(Cv, nu_cSt, ...)
        FR = chooseFRbranch(Cv, Re, D, FL)
        Cv_new = Cv0 / FR
        if converge(Cv_new, Cv): break
        Cv = Cv_new

    # Cv is final required
    # Selection: for each candidate valve in manufacturer_db:
    #   for each size:
    #     Cv_at_design_opening = interp(manufacturer.Cv_table, desired_opening)
    #     if Cv_at_design_opening >= Cv * safety_margin and rangeability_ok:
    #         push candidate into picks
    # sort picks by size asc and present best matches
*******/

//------------------------------------------------------------------------------------
// SECTION 0: DEBUG SCRIPT TOOLS
//------------------------------------------------------------------------------------
/**
 * @const {boolean} DEBUG - Global switch to turn console logs on/off.
 */
// Debug Tools Enable/Disable
const DEBUG = true;

/**
 * Custom logger that respects the DEBUG flag and writes to both
 * the browser console and the on-screen HTML log.
 * @param {string} message - The log message.
 * @param {...any} args - Additional values to log.
 */
function debugLog(message, ...args) {
    // --- 1. Log to the browser console (F12) ---
    // This is still useful for inspecting objects
    if (DEBUG) {
        if (message.startsWith("FAIL")) {
            console.error(`[DEBUG] ${message}`, ...args);
        } else if (message.startsWith("---")) {
            console.warn(message, ...args); // Use warn for major steps
        } else {
            console.log(`[DEBUG] ${message}`, ...args);
        }
    }

    // --- 2. Log to the on-screen HTML element ---
    const logOutputEl = document.getElementById("debugLogOutput");
    if (DEBUG && logOutputEl) {
        // Format the message
        let fullMessage = `[DEBUG] ${message}`;

        // Handle the ...args
        if (args.length > 0) {
            // Convert objects/arrays to a readable JSON string
            const formattedArgs = args.map(arg => {
                if (typeof arg === 'object' && arg !== null) {
                    try {
                        // Use JSON.stringify to show object contents
                        return JSON.stringify(arg, null, 2); 
                    } catch (e) {
                        return "[Circular Object]";
                    }
                }
                return String(arg); // Convert numbers, booleans, etc.
            }).join(' ');
            
            fullMessage += ` ${formattedArgs}`;
        }
        
        // Add a newline
        fullMessage += "\n";

        // Append the new message to the <pre> tag's text
        logOutputEl.textContent += fullMessage;

        // Auto-scroll to the bottom so the latest message is always visible
        logOutputEl.scrollTop = logOutputEl.scrollHeight;
    }
}

//------------------------------------------------------------------------------------
// SECTION 1: GLOBAL CONSTANTS & DATA DEFINITIONS
//------------------------------------------------------------------------------------
//Para outros fluidos, pedir informação: PC(Ponto Crítico), Pressão de Vapor, Densidade Relativa(ou densidade do fluido dividido pela da água)
/**
 * @const {object} numericalConstants
 * @description ANSI/ISA-75.01.01 standard numerical constants for Imperial units.
 */
const P_ATMOSPHERIC_PSI = 14.7;   // Atmospheric pressure in psi
const D_REFERENCE_WATER = 1000.0; // Reference density in kg/m^3

/**
 * @const {object} conversionConstants
 * @description Conversion factors to standardized internal units.
 */
const conversionConstants = {
  KPA_TO_PSI: 0.145038,
  M3H_TO_GPM: 4.40287,
  MM_TO_IN: 0.0393701,
  M2S_TO_CST: 1000000, // 1 m²/s = 1,000,000 cSt
};

const numericalConstants = {
    N1:  1,     // Units conversion factor for Cv (US units)        
    N2:  890,   // Velocity head factor constant
    N4:  17300, // Reynolds number constant
    N18: 645,   // Sizing factor constant
    N32: 17,    // Sizing factor constant
};

/**
 * @const {object} liquidWater
 * @description Default thermophysical properties for water.
 */
const liquidWater = {
    isSelected: false,
    Pc: 3199.308,       // Thermodynamic critical pressure (PSIa)
    Pv: 0,              // Vapor pressure at inlet temperature (PSIa) - will be calculated
    relativeDensity: 0, // Relative density (Specific Gravity) - will be calculated
    viscosity: 0,       // kinematic viscosity (cSt) - will be calculated
};

/**
 * @const {object} liquidCustom
 * @description User-defined fluid for custom entry i am of thermophysical properties.
 */
const liquidCustom = {
    isSelected: false,
    Pc: 0,              // Critical pressure (PSIa)
    Pv: 0,              // Vapor pressure (PSIa)
    relativeDensity: 0, // Relative density (Specific Gravity)
    viscosity: 0,       // kinematic viscosity input by user (cSt) - will be calculated
};

// --- Valve Style Specific Factors ---
const valveGlobesSinglePortContouredPlugOpen = {
    Fl: 0.9,    // Liquid pressure recovery factor        
    Fd: 0.46,   // Valve style modifier
    Cv: 0,      
    Xt: 0.72,   // Differential pressure ratio factor (for gases)
};

const valveGlobesSinglePortContouredPlugClose = {
    Fl: 0.8,
    Fd: 1.0,
    Cv: 0,
    Xt: 0.55,
};

const valveButterfly = {
    Fl: 0.62,
    Fd: 0.57,
    Cv: 0,
    Xt: 0.35,
};

const valveBallFull = {
    Fl: 0.74,
    Fd: 0.99,
    Cv: 0,
    Xt: 0.42,
};

//------------------------------------------------------------------------------------
// SECTION 2: STATE MANAGEMENT MODULE
//------------------------------------------------------------------------------------

/**
 * @module UserInputs
 * @description Manages the application's state (user inputs, fluid/valve data)
 * using an IIFE (Immediately Invoked Function Expression) to create a private scope.
 */
const UserInputs = (() => 
{
    // --- Private State ---
    let inputs = {
        inletPressure: 0,
        outletPressure: 0,
        temperature: 20,   // default 20°C
        flowRate: 0,
        pipeDiameter: 0,
        fluidType: 0,
        valveType: 0,
    };

    const Fluids = Object.freeze({ WATER: 0, CUSTOM: 1 });
    const Valves = Object.freeze({ GLOBEOPEN: 0, GLOBECLOSE: 1, BUTTERFLY: 2, BALL: 3 });

    const fluidData = {
        [Fluids.WATER]: { ...liquidWater, isSelected: true },
        [Fluids.CUSTOM]: { ...liquidCustom, isSelected: false },
    };

    const valveData = {
        [Valves.GLOBEOPEN]: { ...valveGlobesSinglePortContouredPlugOpen },
        [Valves.GLOBECLOSE]: { ...valveGlobesSinglePortContouredPlugClose },
        [Valves.BUTTERFLY]: { ...valveButterfly },
        [Valves.BALL]: { ...valveBallFull },
    };

    /**
     * Safely retrieves a nested property from an object using a string path.
     * @param {object} obj - The object to query.
     * @param {string} path - The dot-separated path to the property (e.g., "getSelectedFluid.Pv").
     * @returns {*} The value of the property, or undefined if not found.
     */    
    const getNestedProperty = (obj, path) => {
        try {
            return path.split(".").reduce((acc, key) => {
                if (typeof acc === "function") acc = acc(); // call methods like getSelectedFluid
                if (acc == null) return undefined;
                return acc[key];
            }, obj);
        } catch (e) {
            return undefined;
        }
    };

    // --- Public Interface ---
    return {
        getInputs: () => ({ ...inputs }),
        getSelectedFluid: () => ({ ...fluidData[inputs.fluidType] }),
        getSelectedValve: () => ({ ...valveData[inputs.valveType] }),

        /**
         * Sets a value in the main user inputs state.
         * @param {string} key - The property name in the inputs object.
         * @param {*} value - The value to set.
         */
        setInput: (key, value) => {
            if (!(key in inputs)) {
                console.warn(`Invalid key: ${key}`);
                return;
            }
            inputs[key] = value;
            if (key === "fluidType" && fluidData[value] !== undefined) {
                Object.keys(fluidData).forEach(f => (fluidData[f].isSelected = false));
                fluidData[value].isSelected = true;
            }
            if (key === "valveType" && valveData[value] === undefined) {
                console.warn("Invalid valve type index");
            }
        },

        /**
         * Sets a property on the currently selected fluid object.
         * @param {string} key - The property name (e.g., "viscosity").
         * @param {*} value - The value to set.
         */
        setFluidProperty: (key, value) => {
            const selectedFluid = fluidData[inputs.fluidType];
            if (selectedFluid && key in selectedFluid) {
                selectedFluid[key] = value;
            } else {
                console.warn(`Invalid fluid property: ${key}`);
            }
        },

        setValveProperty: (key, value) => {
            const selectedValve = valveData[inputs.valveType];
            if (selectedValve && key in selectedValve) {
                selectedValve[key] = value;
            } else {
                console.warn(`Invalid valve property: ${key}`);
            }
        },

        resetInputs: () => {
            inputs = {
                inletPressure: 0,
                outletPressure: 0,
                temperature: 20,
                flowRate: 0,
                pipeDiameter: 0,
                fluidType: Fluids.WATER,
                valveType: Valves.GLOBEOPEN,
            };
            Object.keys(fluidData).forEach(f => (fluidData[f].isSelected = parseInt(f) === Fluids.WATER));
        },

        getProperty: (path) => getNestedProperty(UserInputs, path),

        Fluids,
        Valves,
    };
})();

//------------------------------------------------------------------------------------
// SECTION 3: FLUID PROPERTY CALCULATIONS
//------------------------------------------------------------------------------------
// ---------- Density & specific gravity for water ----------
/**
 * Calculates the density of water based on temperature, pressure, and salinity.
 * Uses a polynomial approximation for temperature and simplified adjustments for pressure and salinity.
 * @param {number} [S_percent=0] - Salinity in percent (%).
 * @returns {number|null} The final water density in kg/m³, or null on error.
 */
/**
 * Calculates the density of water from:
 *   • temperature (°C) — taken automatically from UserInputs
 *   • pressure (psi)  — taken automatically from UserInputs
 *   • salinity (%)    — provided as argument
 *
 * Based on:
 *   • Temperature polynomial for pure water
 *   • Linearized pressure correction
 *   • Practical empirical salinity correction
 *
 * RETURNS:
 *   density in kg/m³
 */
function waterDensityCalc(S_percent = 0) {
    // Temperature (°C) and pressure (psi) come from user inputs
    const pressure = UserInputs.getProperty("getInputs.inletPressure");
    const temperature = UserInputs.getProperty("getInputs.temperature");

    // ---- 1) Pure water density polynomial in °C ----
    const rho_0 = 999.83311;
    const a1 = 0.0752;
    const a2 = 0.0089;
    const a3 = 7.36413e-5;
    const a4 = 4.74639e-7;
    const a5 = 1.34888e-9;

    let rho_T = rho_0 + (a1 * temperature) - (a2 * Math.pow(temperature, 2)) + (a3 * Math.pow(temperature, 3))
                - (a4 * Math.pow(temperature,4)) + (a5 * Math.pow(temperature,5));

    // ---- 2) Pressure correction ----
    // Convert psi → dbar (approx). Water density rises ~4.5e-3 kg/m³ per dbar.
    let P_dbar = 0;
    if (typeof pressure === "number" && !Number.isNaN(pressure)) {
        P_dbar = pressure * 0.689476;
    }

    // ---- 3) Salinity correction ----
    // Convert % → ‰ (practical oceanographic convention)
    let S = S_percent * 10; // % -> ‰
    let pressure_effect = 4.5e-3 * P_dbar;

    let salinity_adjustment = 0;
    if (S > 0) {
        salinity_adjustment = 0.824493 * S - 0.0040899 * S * temperature + 0.000076438 * S * Math.pow(temperature, 2)
                              + 0.0000082467 * Math.pow(S, 1.5) + 0.0000005866 * Math.pow(S,2) + 0.0000053875 * P_dbar * S;
    }

    // ---- Final density ----
    const rho_final = rho_T + pressure_effect + salinity_adjustment;
    if (!isFinite(rho_final) || rho_final <= 0) {
        console.warn("Computed density invalid.");
        return null;
    }
    return rho_final;
}

// ---------- Utility / thermophysical estimators ----------
/**
 * Estimates the dynamic viscosity of water using a Vogel-like empirical formula.
 * This is a practical estimator for water between approximately 0-100 °C.
 * @param {number} T - Temperature in Celsius (°C).
 * @returns {number} The estimated dynamic viscosity in centipoise (cP).
 * 
 *  // ---------- Dynamic viscosity of water (μ) in cP ----------
    /**
     * Estimates the dynamic viscosity of water using a Vogel/Andrade-type
     * empirical correlation. Reasonably accurate for ~0–100 °C.
     *
     * INPUT:
     *   T  = temperature in Celsius (from user)
     *
     * OUTPUT:
     *   dynamic viscosity in centipoise (cP)
     *
     * NOTES:
     *   • Temperature is first converted to Kelvin (Tk)
     *   • Uses a logarithmic viscosity correlation calibrated for water
     *   • Clamps output to reasonable physical bounds to avoid absurd values
     *
*/
function estimateWaterViscosity_cP(T) {
    if (typeof T !== "number" || Number.isNaN(T)) return 1.0;
    // Convert Celsius → Kelvin (required by the Vogel/Andrade expression)
    const Tk = T + 273.15;

    // Vogel/Andrade constants for liquid water
    const A = 2.414e-5; // Pa·s
    const B = 247.8;
    const C = 140;

    // Compute the formula 
    let mu_Pa_s = A * Math.pow(10, B / (Tk - C));

    // Convert Pa·s → cP (1 cP = 0.001 Pa·s)
    let mu_cP = mu_Pa_s * 1000.0;

    // Safety clamps — prevents numerical blowups at edge temperatures
    if (mu_cP < 0.25) mu_cP = 0.25;
    if (mu_cP > 1000) mu_cP = 1000;

    return mu_cP;
}

// ---------- Kinematic viscosity (ν) in cSt ----------
/**
 * Estimates water kinematic viscosity in centistokes (cSt).
 *
 * ν = μ / ρ
 *
 * INPUT:
 *   T_C = temperature in Celsius
 *
 * OUTPUT:
 *   kinematic viscosity in cSt
 *
 * NOTES:
 *   • Dynamic viscosity μ is obtained from estimateWaterViscosity_cP()
 *   • Density ρ is obtained from waterDensityCalc()
 *   • Density function already reads temperature from UserInputs internally
 *   • Water salinity is fixed at 0% here (fresh water)
 */
function estimateWaterKinematicViscosity_cSt(T_C) {
    // Step 1: Calculate dynamic viscosity (cP)
    const mu_cP = estimateWaterViscosity_cP(T_C); // dynamic in cP

    // Step 2: density (kg/m³) — salinity input only, temperature taken internally
    const rho_kg_m3 = waterDensityCalc(0); // fresh water, no salinity
    if (!rho_kg_m3) return null;

    // Convert density to g/cm³ (1 g/cm³ = 1000 kg/m³)
    const rho_g_cm3 = rho_kg_m3 / 1000.0;
    
    // Engineering approximation: ν(cSt) ≈ μ(cP) / ρ(g/cm³)
    return mu_cP / rho_g_cm3;
}

/**
 * Calculates the specific gravity (relative density) of water and updates the state.
 * @returns {number|null} The calculated specific gravity, or null on error.
 */
// ---------- Specific Gravity (SG = ρ_water / ρ_ref) ----------
/**
 * Computes the specific gravity of water from the internally computed density.
 *
 * SG = ρ_water / ρ_reference
 *
 * NOTES:
 *   • ρ_water obtained from waterDensityCalc()
 *   • Salinity is set to 0% (fresh water)
 *   • Temperature and inlet pressure taken automatically inside waterDensityCalc()
 */
function calculateSpecificGravity() {
    const salinityPercent = 0;
    // Step 1: compute water density (kg/m³)
    const waterDensity = waterDensityCalc(salinityPercent);
    if (!waterDensity) return null;

    // Step 2: reference density (usually 999 kg/m³ at 4 °C)
    const referenceDensity = D_REFERENCE_WATER;

    // Step 3: SG = ρ / ρ_ref
    const specificGravity = waterDensity / referenceDensity;

    // Save into the fluid properties
    UserInputs.setFluidProperty("relativeDensity", specificGravity);

    return specificGravity;
}

/**
 * Estimates the vapor pressure of water in PSI using the Antoine equation.
 * Valid for temperatures between 1°C and 100°C.
 * @param {number} T - Temperature in Celsius (°C).
 * @returns {number} The estimated vapor pressure in PSI.
 */
function estimateWaterVaporPressure(T) {
    if (T <= 0) return 0.0886; // Pv at 0°C
    if (T >= 100) return 14.696; // Pv at 100°C

    // Antoine constants for water (T in °C, P in mmHg)
    const A = 8.07131;
    const B = 1730.63;
    const C = 233.426;

    // P_mmHg = 10^(A - (B / (C + T)))
    const P_mmHg = Math.pow(10, A - (B / (C + T)));
    
    // Convert mmHg to PSI (1 PSI = 51.715 mmHg)
    const P_PSI = P_mmHg / 51.715;
    
    return P_PSI;
}

//------------------------------------------------------------------------------------
// SECTION 4: CORE SIZING LOGIC
//------------------------------------------------------------------------------------

/**
 * Calculates the liquid pressure recovery factor (FF) for flashing service.
 * Based on the formula: FF = 0.96 - 0.28 * sqrt(Pv/Pc).
 * @returns {number|null} The calculated FF value, or null if inputs are invalid.
 */
/**
 * Calculates the liquid pressure recovery factor (FF) for flashing service.
 *
 * Purpose:
 *   FF is an empirical factor used in the ISA liquid sizing equations to account
 *   for the vaporization tendency of the liquid (flashing). It reduces the
 *   effective available pressure head when vapor pressure is significant.
 *
 * Formula & notes:
 *   FF = 0.96 - 0.28 * sqrt(Pv / Pc)
 *   - Pv = vapor pressure of the liquid (psi)
 *   - Pc = critical pressure of the fluid (psi)
 *   - The square-root term quantifies how close the fluid is to its critical state.
 *   - FF is clamped to a small positive lower bound (0.05) to avoid numerical issues.
 *
 * Returns:
 *   number (FF) or null on invalid inputs.
 */
function calculateFF() {
    const vaporPressure = UserInputs.getProperty("getSelectedFluid.Pv");
    const criticalPressure = UserInputs.getProperty("getSelectedFluid.Pc");

    // Validate input presence and types
    if (typeof vaporPressure !== "number" || typeof criticalPressure !== "number") {
        console.warn("Pv or Pc is not a valid number.");
        return null;
    }
    if (criticalPressure <= 0) {
        console.warn("Critical pressure must be > 0.");
        return null;
    }
    const ratio = vaporPressure / criticalPressure;
    if (ratio < 0 || ratio > 1) {
        // Not fatal — warn in case of user data entry error or unusual fluid
        console.warn("Pv/Pc ratio outside 0..1 — check inputs.");
    }

    // Core formula
    let FF = 0.96 - 0.28 * Math.sqrt(Math.max(0, ratio));

     // Practical lower bound: prevents negative FF and maintains numeric stability
    if (FF < 0.05) FF = 0.05;

    return FF;
}

/**
 * Determines if the flow through the valve is choked (critical).
 * Flow is choked if the actual pressure drop (ΔP) meets or exceeds the critical pressure drop limit.
 * @param {number} FF - The liquid pressure recovery factor calculated by `calculateFF`.
 * @returns {boolean|null} Returns `true` for choked flow, `false` for subcritical (non-choked) flow, or null on error.
 */
/**
 * Determines whether the flow through the valve is choked (critical) for liquid.
 *
 * Purpose:
 *   Uses the ISA/ANSI choked-liquid criterion:
 *     ΔP_crit = F_L^2 * (P1 - F_F * P_v)
 *   Flow is choked if actual ΔP >= ΔP_crit.
 *
 * Inputs read from UserInputs:
 *   - FL (liquid FL factor)
 *   - Pv (vapor pressure)
 *   - inletPressure (P1)
 *   - outletPressure (P2)
 *   - FF (passed in or previously computed)
 *
 * Returns:
 *   true  -> choked/critical flow
 *   false -> subcritical (non-choked)
 *   null  -> error / invalid inputs
 */
function isChokedFlow(FF) {
    const FL = UserInputs.getProperty("getSelectedValve.Fl");
    const vaporPressure = UserInputs.getProperty("getSelectedFluid.Pv");
    const inP = UserInputs.getProperty("getInputs.inletPressure");
    const outP = UserInputs.getProperty("getInputs.outletPressure");

    // Validate presence and numeric quality of inputs
    if ([FL, vaporPressure, inP, outP, FF].some(v => typeof v !== "number" || Number.isNaN(v))) {
        console.warn("Invalid input(s) in choked flow check.");
        return null;
    }
    if (inP <= outP) {
        console.warn("Inlet pressure must be > outlet pressure.");
        return null;
    }

    // Compute the ISA choked threshold (ΔP_crit)
    const chokedThreshold = Math.pow(FL, 2) * (inP - (FF * vaporPressure));
    if (!isFinite(chokedThreshold) || chokedThreshold < 0) {
        console.warn("Invalid choked threshold calculation.");
        return null;
    }

    const actualDP = inP - outP;
    // If actual ΔP ≥ threshold, the liquid flow is choked
    return actualDP >= chokedThreshold;
}

/**
 * Calculates the required Cv for choked (critical) flow conditions.
 * @returns {number|null} The calculated Cv, or null if inputs are invalid.
 */
/**
 * Calculates required Cv for choked (critical) liquid flow.
 *
 * Purpose & formula:
 *   For choked liquid flow ISA gives:
 *     C_v = (Q / (N1 * F_L)) * sqrt( G_f / (P1 - F_F * P_v) )
 *   where:
 *     - Q is flowrate (gpm)
 *     - N1 is unit conversion constant (numericalConstants.N1)
 *     - G_f is relative density (specific gravity)
 *     - denominator is the effective driving pressure for flashing case
 *
 * Notes:
 *   - This function reads inputs from UserInputs. It currently recomputes FF
 *     internally (matching your original structure).
 *
 * Returns:
 *   Cv (number) or null on error.
 */
function calculateC_choked(FF) {
    const Q = UserInputs.getProperty("getInputs.flowRate");
    const inP = UserInputs.getProperty("getInputs.inletPressure");
    const FL = UserInputs.getProperty("getSelectedValve.Fl");
    const vaporPressure = UserInputs.getProperty("getSelectedFluid.Pv");
    const relDensity = UserInputs.getProperty("getSelectedFluid.relativeDensity");

    if ([Q, inP, FL, vaporPressure, relDensity].some(v => typeof v !== "number" || Number.isNaN(v))) {
        console.warn("Missing inputs for choked C calculation.");
        return null;
    }

    const denom = (inP - (FF * vaporPressure));
    if (denom <= 0) {
        console.warn("Denominator <= 0 in choked C calc.");
        return null;
    }

    // Final choked Cv expression (N1 accounts for unit system)
    const Ci = (Q / (numericalConstants.N1 * FL)) * Math.sqrt((relDensity) / denom);
    return Ci;
}

/**
 * Calculates the required Cv for subcritical (non-choked) flow conditions.
 * @returns {number|null} The calculated Cv, or null if inputs are invalid.
 */
/**
 * Calculates required Cv for subcritical (non-choked) liquid flow.
 *
 * Purpose & formula:
 *   For non-choked (regular) liquid flow:
 *     C_v = (Q / N1) * sqrt( G_f / ΔP )
 *   where ΔP = P1 - P2.
 *
 * Returns:
 *   Cv (number) or null on error.
 */
function calculateC_subcritical() {
    const Q = UserInputs.getProperty("getInputs.flowRate");
    const inP = UserInputs.getProperty("getInputs.inletPressure");
    const outP = UserInputs.getProperty("getInputs.outletPressure");
    const relDensity = UserInputs.getProperty("getSelectedFluid.relativeDensity");
    const deltaP = inP - outP;

    if ([Q, inP, outP, relDensity].some(v => typeof v !== "number" || Number.isNaN(v))) {
        console.warn("Missing inputs for subcritical C calculation.");
        return null;
    }
    if (deltaP <= 0) {
        console.warn("deltaP must be > 0 for subcritical C calculation.");
        return null;
    }

    // Standard ISA subcritical Cv expression
    const Ci = (Q / numericalConstants.N1) * Math.sqrt((relDensity) / deltaP);
    return Ci;
}

/**
 * Calculates the valve Reynolds number (Rev).
 * This dimensionless number helps determine if the flow is laminar, transitional, or turbulent.
 * @param {number} Ci - An estimate of the valve's flow coefficient (Cv).
 * @returns {number|null} The calculated Reynolds number, or null on error.
 */
/**
 * Calculates the valve Reynolds number (Rev) used to determine viscosity correction needs.
 *
 * Purpose:
 *   Rev is an empirical valve-specific Reynolds-like number defined in ISA to decide
 *   whether the flow is turbulent or whether viscosity corrections (FR) must be applied.
 *
 * Formula:
 *   Rev = (N4 * F_D * Q) / (v * sqrt(Cv * F_L)) * [ 1 + (F_L^2 * Cv^2)/(N2 * D^4) ]^(1/4)
 *
 * Inputs:
 *   - Ci: valve Cv estimate (number)
 *   - Reads FD, Q, FL, D, and viscosity from UserInputs
 *
 * Requirements:
 *   - D (pipe diameter) must be > 0 (inches for the constant set)
 *   - - viscosity (v) should be in KINEMATIC units (centistokes, cSt) to match N4 constant
 *
 * Returns:
 *   Rev (number) or null on error.
 */
function calculateValveReynoldsRev(Ci, FLP, d) {
    if (typeof Ci !== "number" || Number.isNaN(Ci) || Ci <= 0) {
        console.warn("Invalid Ci for Rev calculation.");
        return null;
    }

    // Get the needed input data for the calculation
    const FD = UserInputs.getProperty("getSelectedValve.Fd");
    const Q = UserInputs.getProperty("getInputs.flowRate");

    // Get the variable parameters
    const FL = FLP;
    const D = d;

    // Pick the viscosity from fluid data (expected as cP)
    let v = UserInputs.getProperty("getSelectedFluid.viscosity");

    // Validate geometry & fluid inputs
    if (!D || D <= 0) {
        console.warn("Invalid pipe diameter for Rev calculation.");
        return null;
    }
    if (!v || v <= 0) {
        console.warn("Invalid viscosity for Rev calculation.");
        return null;
    }

    // Core equation broken down for clarity:
    // eq1 = (N4 * FD * Q) / (v * sqrt(Ci * FL))
    const eq1 = (numericalConstants.N4 * FD * Q) / (v * Math.sqrt(Ci * FL));

    // inside = 1 + (FL^2 * Ci^2) / (N2 * D^4)
    const inside = (Math.pow(FL, 2) * Math.pow(Ci, 2)) / (numericalConstants.N2 * Math.pow(D, 4)) + 1;
    const eq2 = Math.pow(inside, 0.25);

    const Rev = eq1 * eq2;

    // Check if the calculated Rev is valid 
    if (!isFinite(Rev) || Rev < 0) {
        console.warn("Invalid Rev computed.");
        return null;
    }

    return Rev;
}

/**
 * Calculates the Fp (Piping Geometry) and FLP (Installed) factors
 * for a valve with reducers.
 * @param {number} Cv - The current calculated Cv (or guess).
 * @param {number} valveSize - The nominal valve diameter (d).
 *S @param {number} pipeSize - The nominal pipe diameter (D).
 * @param {number} FL - The valve's intrinsic FL factor.
 * @returns {object} { Fp, FLP }
 */
function calculateFpFactors(Cv, valveSize, pipeSize, FL) {
    // If valve is line-sized, factors are 1.0 and FL.
    if (Math.abs(valveSize - pipeSize) < 0.001) {
        return { Fp: 1.0, FLP: FL };
    }

    const N2 = numericalConstants.N2; // 890
    const d = valveSize;
    const D = pipeSize;

    // Loss coefficients for standard concentric reducers
    const ratioSq = (d * d) / (D * D);
    const K1 = 0.5 * Math.pow(1 - ratioSq, 2);  // Inlet reducer
    const K2 = 1.0 * Math.pow(1 - ratioSq, 2);  // Outlet expander
    const KB1 = 1 - Math.pow(ratioSq, 2);     // Bernoulli (inlet)
    
    const SigmaK = K1 + K2;
    const Ki = K1 + KB1;

    // This term is used repeatedly
    const Cvd2_term = Math.pow(Cv / (d * d), 2);

    // Calculate Fp (Piping Geometry Factor)
    const Fp = Math.pow(1 + (SigmaK / N2) * Cvd2_term, -0.5);

    // Calculate FLP (Installed Pressure Recovery Factor)
    const FLP = Math.pow(
        (1 / (FL * FL)) + (Ki / N2) * Cvd2_term,
        -0.5
    );

    return { Fp, FLP };
}

//------------------------------------------------------------------------------------
// SECTION 5: VISCOSITY CORRECTION FACTOR (FR) LOGIC
//------------------------------------------------------------------------------------

/**
 * Determines which set of viscosity correction (FR) formulas to use based on valve size.
 * This is a check specified by the ANSI/ISA standard.
 * @param {number} C - The calculated flow coefficient (Cv).
 * @returns {boolean} `true` to use one set of formulas, `false` for the other.
 */
/**
 * Decide which FR formula branch to use (full-range vs restricted-range).
 *
 * Purpose:
 *   ISA uses a boundary (Cv / D^2 comparison) to choose between two FR families.
 *   This helper implements that check using your constant scaling.
 *
 * Returns:
 *   true  -> use "Yes" branch (full-range formulas)
 *   false -> use "No" branch (restricted-range formulas)
 */
function verifyFR(C, d) {
    const D = d;
    if (!D || D <= 0) return true; // safest default is full-range

    // boundary comparison: C / D^2  >  0.016 * N18  (your chosen scaling)
    const div = C / Math.pow(D, 2);
    const mul = 0.016 * numericalConstants.N18;

    return div > mul;
}

/**
 * Calculates the viscosity correction factor (FR) for one of the conditions from `verifyFR`.
 * @param {number} Ci - The estimated flow coefficient (Cv).
 * @param {number} Rev - The calculated valve Reynolds number.
 * @returns {number|null} The viscosity correction factor FR (always <= 1), or null on error.
 */
/**
 * Viscosity correction factor FR — restricted-range branch (verifyFR == false)
 *
 * Purpose:
 *   Compute FR (≤ 1) that reduces Cv for viscous (low-Re) flows.
 *   Two limiting forms are used and the smallest physically-valid value is chosen.
 *
 * Inputs:
 *   Ci  - current Cv estimate
 *   Rev - valve Reynolds number
 *
 * Returns:
 *   FR (0 < FR ≤ 1) or null on error
 */
function calculateFRifNo(Ci, Rev, d_pipe, FLP) {
    if (Ci == null || Rev == null) {
        console.warn("Missing inputs to calculate FR (No).");
        return null;
    }

    // Get the variable parameters
    let FR = 0;
    const FL = FLP;
    const d = d_pipe;

    // empirical combination factor n2
    const n2 = 1 + numericalConstants.N32 * Math.pow((Ci / (d * d)), 2/3);

    if (Rev < 10) {
        // very low Re asymptote (laminar-like)
        FR = (0.026 / FL) * Math.sqrt(n2 * Rev);
        FR = Math.min(FR, 1);
    } else {
        // general case: choose min(F3A, F4, 1)
        const parte1 = (0.33 * Math.sqrt(FL)) / Math.pow(n2, 0.25);
        const parte2 = Math.log10(Rev / 10000);
        const FR_F3A = 1 + parte1 * parte2;

        const FR_F4 = (0.026 / FL) * Math.sqrt(n2 * Rev);

        // pick the smallest physically-restrictive correction (and cap at 1)
        FR = Math.min(FR_F3A, FR_F4, 1);
    }

    return FR;
}

/**
 * Calculates the viscosity correction factor (FR) for the other condition from `verifyFR`.
 * @param {number} Ci - The estimated flow coefficient (Cv).
 * @param {number} Rev - The calculated valve Reynolds number.
 * @returns {number|null} The viscosity correction factor FR (always <= 1), or null on error.
 */
/**
 * Viscosity correction factor FR — full-range branch (verifyFR == true)
 *
 * Purpose:
 *   The alternate FR formulation used for valves with different Cv/D^2 sizing.
 *   Same selection logic as restricted-range: choose the limiting expression.
 *
 * Inputs:
 *   Ci  - current Cv estimate
 *   Rev - valve Reynolds number
 *
 * Returns:
 *   FR (0 < FR ≤ 1) or null on error
 */
function calculateFRifYes(Ci, Rev, d_pipe, FLP) {
    if (Ci == null || Rev == null) {
        console.warn("Missing inputs to calculate FR (Yes).");
        return null;
    }

    // Get the variable parameters
    let FR = 0;
    const FL = FLP;
    const d = d_pipe;

    // empirical parameter n1 used in the full-range formulas
    const n1 = numericalConstants.N2 * (Ci / (d * d));

    if (Rev < 10) {
        FR = (0.026 / FL) * Math.sqrt(n1 * Rev);
        FR = Math.min(FR, 1);
    } else {
        const parte1 = (0.33 * Math.sqrt(FL)) / Math.pow(n1, 0.25);
        const parte2 = Math.log10(Rev / 10000);
        const FR_F1A = 1 + parte1 * parte2;
        const FR_F2 = (0.026 / FL) * Math.sqrt(n1 * Rev);
        FR = Math.min(FR_F1A, FR_F2, 1);
    }

    return FR;
}

//------------------------------------------------------------------------------------
// SECTION 6: MAIN SIZING ORCHESTRATOR
//------------------------------------------------------------------------------------

/**
 * Main valve sizing function that orchestrates the calculation process.
 * Now includes Fp/FLP for reducers and is iterative.
 * @param {number} valveSize - The nominal size (d) of the valve being tested.
 * @returns {object|null} { Cv, Re, iterations, choked } or null on error
 */
function valveSizing(valveSize) 
{
    const EPS = 0.0001;      // Convergence tolerance
    const MAX_ITERS = 20;   // Max iterations for Cv and FR

    const pipeSize = UserInputs.getProperty("getInputs.pipeDiameter");
    const FL = UserInputs.getProperty("getSelectedValve.Fl"); // Generic or specific FL
    const Q = UserInputs.getProperty("getInputs.flowRate");
    const SG = UserInputs.getProperty("getSelectedFluid.relativeDensity");
    const N1 = numericalConstants.N1;

    // Get pressures
    const inP = UserInputs.getProperty("getInputs.inletPressure");
    const outP = UserInputs.getProperty("getInputs.outletPressure");
    const actualDP = inP - outP;
    if (actualDP <= 0) {
        console.warn("deltaP must be > 0.");
        return null;
    }
    
    // --- Step 1: Calculate FF (as before) ---
    const FF = calculateFF();
    if (FF == null) return null;

    // --- Step 2: Iteratively find Initial (Turbulent) Cv ---
    let Cv_initial = Q * Math.sqrt(SG / actualDP); // First guess
    let Fp = 1.0;
    let FLP = FL;

    for (let i = 0; i < MAX_ITERS; i++) {
        // Get new Fp/FLP based on the last Cv guess
        const factors = calculateFpFactors(Cv_initial, valveSize, pipeSize, FL);
        Fp = factors.Fp;
        FLP = factors.FLP;

        // Check for choked flow using *installed* FLP
        const chokedThreshold = Math.pow(FLP, 2) * (inP - (FF * UserInputs.getProperty("getSelectedFluid.Pv")));
        const choked = (actualDP >= chokedThreshold);

        let Cv_new = 0;
        if (choked) {
            // Use choked flow formula with FLP
            const denom = (inP - (FF * UserInputs.getProperty("getSelectedFluid.Pv")));
            if (denom <= 0) return null; // Error
            Cv_new = (Q / (N1 * FLP)) * Math.sqrt(SG / denom);
        } else {
            // Use subcritical formula with Fp
            Cv_new = (Q / (N1 * Fp)) * Math.sqrt(SG / actualDP);
        }

        // Check for convergence
        if (Math.abs(Cv_new - Cv_initial) / Cv_new < EPS) {
            Cv_initial = Cv_new; // Converged
            break;
        }
        Cv_initial = Cv_new; // Update guess
        if (i === MAX_ITERS - 1) console.warn("Cv(Fp) iteration did not converge.");
    }
    
    if (Cv_initial == null || Cv_initial <= 0) {
        console.warn("Computed initial turbulent Cv is invalid.");
        return null;
    }

    // --- Step 3: Reynolds number with initial Cv ---
    // We must pass FLP and valveSize to the Reynolds number calc
    let Re = calculateValveReynoldsRev(Cv_initial, FLP, valveSize);
    debugLog("Reynolds Number(Initial CV):", { Re });
    if (Re == null) return null;

    if (Re > 10000) {
        return { Cv: Cv_initial, Re, iterations: 0, choked: (actualDP >= (Math.pow(FLP, 2) * (inP - (FF * UserInputs.getProperty("getSelectedFluid.Pv"))))) };
    }

    // --- Step 4: Iterative correction for low Reynolds number ---
    let Cv_corrected = Cv_initial;

    for (let i = 0; i < MAX_ITERS; i++) {
        Re = calculateValveReynoldsRev(Cv_corrected, FLP, valveSize);
        if (Re == null) return null;

        // Note: verifyFR and FR calcs must also use valveSize
        const fullRange = verifyFR(Cv_corrected, valveSize);
        const FR = fullRange
            ? calculateFRifYes(Cv_corrected, Re, valveSize, FLP)
            : calculateFRifNo(Cv_corrected, Re, valveSize, FLP);

        if (FR == null || FR <= 0) return null;

        const new_Cv = Cv_initial / FR; // Correct the *initial* turbulent Cv

        if(Math.abs(new_Cv - Cv_corrected) / Cv_corrected < EPS) {
            const Re_final = calculateValveReynoldsRev(new_Cv, FLP, valveSize);
            return { Cv: new_Cv, Re: Re_final, iterations: i + 1, choked: (actualDP >= (Math.pow(FLP, 2) * (inP - (FF * UserInputs.getProperty("getSelectedFluid.Pv"))))) };
        }
        Cv_corrected = new_Cv;
    }

    console.warn("Viscosity correction iteration did not converge.");
    const Re_final = calculateValveReynoldsRev(Cv_corrected, FLP, valveSize);
    return { Cv: Cv_corrected, Re: Re_final, iterations: MAX_ITERS, choked: (actualDP >= (Math.pow(FLP, 2) * (inP - (FF * UserInputs.getProperty("getSelectedFluid.Pv"))))) };
}

//------------------------------------------------------------------------------------
// SECTION 7: CALCULATION AND RECOMMENDATION ORCHESTRATOR
//------------------------------------------------------------------------------------

/**
 * Orchestrates the full 2-stage calculation based on new requirements:
 * 1. Gets an initial Cv with generic valve data.
 * 2. Finds matching valves (by type) and all sizes.
 * 3. Validates each and ranks them by "Fit" (Ideal, Acceptable).
 * @returns {Array} An array of recommendation objects.
 */
/**
 * Orchestrates the full 2-stage calculation.
 * @returns {object} An object: { recommendations, initialResult }
 */
function calculateAndRecommend() 
{
    // --- STEP 1: Initial Sizing (using generic Fl/Fd) ---
    const genericFl = UserInputs.getProperty("getSelectedValve.Fl");
    const genericFd = UserInputs.getProperty("getSelectedValve.Fd");
    const pipeSize  = UserInputs.getProperty("getInputs.pipeDiameter");
    const initialResult = valveSizing(pipeSize);
    
    // Handle initial calculation failure
    if (initialResult == null) {
        return { recommendations: [], initialResult: null };
    }
    
    const initialCv = initialResult.Cv;
    const inputs = UserInputs.getInputs();

    // --- STEP 2: Filter Database (by Type Only) ---
    const candidates = valveDatabase.filter(valve => 
        valve.type === inputs.valveType
    );

    if (candidates.length === 0) {
        // No candidates of this type, but we still have the initial Cv
        return { recommendations: [], initialResult: initialResult };
    }

    // --- STEP 3: Validate Each Candidate ---
    let recommendations = [];
    
    for (const candidate of candidates) {
        // ... (your existing loop to validate candidates is perfect) ...
        UserInputs.setValveProperty("Fl", candidate.Fl);
        UserInputs.setValveProperty("Fd", candidate.Fd);

        const validatedResult = valveSizing(candidate.size);
        
        if (validatedResult == null) continue; 
        const finalCv = validatedResult.Cv;
        
        if (candidate.ratedCv < finalCv) continue; 
        
        const opening = (finalCv / candidate.ratedCv) * 100;

        let fit;
        if (opening >= 30 && opening <= 70) {
            fit = "Ideal";
        } else if (opening >= 20 && opening <= 85) {
            fit = "Acceptable";
        } else {
            fit = "Unsuitable"; 
        }

        recommendations.push({
            name: `${candidate.manufacturer} ${candidate.model}`,
            size: candidate.size,
            validatedCv: finalCv,
            regime: validatedResult.choked ? 'Choked' : 'Subcritical',
            openingPercent: opening,
            fit: fit,
            ratedCv: candidate.ratedCv
        });
    }

    // Restore generic Fl/Fd
    UserInputs.setValveProperty("Fl", genericFl);
    UserInputs.setValveProperty("Fd", genericFd);

    // --- STEP 4: Sort Recommendations ---
    // ... (your existing sorting logic is perfect) ...
    const getFitScore = (fit) => {
        if (fit === "Ideal") return 1;
        if (fit === "Acceptable") return 2;
        return 3;
    };
    recommendations.sort((a, b) => {
        const fitScoreA = getFitScore(a.fit);
        const fitScoreB = getFitScore(b.fit);
        if (fitScoreA !== fitScoreB) {
            return fitScoreA - fitScoreB; 
        }
        // If fit is the same, sort by closest to 70% open
        const proximityA = Math.abs(a.openingPercent - 70);
        const proximityB = Math.abs(b.openingPercent - 70);
        return proximityA - proximityB;
    });

    // Return all data
    return { recommendations, initialResult };
}

//------------------------------------------------------------------------------------
// SECTION 8: UI INTEGRATION & EVENT LISTENERS
//------------------------------------------------------------------------------------
/**
 * Reads all input values from the HTML document (DOM) and updates the application state.
 * It also triggers necessary fluid property calculations (e.g., density, viscosity).
 */
function parameterUpdate() {
  // --- 1. Read DOM Elements (Values) ---
  const inletPressureEl = document.getElementById("inletPressure");
  const outletPressureEl = document.getElementById("outletPressure");
  const temperatureEl = document.getElementById("temperature");
  const flowRateEl = document.getElementById("flowRate");
  const pipeDiameterEl = document.getElementById("pipeDiameter");
  const fluidTypeEl = document.getElementById("fluidType");
  const valveTypeEl = document.getElementById("valveType");

  // --- 2. Read DOM Elements (Units) ---
  const flowRateUnit = document.getElementById("flowRateUnit").value;
  const inletPressureUnit = document.getElementById("inletPressureUnit").value;
  const outletPressureUnit = document.getElementById("outletPressureUnit").value;
  const pipeDiameterUnit = document.getElementById("pipeDiameterUnit").value;
  const temperatureUnit = document.getElementById("temperatureUnit").value;

  // --- 3. Parse and Convert Main Inputs ---

  // Get Inlet Pressure (convert to psig)
  let inletPressure = parseFloat(inletPressureEl.value) || 0;
  if (inletPressureUnit === "kpa_g") {
    inletPressure *= conversionConstants.KPA_TO_PSI; // kPa(g) -> psig
  }
  // Convert from (psig) to (PSIa)
  inletPressure += P_ATMOSPHERIC_PSI;

  // Get Outlet Pressure (convert to psig)
  let outletPressure = parseFloat(outletPressureEl.value) || 0;
  if (outletPressureUnit === "kpa_g") {
    outletPressure *= conversionConstants.KPA_TO_PSI; // kPa(g) -> psig
  }
  // Convert from (psig) to (PSIa)
  outletPressure += P_ATMOSPHERIC_PSI;

  // Get Temperature (convert to °C)
  let temperature = parseFloat(temperatureEl.value) || 0;
  if (temperatureUnit === "f") {
    temperature = (temperature - 32) * 5 / 9; // F -> C
  } else if (temperatureUnit === "k") {
    temperature = temperature - 273.15; // K -> C
  }

  // Get Flow Rate (convert to GPM)
  let flowRate = parseFloat(flowRateEl.value) || 0;
  if (flowRateUnit === "m3h") {
    flowRate *= conversionConstants.M3H_TO_GPM; // m³/h -> GPM
  }

  // Get Pipe Diameter (convert to Inches)
  let pipeDiameter = parseFloat(pipeDiameterEl.value) || 0;
  if (pipeDiameterUnit === "mm") {
    pipeDiameter *= conversionConstants.MM_TO_IN; // mm -> in
  }

  // Get types
  const fluidType =
    (fluidTypeEl
      ? parseInt(fluidTypeEl.value, 10)
      : UserInputs.Fluids.WATER) || UserInputs.Fluids.WATER;
  const valveType = (valveTypeEl ? parseInt(valveTypeEl.value, 10) : 0) || 0;

  // --- 4. Update State (Main Inputs) ---
  UserInputs.setInput("inletPressure", inletPressure); // (PSIa)
  UserInputs.setInput("outletPressure", outletPressure); // (PSIa)
  UserInputs.setInput("temperature", temperature); // (C°)
  UserInputs.setInput("flowRate", flowRate); // (GPM)
  UserInputs.setInput("pipeDiameter", pipeDiameter); // (Inches)
  UserInputs.setInput("fluidType", fluidType);
  UserInputs.setInput("valveType", valveType);

  debugLog("Inputs (Converted to Standard Units):", {
    inletPressure: inletPressure.toFixed(2),
    outletPressure: outletPressure.toFixed(2),
    temperature: temperature.toFixed(2),
    flowRate: flowRate.toFixed(2),
    pipeDiameter: pipeDiameter.toFixed(2),
    fluidType,
    valveType,
  });

  // --- 5. Handle Fluid-Specific Properties ---
  if (fluidType === UserInputs.Fluids.WATER) {
    // Estimate the water vapor pressure according to the process temperature
    const waterVaporPressure = estimateWaterVaporPressure(temperature);
    UserInputs.setFluidProperty("Pv", waterVaporPressure);

    // Calculate the water specific gravity
    calculateSpecificGravity();

    // Estimate the water kinematic viscosity at process temperature
    const estViscosity = estimateWaterKinematicViscosity_cSt(temperature);
    UserInputs.setFluidProperty("viscosity", estViscosity);

    debugLog("Fluid - Water (Calculated):", {
      Pv: waterVaporPressure.toFixed(3),
      SG: UserInputs.getProperty("getSelectedFluid.relativeDensity").toFixed(3),
      visc_cSt: estViscosity.toFixed(3),
    });

  } else if (fluidType === UserInputs.Fluids.CUSTOM) {
    // Get custom fluid elements
    const PcEl = document.getElementById("Pc");
    const PvEl = document.getElementById("Pv");
    const rhoEl = document.getElementById("density");
    const viscEl = document.getElementById("viscosity");

    // Get custom fluid units
    const pcUnit = document.getElementById("pcUnit").value;
    const pvUnit = document.getElementById("pvUnit").value;
    const densityUnit = document.getElementById("densityUnit").value;
    const viscosityUnit = document.getElementById("viscosityUnit").value;

    // Parse and convert custom inputs
    let Pc = PcEl ? parseFloat(PcEl.value) || 0 : 0;
    if (pcUnit === "kpa_a") {
      Pc *= conversionConstants.KPA_TO_PSI; // kPa(a) -> psia
    }

    let Pv = PvEl ? parseFloat(PvEl.value) || 0 : 0;
    if (pvUnit === "kpa_a") {
      Pv *= conversionConstants.KPA_TO_PSI; // kPa(a) -> psia
    }

    let visc = viscEl ? parseFloat(viscEl.value) || 0 : 0;
    if (viscosityUnit === "m2s") {
      visc *= conversionConstants.M2S_TO_CST; // m²/s -> cSt
    }

    // Handle Density/SG
    let rho_input = rhoEl ? parseFloat(rhoEl.value) || 0 : 0;
    let specificGravity = 0;

    if (densityUnit === "sg") {
      specificGravity = rho_input;
    } else {
      // 'kgm3'
      if (rho_input > 0) {
        specificGravity = rho_input / D_REFERENCE_WATER;
      }
    }
    // Safety check
    if (specificGravity <= 0) {
      debugLog(
        "FAIL: Invalid custom density/SG <= 0. Defaulting to 1.0 SG."
      );
      specificGravity = 1.0;
    }

    // Update custom fluid state
    UserInputs.setFluidProperty("Pc", Pc);
    UserInputs.setFluidProperty("Pv", Pv);
    UserInputs.setFluidProperty("viscosity", visc);
    UserInputs.setFluidProperty("relativeDensity", specificGravity);

    // Print the calculated custom fluid inputs for debugging purposes
    debugLog("Fluid - Custom (Converted to Standard Units):", {
      Pc: Pc.toFixed(2),
      Pv: Pv.toFixed(2),
      specificGravity: specificGravity.toFixed(3),
      visc: visc.toFixed(2),
    });
  }

  // Update the console log
  console.log(
    "Inputs Updated:",
    UserInputs.getInputs(),
    "Selected fluid:",
    UserInputs.getSelectedFluid()
  );
}

/**
 * @description Main execution block. Attaches an event listener to the "Calculate" button.
 */
const calcBtn = document.getElementById("calculateBtn");
if (calcBtn) {
  calcBtn.addEventListener("click", () => {

    // Show the debug log container on the first click
    if (DEBUG) {
      const debugContainer = document.getElementById("debugContainer");
      if (debugContainer) {
        debugContainer.style.display = "block";
      }
    }

    // Load all inputs
    parameterUpdate(); 

    const { recommendations, initialResult } = calculateAndRecommend();
    const resultEl = document.getElementById("result");

    // --- 1. Handle Initial Calculation ---
    if (initialResult == null) {
      if (resultEl)
        resultEl.innerHTML = `<div class="result-summary"><strong>Calculation failed.</strong><br>Check inputs and console for details.</div>`;
      return;
    }

    const initialCv = initialResult.Cv.toFixed(2);
    const initialRegime = initialResult.choked ? "Choked" : "Subcritical";

    // NEW: Build result HTML with new structure
    let resultHTML = `
      <div class="result-summary">
        <h3>Initial Required Cv: ${initialCv}</h3>
        <p>Flow Regime: <strong>${initialRegime}</strong></p>
      </div>
      <h3>Valve Recommendations</h3>
    `;

    // --- 2. Find the ONE best for each category ---
    const bestIdeal = recommendations.find((v) => v.fit === "Ideal");
    const bestAcceptable = recommendations.find(
      (v) => v.fit === "Acceptable"
    );

    if (!bestIdeal && !bestAcceptable) {
      resultHTML +=
        "<p>No suitable valves found in the desired 20-80% operating range.</p>";
    } else {
      // --- 3. Display the Best Ideal (if it exists) ---
      if (bestIdeal) {
        resultHTML += `
          <div class="recommendation-card ideal-fit">
            <span class="card-title">Best Fit (Ideal): ${
              bestIdeal.name
            } (Size: ${bestIdeal.size}")</span>
            Max Cv: ${bestIdeal.ratedCv.toFixed(
              0
            )} | Required Cv: <strong>${bestIdeal.validatedCv.toFixed(
          2
        )}</strong><br>
            Operating Point: <strong>${bestIdeal.openingPercent.toFixed(
              1
            )}% Open</strong><br>
            Flow: ${bestIdeal.regime}
          </div>
        `;
      }

      // --- 4. Display the Best Acceptable (if it exists) ---
      if (bestAcceptable) {
        resultHTML += `
          <div class="recommendation-card acceptable-fit">
            <span class="card-title">Best Fit (Acceptable): ${
              bestAcceptable.name
            } (Size: ${bestAcceptable.size}")</span>
            Max Cv: ${bestAcceptable.ratedCv.toFixed(
              0
            )} | Required Cv: <strong>${bestAcceptable.validatedCv.toFixed(
          2
        )}</strong><br>
            Operating Point: <strong>${bestAcceptable.openingPercent.toFixed(
              1
            )}% Open</strong><br>
            Flow: ${bestAcceptable.regime}
          </div>
        `;
      }
    }

    if (resultEl) resultEl.innerHTML = resultHTML;
  });
}

//------------------------------------------------------------------------------------
// SECTION 9: AUXILIARY UI INTEGRATION & EVENT LISTENERS
//------------------------------------------------------------------------------------
const logEl = document.getElementById('log');
function logToConsole(msg) {
    if (!logEl) return;
    const ts = new Date().toISOString();
    logEl.textContent += `[${ts}] ${msg}\n`;
    logEl.scrollTop = logEl.scrollHeight;
}


document.addEventListener("DOMContentLoaded", () => {
    const fluidTypeEl = document.getElementById("fluidType");
    const customFluidDiv = document.getElementById("customFluidInputs");

    if(fluidTypeEl && customFluidDiv) 
    {
        fluidTypeEl.addEventListener("change", (e) => {
        const selected = parseInt(e.target.value, 10);
        customFluidDiv.style.display = selected === 1 ? "block" : "none";
        });
    }
});

document.getElementById("resetBtn").addEventListener("click", () => {
  // Reset values
  document.getElementById("flowRate").value = 530;
  document.getElementById("inletPressure").value = 85;
  document.getElementById("outletPressure").value = 30;
  document.getElementById("pipeDiameter").value = 4;
  document.getElementById("temperature").value = 77;
  document.getElementById("fluidType").value = 0;
  document.getElementById("valveType").value = 0;

  // Reset units
  document.getElementById("flowRateUnit").value = "gpm";
  document.getElementById("inletPressureUnit").value = "psig";
  document.getElementById("outletPressureUnit").value = "psig";
  document.getElementById("pipeDiameterUnit").value = "in";
  document.getElementById("temperatureUnit").value = "f";

  // Reset custom fluid inputs
  document.getElementById("Pc").value = "";
  document.getElementById("Pv").value = "";
  document.getElementById("density").value = "";
  document.getElementById("viscosity").value = "";
  document.getElementById("pcUnit").value = "psia";
  document.getElementById("pvUnit").value = "psia";
  document.getElementById("densityUnit").value = "kgm3";
  document.getElementById("viscosityUnit").value = "cst";
  
  // Trigger UI update for custom panel
  document.getElementById("customFluidInputs").style.display = "none";

  // 4. Hide the Debug Log
  const debugContainer = document.getElementById("debugContainer");

  if (debugContainer) {
      debugContainer.style.display = "none"; // Hide the whole container
  }

  document.getElementById("result").innerHTML = "Reset to defaults.";
  logToConsole("Inputs reset to defaults.");
});

/* Quick tests */
/* Quick tests */
function runPreset(preset) {
  // Ensure default Imperial/US units are selected for these tests
  document.getElementById("flowRateUnit").value = "gpm";
  document.getElementById("inletPressureUnit").value = "psig";
  document.getElementById("outletPressureUnit").value = "psig";
  document.getElementById("pipeDiameterUnit").value = "in";
  document.getElementById("temperatureUnit").value = "f";
  document.getElementById("fluidType").value = 0; // Water

  if (preset === 1) {
    document.getElementById("flowRate").value = 5;
    document.getElementById("inletPressure").value = 50;
    document.getElementById("outletPressure").value = 49; // dp = 1 psi
    document.getElementById("pipeDiameter").value = 0.5;
    document.getElementById("temperature").value = 68; // 20C
  } else if (preset === 2) {
    document.getElementById("flowRate").value = 10;
    document.getElementById("inletPressure").value = 100;
    document.getElementById("outletPressure").value = 90; // dp = 10 psi
    document.getElementById("pipeDiameter").value = 1;
    document.getElementById("temperature").value = 68; // 20C
  } else {
    document.getElementById("flowRate").value = 50;
    document.getElementById("inletPressure").value = 200;
    document.getElementById("outletPressure").value = 100; // dp = 100 psi
    document.getElementById("pipeDiameter").value = 2;
    document.getElementById("temperature").value = 68; // 20C
  }
  // parameterUpdate(); // This is called by the click()
  document.getElementById("calculateBtn").click();
}

document.getElementById('test1').addEventListener('click', () => runPreset(1));
document.getElementById('test2').addEventListener('click', () => runPreset(2));
document.getElementById('test3').addEventListener('click', () => runPreset(3));
