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
 * @description Standard numerical constants for Imperial units.
 */
const P_ATMOSPHERIC_PSI = 14.7;   // Atmospheric pressure in psi
const D_REFERENCE_WATER = 1000.0; // Reference density in kg/m^3

/**
 * @const {object} conversionConstants
 * @description Conversion factors to standardized internal units.
 */
const conversionConstants = {
    KPA_TO_PSI: 0.145038,
    BAR_TO_PSI: 14.5038,
    M3H_TO_GPM: 4.40287,
    MM_TO_IN:   0.0393701,
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
        // NEW: Advanced Sizing Inputs
        flowRateMin: 0,
        inletPressureMin: 0,
        outletPressureMin: 0,
        flowRateMax: 0,
        inletPressureMax: 0,
        outletPressureMax: 0,
        // NEW: User Valve Inputs
        customValveSize: 4,
        customFl: 0.9,
        customFd: 0.46
    };

    const Fluids = Object.freeze({ WATER: 0, CUSTOM: 1 });
    const Valves = Object.freeze({ GLOBEOPEN: 0, GLOBECLOSE: 1, BUTTERFLY: 2, BALL: 3, CUSTOM: 99 });

    const fluidData = {
        [Fluids.WATER]: { ...liquidWater, isSelected: true },
        [Fluids.CUSTOM]: { ...liquidCustom, isSelected: false },
    };

    const valveData = {
        [Valves.GLOBEOPEN]: { ...valveGlobesSinglePortContouredPlugOpen },
        [Valves.GLOBECLOSE]: { ...valveGlobesSinglePortContouredPlugClose },
        [Valves.BUTTERFLY]: { ...valveButterfly },
        [Valves.BALL]: { ...valveBallFull },
        [Valves.CUSTOM]: { Fl: 0.9, Fd: 0.46, Cv: 0, Xt: 0.72 }
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
    if (mu_cP < 0.1) mu_cP = 0.1;
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
 *   FF is an empirical factor used in the liquid sizing equations to account
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
 *   Uses the choked-liquid criterion:
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

    // Compute the choked threshold (ΔP_crit)
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
 *   For choked liquid flow gives:
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

    // Standard subcritical Cv expression
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
 *   Rev is an empirical valve-specific Reynolds-like number defined to decide
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

    // Pick the viscosity from fluid data (expected as cP)
    let v = UserInputs.getProperty("getSelectedFluid.viscosity");

    // Validate geometry & fluid inputs
    if (!d || d <= 0) {
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
    const inside = (Math.pow(FL, 2) * Math.pow(Ci, 2)) / (numericalConstants.N2 * Math.pow(d, 4)) + 1;
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
 * @param {number} C - The calculated flow coefficient (Cv).
 * @returns {boolean} `true` to use one set of formulas, `false` for the other.
 */
/**
 * Decide which FR formula branch to use (full-range vs restricted-range).
 *
 * Purpose:
 *   Uses a boundary (Cv / D^2 comparison) to choose between two FR families.
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
// SECTION 5.5: PRESSURE-TEMPERATURE RATINGS (ASME B16.34 - Group 1.1 WCB)
//------------------------------------------------------------------------------------
/**
 * Pressure-Temperature Ratings for ASTM A216 WCB (Standard Class).
 * Source: ASME B16.34-2017, Table 2-1.1
 * Units: Temp in °F, Pressure in psig
 */
const WCB_Ratings = {
    // Temperature points corresponding to the pressure arrays below
    tempPoints: [-20, 100, 200, 300, 400, 500, 600, 650, 700, 750, 800],
    
    classes: {
        150:  [285, 285, 260, 230, 200, 170, 140, 125, 110, 95, 80],
        300:  [740, 740, 675, 655, 635, 600, 550, 535, 535, 505, 410],
        600:  [1480, 1480, 1350, 1315, 1270, 1200, 1095, 1075, 1065, 1010, 825],
        
        // High Pressure Classes
        900:  [2220, 2220, 2025, 1970, 1900, 1795, 1640, 1610, 1600, 1510, 1235],
        1500: [3705, 3705, 3375, 3280, 3170, 2995, 2735, 2685, 2665, 2520, 2060],
        2500: [6170, 6170, 5625, 5470, 5280, 4990, 4560, 4475, 4440, 4200, 3430],
        4500: [11110, 11110, 10120, 9845, 9505, 8980, 8210, 8055, 7990, 7560, 6170]
    }
};

/**
 * Calculates the Maximum Allowable Working Pressure (MAWP) for WCB steel
 * at a specific temperature and pressure class using linear interpolation.
 * @param {number} tempF - Temperature in Fahrenheit
 * @param {number} classRating - ANSI Class (150, 300, 600...)
 * @returns {number} MAWP in psig
 */
function getMaxPressureWCB(tempF, classRating) {
    if (!WCB_Ratings.classes[classRating]) return 99999; // Unknown class, assume safe

    // Clamp temp to data limits
    const temps = WCB_Ratings.tempPoints;
    const pressures = WCB_Ratings.classes[classRating];
    
    // Case 1: Below min temp
    if (tempF <= temps[0]) return pressures[0];

    // Case 2: Above max temp (Use last known - unsafe to extrapolate up)
    if (tempF >= temps[temps.length - 1]) return pressures[pressures.length - 1];

    // Case 3: Interpolate
    for (let i = 0; i < temps.length - 1; i++) {
        if (tempF >= temps[i] && tempF <= temps[i+1]) {
            const t1 = temps[i];
            const t2 = temps[i+1];
            const p1 = pressures[i];
            const p2 = pressures[i+1];
            
            // Linear Interpolation: P = P1 + (T - T1) * (P2 - P1) / (T2 - T1)
            return p1 + (tempF - t1) * (p2 - p1) / (t2 - t1);
        }
    }
    return pressures[0];
}

/**
 * Validates if the valve class is sufficient for the process pressure
 * considering a safety factor (derating).
 */
function checkPressureClass(valve, maxInletPressure, tempF) {
    // If we're passing a specific class (number) or object with pressureClass property
    const rating = typeof valve === 'number' ? valve : valve.pressureClass;
    
    if (!rating) return { pass: true, limit: 0, rated: 0 }; // Pass if no class defined (legacy data)

    // 1. Get Base ASME Rating
    const asmeLimit = getMaxPressureWCB(tempF, rating);

    // 2. Apply 75% Safety Factor
    const safeLimit = asmeLimit * 0.75;

    // 3. Compare
    const pass = maxInletPressure <= safeLimit;

    return { pass, limit: safeLimit, rated: asmeLimit };
}

/**
 * SENAI-Derived Flow Characteristic Recommender
 * Strictly follows the logic table provided in the prompt.
 * * @param {object} process - { qMin, qMax, dpMin, dpMax }
 * @param {boolean} isChoked - Choked status from sizing engine
 * @returns {object} { type: string, reason: string }
 */
function recommendSenaiCharacteristic(process, isChoked) {
    
    // --- 1. Validation ---
    if (!process.qMin || !process.qMax || !process.dpMin || !process.dpMax || process.dpMin <= 0) {
        return {
            type: "Equal Percentage",
            reason: "Insufficient process data. Defaulting to Equal Percentage (Safety)."
        };
    }

    // --- 2. Calculate Indicators ---
    const dpRatio = process.dpMax / process.dpMin;
    const flowRatio = process.qMax / process.qMin;

    // --- 3. PRIORITY 1: CRITICAL OVERRIDE ---
    if (isChoked) {
        return {
            type: "Equal Percentage",
            reason: "Critical Override: Choked flow detected. Equal Percentage selected for stability."
        };
    }

    // --- 4. PRIORITY 2: STRONGLY INCREASING DP (> 2.0) ---
    if (dpRatio > 2.0) {
        if (flowRatio > 3.0) {
            return {
                type: "Quick Opening",
                reason: `Strongly increasing ΔP (Ratio: ${dpRatio.toFixed(2)}) with large flow range. Quick Opening matches gain.`
            };
        } else {
            return {
                type: "Equal Percentage",
                reason: `Strongly increasing ΔP (Ratio: ${dpRatio.toFixed(2)}) but small flow range. Defaulting to Eq%.`
            };
        }
    }

    // --- 5. PRIORITY 3: APPROXIMATELY CONSTANT DP (0.8 to 1.25) ---
    if (dpRatio >= 0.8 && dpRatio <= 1.25) {
        return {
            type: "Linear",
            reason: `Approximately constant ΔP (Ratio: ${dpRatio.toFixed(2)}). Linear trim matches constant system gain.`
        };
    }

    // --- 6. PRIORITY 4: VARIABLE DP (Decreasing < 0.8 OR Increasing 1.25-2.0) ---
    let trend = dpRatio < 0.8 ? "Decreasing" : "Increasing";
    return {
        type: "Equal Percentage",
        reason: `Variable ΔP (${trend}, Ratio: ${dpRatio.toFixed(2)}). Equal Percentage compensates for gain changes.`
    };
}

/**
 * Calculates Hydrodynamic Noise using the SENAI/Masoneilan Method.
 * Formula: SPL = SPL_dP + SPL_delta - SPL_C
 * * @param {number} Cv - Valve Cv
 * @param {number} P1 - Inlet Pressure (psia)
 * @param {number} P2 - Outlet Pressure (psia)
 * @param {number} Pv - Vapor Pressure (psia)
 * @param {number} FL - Liquid Pressure Recovery Factor
 * @param {number} d_pipe - Pipe Diameter (inches)
 * @returns {number} Noise in dBA
 */
function calculateNoiseSenai(Cv, P1, P2, Pv, FL, d_pipe) {
    // 0. Safety Checks
    if (Cv <= 0 || P1 <= P2) return 0;
    
    // 1. Define Variables
    const dP = P1 - P2;
    const P1_Pv = P1 - Pv;
    const P2_Pv = P2 - Pv;
    
    // Critical pressure drop where cavitation starts (approximate)
    // In this method, dP_cav is roughly related to FL^2
    const dP_choked = FL * FL * (P1 - Pv);
    
    // Kc (Cavitation Index for this method)
    // Low Kc means high cavitation risk. High Kc means safe.
    // Kc = (P1 - Pv) / dP  <-- Standard ISA definition
    // Note: The prompt implies a comparison. Let's calculate the ratio xF.
    const xF = dP / (P1 - Pv); // Operating pressure ratio
    const xF_limit = 0.9 * (FL * FL); // The limit mentioned in the prompt (0.9 FL^2)

    // --- 2. Calculate SPL_delta (Fig. 56) ---
    // Represents Base Energy function of (P2 - Pv) and Cv
    // Approx: 10*log(Cv) + 20*log(P2-Pv) + Constant
    // We assume P2-Pv must be > 0.1 to avoid log(0) errors
    const safe_P2_Pv = Math.max(0.1, P2_Pv);
    let SPL_delta = 10 * Math.log10(Cv) + 20 * Math.log10(safe_P2_Pv) + 38;

    // --- 3. Calculate SPL_dP (Fig. 55) ---
    // Represents Cavitation Intensity. 
    // It rises as dP approaches dP_choked.
    // We model this curve as a ratio of dP/dP_choked.
    let ratio = dP / dP_choked;
    let SPL_dP = 0;
    
    if (ratio < 0.5) {
        // Linear rise in low turbulence
        SPL_dP = 10 * ratio; 
    } else if (ratio < 1.0) {
        // Exponential rise as we hit cavitation
        SPL_dP = 10 + 40 * Math.pow((ratio - 0.5) * 2, 2);
    } else {
        // Choked / Full Cavitation (Saturation plateau)
        SPL_dP = 50; 
    }

    // --- 4. Calculate SPL_C (Fig. 57) ---
    // "O valor de SPLC pode ser desprezado caso KC seja superior a 0,9 FL^2"
    // Note: Inverted logic for xF. 
    // xF is (dP / P1-Pv). 0.9 FL^2 is the limit. 
    // If xF < 0.9 FL^2, we are SAFE (No Cavitation). 
    // If xF > 0.9 FL^2, we are in CAVITATION.
    
    let SPL_C = 0;
    
    // If we are in the "Safe Zone" (dP is small compared to limit), apply attenuation
    if (xF < xF_limit) {
         // Pipe wall attenuation estimate (Schedule 40)
         // Thicker pipes reduce noise.
         // Approx: 10dB reduction per inch of thickness (simplified)
         // We use a simplified curve: 30dB baseline - size factor
         SPL_C = 30 - (10 * Math.log10(d_pipe));
    } else {
         // In cavitation (High Noise), SPL_C is ignored (0), meaning full noise transmission
         SPL_C = 0;
    }

    // --- 5. Final Sum ---
    // SPL = SPL_dP + SPL_delta - SPL_C
    let SPL = SPL_dP + SPL_delta - SPL_C;

    // Clamp to realistic bounds
    if (SPL < 30) SPL = 30;
    if (SPL > 130) SPL = 130;

    return parseFloat(SPL.toFixed(1));
}

/**
 * Calculates ISA-RP75.23 Cavitation Index (Sigma)
 * and determines consequences based on specific user thresholds.
 * Formula: Sigma = (P1 - Pv) / (P1 - P2)
 */
function calculateCavitationSigma(P1, P2, Pv) {
    const dP = P1 - P2;
    // Safety check for invalid dP
    if (dP <= 0) return { sigma: 999, severity: "Safe", consequence: "No Risk", color: "#34d399" };

    const sigma = (P1 - Pv) / dP;
    const sigmaVal = parseFloat(sigma.toFixed(2));

    let result = {
        sigma: sigmaVal,
        severity: "Safe",
        consequence: "No Risk of Cavitation",
        color: "#34d399" // Green
    };

    if (sigmaVal <= 1.0) {
        result.severity = "Flashing";
        result.consequence = "Flashing is occurring.";
        result.color = "#7f1d1d"; // Dark Red / Purple
    } 
    else if (sigmaVal <= 1.5) {
        result.severity = "Severe";
        result.consequence = "Potential for severe cavitation.";
        result.color = "#ef4444"; // Red
    } 
    else if (sigmaVal <= 1.7) {
        result.severity = "Moderate";
        result.consequence = "Some cavitation control required.";
        result.color = "#f97316"; // Orange
    } 
    else if (sigmaVal <= 2.0) {
        result.severity = "Incipient";
        result.consequence = "No cavitation control required (Hardened trim recommended).";
        result.color = "#eab308"; // Yellow
    }

    return result;
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
function valveSizing(valveSize, overrideQ, overrideP1, overrideP2) 
{
    const EPS = 0.0001;      // Convergence tolerance
    const MAX_ITERS = 20;   // Max iterations for Cv and FR

    const pipeSize = UserInputs.getProperty("getInputs.pipeDiameter");
    const FL = UserInputs.getProperty("getSelectedValve.Fl"); // Generic or specific FL
    const SG = UserInputs.getProperty("getSelectedFluid.relativeDensity");
    const N1 = numericalConstants.N1;

    // Use overrides if provided, otherwise standard inputs
    const Q = (overrideQ !== undefined) ? overrideQ : UserInputs.getProperty("getInputs.flowRate");
    const inP = (overrideP1 !== undefined) ? overrideP1 : UserInputs.getProperty("getInputs.inletPressure");
    const outP = (overrideP2 !== undefined) ? overrideP2 : UserInputs.getProperty("getInputs.outletPressure");

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

    if (Cv_initial >= 10000) {
        console.warn("Pipe undersized.");
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
 * Calculates Valve Opening % based on Required Cv using Piecewise Linear Interpolation.
 */
function calculateOpeningFromCv(requiredCv, cvCurve) {
    if (requiredCv <= 0) return 0;
    const maxCv = cvCurve[cvCurve.length - 1];
    
    // Extrapolate if over max (to show >100%)
    if (requiredCv >= maxCv) {
        return 100 * (requiredCv / maxCv);
    }

    const fullCvCurve = [0, ...cvCurve]; 
    const fullOpeningCurve = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];

    for (let i = 0; i < fullCvCurve.length - 1; i++) {
        const cvLow = fullCvCurve[i];
        const cvHigh = fullCvCurve[i+1];

        if (requiredCv >= cvLow && requiredCv <= cvHigh) {
            const percentLow = fullOpeningCurve[i];
            const percentHigh = fullOpeningCurve[i+1];
            const fraction = (requiredCv - cvLow) / (cvHigh - cvLow);
            return percentLow + (fraction * (percentHigh - percentLow));
        }
    }
    return 100;
}

/**
 * Orchestrates the full 3-point calculation (Min, Normal, Max)
 * and determines process requirements.
 */
function calculateAndRecommend() 
{
    const inputs = UserInputs.getInputs();
    const SG = UserInputs.getProperty("getSelectedFluid.relativeDensity");
    const Pv = UserInputs.getProperty("getSelectedFluid.Pv"); 
    const N1 = numericalConstants.N1;

    // --- 0. SAFETY CHECK: INLET FLASHING ---
    if (inputs.inletPressure <= Pv) {
        return { 
            recommendations: [], 
            processInfo: {
                error: true,
                errorMessage: "CRITICAL: Inlet Pressure is below Vapor Pressure. Fluid is flashing/boiling at the inlet.",
                baseCvOp: 0, baseCvMin: 0, baseCvMax: 0,
                preferredChar: "Error", charReason: "Check Process Conditions",
                minProcessClass: "N/A", sigma: 0, sigmaSeverity: "N/A", sigmaColor: "#ef4444", sigmaLabel: "N/A"
            } 
        };
    }

    // --- 1. REFERENCE CALCULATION ---
    const calcSimpleCv = (Q, P1, P2) => {
        const dP = P1 - P2;
        if (dP <= 0 || Q <= 0) return 0;
        return (Q / N1) * Math.sqrt(SG / dP);
    };

    const baseCvOp = calcSimpleCv(inputs.flowRate, inputs.inletPressure, inputs.outletPressure);
    const baseCvMin = calcSimpleCv(inputs.flowRateMin, inputs.inletPressureMin, inputs.outletPressureMin);
    const baseCvMax = calcSimpleCv(inputs.flowRateMax, inputs.inletPressureMax, inputs.outletPressureMax);

    // --- 2. PROCESS ANALYSIS ---
    const processSigmaData = calculateCavitationSigma(inputs.inletPressure, inputs.outletPressure, Pv);

    const standardClasses = [150, 300, 600, 900, 1500, 2500, 4500];
    let minProcessClass = 4500; // Default high
    const tempInC = UserInputs.getProperty("getInputs.temperature");
    const tempF = (tempInC * 9/5) + 32;
    const pMaxProcess = inputs.inletPressureMax;

    for (const cls of standardClasses) {
        const check = checkPressureClass(cls, pMaxProcess, tempF);
        if (check.pass) { minProcessClass = cls; break; }
    }

    // --- 3. SENAI LOGIC ---
    const resBaseMax = valveSizing(inputs.pipeDiameter, inputs.flowRateMax, inputs.inletPressureMax, inputs.outletPressureMax);
    const isChokedSafety = resBaseMax ? resBaseMax.choked : false;

    const processData = {
        qMin: inputs.flowRateMin, qMax: inputs.flowRateMax,
        dpMin: inputs.inletPressureMin - inputs.outletPressureMin, 
        dpMax: inputs.inletPressureMax - inputs.outletPressureMax  
    };
    const senaiRec = recommendSenaiCharacteristic(processData, isChokedSafety);
    let preferredChar = senaiRec.type;
    let charReason = senaiRec.reason;
    
    // --- SCORING HELPER FUNCTION ---
    const calculateScore = (valve, opData) => {
        let score = 0;

        // 1. Operating Point (40 pts) - Target: 60-75%
        const open = opData.openOp;
        let dist = 0;
        if (open < 60) dist = 60 - open;
        else if (open > 75) dist = open - 75;
        
        if (dist === 0) score += 40;
        else {
            // Deduct 0.8 pts per % deviation. 
            // 10% open (dist 50) -> Score 0. 85% open (dist 10) -> Score 32.
            score += Math.max(0, 40 - (dist * 0.8));
        }

        // 2. Characteristic (20 pts)
        if (opData.isPreferredChar) score += 20;
        else score += 10; // Fallback (e.g. Eq% used for Linear)

        // 3. Noise (15 pts) - Target < 75 dBA
        const noise = opData.noise;
        if (noise <= 75) score += 15;
        else if (noise >= 85) score += 0;
        else {
            // Linear drop between 75 and 85
            score += 15 - ((noise - 75) * 1.5);
        }

        // 4. Rangeability (10 pts) - Target > 50:1
        const ratio = opData.rangeability;
        if (ratio >= 50) score += 10;
        else if (ratio <= 10) score += 0;
        else {
            score += ((ratio - 10) / 40) * 10;
        }

        // 5. Class Efficiency (10 pts)
        const classes = [150, 300, 600, 900, 1500, 2500, 4500];
        const reqIdx = classes.indexOf(minProcessClass);
        const valIdx = classes.indexOf(valve.pressureClass);
        
        if (valIdx !== -1 && reqIdx !== -1) {
            const diff = valIdx - reqIdx;
            if (diff === 0) score += 10;
            else if (diff === 1) score += 8;
            else score += 2; // Massive over-rating
        } else {
            score += 10; // Custom/Unknown
        }

        // 6. Size Match (5 pts)
        // Standard pipe sizes to determine steps
        const sizes = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24];
        const pIdx = sizes.indexOf(inputs.pipeDiameter);
        const vIdx = sizes.indexOf(valve.size);
        
        if (pIdx !== -1 && vIdx !== -1) {
            const diff = pIdx - vIdx;
            // Exact match or 1 size smaller is ideal (Max points)
            if (diff === 0 || diff === 1) score += 5;
            else score += 0; // 2+ sizes smaller is penalized
        } else {
            // Fallback for non-standard sizes
            if (valve.size >= inputs.pipeDiameter * 0.6) score += 5;
        }

        return score;
    };

    let recommendations = [];

    // ============================================================
    // BRANCH A: CUSTOM VALVE ANALYSIS
    // ============================================================
    if (inputs.valveType === 99) {
        // [Same logic as before, just added score placeholder]
        const customCandidate = {
            manufacturer: "User Defined", model: "Custom Analysis",
            size: inputs.customValveSize || inputs.pipeDiameter, 
            Fl: inputs.customFl, Fd: inputs.customFd,            
            pressureClass: "N/A", characteristic: "User Defined"
        };

        UserInputs.setValveProperty("Fl", customCandidate.Fl);
        UserInputs.setValveProperty("Fd", customCandidate.Fd);

        const resOp = valveSizing(customCandidate.size, inputs.flowRate, inputs.inletPressure, inputs.outletPressure);
        let resMin = null, resMax = null;
        try {
            resMin = valveSizing(customCandidate.size, inputs.flowRateMin, inputs.inletPressureMin, inputs.outletPressureMin);
            resMax = valveSizing(customCandidate.size, inputs.flowRateMax, inputs.inletPressureMax, inputs.outletPressureMax);
        } catch(e) {}

        if (resOp) {
             const noiseVal = calculateNoiseSenai(resOp.Cv, inputs.inletPressure, inputs.outletPressure, Pv, customCandidate.Fl, customCandidate.size);
            let processTurndown = (resMin && resMax && resMin.Cv > 0) ? (resMax.Cv / resMin.Cv).toFixed(2) : "N/A";
            
            recommendations.push({
                name: "Custom Valve Analysis",
                size: customCandidate.size,
                ratedCv: 0, characteristic: "N/A", pressureClass: "N/A",
                noise: noiseVal,
                sigma: processSigmaData.sigma, cavSeverity: processSigmaData.severity, cavConsequence: processSigmaData.consequence, cavColor: processSigmaData.color,
                cvOp: resOp.Cv, cvMin: resMin ? resMin.Cv : 0, cvMax: resMax ? resMax.Cv : 0,
                openOp: 0, openMin: 0, openMax: 0, 
                fit: "Custom", score: 100, // Custom always 100
                isPreferredChar: false, turndownPass: true,
                valveRangeability: 0, processTurndown: processTurndown,
                regime: resOp.choked ? "Choked" : "Non-Choked"
            });
        }
    }
    // ============================================================
    // BRANCH B: STANDARD DATABASE SEARCH
    // ============================================================
    else {
        const originalFl = UserInputs.getProperty("getSelectedValve.Fl");
        const originalFd = UserInputs.getProperty("getSelectedValve.Fd");
        const candidates = typeof valveDatabase !== 'undefined' ? valveDatabase.filter(valve => valve.type === inputs.valveType) : [];

        for (const originalCandidate of candidates) {
            let candidate = { ...originalCandidate };

            // A. Basic Geometry Filters
            if (candidate.size > inputs.pipeDiameter) continue; 
            if (candidate.ratedCv < (baseCvMax * 0.90)) continue; 

            // B. Flow Characteristic Filtering
            const valveChar = candidate.characteristic || "Linear"; 
            if (preferredChar === "Equal Percentage") {
                if (valveChar !== "Equal Percentage") continue;
            } else {
                if (valveChar !== preferredChar && valveChar !== "Equal Percentage") continue;
            }

            // C. Pressure Class Filter
            let selectedClass = null;
            let availableClasses = [];
            if (candidate.pressureClasses && Array.isArray(candidate.pressureClasses)) availableClasses = candidate.pressureClasses;
            else if (candidate.pressureClass) availableClasses = [candidate.pressureClass];
            
            availableClasses.sort((a, b) => a - b);
            for (const cls of availableClasses) {
                const check = checkPressureClass(cls, inputs.inletPressureMax, tempF);
                if (check.pass) { selectedClass = cls; break; }
            }
            if (selectedClass === null) continue; 
            candidate.pressureClass = selectedClass;

            // D. Run Sizing
            UserInputs.setValveProperty("Fl", candidate.Fl);
            UserInputs.setValveProperty("Fd", candidate.Fd);

            const resOp = valveSizing(candidate.size, inputs.flowRate, inputs.inletPressure, inputs.outletPressure);
            if (!resOp) continue;
            const resMin = valveSizing(candidate.size, inputs.flowRateMin, inputs.inletPressureMin, inputs.outletPressureMin);
            if (!resMin) continue;
            const resMax = valveSizing(candidate.size, inputs.flowRateMax, inputs.inletPressureMax, inputs.outletPressureMax);
            if (!resMax) continue;

            const noiseVal = calculateNoiseSenai(resOp.Cv, inputs.inletPressure, inputs.outletPressure, Pv, candidate.Fl, candidate.size);

            const getOpen = (cv, val) => {
                if (val.cvCurve && val.cvCurve.length === 10) return calculateOpeningFromCv(cv, val.cvCurve);
                return (cv / val.ratedCv) * 100;
            };
            const openOp = getOpen(resOp.Cv, candidate);
            const openMin = getOpen(resMin.Cv, candidate);
            const openMax = getOpen(resMax.Cv, candidate);

            // E. Fit Filtering
            let fit = "Marginal";
            if (openOp > 0 && openMax <= 100) {
                if (openOp >= 50 && openOp <= 80 && openMax < 90 && openMin > 10) fit = "Ideal";
                else if (openOp >= 10 && openOp <= 90 && openMax <= 95) fit = "Acceptable";
            } else {
                fit = "Invalid"; 
            }
            if (fit === "Invalid") continue;

            // F. Process Data
            const processTurndown = resMax.Cv / resMin.Cv;
            let cvAt90 = 0, cvAt10 = 0;
            if (candidate.cvCurve && candidate.cvCurve.length === 10) {
                cvAt10 = candidate.cvCurve[0]; cvAt90 = candidate.cvCurve[8]; 
            } else {
                cvAt10 = candidate.ratedCv * 0.1; cvAt90 = candidate.ratedCv * 0.9;
            }
            const valveRangeability = cvAt90 / cvAt10;
            const turndownPass = (valveRangeability >= processTurndown);
            const isPreferredChar = (candidate.characteristic === preferredChar);

            // G. Calculate Score
            const opData = {
                openOp, noise: noiseVal, 
                rangeability: valveRangeability, 
                isPreferredChar
            };
            const score = calculateScore(candidate, opData);

            recommendations.push({
                name: `${candidate.manufacturer} ${candidate.model}`,
                size: candidate.size, ratedCv: candidate.ratedCv,
                characteristic: candidate.characteristic || "Linear",
                pressureClass: candidate.pressureClass,
                noise: noiseVal,
                sigma: processSigmaData.sigma, cavSeverity: processSigmaData.severity, cavConsequence: processSigmaData.consequence, cavColor: processSigmaData.color,
                cvOp: resOp.Cv, cvMin: resMin.Cv, cvMax: resMax.Cv,
                openOp: openOp, openMin: openMin, openMax: openMax,
                fit: fit, score: Math.round(score), // Save Score
                isPreferredChar: isPreferredChar, turndownPass: turndownPass,
                valveRangeability: valveRangeability, processTurndown: processTurndown.toFixed(2),
                regime: resOp.choked ? "Choked" : "Non-Choked"
            });
        }

        UserInputs.setValveProperty("Fl", originalFl);
        UserInputs.setValveProperty("Fd", originalFd);

        // H. Sort by Score
        recommendations.sort((a, b) => {
            // Safety first: Turndown Pass is mandatory
            if (a.turndownPass && !b.turndownPass) return -1;
            if (!a.turndownPass && b.turndownPass) return 1;
            
            // Then purely by Weighted Score
            return b.score - a.score;
        });
    }

    return { 
        recommendations, 
        processInfo: { 
            preferredChar, charReason, baseCvOp, baseCvMin, baseCvMax, minProcessClass,
            sigma: processSigmaData.sigma, sigmaSeverity: processSigmaData.severity, sigmaColor: processSigmaData.color, sigmaLabel: processSigmaData.consequence
        } 
    };
}

//------------------------------------------------------------------------------------
// SECTION 8: UI INTEGRATION & EVENT LISTENERS
//------------------------------------------------------------------------------------
/**
 * Reads all input values from the HTML document (DOM) and updates the application state.
 * It also triggers necessary fluid property calculations (e.g., density, viscosity).
 */
function parameterUpdate() {
  // --- 1. Read Standard DOM Elements ---
  const inletPressureEl = document.getElementById("inletPressure");
  const outletPressureEl = document.getElementById("outletPressure");
  const temperatureEl = document.getElementById("temperature");
  const flowRateEl = document.getElementById("flowRate");
  const pipeDiameterEl = document.getElementById("pipeDiameter");
  const fluidTypeEl = document.getElementById("fluidType");
  const valveTypeEl = document.getElementById("valveType");

  // Main Units
  const flowRateUnit = document.getElementById("flowRateUnit").value;
  const inletPressureUnit = document.getElementById("inletPressureUnit").value;
  const outletPressureUnit = document.getElementById("outletPressureUnit").value;
  const pipeDiameterUnit = document.getElementById("pipeDiameterUnit").value;
  const temperatureUnit = document.getElementById("temperatureUnit").value;

  // --- 2. Convert Main Inputs ---
  // Pressure: Convert to PSI Gauge, then add Atmospheric to get PSIA
  let inletPressure = parseFloat(inletPressureEl.value) || 0;
  if (inletPressureUnit === "kpa_g") inletPressure *= conversionConstants.KPA_TO_PSI;
  else if (inletPressureUnit === "bar_g") inletPressure *= conversionConstants.BAR_TO_PSI;
  inletPressure += P_ATMOSPHERIC_PSI;

  let outletPressure = parseFloat(outletPressureEl.value) || 0;
  if (outletPressureUnit === "kpa_g") outletPressure *= conversionConstants.KPA_TO_PSI;
  else if (outletPressureUnit === "bar_g") outletPressure *= conversionConstants.BAR_TO_PSI;
  outletPressure += P_ATMOSPHERIC_PSI;

  // Temperature: Internal calculations use Celsius
  let temperature = parseFloat(temperatureEl.value) || 0;
  if (temperatureUnit === "f") temperature = (temperature - 32) * 5 / 9;
  else if (temperatureUnit === "k") temperature = temperature - 273.15;
  // If 'c', do nothing (already correct)

  let flowRate = parseFloat(flowRateEl.value) || 0;
  if (flowRateUnit === "m3h") flowRate *= conversionConstants.M3H_TO_GPM;

  let pipeDiameter = parseFloat(pipeDiameterEl.value) || 0;
  if (pipeDiameterUnit === "mm") pipeDiameter *= conversionConstants.MM_TO_IN;

  const fluidType = (fluidTypeEl ? parseInt(fluidTypeEl.value, 10) : 0) || 0;
  const valveType = (valveTypeEl ? parseInt(valveTypeEl.value, 10) : 0) || 0;

  UserInputs.setInput("inletPressure", inletPressure);
  UserInputs.setInput("outletPressure", outletPressure);
  UserInputs.setInput("temperature", temperature);
  UserInputs.setInput("flowRate", flowRate);
  UserInputs.setInput("pipeDiameter", pipeDiameter);
  UserInputs.setInput("fluidType", fluidType);
  UserInputs.setInput("valveType", valveType);

  // ======================================================
  // NEW: READ CUSTOM VALVE INPUTS (Direct Conversion)
  // ======================================================
  
  const custSizeInput = document.getElementById("customValveSize");
  const custSizeUnit = document.getElementById("customValveSizeUnit")?.value || "in";
  let custSize;

  // Logic: Only apply conversion if user actually typed a number.
  if (custSizeInput && custSizeInput.value !== "") {
      let rawVal = parseFloat(custSizeInput.value);
      if (isNaN(rawVal)) rawVal = 0;

      if (custSizeUnit === "mm") {
          // Direct mathematical conversion (e.g. 100mm -> 3.937")
          custSize = rawVal * conversionConstants.MM_TO_IN;
      } else {
          custSize = rawVal;
      }
  } else {
      // Fallback: Use Pipe Diameter (Already inches)
      custSize = pipeDiameter;
  }

  // --- UPDATED SECTION START: Read FL and FD directly from Inputs --- 
  // We removed the switch statement here because the inputs are already 
  // pre-filled by the event listener in DOMContentLoaded.
  
  const custFl = parseFloat(document.getElementById("customFl")?.value) || 0.9;
  const custFd = parseFloat(document.getElementById("customFd")?.value) || 0.46; 

  UserInputs.setInput("customValveSize", custSize);
  UserInputs.setInput("customFl", custFl);
  UserInputs.setInput("customFd", custFd);
  // --- UPDATED SECTION END ---

  // --- 3. Handle Advanced Inputs (Min/Max) ---
  const useAdvanced = document.getElementById("useAdvancedInputs")?.checked;
  let Q_min, P1_min, Q_max, P1_max, P2_min, P2_max;

  if (useAdvanced) {
      const qMinUnit = document.getElementById("flowRateMinUnit").value;
      let rawQmin = parseFloat(document.getElementById("flowRateMin").value);
      Q_min = (!isNaN(rawQmin)) ? ((qMinUnit === "m3h") ? rawQmin * conversionConstants.M3H_TO_GPM : rawQmin) : flowRate * 0.9;
      
      const pMinUnit = document.getElementById("inletPressureMinUnit").value;
      let rawP1min = parseFloat(document.getElementById("inletPressureMin").value);
      if (!isNaN(rawP1min)) {
          if (pMinUnit === "kpa_g") rawP1min *= conversionConstants.KPA_TO_PSI;
          else if (pMinUnit === "bar_g") rawP1min *= conversionConstants.BAR_TO_PSI;
          P1_min = rawP1min + P_ATMOSPHERIC_PSI;
      } else P1_min = inletPressure;

      const p2MinUnit = document.getElementById("outletPressureMinUnit").value;
      let rawP2min = parseFloat(document.getElementById("outletPressureMin").value);
      if (!isNaN(rawP2min)) {
          if (p2MinUnit === "kpa_g") rawP2min *= conversionConstants.KPA_TO_PSI;
          else if (p2MinUnit === "bar_g") rawP2min *= conversionConstants.BAR_TO_PSI;
          P2_min = rawP2min + P_ATMOSPHERIC_PSI;
      } else P2_min = outletPressure;

      const qMaxUnit = document.getElementById("flowRateMaxUnit").value;
      let rawQmax = parseFloat(document.getElementById("flowRateMax").value);
      Q_max = (!isNaN(rawQmax)) ? ((qMaxUnit === "m3h") ? rawQmax * conversionConstants.M3H_TO_GPM : rawQmax) : flowRate * 1.1;

      const pMaxUnit = document.getElementById("inletPressureMaxUnit").value;
      let rawP1max = parseFloat(document.getElementById("inletPressureMax").value);
      if (!isNaN(rawP1max)) {
          if (pMaxUnit === "kpa_g") rawP1max *= conversionConstants.KPA_TO_PSI;
          else if (pMaxUnit === "bar_g") rawP1max *= conversionConstants.BAR_TO_PSI;
          P1_max = rawP1max + P_ATMOSPHERIC_PSI;
      } else P1_max = inletPressure;

      const p2MaxUnit = document.getElementById("outletPressureMaxUnit").value;
      let rawP2max = parseFloat(document.getElementById("outletPressureMax").value);
      if (!isNaN(rawP2max)) {
          if (p2MaxUnit === "kpa_g") rawP2max *= conversionConstants.KPA_TO_PSI;
          else if (p2MaxUnit === "bar_g") rawP2max *= conversionConstants.BAR_TO_PSI;
          P2_max = rawP2max + P_ATMOSPHERIC_PSI;
      } else P2_max = outletPressure;

  } else {
      Q_min = flowRate * 0.9; P1_min = inletPressure; P2_min = outletPressure;
      Q_max = flowRate * 1.1; P1_max = inletPressure; P2_max = outletPressure;
  }

  UserInputs.setInput("flowRateMin", Q_min);
  UserInputs.setInput("inletPressureMin", P1_min);
  UserInputs.setInput("outletPressureMin", P2_min);
  UserInputs.setInput("flowRateMax", Q_max);
  UserInputs.setInput("inletPressureMax", P1_max);
  UserInputs.setInput("outletPressureMax", P2_max);

  // --- 4. Fluid Properties ---
  if (fluidType === UserInputs.Fluids.WATER) {
    const waterVaporPressure = estimateWaterVaporPressure(temperature);
    UserInputs.setFluidProperty("Pv", waterVaporPressure);
    calculateSpecificGravity();
    const estViscosity = estimateWaterKinematicViscosity_cSt(temperature);
    UserInputs.setFluidProperty("viscosity", estViscosity);
  } else if (fluidType === UserInputs.Fluids.CUSTOM) {
      const PcEl = document.getElementById("Pc");
      const PvEl = document.getElementById("Pv");
      const rhoEl = document.getElementById("density");
      const viscEl = document.getElementById("viscosity");
      const pcUnit = document.getElementById("pcUnit").value;
      const pvUnit = document.getElementById("pvUnit").value;
      const densityUnit = document.getElementById("densityUnit").value;
      const viscosityUnit = document.getElementById("viscosityUnit").value;

      let Pc = PcEl ? parseFloat(PcEl.value) || 0 : 0;
      if (pcUnit === "kpa_a") Pc *= conversionConstants.KPA_TO_PSI;
      else if (pcUnit === "bar_a") Pc *= conversionConstants.BAR_TO_PSI;

      let Pv = PvEl ? parseFloat(PvEl.value) || 0 : 0;
      if (pvUnit === "kpa_a") Pv *= conversionConstants.KPA_TO_PSI;
      else if (pvUnit === "bar_a") Pv *= conversionConstants.BAR_TO_PSI;

      let visc = viscEl ? parseFloat(viscEl.value) || 0 : 0;
      if (viscosityUnit === "m2s") visc *= conversionConstants.M2S_TO_CST;

      let rho_input = rhoEl ? parseFloat(rhoEl.value) || 0 : 0;
      let specificGravity = 1.0;
      if (densityUnit === "sg") specificGravity = rho_input;
      else if (rho_input > 0) specificGravity = rho_input / D_REFERENCE_WATER;

      if (specificGravity <= 0) specificGravity = 1.0;

      UserInputs.setFluidProperty("Pc", Pc);
      UserInputs.setFluidProperty("Pv", Pv);
      UserInputs.setFluidProperty("viscosity", visc);
      UserInputs.setFluidProperty("relativeDensity", specificGravity);
  }
}

/**
 * @description Main execution block. Attaches an event listener to the "Calculate" button.
 */
const calcBtn = document.getElementById("calculateBtn");
if (calcBtn) {
  calcBtn.addEventListener("click", () => {
    // Show debug log container if debug is enabled
    if (typeof DEBUG !== 'undefined' && DEBUG) {
      const debugContainer = document.getElementById("debugContainer");
      if (debugContainer) debugContainer.style.display = "block";
    }

    // 1. Update Parameters & Run Calculation
    parameterUpdate(); 
    const { recommendations, processInfo } = calculateAndRecommend();
    
    // 2. Get DOM Elements
    const resultEl = document.getElementById("result");
    const useAdvanced = document.getElementById("useAdvancedInputs")?.checked;

    // 3. Start building Result HTML
    let resultHTML = "";
    
    // --- PART 1: PROCESS ANALYSIS ---
    if (processInfo) {
        if (processInfo.error) {
             resultEl.innerHTML = `
                <div style="background:rgba(239, 68, 68, 0.1); border-left:3px solid var(--accent-red); padding:12px; border-radius:4px; font-size:14px; color:#fca5a5;">
                    <strong>${processInfo.errorMessage}</strong>
                </div>`;
             return;
        }

        resultHTML += `
          <div class="result-summary" style="display:flex; flex-direction:column; gap:16px;">
            <h3 style="margin:0; padding-bottom:8px; border-bottom:1px solid rgba(255,255,255,0.1);">Process Analysis</h3>

            <div>
                <div style="color:var(--muted); font-size:11px; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:4px;">
                    Required Cv (Operating)
                </div>
                <div style="font-size:18px; font-weight:700; color:var(--text-light); line-height:1.2;">
                    ${processInfo.baseCvOp.toFixed(2)}
                </div>
                ${useAdvanced ? `
                    <div style="display:flex; gap:12px; margin-top:4px; font-size:11px; align-items:baseline;">
                        <div>
                            <span style="text-transform:uppercase; color:var(--muted);">Min</span> 
                            <span style="color:var(--text-light); font-weight:600; margin-left:4px;">${processInfo.baseCvMin.toFixed(2)}</span>
                        </div>
                        <div style="color:var(--muted);">·</div>
                        <div>
                            <span style="text-transform:uppercase; color:var(--muted);">Max</span> 
                            <span style="color:var(--text-light); font-weight:600; margin-left:4px;">${processInfo.baseCvMax.toFixed(2)}</span>
                        </div>
                    </div>
                ` : ``}
            </div>

            <div>
                <div style="color:var(--muted); font-size:11px; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:4px;">
                    Recommended Characteristic
                </div>
                <div style="font-size:16px; font-weight:600; color:var(--accent);">
                    ${processInfo.preferredChar}
                </div>
                <div style="font-size:12px; color:var(--muted); margin-top:2px;">
                    Reason: ${processInfo.charReason}
                </div>
            </div>

            <div>
                <div style="font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:0.5px; margin-bottom:4px;">
                    Min. Required Body Class
                </div>
                <div style="font-size:16px; font-weight:700; color:var(--accent-yellow);">
                    Class ${processInfo.minProcessClass}
                </div>
                <div style="font-size:12px; color:var(--muted); margin-top:2px;">
                    Based on the ASTM A216 WCB Pressure-Temperature ratings, applied with a 25% safety margin (0.75 factor).
                </div>
            </div>

            <div>
                <div style="color:var(--muted); font-size:11px; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:4px;">
                    Cavitation Index <span style="text-transform:none; font-family:serif;">(σ)</span>
                </div>
                <div style="font-size:18px; font-weight:700; color:${processInfo.sigmaColor};">
                    ${processInfo.sigma}
                </div>
                <div style="font-size:12px; color:var(--muted); margin-top:2px;">
                    ${processInfo.sigmaLabel}
                </div>
            </div>
          </div>
        `;
        resultHTML += `</div>`; 
        resultHTML += `<h3 style="margin-top:24px;">Valve Recommendations</h3>`;
    }

    // --- PART 2: VALVE CARDS ---
    let finalPicks = [];
    
    // Logic for Custom vs Standard
    const customPick = recommendations.find(r => r.fit === "Custom");
    if (customPick) {
        finalPicks = [customPick]; 
    } else {
        const goodCandidates = recommendations.filter(r => r.fit === "Ideal" || r.fit === "Acceptable");
        const marginalCandidates = recommendations.filter(r => r.fit === "Marginal");
        
        if (goodCandidates.length > 0) finalPicks = goodCandidates.slice(0, 3);
        else if (marginalCandidates.length > 0) finalPicks = marginalCandidates.slice(0, 3);
    }

    // --- DEBUG LOG: SCORING LEADERBOARD ---
    if (typeof DEBUG !== 'undefined' && DEBUG && finalPicks.length > 0) {
        debugLog("--- 🏆 VALVE SCORING LEADERBOARD ---");
        finalPicks.forEach((v, index) => {
            debugLog(`RANK #${index + 1}: ${v.name}`, {
                "Total Score": v.score + "/100",
                "Operating Open": v.openOp.toFixed(1) + "%",
                "Noise": v.noise + " dBA",
                "Regime": v.regime
            });
        });
        debugLog("------------------------------------");
    }
    // --------------------------------------
    
    if (finalPicks.length === 0) {
        resultHTML += `
        <div style="background:rgba(239, 68, 68, 0.1); border-left:3px solid var(--accent-red); padding:12px; border-radius:4px; font-size:14px;">
            <strong>No suitable valves found.</strong><br>
            <span style="font-size:13px; color:var(--muted);">All candidates failed safety checks. Try increasing pipe size.</span>
        </div>`;
    } else {
        finalPicks.forEach(valve => {
            let cssClass = "";
            let cardStyle = "";
            let titleStyle = "";
            let gridHtml = "";
            let footerHtml = "";
            let charHtml = "";

            // --- BUILD CHARACTERISTIC BADGE (RIGHT SIDE - TOP) ---
            // Styling: White (#f3f4f6), Bold (700), Pill background
            if (valve.fit !== "Custom") {
                charHtml = `
                    <span style="margin-bottom:6px; font-size:10px; color:#f3f4f6; font-weight:700; background:rgba(255,255,255,0.05); padding:2px 6px; border-radius:4px; text-transform:uppercase; letter-spacing:0.5px;">
                        ${valve.characteristic}
                    </span>`;
            }

            if (valve.fit === "Custom") {
                // --- CUSTOM CARD STYLE ---
                cardStyle = `
                    background: rgba(255, 255, 255, 0.07); 
                    border: 1px solid rgba(255, 255, 255, 0.2);
                    border-left: 3px solid rgba(255, 255, 255, 0.5);
                `;
                titleStyle = "color: #f3f4f6; font-weight:600;";
                
                // Grid (Simple: Cv + Process TD only)
                gridHtml = `
                    <div style="background:rgba(0,0,0,0.2); padding:8px; border-radius:4px; grid-column: 1 / -1;">
                        <div style="color:var(--muted); font-size:10px; text-transform:uppercase; margin-bottom:2px;">Calculated Requirement</div>
                        <div style="display:flex; justify-content:space-between; align-items:baseline;">
                             <strong style="font-size:15px; color:#f3f4f6;">Cv: ${valve.cvOp.toFixed(2)}</strong>
                             <div style="color:var(--muted); font-size:11px;">Process TD: ${valve.processTurndown}</div>
                        </div>
                    </div>
                `;

                // Footer (Regime only, Right aligned)
                footerHtml = `
                    <div style="margin-top:8px; font-size:11px; color:var(--muted); text-align:right; font-style:italic;">
                        Flow Regime: ${valve.regime}
                    </div>
                `;

            } else {
                // --- STANDARD CARD STYLE ---
                if (valve.fit === "Ideal") cssClass = "ideal-fit";
                else if (valve.fit === "Acceptable") cssClass = "acceptable-fit";
                else cssClass = "marginal-fit";
                
                // Grid (Standard: Shows Min-Max Cv)
                gridHtml = `
                    <div style="background:rgba(0,0,0,0.2); padding:8px; border-radius:4px;">
                        <div style="color:var(--muted); font-size:10px; text-transform:uppercase; margin-bottom:2px;">Operating Point</div>
                        <strong style="font-size:15px; color:#f3f4f6;">${valve.openOp.toFixed(1)}% Open</strong>
                        <div style="color:var(--muted); font-size:11px; margin-top:2px;">Cv: ${valve.cvOp.toFixed(2)}</div>
                    </div>
                    <div style="background:rgba(0,0,0,0.2); padding:8px; border-radius:4px;">
                        <div style="color:var(--muted); font-size:10px; text-transform:uppercase; margin-bottom:2px;">Range (Min-Max)</div>
                        <strong style="font-size:15px; color:#f3f4f6;">${valve.openMin.toFixed(0)}% - ${valve.openMax.toFixed(0)}%</strong>
                        <div style="color:var(--muted); font-size:11px; margin-top:2px;">Cv: ${valve.cvMin.toFixed(2)} - ${valve.cvMax.toFixed(2)}</div>
                    </div>
                `;

                // Footer (Rangeability + Regime)
                const rangeText = valve.valveRangeability > 0 ? valve.valveRangeability.toFixed(0) + ":1" : "N/A";
                footerHtml = `
                    <div style="margin-top:8px; font-size:11px; color:var(--muted); display:flex; justify-content:space-between; align-items:center;">
                        <span style="font-style:italic;" title="Ratio of Cv at 90% to Cv at 10% open">
                            Turndown: ${rangeText}
                        </span>
                        <span style="font-style:italic;">Flow Regime: ${valve.regime}</span>
                    </div>
                `;
            }

            let noiseColor = "#34d399"; 
            let noiseBorder = "rgba(52, 211, 153, 0.3)";
            
            if (valve.noise >= 85) {
                noiseColor = "#ef4444"; 
                noiseBorder = "rgba(239, 68, 68, 0.3)";
            } else if (valve.noise >= 80) {
                noiseColor = "#f59e0b"; 
                noiseBorder = "rgba(245, 158, 11, 0.3)";
            }

            resultHTML += `
                <div class="recommendation-card ${cssClass}" style="${cardStyle}">
                    <div style="display:flex; justify-content:space-between; align-items:start; margin-bottom:12px;">
                        
                        <div>
                            <div class="card-title" style="margin:0; font-size:15px; margin-bottom:4px; ${titleStyle}">${valve.name}</div>
                            <span style="font-size:10px; color:#f3f4f6; font-weight:700; background:rgba(255,255,255,0.05); padding:2px 6px; border-radius:4px; text-transform:uppercase; letter-spacing:0.5px;">
                                ${valve.fit === "Custom" ? "Manual Check" : `Class ${valve.pressureClass}`}
                            </span>
                        </div>

                        <div style="display:flex; flex-direction:column; align-items:flex-end;">
                            ${charHtml}
                            <span style="font-size:10px; padding:2px 6px; border-radius:4px; border:1px solid ${noiseBorder}; color:${noiseColor}; font-weight:700;" title="Hydrodynamic Noise">
                                ${valve.noise} dBA
                            </span>
                        </div>
                    </div>
                    
                    <div style="display:grid; grid-template-columns: 1fr 1fr; gap:8px; font-size:13px;">
                        ${gridHtml}
                    </div>
                    
                    ${footerHtml}
                </div>
            `;
        });
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
  
  // --- 1. Initialize Tabs ---
  const tabStandard = document.getElementById("tabStandard");
  const tabCustom = document.getElementById("tabCustom");
  const panelStandard = document.getElementById("panelStandard");
  const panelCustom = document.getElementById("panelCustom");
  const valveTypeInput = document.getElementById("valveType");

  const switchTab = (isCustom) => {
    if (tabStandard) tabStandard.classList.toggle("active", !isCustom);
    if (tabCustom) tabCustom.classList.toggle("active", isCustom);
    if (panelStandard) panelStandard.style.display = isCustom ? "none" : "block";
    if (panelCustom) panelCustom.style.display = isCustom ? "block" : "none";
    if (valveTypeInput) valveTypeInput.value = isCustom ? "99" : "0";
  };

  if (tabStandard) tabStandard.addEventListener("click", () => switchTab(false));
  if (tabCustom) tabCustom.addEventListener("click", () => switchTab(true));

  // --- 2. Advanced Inputs Toggle ---
  const advCheckbox = document.getElementById("useAdvancedInputs");
  const advPanel = document.getElementById("advancedInputs");
  if (advCheckbox && advPanel) {
    advCheckbox.addEventListener("change", (e) => {
      advPanel.style.display = e.target.checked ? "grid" : "none";
    });
  }

  // --- 3. Fluid & Valve Type Visibility ---
  const fluidTypeEl = document.getElementById("fluidType");
  const customFluidDiv = document.getElementById("customFluidInputs");
  if(fluidTypeEl && customFluidDiv) {
      fluidTypeEl.addEventListener("change", (e) => {
          customFluidDiv.style.display = parseInt(e.target.value) === 1 ? "block" : "none";
      });
  }

  const valveTypeEl = document.getElementById("valveType");
  const customValveDiv = document.getElementById("customValveInputs");
  if (valveTypeEl && customValveDiv) {
      valveTypeEl.addEventListener("change", (e) => {
          customValveDiv.style.display = parseInt(e.target.value) === 99 ? "block" : "none";
      });
  }

  // --- 4. Custom Valve Style: Sync FL and FD (ROBUST) ---
  const styleSelect = document.getElementById("customValveStyle");
  
  if (styleSelect) {
    styleSelect.addEventListener("change", () => {
        const val = parseInt(styleSelect.value);
        let newFl = 0.90; 
        let newFd = 0.46;

        // Hardcoded Physics Defaults
        switch(val) {
            case 0: newFl = 0.90; newFd = 0.46; break; // Globe Open
            case 1: newFl = 0.80; newFd = 1.00; break; // Globe Close
            case 2: newFl = 0.60; newFd = 0.57; break; // Butterfly
            case 3: newFl = 0.60; newFd = 0.99; break; // Ball
        }

        // Safely update inputs if they exist
        const flInput = document.getElementById("customFl");
        if (flInput) {
            flInput.value = newFl;
            flashInput(flInput);
        }

        const fdInput = document.getElementById("customFd");
        if (fdInput) {
            fdInput.value = newFd;
            flashInput(fdInput);
        }
    });
  }

  // Helper for visual feedback
  function flashInput(el) {
      el.style.transition = "background-color 0.4s ease";
      el.style.backgroundColor = "rgba(255, 255, 255, 0.6)";
      setTimeout(() => el.style.backgroundColor = "", 400);
  }

  // --- 5. Populate Fluids Dropdown ---
  const fluidSelect = document.getElementById("fluidSelect");
  if (fluidSelect && typeof fluidsDatabase !== 'undefined') {
    fluidSelect.innerHTML = "";
    fluidsDatabase.forEach((fluid, index) => {
      const opt = document.createElement("option");
      opt.value = index;
      opt.textContent = fluid.name;
      fluidSelect.appendChild(opt);
    });

    fluidSelect.addEventListener("change", (e) => {
      const fluid = fluidsDatabase[e.target.value];
      if (fluid) {
        const setVal = (id, val) => { const el = document.getElementById(id); if(el) el.value = val; };
        setVal("sg", fluid.relativeDensity);
        setVal("pv", fluid.Pv);
        setVal("pc", fluid.Pc);
        setVal("viscosity", fluid.dynamicViscosity);
      }
    });
    // Initial Load
    fluidSelect.dispatchEvent(new Event("change"));
  }
});

document.getElementById("resetBtn").addEventListener("click", () => {
  // 1. Reset Main Process values (Flow, Pressure, Temp, Pipe)
  document.getElementById("flowRate").value = 530;
  document.getElementById("inletPressure").value = 85;
  document.getElementById("outletPressure").value = 30;
  document.getElementById("pipeDiameter").value = 4;
  document.getElementById("temperature").value = 77;
  
  // Reset Selectors (Back to defaults)
  document.getElementById("fluidType").value = 0; // Reset selection to Water
  document.getElementById("valveType").value = 0; // Reset selection to Globe
  document.getElementById("customValveStyle").value = 0; // Reset Custom Style

  // Reset Main units
  document.getElementById("flowRateUnit").value = "gpm";
  document.getElementById("inletPressureUnit").value = "psig";
  document.getElementById("outletPressureUnit").value = "psig";
  document.getElementById("pipeDiameterUnit").value = "in";
  document.getElementById("temperatureUnit").value = "f";

  // --- NEW: Reset Advanced Min/Max Checkbox ---
  const advCheckbox = document.getElementById("useAdvancedInputs");
  const advPanel = document.getElementById("advancedInputs");
  if (advCheckbox && advPanel) {
      advCheckbox.checked = false; // Uncheck box
      advPanel.style.display = "none"; // Hide panel immediately
  }
  
  // Reset Custom Valve inputs (These are usually specific to the attempt, so we reset them)
  document.getElementById("customValveSize").value = 4;
  document.getElementById("customValveSizeUnit").value = "in";
  document.getElementById("customFl").value = 0.90;

  // 2. Hide Custom Panels
  // The data is still there, just hidden.
  document.getElementById("customFluidInputs").style.display = "none";
  document.getElementById("customValveInputs").style.display = "none";

  // 3. Hide Debug Log
  const debugContainer = document.getElementById("debugContainer");
  if (debugContainer) {
      debugContainer.style.display = "none"; 
  }

  // 4. Reset Result Text
  document.getElementById("result").innerHTML = "Reset to defaults.";
  
  // 5. Trigger Update
  parameterUpdate();
  logToConsole("Inputs reset (Custom Fluid data preserved).");
});

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
    document.getElementById("flowRate").value = 530;
    document.getElementById("inletPressure").value = 85;
    document.getElementById("outletPressure").value = 30; // dp = 1 psi
    document.getElementById("pipeDiameter").value = 4;
    document.getElementById("temperature").value = 77; // 20C
  } else if (preset === 2) {
    document.getElementById("flowRate").value = 10;
    document.getElementById("inletPressure").value = 100;
    document.getElementById("outletPressure").value = 90; // dp = 10 psi
    document.getElementById("pipeDiameter").value = 1;
    document.getElementById("temperature").value = 68; // 20C
  } else {
    document.getElementById("flowRate").value = 390;
    document.getElementById("inletPressure").value = 20;
    document.getElementById("outletPressure").value = 10; // dp = 100 psi
    document.getElementById("pipeDiameter").value = 4;
    document.getElementById("temperature").value = 68; // 20C
  }
  // parameterUpdate(); // This is called by the click()
  document.getElementById("calculateBtn").click();
}

document.getElementById('test1').addEventListener('click', () => runPreset(1));
document.getElementById('test2').addEventListener('click', () => runPreset(2));
document.getElementById('test3').addEventListener('click', () => runPreset(3));
