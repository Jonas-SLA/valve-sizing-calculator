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

    FF = LTS_ValveSizing_CalculateFF(Pv, Pc)
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
// --- Updated Valve Style Specific Factors ---
const valveGlobesSinglePortContouredPlugOpen = {
    Fl: 0.9,
    Fd: 0.46,
    Kc: 0.65, // Typical for contoured globe
    Cv: 0,
    Xt: 0.72,
};

const valveGlobesSinglePortContouredPlugClose = {
    Fl: 0.8,
    Fd: 1.0,
    Kc: 0.68, 
    Cv: 0,
    Xt: 0.55,
};

const valveButterfly = {
    Fl: 0.62,
    Fd: 0.57,
    Kc: 0.38, // Lower recovery index
    Cv: 0,
    Xt: 0.35,
};

const valveBallFull = {
    Fl: 0.74,
    Fd: 0.99,
    Kc: 0.28, 
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
        customFd: 0.46,
        // NEW: User Controlled Variable Inputs for Flow Characteristics Recommendation
        controlledVariable: 0
    };

    const Fluids = Object.freeze({ WATER: 0, CUSTOM: 1 });
    const Valves = Object.freeze({ GLOBEOPEN: 0, GLOBECLOSE: 1, BUTTERFLY: 2, BALL: 3, CUSTOM: 99 });
    const ControlledVariables = Object.freeze({ NOT_SPECIFIED: 0, LEVEL: 1, PRESSURE: 2, FLOW: 3, FLOW_SPECIFIC_CASE: 4});

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
        ControlledVariables
    };
})();

//------------------------------------------------------------------------------------
// SECTION 3: FLUID PROPERTY CALCULATIONS
//------------------------------------------------------------------------------------
/**
 * LTS_ValveSizing_CalculateWaterDensity
 * Calculates the density of water based on temperature, pressure, and salinity.
 * Uses a polynomial approximation for temperature and simplified adjustments for pressure/salinity.
 * @param {number} [S_percent=0] - Salinity in percent (%).
 * @returns {number|null} The final water density in kg/m³, or null on error.
 * * Purpose:
 * Determines the fluid density required for Reynolds number and Specific Gravity
 * calculations. It accounts for compressibility (pressure) and salinity effects
 * on top of the standard temperature curve.
 *
 * Formula & notes:
 * 1. Base (1 atm): 5th-order polynomial of T(°C).
 * rho_T = rho_0 + a1*T - a2*T^2 + ...
 * 2. Pressure Correction:
 * Linear adjustment based on bulk modulus approx (4.5e-3 kg/m³ per dbar).
 * 3. Salinity Correction:
 * Empirical polynomial adjustment for dissolved salts.
 *
 * Variables:
 * - T = Temperature (°C)
 * - P = Pressure (psi -> converted to dbar)
 * - S = Salinity (converted to parts per thousand)
 *
 * Returns:
 * number (Density in kg/m³).
 */
function LTS_ValveSizing_CalculateWaterDensity(S_percent = 0) {
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
    // Convert psi -> dbar (approx). Water density rises ~4.5e-3 kg/m³ per dbar.
    let P_dbar = 0;
    if (typeof pressure === "number" && !Number.isNaN(pressure)) {
        P_dbar = pressure * 0.689476;
    }

    // ---- 3) Salinity correction ----
    // Convert % -> ‰ (practical oceanographic convention)
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

/**
 * LTS_ValveSizing_EstimateWaterViscosityCp
 * Estimates the dynamic viscosity of water using a Vogel-like empirical formula.
 * This is a practical estimator for water between approximately 0-100 °C.
 * @param {number} T - Temperature in Celsius (°C).
 * @returns {number} The estimated dynamic viscosity in centipoise (cP).
 * * Purpose:
 * Calculates the fluid's resistance to flow (dynamic viscosity). 
 * Used as an intermediate step to find Kinematic Viscosity (v).
 *
 * Formula & notes:
 * Uses the Vogel/Andrade equation:
 * mu (Pa·s) = A * 10^(B / (Tk - C))
 * - Result is multiplied by 1000 to convert Pa·s to cP.
 * - Clamped between 0.1 and 1000 cP for safety.
 *
 * Variables:
 * - Tk = Temperature in Kelvin (T + 273.15)
 * - A, B, C = Empirical constants for water
 *
 * Returns:
 * number (Dynamic Viscosity in cP).
 */
function LTS_ValveSizing_EstimateWaterViscosityCp(T) {
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
 * LTS_ValveSizing_EstimateWaterKinematicViscosityCSt
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
 *   • Dynamic viscosity μ is obtained from LTS_ValveSizing_EstimateWaterViscosityCp()
 *   • Density ρ is obtained from LTS_ValveSizing_CalculateWaterDensity()
 *   • Density function already reads temperature from UserInputs internally
 *   • Water salinity is fixed at 0% here (fresh water)
 */
function LTS_ValveSizing_EstimateWaterKinematicViscosityCSt(T_C) {
    // Step 1: Calculate dynamic viscosity (cP)
    const mu_cP = LTS_ValveSizing_EstimateWaterViscosityCp(T_C); // dynamic in cP

    // Step 2: density (kg/m³) — salinity input only, temperature taken internally
    const rho_kg_m3 = LTS_ValveSizing_CalculateWaterDensity(0); // fresh water, no salinity
    if (!rho_kg_m3) return null;

    // Convert density to g/cm³ (1 g/cm³ = 1000 kg/m³)
    const rho_g_cm3 = rho_kg_m3 / 1000.0;
    
    // Engineering approximation: ν(cSt) ≈ μ(cP) / ρ(g/cm³)
    return mu_cP / rho_g_cm3;
}

/**
 * LTS_ValveSizing_CalculateSpecificGravity
 * Calculates the specific gravity (relative density) of water and updates the state.
 * @returns {number|null} The calculated specific gravity, or null on error.
 * * Purpose:
 * SG is the ratio of fluid density to the density of a reference fluid (water at 4°C).
 * It appears in almost all Cv equations (square root term).
 *
 * Formula & notes:
 * SG = rho_water / rho_ref
 * - rho_water: calculated at current T and P.
 * - rho_ref: 1000 kg/m³ (Standard water density).
 *
 * Variables:
 * - rho_water = Density from LTS_ValveSizing_CalculateWaterDensity
 * - D_REFERENCE_WATER = 1000.0
 *
 * Returns:
 * number (Specific Gravity, dimensionless).
 */
function LTS_ValveSizing_CalculateSpecificGravity() {
    const salinityPercent = 0;
    // Step 1: compute water density (kg/m³)
    const waterDensity = LTS_ValveSizing_CalculateWaterDensity(salinityPercent);
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
 * LTS_ValveSizing_EstimateWaterVaporPressure
 * Estimates the vapor pressure of water in PSI using the Antoine equation.
 * Valid range extended for pressurized water.
 * @param {number} T - Temperature in Celsius (°C).
 * @returns {number} The estimated vapor pressure in PSI.
 * * Purpose:
 * Determines the pressure at which the liquid will boil (flash).
 * Critical for calculating the Choked Flow threshold and FF factor.
 *
 * Formula & notes:
 * Antoine Equation: log10(P_mmHg) = A - (B / (C + T))
 * - Result converted from mmHg to PSI.
 * - Constants A, B, C are specific to water.
 *
 * Variables:
 * - T = Temperature (°C)
 * - A=8.07131, B=1730.63, C=233.426
 *
 * Returns:
 * number (Vapor Pressure in psia).
 */
function LTS_ValveSizing_EstimateWaterVaporPressure(T) {
    // Physical lower bound check (Freezing point approximation)
    if (T <= 0) return 0.0886; 

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
 * LTS_ValveSizing_CalculateFF
 * Calculates the liquid pressure recovery factor (FF) for per ISA-75.01.01 standard.
 * Based on the formula: FF = 0.96 - 0.28 * sqrt(Pv/Pc).
 * @param {number} vaporPressure - Vapor pressure (psia)
 * @param {number} criticalPressure - Critical pressure (psia)
 * @returns {number|null} FF factor
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
function LTS_ValveSizing_CalculateFF(vaporPressure, criticalPressure) {
    // Validate inputs are numbers and Pc is positive to prevent division by zero errors
    if (typeof vaporPressure !== "number" || typeof criticalPressure !== "number") {
        console.warn("Pv or Pc is not a valid number.");
        return null;
    }
    if (criticalPressure <= 0) {
        console.warn("Critical pressure must be > 0.");
        return null;
    }

    // Calculate pressure ratio (Pv/Pc); warn if thermodynamic state is invalid (e.g. Pv > Pc)
    const ratio = vaporPressure / criticalPressure;
    if (ratio < 0 || ratio > 1) {
        console.warn("Pv/Pc ratio outside 0..1 — check inputs.");
    }

    // Apply standard formula: FF = 0.96 - 0.28 * sqrt(Pv/Pc), ensuring non-negative root
    let FF = 0.96 - 0.28 * Math.sqrt(Math.max(0, ratio));

    // Clamp FF to a minimum of 0.05 to prevent numerical instability (divide-by-zero) downstream
    if (FF < 0.05) FF = 0.05;

    return FF;
}

/**
 * LTS_ValveSizing_IsChokedFlow
 * Determines if the liquid flow regime is choked (critical) per ISA-75.01.01.
 * based on the comparison: Actual ΔP vs. Maximum Allowable (Choked) ΔP.
 * @param {number} FF - Liquid critical pressure ratio factor (calculated via Pv/Pc)
 * @returns {boolean|null} True if choked, False if not, Null on error
 * * Purpose:
 * Identifies if the valve has reached its maximum flow capacity due to
 * vaporization (cavitation/flashing) at the vena contracta.
 * - Subcritical: Flow increases as outlet pressure drops.
 * - Choked: Flow is limited; lowering outlet pressure further adds no flow.
 *
 * Formula & notes:
 * Condition: ΔP_actual >= ΔP_choked
 * - ΔP_actual = P1 - P2
 * - ΔP_choked = FL^2 * (P1 - FF * Pv)
 * * Variables:
 * - P1 = Inlet Pressure, P2 = Outlet Pressure
 * - Pv = Vapor Pressure
 * - FL = Liquid Pressure Recovery Factor (Valve geometry constant)
 * - FF = Critical Pressure Ratio Factor (Fluid property)
 *
 * Returns:
 * boolean (true = Choked/Critical, false = Subcritical).
 */
function LTS_ValveSizing_IsChokedFlow(inletPressure, outletPressure, vaporPressure, FL, FF) {
    // Validate that all inputs (including the passed FF) are valid numbers
    if ([FL, vaporPressure, inletPressure, outletPressure, FF].some(v => typeof v !== "number" || Number.isNaN(v))) {
        console.warn("Invalid input(s) in choked flow check.");
        return null;
    }

    // Ensure physical validity: flow requires inlet pressure > outlet pressure
    if (inletPressure <= outletPressure) {
        console.warn("Inlet pressure must be > outlet pressure.");
        return null;
    }

    // Calculate critical ΔP limit using ISA formula: ΔP_crit = FL^2 * (P1 - FF * Pv)
    const chokedThreshold = Math.pow(FL, 2) * (inletPressure - (FF * vaporPressure));

    // Verify the threshold is valid (e.g., P1 must be > FF*Pv for liquid state to exist at inlet)
    if (!isFinite(chokedThreshold) || chokedThreshold < 0) {
        console.warn("Invalid choked threshold calculation.");
        return null;
    }

    // Return true if actual ΔP (P1 - P2) meets or exceeds the critical choked limit
    return (inletPressure - outletPressure) >= chokedThreshold;
}

/**
 * LTS_ValveSizing_CalculateChokedC
 * Calculates the required Valve Coefficient (C) under choked (critical) flow conditions.
 * Used when actual ΔP exceeds the critical ΔP, limiting flow due to flashing/cavitation.
 *
 * @param {number} FF - Critical Pressure Ratio Factor (calculated via Pv/Pc).
 * @returns {number|null} The required Cv value, or null if inputs are invalid.
 *
 * Purpose:
 * In choked flow, increasing the pressure drop (ΔP) further does NOT increase flow rate.
 * Therefore, the standard Cv equation (based on actual ΔP) is invalid.
 * This function uses the "Choked Flow Equation" per ISA-75.01.01, which substitutes
 * actual ΔP with the limiting pressure drop term: (P1 - FF * Pv).
 *
 * Formula:
 * C = Q / (N1 * FL) * sqrt( Gf / (P1 - FF * Pv) )
 *
 * Variables:
 * - Q  = Flow Rate (e.g., gpm)
 * - N1 = Unit conversion constant (numericalConstants.N1)
 * - FL = Liquid Pressure Recovery Factor (valve geometry specific)
 * - Gf = Specific Gravity (relative density) of the fluid
 * - P1 = Inlet Pressure
 * - Pv = Vapor Pressure
 * - FF = Liquid Critical Pressure Ratio Factor
 *
 * Returns:
 * number (Cv) or null on calculation error.
 */
function LTS_ValveSizing_CalculateChokedC(flowRate, inletPressure, vaporPressure, specificGravity, FL, FF) {
    // Validate that all required physical properties are valid numbers
    if ([flowRate, inletPressure, FL, vaporPressure, specificGravity].some(v => typeof v !== "number" || Number.isNaN(v))) {
        console.warn("Missing inputs for choked Cv calculation.");
        return null;
    }

    // Calculate effective driving pressure (P1 - FF*Pv); ensures fluid is not already flashing at inlet
    const effectivePressureDrop = (inletPressure - (FF * vaporPressure));
    if (effectivePressureDrop <= 0) {
        console.warn("Denominator <= 0 in choked Cv calc (flashing at inlet).");
        return null;
    }

    // Compute Cv using choked formula where flow is limited by vaporization, not outlet pressure
    // Note: numericalConstants.N1 handles unit conversion (e.g., gpm/psi vs m3h/bar)
    const Cv = (flowRate / (numericalConstants.N1 * FL)) * Math.sqrt(specificGravity / effectivePressureDrop);

    return Cv;
}

/**
 * LTS_ValveSizing_CalculateSubcriticalC
 * Calculates required Cv for subcritical (turbulent, non-choked) liquid flow.
 * Standard: ISA-75.01.01 / IEC 60534-2-1.
 * * Purpose:
 * In the subcritical region, flow velocity is relatively low, and no vaporization 
 * occurs. The flow rate is directly proportional to the square root of the 
 * pressure differential (ΔP).
 *
 * Formula:
 * Cv = (Q / N1) * sqrt( Gf / (P1 - P2) )
 *
 * Variables:
 * - Q  = Flow Rate (e.g., gpm)
 * - N1 = Unit conversion constant (numericalConstants.N1)
 * - Gf = Specific Gravity (relative density)
 * - P1 = Inlet Pressure, P2 = Outlet Pressure
 *
 * Returns:
 * number (Cv) or null on error.
 */
function LTS_ValveSizing_CalculateSubcriticalC(flowRate, inletPressure, outletPressure, specificGravity, FP) {  
    // Default Fp to 1.0 (no piping losses) if argument is missing
    const _Fp = (typeof FP === "number") ? FP : 1.0;

    // Calculate actual pressure drop (ΔP); basic driving force for subcritical flow
    const deltaP = inletPressure - outletPressure;

    // Validate that inputs are numbers and check for valid pressure drop (> 0)
    if ([flowRate, inletPressure, outletPressure, specificGravity, _Fp].some(v => typeof v !== "number" || Number.isNaN(v))) {
        console.warn("Missing inputs for subcritical C calculation.");
        return null;
    }
    if (deltaP <= 0) {
        console.warn("deltaP must be > 0 (Flow requires P1 > P2).");
        return null;
    }

    // Calculate Cv using standard square-root law: Flow ~ sqrt(ΔP)
    // numericalConstants.N1 adjusts for units (e.g. Imperial vs Metric)
    const Cv = (flowRate / (numericalConstants.N1 * _Fp)) * Math.sqrt(specificGravity / deltaP);

    return Cv;
}

/**
 * LTS_ValveSizing_CalculateValveReynoldsRev
 * Calculates the Valve Reynolds Number (Rev) per ISA-75.01.01 standard.
 * Used to correct the flow coefficient (Cv) for viscous (laminar/transitional) flow.
 *
 * @param {number} Ci  - Estimated Flow Coefficient (Cv) for the current iteration.
 * @param {number} FLP - Liquid Pressure Recovery Factor (FL or FLP if fittings exist).
 * @param {number} d   - Valve inlet internal diameter (in inches).
 * @returns {number|null} The calculated Reynolds number, or null on error.
 *
 * Purpose:
 * Standard Cv equations assume turbulent flow. When viscosity is high or velocity 
 * is low (low Rev), the flow becomes laminar, requiring a correction factor (FR).
 * This function calculates Rev to determine that regime.
 *
 * Formula:
 * Rev = [ (N4 * Fd * Q) / (v * sqrt(FL * Ci)) ] * [ (FL^2 * Ci^2) / (N2 * d^4) + 1 ]^(1/4)
 *
 * Variables:
 * - N4, N2 = Unit conversion constants
 * - Fd = Valve Style Modifier (accounts for hydraulic diameter differences)
 * - v  = Kinematic Viscosity (cSt)
 * - FL = Pressure Recovery Factor
 *
 * Returns:
 * number (Rev) or null if inputs are invalid.
 */
function LTS_ValveSizing_CalculateValveReynoldsRev(Ci, FLP, d, flowRate, v, FD) {
// Basic validation for all inputs
    if ([Ci, FLP, d, flowRate, v, FD].some(x => typeof x !== 'number' || isNaN(x))) {
        console.warn("[LTS_ValveSizing] Invalid inputs for Rev calculation.");
        return null;
    }
    
    // Physical validation (dimensions must be positive)
    if (Ci <= 0 || d <= 0 || v <= 0) {
        console.warn("[LTS_ValveSizing] Ci, d, or v must be > 0.");
        return null;
    }

    // Use the passed recovery factor (FL or FLP)
    const FL = FLP;

    // Term 1: Main Reynolds expression based on flow velocity and hydraulic radius
    // Rev ~ Velocity * Diameter / Viscosity
    const eq1 = (numericalConstants.N4 * FD * flowRate) / (v * Math.sqrt(Ci * FL));

    // Term 2: Geometry correction factor for the valve inlet shape
    // [ 1 + (FL^2 * Cv^2) / (N2 * D^4) ]^(1/4)
    const inside = (Math.pow(FL, 2) * Math.pow(Ci, 2)) / (numericalConstants.N2 * Math.pow(d, 4)) + 1;
    const eq2 = Math.pow(inside, 0.25);

    // Combine terms to get final Valve Reynolds Number
    const Rev = eq1 * eq2;

    // Final Check
    if (!isFinite(Rev) || Rev < 0) {
        console.warn("[LTS_ValveSizing] Invalid Rev computed (Infinite or Negative).");
        return null;
    }

    return Rev;
}

/**
 * LTS_ValveSizing_CalculateFpFactors
 * Calculates the Piping Geometry Factor (Fp) and the Liquid Pressure Recovery Factor 
 * with Piping (FLP) for a valve installed between reducers.
 * * Standard: ISA-75.01.01 / IEC 60534-2-3 (Annex C).
 * * @param {number} Cv - The estimated Flow Coefficient (Cv) for the current sizing iteration.
 * @param {number} valveSize - The internal diameter of the valve (d) in inches.
 * @param {number} pipeSize - The internal diameter of the pipe (D) in inches.
 * @param {number} FL - The valve's intrinsic Liquid Pressure Recovery Factor.
 * @returns {object} { Fp, FLP } - Piping Geometry Factor and Installed Recovery Factor.
 * * Purpose:
 * When a valve is smaller than the line size, reducers and expanders create additional 
 * pressure losses. 
 * - Fp corrects the Cv for these losses (effective capacity decreases).
 * - FLP corrects the choking limit (choking happens earlier due to fitting losses).
 * * Formulas:
 * - Fp = 1 / sqrt( 1 + (ΣK / N2) * (Cv/d^2)^2 )
 * - FLP = FL / sqrt( 1 + FL^2 * (Ki / N2) * (Cv/d^2)^2 )
 * (Simplified in code as FLP = [ (1/FL^2) + (Ki/N2)*(Cv/d^2)^2 ]^-0.5 )
 */
function LTS_ValveSizing_CalculateFpFactors(Cv, valveSize, pipeSize, FL) {
    // ---------------------------------------------------------
    // 1. LINE-SIZED CHECK
    // ---------------------------------------------------------
    // If valve diameter matches pipe diameter, there are no fittings.
    // Fp = 1.0 (no geometric loss) and FLP = FL (intrinsic recovery only).
    if (Math.abs(valveSize - pipeSize) < 0.001) {
        return { Fp: 1.0, FLP: FL };
    }

    // ---------------------------------------------------------
    // 2. LOSS COEFFICIENTS (K-FACTORS)
    // ---------------------------------------------------------
    const N2 = numericalConstants.N2; // Unit constant (e.g., 890 for inch units)
    const d = valveSize;
    const D = pipeSize;

    // Calculate diameter ratio squared (d/D)^2 for standard concentric reducers
    const ratioSq = (d * d) / (D * D);

    // K1: Resistance coefficient for inlet reducer (contraction)
    // Formula: K1 = 0.5 * (1 - (d/D)^2)^2
    const K1 = 0.5 * Math.pow(1 - ratioSq, 2);

    // K2: Resistance coefficient for outlet expander (expansion)
    // Formula: K2 = 1.0 * (1 - (d/D)^2)^2
    const K2 = 1.0 * Math.pow(1 - ratioSq, 2);

    // KB1: Bernoulli coefficient for inlet head loss
    // Formula: KB1 = 1 - (d/D)^4
    const KB1 = 1 - Math.pow(ratioSq, 2);
    
    // Sum of loss coefficients (Inlet + Outlet)
    const SigmaK = K1 + K2;
    
    // Inlet Loss Factor for Choking Calculation (K1 + KB1)
    const Ki = K1 + KB1;

    // ---------------------------------------------------------
    // 3. FACTOR CALCULATIONS
    // ---------------------------------------------------------
    // Pre-calculate the velocity head term: (Cv / d^2)^2
    const Cvd2_term = Math.pow(Cv / (d * d), 2);

    // Calculate Fp (Piping Geometry Factor)
    // Measures the flow reduction caused by the fittings.
    const Fp = Math.pow(1 + (SigmaK / N2) * Cvd2_term, -0.5);

    // Calculate FLP (Installed Pressure Recovery Factor)
    // Measures the combined recovery of the valve + fittings.
    // Used instead of FL in choked flow equations when fittings are present.
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
 * LTS_ValveSizing_VerifyFR
 * Determines which Reynolds Number Factor (FR) formula branch to use.
 * * Standard: ISA-75.01.01 / IEC 60534-2-1.
 * * @param {number} C - The Flow Coefficient (Cv) being evaluated.
 * @param {number} d - The valve internal diameter (in inches).
 * @returns {boolean} True if "Full Range" FR formulas apply; False if "Restricted" (Small Flow) formulas apply.
 * * Purpose:
 * ISA standards distinguish between "standard" valve trim and "low-flow" trim.
 * This check compares the capacity (Cv) to the valve area (d^2) to decide which regime applies.
 *
 * Formula:
 * IF (Cv / d^2) > (0.016 * N18) THEN Use Standard FR
 * ELSE Use Low-Flow FR
 *
 * Variables:
 * - N18 = Unit conversion constant (numericalConstants.N18)
 * - 0.016 = Empirical boundary constant for low-flow trim
 */
function LTS_ValveSizing_VerifyFR(C, d) {
    const D = d;
    
    // Safety check: if diameter is invalid, default to standard formulas (True)
    if (!D || D <= 0) return true; 

    // ---------------------------------------------------------
    // 1. CAPACITY DENSITY CALCULATION
    // ---------------------------------------------------------
    // Calculate Capacity Density: Cv per square inch of valve area.
    // Low values imply a very restrictive trim (e.g., needle valve or micro-flow).
    const capacityDensity = C / Math.pow(D, 2);

    // ---------------------------------------------------------
    // 2. BOUNDARY THRESHOLD
    // ---------------------------------------------------------
    // Compare against the standard low-flow boundary limit.
    // N18 scales the 0.016 constant based on the unit system (e.g., metric vs imperial).
    const threshold = 0.016 * numericalConstants.N18;

    // Returns TRUE if capacity is large enough for standard equations.
    // Returns FALSE if capacity is small (Low-Flow regime).
    return capacityDensity > threshold;
}

/**
 * Calculates the viscosity correction factor (FR) for the restricted-range condition.
 * @param {number} Ci - The estimated flow coefficient (Cv).
 * @param {number} Rev - The calculated valve Reynolds number.
 * @param {number} d_pipe - The valve trim diameter.
 * @param {number} FLP - The liquid pressure recovery factor.
 * @returns {number|null} The viscosity correction factor FR (always <= 1), or null on error.
 *
 * Viscosity correction factor FR — restricted-range branch (verifyFR == false).
 * Used when valve capacity is small relative to port size (Cv/d^2 is low).
 *
 * Purpose:
 * Compute FR (≤ 1) for viscous (low-Re) flows in small/restrictive trims.
 * Calculates two curve fits (Laminar & Transitional) and picks the most restrictive.
 *
 * Inputs:
 * - Ci (Cv estimate)
 * - Rev (Reynolds Number)
 * - d_pipe (Valve diameter)
 * - FLP (Recovery factor)
 *
 * Returns:
 * FR (0 < FR ≤ 1) or null on error.
 */
function LTS_ValveSizing_CalculateFRifNo(Ci, Rev, d_pipe, FLP) {
    if (Ci == null || Rev == null) {
        console.warn("Missing inputs to calculate FR (No).");
        return null;
    }

    // Map inputs: FL (Recovery Factor) and d (Valve/Trim Diameter)
    let FR = 0;
    const FL = FLP;
    const d = d_pipe;

    // Calculate 'n2': Empirical geometry factor for low-flow trims based on Cv/d^2
    const n2 = 1 + numericalConstants.N32 * Math.pow((Ci / (d * d)), 2/3);

    if (Rev < 10) {
        // Laminar Regime (Rev < 10): Flow is purely viscous; FR follows linear slope
        FR = (0.026 / FL) * Math.sqrt(n2 * Rev);
        FR = Math.min(FR, 1);
    } else {
        // Transitional Regime: Flow is between laminar and turbulent
        // We calculate two curve fits (F3A = Turbulent-transition, F4 = Laminar-extension)
        
        // F3A: Logarithmic curve approaching turbulent limit (1.0)
        const parte1 = (0.33 * Math.sqrt(FL)) / Math.pow(n2, 0.25);
        const parte2 = Math.log10(Rev / 10000);
        const FR_F3A = 1 + parte1 * parte2;

        // F4: Continuation of the laminar slope
        const FR_F4 = (0.026 / FL) * Math.sqrt(n2 * Rev);

        // Physical Reality: Flow follows the most restrictive curve (Minimum of curves & 1.0)
        FR = Math.min(FR_F3A, FR_F4, 1);
    }

    return FR;
}

/**
 * Calculates the viscosity correction factor (FR) for the full-range condition.
 * @param {number} Ci - The estimated flow coefficient (Cv).
 * @param {number} Rev - The calculated valve Reynolds number.
 * @param {number} d_pipe - The valve trim diameter.
 * @param {number} FLP - The liquid pressure recovery factor.
 * @returns {number|null} The viscosity correction factor FR (always <= 1), or null on error.
 *
 * Viscosity correction factor FR — full-range branch (verifyFR == true).
 * Used when valve capacity is normal relative to port size (Cv/d^2 is high).
 *
 * Purpose:
 * Compute FR (≤ 1) for viscous (low-Re) flows in standard control valves.
 * Calculates two limiting curve fits and selects the most restrictive one.
 *
 * Inputs:
 * - Ci (Cv estimate)
 * - Rev (Reynolds Number)
 * - d_pipe (Valve diameter)
 * - FLP (Recovery factor)
 *
 * Returns:
 * FR (0 < FR ≤ 1) or null on error.
 */
function LTS_ValveSizing_CalculateFRifYes(Ci, Rev, d_pipe, FLP) {
    // Validate inputs; Rev drives the regime, Ci drives the geometry factor
    if (Ci == null || Rev == null) {
        console.warn("Missing inputs to calculate FR (Yes).");
        return null;
    }

    // Map inputs: FL (Recovery Factor) and d (Valve/Trim Diameter)
    let FR = 0;
    const FL = FLP;
    const d = d_pipe;

    // Calculate 'n1': Empirical geometry factor for standard trims based on Capacity Density
    const n1 = numericalConstants.N2 * (Ci / (d * d));

    if (Rev < 10) {
        // Laminar Regime (Rev < 10): Flow is purely viscous; FR follows linear slope
        FR = (0.026 / FL) * Math.sqrt(n1 * Rev);
        // Clamp result to 1.0 (physics dictates viscosity cannot increase capacity)
        FR = Math.min(FR, 1);
    } else {
        // Transitional Regime: Flow is between laminar and turbulent
        // We calculate two curve fits (F1A = Turbulent-transition, F2 = Laminar-extension)

        // F1A: Logarithmic curve approaching turbulent limit (1.0)
        const parte1 = (0.33 * Math.sqrt(FL)) / Math.pow(n1, 0.25);
        const parte2 = Math.log10(Rev / 10000);
        const FR_F1A = 1 + parte1 * parte2;

        // F2: Continuation of the laminar slope
        const FR_F2 = (0.026 / FL) * Math.sqrt(n1 * Rev);

        // Physical Reality: Flow follows the most restrictive curve (Minimum of curves & 1.0)
        FR = Math.min(FR_F1A, FR_F2, 1);
    }

    return FR;
}

//------------------------------------------------------------------------------------
// SECTION 6: PRESSURE-TEMPERATURE RATINGS (ASME B16.34 - Group 1.1 WCB)
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
 * LTS_ValveSizing_GetMaxPressureWCB
 * Calculates the Maximum Allowable Working Pressure (MAWP) for WCB steel.
 * Interpolates between standard ASME B16.34 data points.
 * @param {number} tempF - Temperature in Fahrenheit.
 * @param {number} classRating - ANSI Class (150, 300, 600...).
 * @returns {number} MAWP in psig.
 * * Purpose:
 * Determines the absolute physical pressure limit of the valve body material
 * at the specific operating temperature provided by the user.
 *
 * Formula & notes:
 * Uses Linear Interpolation:
 * P = P1 + (T - T1) * (P2 - P1) / (T2 - T1)
 * - T = Target Temp, P = Target Pressure
 * - T1, T2 = Nearest temperature points in WCB_Ratings
 * - P1, P2 = Corresponding pressures
 *
 * Variables:
 * - tempF = Operating Temperature (°F)
 * - classRating = ANSI Pressure Class
 *
 * Returns:
 * number (Pressure in psig).
 */
function LTS_ValveSizing_GetMaxPressureWCB(tempF, classRating) {
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
 * LTS_ValveSizing_CheckPressureClass
 * Validates if the selected valve class is sufficient for the process pressure.
 * Applies a safety factor (derating) to the standard ASME rating.
 * @param {object|number} valve - Valve object (with .pressureClass) or raw class number.
 * @param {number} maxInletPressure - Maximum inlet pressure (psig).
 * @param {number} tempF - Temperature in Fahrenheit.
 * @returns {object} { pass: boolean, limit: number, rated: number }
 * * Purpose:
 * Safety check to ensure the valve body will not fail under pressure.
 * We enforce a 75% utilization limit (25% safety margin) on top of the code standard.
 *
 * Formula & notes:
 * Limit = MAWP_ASME * 0.75
 * Pass = InletPressure <= Limit
 *
 * Variables:
 * - MAWP_ASME = Result from LTS_ValveSizing_GetMaxPressureWCB
 * - 0.75 = Internal Safety Factor
 *
 * Returns:
 * Object containing pass/fail status and the specific limits calculated.
 */
function LTS_ValveSizing_CheckPressureClass(valve, maxInletPressure, tempF) {
    // If we're passing a specific class (number) or object with pressureClass property
    const rating = typeof valve === 'number' ? valve : valve.pressureClass;
    
    if (!rating) return { pass: true, limit: 0, rated: 0 }; // Pass if no class defined (legacy data)

    // 1. Get Base ASME Rating
    const asmeLimit = LTS_ValveSizing_GetMaxPressureWCB(tempF, rating);

    // 2. Apply 75% Safety Factor
    const safeLimit = asmeLimit * 0.75;

    // 3. Compare
    const pass = maxInletPressure <= safeLimit;

    return { pass, limit: safeLimit, rated: asmeLimit };
}

/**
 * LTS_ValveSizing_RecommendFlowCharacteristic
 * Recommends the optimal valve trim characteristic based on the SENAI standard.
 * @param {object} process - { qMin, qMax, dpMin, dpMax }
 * @param {number} controlledVar - Selected from UserInputs.ControlledVariables
 * @returns {object} { type: string, reason: string }
 * * Purpose:
 * To maintain a constant control loop gain, the valve characteristic should complement
 * the process characteristic.
 * - Constant dP -> Linear Valve.
 * - Variable dP (High drop at low flow) -> Equal Percentage Valve.
 *
 * Formula & notes:
 * 1. Ratio = dpMax / dpMin
 * 2. If Choked -> Equal Percentage (Critical override)
 * 3. If Ratio > 2.0 -> Quick Opening (if large flow range) or Eq% (default)
 * 4. If Ratio ~ 1.0 (0.8 to 1.25) -> Linear
 * 5. Otherwise -> Equal Percentage
 *
 * Variables:
 * - dpMax = Pressure drop at maximum flow (maximum drop)
 * - dpMin = Pressure drop at minimum flow (minimum drop)
 *
 * Returns:
 * Object { type: "Linear"|"Equal Percentage"|"Quick Opening", reason: "..." }
 */
function LTS_ValveSizing_RecommendFlowCharacteristic(process, controlledVar) {
    
    // --- 1. Safety Validation ---
    if (!process.qMin || !process.qMax || !process.dpMin || !process.dpMax) {
        return { type: "Equal Percentage", reason: "Insufficient data. Defaulting to Equal Percentage for safety." };
    }

    // --- 2. Physics Translation to SENAI Terminology ---
    // The lowest pressure drop (dpMin) occurs at maximum flow (due to piping friction).
    const dpAtMaxFlow = process.dpMin; 
    const dpAtMinFlow = process.dpMax; 

    // --- 3. SENAI Guidelines Application ---
    switch (controlledVar) {
        
        case UserInputs.ControlledVariables.LEVEL:
            // Rules for Liquid Level (SENAI)
            // Using a 1% tolerance to define "Constant"
            if (Math.abs(dpAtMaxFlow - dpAtMinFlow) / dpAtMinFlow < 0.01) {
                return { type: "Linear", reason: "Level control with Constant pressure drop. Linear response recommended." };
            }

            if (dpAtMaxFlow < dpAtMinFlow) {
                // Decreasing pressure drop as flow increases
                if (dpAtMaxFlow > (0.20 * dpAtMinFlow)) {
                    return { type: "Linear", reason: "Level control with decreasing ΔP, and ΔP at max flow is > 20% of ΔP at min flow. Linear recommended." };
                } else {
                    return { type: "Equal Percentage", reason: "Level control with decreasing ΔP, and ΔP at max flow is < 20% of ΔP at min flow. Equal Percentage recommended." };
                }
            } else {
                // Increasing pressure drop as flow increases
                if (dpAtMaxFlow > (2.00 * dpAtMinFlow)) {
                    return { type: "Linear", reason: "Level control with increasing ΔP, and ΔP at max flow is > 200% of ΔP at min flow. Linear recommended." };
                } else {
                    return { type: "Quick Opening", reason: "Level control with increasing ΔP, and ΔP at max flow is < 200% of ΔP at min flow. Quick Opening recommended." };
                }
            }

        case UserInputs.ControlledVariables.PRESSURE:
            // Rule for Pressure (Liquids)
            return { type: "Equal Percentage", reason: "For liquid systems controlling Pressure, Equal Percentage is recommended." };

        case UserInputs.ControlledVariables.FLOW:
            // Regra SENAI: Controle de Vazão (Série / Proporcional ao fluxo)
            // Calculamos a razão de variação do fluxo (Turndown)
            const flowRatio = process.qMax / process.qMin;

            if (flowRatio > 1.5) {
                // "Grandes variações de fluxo" -> Linear
                return { type: "Linear", reason: "Flow control with Large flow variations. Linear response recommended." };
            } else {
                // "Pequenas variações de fluxo..." -> Equal Percentage
                return { type: "Equal Percentage", reason: "Flow control but with Small flow variations. Equal Percentage recommended." };
            }

        case UserInputs.ControlledVariables.FLOW_SPECIFIC_CASE:
            // Rule for Flow (Bypass / Proportional to flow squared)
            return { type: "Equal Percentage", reason: "Flow (Bypass/Squared Signal). Equal Percentage recommended." };

        case UserInputs.ControlledVariables.NOT_SPECIFIED:
        default:
            return { type: "Equal Percentage", reason: "Controlled variable not specified. Defaulting to Equal Percentage for safety." };
    }
}

/**
 * Data tables for Noise Calculation (Figures 55, 56, 57).
 * Source: Masoneilan Control Valve Sizing Handbook / SENAI standard.
 */
const valveNoiseChart = {
    // X-Axis: Pressure Drop (psi)
    dp: [1, 2, 3, 5, 7, 10, 20, 50, 100, 200, 500, 1000],

    // Y-Axis: Sound Pressure Level (dBA) for specific Cv values
    curves: [
        { Cv: 1,    spl: [33,35,36,38,39,41,43,46,48,51,54,56] },
        { Cv: 2,    spl: [36,38,39,41,42,44,46,49,51,54,57,59] },
        { Cv: 5,    spl: [40,42,43,45,46,48,50,53,55,58,61,63] },
        { Cv: 10,   spl: [43,45,46,48,49,51,53,56,58,61,64,66] },
        { Cv: 20,   spl: [46,48,49,51,52,54,56,59,61,64,67,69] },
        { Cv: 50,   spl: [50,52,53,54,56,58,60,63,65,68,71,73] },
        { Cv: 100,  spl: [53,55,56,57,59,61,63,66,68,71,74,76] },
        { Cv: 200,  spl: [56,58,59,60,62,64,66,69,71,74,77,79] },
        { Cv: 500,  spl: [60,62,63,64,66,68,70,73,75,78,81,83] },
        { Cv: 1000, spl: [64,66,67,68,70,72,74,77,79,82,85,87] },
        { Cv: 2000, spl: [66,68,69,70,72,74,76,79,81,84,87,89] },
        { Cv: 5000, spl: [70,72,73,74,76,78,80,83,85,88,91,93] },
        { Cv: 10000,spl: [73,75,76,77,79,81,83,86,88,91,94,97] }
    ]
};

const spldChart = {
    // X-Axis: Pressure Drop (psi)
    dp: [1, 2, 3, 5, 7, 10, 20, 50, 100, 200, 500, 1000],
    // Y-Axis: Correction Factor (dBA)
    spld: [0, 3, 8, 13, 17, 20, 26, 33, 40, 45, 54, 61]
};

// Pre-calculate Log values for faster interpolation
const preparedSPL = {
  dpLog: valveNoiseChart.dp.map(Math.log10),
  curves: valveNoiseChart.curves.map(c => ({ 
    ...c, 
    cvLog: Math.log10(c.Cv) 
  }))
};

const preparedSPLd = {
  dpLog: spldChart.dp.map(Math.log10),
  spld: spldChart.spld
};

// Linear Interpolation Utility
const lerp = (x, x0, x1, y0, y1) => {
  if (Math.abs(x1 - x0) < 1e-9) return y0; 
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
};

/**
 * LTS_ValveSizing_CalculateSPL
 * Calculates the Base Sound Pressure Level (SPL) from Cv and Delta P.
 * Interpolates from the standard "Figure 55" curves.
 * @param {number} deltaP - Pressure drop across the valve (psi).
 * @param {number} Cv - Valve flow coefficient.
 * @param {object} chart - The prepared noise chart data.
 * @returns {number} Base SPL in dBA.
 * * Purpose:
 * Determines the baseline noise generated by fluid turbulence before corrections.
 *
 * Formula & notes:
 * Uses Bilinear Interpolation on Log-Log scale:
 * 1. Interpolate SPL for the specific Delta P on each Cv curve.
 * 2. Interpolate between Cv curves to find the final SPL.
 */
function LTS_ValveSizing_CalculateSPL(deltaP, Cv, chart) {
    if (deltaP <= 0 || Cv <= 0) throw new Error("deltaP and Cv must be > 0");

    const targetLogDp = Math.log10(deltaP);
    const targetLogCv = Math.log10(Cv);
    const dpGrid = chart.dpLog;

    // 1. Interpolate SPL vs dP for every Cv curve in the library
    const splAtDp = chart.curves.map(curve => {
        const spl = curve.spl;
        // Clamp lower bound
        if (targetLogDp <= dpGrid[0]) return { cvLog: curve.cvLog, spl: spl[0] };
        // Clamp upper bound
        if (targetLogDp >= dpGrid[dpGrid.length - 1]) return { cvLog: curve.cvLog, spl: spl[spl.length - 1] };

        // Find interval
        for (let i = 0; i < dpGrid.length - 1; i++) {
            if (targetLogDp <= dpGrid[i + 1]) {
                return {
                    cvLog: curve.cvLog,
                    spl: lerp(targetLogDp, dpGrid[i], dpGrid[i + 1], spl[i], spl[i + 1])
                };
            }
        }
    });

    // 2. Interpolate SPL vs Cv using the results from Step 1
    // Clamp lower bound
    if (targetLogCv <= splAtDp[0].cvLog) return splAtDp[0].spl;
    // Clamp upper bound
    if (targetLogCv >= splAtDp[splAtDp.length - 1].cvLog) return splAtDp[splAtDp.length - 1].spl;

    // Find interval
    for (let i = 0; i < splAtDp.length - 1; i++) {
        if (targetLogCv <= splAtDp[i + 1].cvLog) {
            return lerp(targetLogCv, splAtDp[i].cvLog, splAtDp[i + 1].cvLog, splAtDp[i].spl, splAtDp[i + 1].spl);
        }
    }
    return 0; // Should not reach here
}

/**
 * LTS_ValveSizing_CalculateSPLd
 * Calculates the Outlet Pressure Correction (SPL_delta) or Downstream correction.
 * Interpolates from "Figure 56".
 * @param {number} deltaP - Specifically (Outlet Pressure - Vapor Pressure).
 * @param {object} chart - The prepared SPLd chart data.
 * @returns {number} Correction factor in dBA.
 * * Purpose:
 * Adjusts the noise level based on the downstream static pressure condition.
 * Higher outlet pressure suppresses cavitation noise formation.
 *
 * Formula & notes:
 * Uses Linear Interpolation on Log scale for dP.
 */
function LTS_ValveSizing_CalculateSPLd(deltaP, chart) {
  if (deltaP <= 0) throw new Error("deltaP must be > 0");

  const targetLogDp = Math.log10(deltaP);
  const dpGrid = chart.dpLog;

  // Clamp bounds
  if (targetLogDp <= dpGrid[0]) return chart.spld[0];
  if (targetLogDp >= dpGrid[11]) return chart.spld[11];

  // Find interval
  for (let i = 0; i < 11; i++) {
    if (targetLogDp <= dpGrid[i + 1]) {
      return lerp(targetLogDp, dpGrid[i], dpGrid[i + 1], chart.spld[i], chart.spld[i + 1]);
    }
  }
}

/**
 * LTS_ValveSizing_CalculateSPLc
 * Calculates the Cavitation Correction Factor (SPL_c).
 * Based on "Figure 57" geometry formulas.
 * @param {number} InletPressure - P1 (psia)
 * @param {number} OutletPressure - P2 (psia)
 * @param {number} Pv - Vapor Pressure (psia)
 * @param {number} Fl - Pressure Recovery Factor
 * @param {number} Kc - Cavitation Index
 * @returns {number} Correction factor in dBA (0 to 10).
 * * Purpose:
 * Reduces the predicted noise as the valve approaches full cavitation (super-cavitation).
 * When flow is fully choked, the noise stops increasing or dampens due to vapor bubbles.
 *
 * Formula & notes:
 * X = (Actual_dP - Kc(P1-Pv)) / [ (P1-Pv)(FL^2 - Kc) ]
 * SPL_c = 10 * (1 - X)
 * - Result is clamped between 0 and 10 dBA.
 *
 * Variables:
 * - Kc: Incipient cavitation index (usually 0.65 * FL^2 for globes)
 */
function LTS_ValveSizing_CalculateSPLc(InletPressure, OutletPressure, Pv, Fl, Kc) {
    // Threshold check: If Kc is very close to FL^2, the formula is unstable or N/A
    const threshold = 0.9 * (Fl * Fl);
    if (Kc > threshold) {
        return 0;
    }

    const actualDeltaP = InletPressure - OutletPressure;
    const pMinusPv = InletPressure - Pv;

    // Denominator: Operating range between Incipient (Kc) and Choked (FL^2)
    const denom = pMinusPv * (Fl * Fl - Kc);

    // Safety check for divide-by-zero
    if (Math.abs(denom) < 1e-9 || denom < 0) return 10;

    // Calculate ratio X (0 = Incipient, 1 = Fully Choked)
    let X = (actualDeltaP - Kc * pMinusPv) / denom;

    // Clamp X to 0..1 range
    X = Math.max(0, Math.min(1, X));

    // Calculate attenuation: 10 dBA max reduction at full choke
    return 10 - (10 * X);
}

/**
 * LTS_ValveSizing_CalculateFinalSPL
 * Master function to compute total Hydrodynamic Noise using the SENAI/Masoneilan method.
 * Combines Base SPL, Outlet Correction, and Cavitation Attenuation.
 * @param {object} params - Input object { InletPressure, OutletPressure, Pv, Fl, Kc, Cv ... }
 * @returns {number} Final Noise Level in dBA.
 * * Purpose:
 * Predicts the aerodynamic/hydrodynamic noise level at 1 meter from the valve.
 * Used to check against safety limits (typically 85 dBA).
 *
 * Formula & notes:
 * SPL_Total = SPL_Base + SPL_OutletCorrection - SPL_CavitationCorrection
 * - SPL_Base: Function of Cv and DeltaP.
 * - SPL_Outlet: Correction for downstream pressure (P2).
 * - SPL_Cavitation: Attenuation factor as flow becomes fully choked.
 *
 * Variables:
 * - actualDP = P1 - P2
 * - chokedDP = FL^2 * (P1 - Pv)
 *
 * Returns:
 * number (dBA).
 */
function LTS_ValveSizing_CalculateFinalSPL({InletPressure, OutletPressure, Pv, Fl, Kc, Cv, splChart = preparedSPL, spldChart = preparedSPLd}) {
    // Calculate physical pressure drops
    const actualDP = InletPressure - OutletPressure;
    const chokedDP = (Fl * Fl) * (InletPressure - Pv);

    // Downstream pressure relative to vapor pressure
    const deltaP_SPL = Math.max(0.001, actualDP - chokedDP);
    const deltaP_SPLd = Math.max(0.001, OutletPressure - Pv);

    // 1. Base SPL (Function of Energy Conversion)
    const SPL  = LTS_ValveSizing_CalculateSPL(deltaP_SPL, Cv, splChart);

    // 2. Outlet Pressure Correction (Suppression of bubble formation)
    const SPLd = LTS_ValveSizing_CalculateSPLd(deltaP_SPLd, spldChart);

    // 3. Cavitation Attenuation (Super-cavitation masking)
    const SPLc = LTS_ValveSizing_CalculateSPLc(InletPressure, OutletPressure, Pv, Fl, Kc);

    // Final Summation
    return SPL + SPLd - SPLc;
}


/**
 * LTS_ValveSizing_CalculateCavitationSigma
 * Calculates the ISA-RP75.23 Cavitation Index (Sigma) and determines severity.
 * @param {number} P1 - Inlet Pressure (psia).
 * @param {number} P2 - Outlet Pressure (psia).
 * @param {number} Pv - Vapor Pressure (psia).
 * @returns {object} Result object { sigma, severity, consequence, color }.
 * * Purpose:
 * Sigma (σ) is a dimensionless number describing how close the fluid is to flashing.
 * Lower Sigma = Higher potential for cavitation damage.
 *
 * Formula & notes:
 * Sigma = (P1 - Pv) / (P1 - P2)
 *
 * Thresholds:
 * - σ > 2.0: Safe / Incipient
 * - 1.7 < σ < 2.0: Moderate
 * - 1.0 < σ < 1.5: Severe
 * - σ < 1.0: Flashing (Fluid is boiling at outlet)
 *
 * Variables:
 * - P1, P2, Pv must be in absolute units (psia).
 *
 * Returns:
 * Object with classification strings and UI color codes.
 */
function LTS_ValveSizing_CalculateCavitationSigma(P1, P2, Pv) {
    const dP = P1 - P2;

    // Safety check for invalid dP (reverse flow or zero drop)
    if (dP <= 0) return { sigma: 999, severity: "Safe", consequence: "No Risk", color: "#34d399" };

    // Calculate Sigma
    const sigma = (P1 - Pv) / dP;
    const sigmaVal = parseFloat(sigma.toFixed(2));

    // Default Result (Safe)
    let result = {
        sigma: sigmaVal,
        severity: "Safe",
        consequence: "No Risk of Cavitation",
        color: "#34d399" // Green
    };

    // Fetch the current valve type from state
    const currentValveType = UserInputs.getProperty("getInputs.valveType");
    const isBallValve = (currentValveType === UserInputs.Valves.BALL);

    // Determine Severity based on industry standard ranges
    if (sigmaVal <= 1.0) {
        result.severity = "Flashing";
        result.consequence = "Flashing is occurring.";
        result.color = "#7f1d1d"; // Dark Red / Purple
    }
    else if(isBallValve) {
        if (sigmaVal <= 3.5) {
            result.severity = "Severe";
            result.consequence = "Potential for severe cavitation.";
            result.color = "#ef4444"; // Red
        } 
        else if (sigmaVal <= 4.0) {
            result.severity = "Moderate";
            result.consequence = "Some cavitation control required.";
            result.color = "#f97316"; // Orange
        } 
        else if (sigmaVal < 4.7) {
            result.severity = "Incipient";
            result.consequence = "No cavitation control required (Hardened trim recommended).";
            result.color = "#eab308"; // Yellow
        }
    }
    else {
        if (sigmaVal <= 1.5) {
            result.severity = "Severe";
            result.consequence = "Potential for severe cavitation.";
            result.color = "#ef4444"; // Red
        } 
        else if (sigmaVal <= 1.7) {
            result.severity = "Moderate";
            result.consequence = "Some cavitation control required.";
            result.color = "#f97316"; // Orange
        } 
        else if (sigmaVal < 2.0) {
            result.severity = "Incipient";
            result.consequence = "No cavitation control required (Hardened trim recommended).";
            result.color = "#eab308"; // Yellow
        }
    }

    return result;
}

/**
 * LTS_ValveSizing_InterpolateDp
 * Linearly interpolates Pressure Drop (dP) based on Valve Opening %.
 * Simulates the "System Curve" effect where dP changes as the valve moves.
 * @param {number} targetOpen - Target valve opening (0-100%).
 * @param {number} openMin - Opening at minimum flow condition.
 * @param {number} dpAtMin - Pressure drop at minimum flow.
 * @param {number} openMax - Opening at maximum flow condition.
 * @param {number} dpAtMax - Pressure drop at maximum flow.
 * @returns {number} Interpolated pressure drop (psi).
 * * Purpose:
 * In a real piping system, as the valve opens, flow increases, and line losses increase.
 * This causes the available pressure drop across the valve (dP) to decrease.
 * This function models that relationship linearly between the Min and Max operating points.
 *
 * Formula & notes:
 * dP = dP_min + (Open - Open_min) * Slope
 * - Clamps 'Open' to [Min, Max] to prevent extrapolation errors.
 * - Enforces dP > 0.1 psi to avoid division-by-zero in downstream math.
 */
function LTS_ValveSizing_InterpolateDp(targetOpen, openMin, dpAtMin, openMax, dpAtMax) {
    // 1. Safety Clamp: Ensure we don't extrapolate into non-physical pressures.
    // This keeps interpolation strictly bounded by the Min/Max operating points.
    const clampedOpen = Math.min(Math.max(targetOpen, openMin), openMax);

    // 2. Safety: Avoid division by zero if operating range is negligible
    if (Math.abs(openMax - openMin) < 0.1) return (dpAtMin + dpAtMax) / 2;

    // Linear mapping
    const slope = (dpAtMax - dpAtMin) / (openMax - openMin);
    const dP = dpAtMin + (clampedOpen - openMin) * slope;

    // 3. Final Physics Safety: dP cannot be zero or negative
    return Math.max(0.1, dP);
}

/**
 * LTS_ValveSizing_CalculateInstalledGain
 * Calculates the "Installed Gain" (dQ/dStroke) of the valve in the system.
 * Estimates the change in Flow (Q) for a small change in Stroke.
 * @param {number[]} cvCurve - 10-point Cv curve from manufacturer database.
 * @param {number} currentOpen - The operating opening % (0-100).
 * @param {number} sg - Specific Gravity of fluid.
 * @param {object} sys - System Reference Points { openMin, dpAtMin, openMax, dpAtMax }.
 * @returns {number|null} The installed gain value (slope), or null on error.
 * * Purpose:
 * Determines if the valve is controllable at the operating point.
 * - Gain too high (> 3.0): Unstable control (hunting).
 * - Gain too low (< 0.2): Sluggish control (deadband).
 * - Ideal Gain is ~1.0 (Linear relationship installed).
 *
 * Formula & notes:
 * Gain = (Q_high - Q_low) / (Stroke_high - Stroke_low)
 * - Q = N1 * Cv * sqrt(dP / SG)
 * - Uses Central Difference method (checking points above and below current open).
 *
 * Variables:
 * - sys.dpAtMin/Max: System pressure drop at min/max flow (defines the system curve).
 */
function LTS_ValveSizing_CalculateInstalledGain(cvCurve, currentOpen, sg, sys) {
    if (!cvCurve || cvCurve.length !== 10) return null;
    if (sg <= 0) return null; // Safety against bad SG

    // 0. Sanity Swap: Ensure Min < Max for interpolation logic
    // This handles user error (e.g. Qmin > Qmax) or weird sizing results
    if (sys.openMin > sys.openMax) {
        [sys.openMin, sys.openMax] = [sys.openMax, sys.openMin];
        [sys.dpAtMin, sys.dpAtMax] = [sys.dpAtMax, sys.dpAtMin];
    }

    // 1. Identify Curve Index (0=10%, 9=100%)
    let idx = Math.round(currentOpen / 10) - 1;
    if (idx < 0) idx = 0;
    if (idx > 9) idx = 9;

    let idxLow, idxHigh;
    let strokeLow, strokeHigh;

    // 2. Determine Interval (Forward, Backward, or Central Difference)
    // We prefer Central (looking both ways) but must use Forward/Backward at edges.
    if (idx === 0) { // Forward (at 10% open)
        idxLow = 0; idxHigh = 1;
        strokeLow = 10; strokeHigh = 20;
    } else if (idx === 9) { // Backward (at 100% open)
        idxLow = 8; idxHigh = 9;
        strokeLow = 90; strokeHigh = 100;
    } else { // Central (Middle) - Best Accuracy
        idxLow = idx - 1; idxHigh = idx + 1;
        strokeLow = (idx - 1 + 1) * 10; // e.g. index 3 (40%) -> 40
        strokeHigh = (idx + 1 + 1) * 10; // e.g. index 5 (60%) -> 60
    }

    // 3. Get Cv at neighbors
    const cvLow = cvCurve[idxLow];
    const cvHigh = cvCurve[idxHigh];

    // 4. Estimate dP at neighbors (The "Installed" part)
    // We infer the dP at the neighbor points using the process Min/Max anchors
    const dpLow = LTS_ValveSizing_InterpolateDp(strokeLow, sys.openMin, sys.dpAtMin, sys.openMax, sys.dpAtMax);
    const dpHigh = LTS_ValveSizing_InterpolateDp(strokeHigh, sys.openMin, sys.dpAtMin, sys.openMax, sys.dpAtMax);

    // 5. Calculate Flow at neighbors (Q = Cv * sqrt(dP/SG))
    // Note: N1 constant (1.0) is assumed standard for GPM/PSI
    const qLow = cvLow * Math.sqrt(dpLow / sg);
    const qHigh = cvHigh * Math.sqrt(dpHigh / sg);

    // 6. Calculate Gain (dQ / dStroke)
    const dQ = qHigh - qLow;
    const dStroke = strokeHigh - strokeLow;

    return dQ / dStroke;
}

//------------------------------------------------------------------------------------
// SECTION 7: MAIN SIZING ORCHESTRATOR
//------------------------------------------------------------------------------------
/**
 * LTS_ValveSizing_SizingOrchestrator
 * Main driver function for Liquid Valve Sizing per ISA-75.01.01.
 * Coordinates data gathering, regime determination (Turbulent vs Laminar), 
 * and the iterative solution for piping geometry factors.
 *
 * @param {number} valveSize - Nominal valve internal diameter (d).
 * @param {number} [overrideQ] - Optional Flow Rate override (for graphing/simulation).
 * @param {number} [overrideP1] - Optional Inlet Pressure override.
 * @param {number} [overrideP2] - Optional Outlet Pressure override.
 * @returns {object|null} Result { Cv, Re, iterations, choked, warning } or null on error.
 */
/**
 * Purpose:
 * Solves the non-linear system of equations for valve sizing.
 * * The Problem:
 * 1. To calculate Cv, we need the Piping Geometry Factor (Fp).
 * 2. To calculate Fp, we need the Cv (to know the velocity loss in reducers).
 * * Solution:
 * Uses a fixed-point iteration loop:
 * - Guess Cv (using intrinsic valve FL).
 * - Calculate Fp based on Guess.
 * - Recalculate Cv using new Fp.
 * - Repeat until Cv converges (change < 0.01%).
 */
function LTS_ValveSizing_SizingOrchestrator(valveSize, overrideQ, overrideP1, overrideP2) {
    const EPS = 0.0001;      // Convergence tolerance (0.01%)
    const MAX_ITERS = 20;    // Safety break to prevent infinite loops

    // -------------------------------------------------------
    // 1. GATHER INPUTS
    // -------------------------------------------------------
    // Resolve Process Conditions (Use override if provided, else UserInputs)
    const Q  = (overrideQ !== undefined)  ? overrideQ  : UserInputs.getProperty("getInputs.flowRate");
    const P1 = (overrideP1 !== undefined) ? overrideP1 : UserInputs.getProperty("getInputs.inletPressure");
    const P2 = (overrideP2 !== undefined) ? overrideP2 : UserInputs.getProperty("getInputs.outletPressure");

    // Fluid Properties: Vapor Pressure (Pv), Critical Pressure (Pc), Specific Gravity (SG), Viscosity (v)
    const Pv = UserInputs.getProperty("getSelectedFluid.Pv");
    const Pc = UserInputs.getProperty("getSelectedFluid.Pc");
    const SG = UserInputs.getProperty("getSelectedFluid.relativeDensity");
    const v  = UserInputs.getProperty("getSelectedFluid.viscosity"); 

    // Valve Properties: Pipe Size (D), Recovery Factor (FL), Valve Style Modifier (Fd)
    const pipeSize = UserInputs.getProperty("getInputs.pipeDiameter");
    const FL       = UserInputs.getProperty("getSelectedValve.Fl");
    const Fd       = UserInputs.getProperty("getSelectedValve.Fd");

    // Basic physics sanity checks
    if (P1 <= P2) { console.warn("Sizing Error: P1 must be > P2"); return null; }
    if (!FL || !SG || !Pv || !Pc) { console.warn("Sizing Error: Missing fluid/valve properties"); return null; }


    // -------------------------------------------------------
    // 2. PRE-CALCULATION (Intrinsic State)
    // -------------------------------------------------------
    // --- Step 1: Calculate FF (Critical Pressure Ratio) ---
    // Determines the fluid's tendency to flash/cavitate relative to its critical point
    const FF = LTS_ValveSizing_CalculateFF(Pv, Pc);
    if (FF == null) return null;

    // Step 2: Calculate Base Cv (Initial Guess)
    // Assume no fittings (Fp=1.0) to get a starting point for the iteration
    const isChoked_Base = LTS_ValveSizing_IsChokedFlow(P1, P2, Pv, FL, FF);
    let C_base = 0;

    if (isChoked_Base) {
        // Choked Flow: Capacity limited by FL (Recovery Factor)
        C_base = LTS_ValveSizing_CalculateChokedC(Q, P1, Pv, SG, FL, FF);
    } else {
        // Subcritical Flow: Capacity driven by sqrt(P1 - P2)
        C_base = LTS_ValveSizing_CalculateSubcriticalC(Q, P1, P2, SG, 1.0);
    }

    if (!C_base) return null;

    // Step 3: Calculate Reynolds Number (Rev)
    // Determines if we are in Turbulent (Standard) or Laminar (Viscous) regime
    let Rev = LTS_ValveSizing_CalculateValveReynoldsRev(C_base, FL, valveSize, Q, v, Fd);
    if (Rev == null) return null;


    // -------------------------------------------------------
    // 3. REGIME SELECTION & ITERATION
    // -------------------------------------------------------
    // === BRANCH A: TURBULENT FLOW (Rev > 10,000) ===
    // Standard industrial sizing equations apply here.
    if (Rev > 10000) {
        
        // Case A1: Valve Size = Pipe Size (No geometric corrections needed)
        // Return the Base Cv calculated in Step 2.
        if (Math.abs(valveSize - pipeSize) < 0.01) {
            return { Cv: C_base, Re: Rev, iterations: 0, choked: isChoked_Base };
        } 
        
        // Case A2: Valve Size != Pipe Size (Geometry Correction Loop)
        // Fittings reduce capacity (Fp < 1) and alter choking (FLP != FL).
        else {
            let Ci = C_base; 
            let final_choked = false;

            for (let i = 0; i < MAX_ITERS; i++) {
                // 1. Update Geometry Factors (Fp, FLP) based on current Cv guess
                const factors = LTS_ValveSizing_CalculateFpFactors(Ci, valveSize, pipeSize, FL);
                let Fp = factors.Fp || 1.0;
                const FLP = factors.FLP;

                // 2. Determine "Effective FL" for Choking Check
                // Fittings cause pressure drop *before* the valve, lowering the choking threshold
                const effectiveFL = FLP / Fp;
                
                // 3. Re-evaluate Choked Flow condition with new geometry factors
                const isIterChoked = LTS_ValveSizing_IsChokedFlow(P1, P2, Pv, effectiveFL, FF);
                
                let C_next = 0;

                if (isIterChoked) {
                    // Choked: Calculate using FLP (Installed Recovery Factor)
                    C_next = LTS_ValveSizing_CalculateChokedC(Q, P1, Pv, SG, FLP, FF);
                    final_choked = true;
                } else {
                    // Subcritical: Calculate using Fp (Piping Geometry Factor)
                    C_next = LTS_ValveSizing_CalculateSubcriticalC(Q, P1, P2, SG, Fp);
                    final_choked = false;
                }

                // 4. Convergence Check: Stop if Cv stops changing
                if (Math.abs(C_next - Ci) / C_next < EPS) {
                    return { Cv: C_next, Re: Rev, iterations: i+1, choked: final_choked };
                }
                
                // Update guess for next loop iteration
                Ci = C_next;
            }

            // Fallback: Return last calculated value if convergence fails
            return { Cv: Ci, Re: Rev, iterations: MAX_ITERS, choked: final_choked };
        }
    } 
    
    // === BRANCH B: LAMINAR/TRANSITIONAL FLOW (Rev <= 10,000) ===
    //  Non-turbulent flow logic per ISA-75.01.01 flowchart
    else {
        // Establish initial guess: Ci = 1.3 * C_base
        let Ci = 1.3 * C_base; 
        let FR = 1.0;

        for (let i = 0; i < MAX_ITERS; i++) {
            // 1. Recalculate Reynolds Number with current guess (Ci)
            // Note: We use the *intrinsic* FL here as per standard for Re calc
            Rev = LTS_ValveSizing_CalculateValveReynoldsRev(Ci, FL, valveSize, Q, v, Fd);
            if (Rev == null) break; // Safety break

            // 2. Determine Regime (Full Range vs Restricted)
            // Uses pure helper from Section 5
            const useFullRange = LTS_ValveSizing_VerifyFR(Ci, valveSize); 

            // 3. Calculate Viscosity Correction Factor (FR)
            if (useFullRange) {
                // Full Range equations
                FR = LTS_ValveSizing_CalculateFRifYes(Ci, Rev, valveSize, FL);
            } else {
                // Restricted Range equations
                FR = LTS_ValveSizing_CalculateFRifNo(Ci, Rev, valveSize, FL);
            }

            // Safety check for invalid FR
            if (!FR || FR <= 0) FR = 1.0; 

            // 4. Calculate Required Cv 
            // The required Cv is the Turbulent Cv divided by the efficiency factor FR
            const C_required = C_base / FR;

            // 5. Convergence Check
            // We check if our guess (Ci) matches the required capacity (C_required)
            if (Math.abs(C_required - Ci) / C_required < EPS) {
                 return { Cv: C_required, Re: Rev, iterations: i+1, choked: isChoked_Base };
            }

            // Update guess for next loop
            // If we aren't there yet, move Ci towards C_required
            Ci = C_required;
        }
        
        // Return best guess if max iterations reached
        return { Cv: Ci, Re: Rev, iterations: MAX_ITERS, choked: isChoked_Base };
    }
}

//------------------------------------------------------------------------------------
// SECTION 8: CALCULATION AND RECOMMENDATION ORCHESTRATOR
//------------------------------------------------------------------------------------
/**
 * LTS_ValveSizing_CalculateOpeningFromCv
 * Inverse lookup: Finds required Valve Opening % for a given Cv.
 * Uses Piecewise Linear Interpolation on the manufacturer's Cv curve.
 * @param {number} requiredCv - The calculated Cv requirement.
 * @param {number[]} cvCurve - 10-point Cv array (10%, 20%... 100%).
 * @returns {number} Valve opening percentage (0-100+).
 * * Purpose:
 * Maps the calculated requirement (Step 4 of sizing) back to a physical valve position.
 * Used to determine if a valve is too big (open < 10%) or too small (open > 90%).
 *
 * Formula & notes:
 * Iterates through the 10-point curve to find the segment [Cv_low, Cv_high]
 * that contains the required Cv, then interpolates the Opening % linearly:
 * Open = Open_low + ( (Cv - Cv_low) / (Cv_high - Cv_low) ) * (Open_high - Open_low)
 *
 * Variables:
 * - fullCvCurve: Prepend 0 to handle 0-10% range.
 * - fullOpeningCurve: [0, 10, 20... 100].
 */
function LTS_ValveSizing_CalculateOpeningFromCv(requiredCv, cvCurve) {
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
 * LTS_ValveSizing_CalculateAndRecommend
 * Orchestrates the full 3-point calculation (Min, Normal, Max) and determines 
 * process requirements based on the WORST CASE scenario.
 * * * Purpose:
 * 1. Resolves 3 scenarios (Min, Normal, Max) from user inputs.
 * 2. Sorts scenarios by Cv demand (Smallest -> Largest) to find true Min/Max loads.
 * 3. Determines critical process limits:
 * - Max Pressure Class (ASME B16.34)
 * - Cavitation risk (Sigma analysis)
 * - Optimal Flow Characteristic (SENAI logic)
 * 4. Iterates through the valve database to find valid candidates.
 * 5. Scores each candidate based on fit (Opening %, Noise, Gain, etc.).
 * * @returns {object} { recommendations: Array, processInfo: Object }
 */
function LTS_ValveSizing_CalculateAndRecommend() 
{
    const inputs = UserInputs.getInputs();
    const SG = UserInputs.getProperty("getSelectedFluid.relativeDensity");
    const Pv = UserInputs.getProperty("getSelectedFluid.Pv"); 
    const N1 = numericalConstants.N1;

    // --- 1. SCENARIO SORTING (Refined Labels & Logic) ---
    // Helper to get a quick Cv estimate for sorting purposes
    const calcSimpleCv = (Q, P1, P2) => {
        const dP = P1 - P2;
        if (dP <= 0 || Q <= 0) return 0;
        return (Q / N1) * Math.sqrt(SG / dP);
    };

    // Define the 3 input scenarios
    let scenarios = [
        { id: "min", label: "Minimum Conditions Input",          q: inputs.flowRateMin, p1: inputs.inletPressureMin, p2: inputs.outletPressureMin, cv: 0 },
        { id: "nor", label: "Normal Operating Conditions Input", q: inputs.flowRate,    p1: inputs.inletPressure,    p2: inputs.outletPressure,    cv: 0 },
        { id: "max", label: "Maximum Conditions Input",          q: inputs.flowRateMax, p1: inputs.inletPressureMax, p2: inputs.outletPressureMax, cv: 0 }
    ];

    // Calculate preliminary Cv for sorting
    scenarios.forEach(s => s.cv = calcSimpleCv(s.q, s.p1, s.p2));

    // Filter out invalid scenarios (zero flow/pressure)
    const validScenarios = scenarios.filter(s => s.cv > 0);
    
    // Defaults (in case of insufficient data)
    let minCase    = scenarios[0]; 
    let opCase     = scenarios[1]; 
    let designCase = scenarios[2]; 
    let sortWarning = null;

    if (validScenarios.length > 0) {
        // Sort: Smallest Cv -> Largest Cv
        validScenarios.sort((a, b) => a.cv - b.cv);
        
        // 1. TURNDOWN CASE (Smallest Cv)
        minCase = validScenarios[0];

        // 2. DESIGN CASE (Largest Cv)
        designCase = validScenarios[validScenarios.length - 1];

        // 3. OPERATING CASE (Middle Cv)
        if (validScenarios.length >= 3) {
            opCase = validScenarios[1]; 
        } else if (validScenarios.length === 2) {
            opCase = validScenarios[1]; 
        } else {
            opCase = validScenarios[0];
        }

        // --- DETECT SWAPS (Clearer Third-Person Language) ---
        let swaps = [];
        
        // Check if Sizing (Max Cv) is coming from somewhere unexpected
        if (designCase.id !== "max") {
            swaps.push(`the <strong>${designCase.label}</strong> required a larger valve size (Cv ${designCase.cv.toFixed(2)}) than the 'Maximum Conditions Input'`);
        }
        
        // Check if Turndown (Min Cv) is coming from somewhere unexpected
        if (minCase.id !== "min" && validScenarios.length > 1) {
            swaps.push(`the <strong>${minCase.label}</strong> resulted in a lower flow requirement (Cv ${minCase.cv.toFixed(2)}) than the 'Minimum Conditions Input'`);
        }

        // Build the final sentence (Third Person)
        if (swaps.length > 0) {
            const detail = swaps.join(", and ");
            sortWarning = `<strong>Automatic Adjustment:</strong> The sizing conditions were switched because ${detail}.`;
        }
    }

    // Assign corrected variables for use downstream
    const baseCvOp  = opCase.cv;
    const baseCvMin = minCase.cv;
    const baseCvMax = designCase.cv; 

    // --- 2. PROCESS ANALYSIS ---
    // Cavitation Check (using Operating Conditions)
    const processSigmaData = LTS_ValveSizing_CalculateCavitationSigma(inputs.inletPressure, inputs.outletPressure, Pv);

    // Pressure Class Check (ASME B16.34)
    const standardClasses = [150, 300, 600, 900, 1500, 2500, 4500];
    let minProcessClass = 4500; 
    const tempInC = UserInputs.getProperty("getInputs.temperature");
    const tempF = (tempInC * 9/5) + 32;

    // Find highest pressure across all scenarios
    const pMaxProcess = Math.max(inputs.inletPressure, inputs.inletPressureMax, inputs.inletPressureMin);
    for (const cls of standardClasses) {
        const check = LTS_ValveSizing_CheckPressureClass(cls, pMaxProcess, tempF);
        if (check.pass) { minProcessClass = cls; break; }
    }

    // --- 3. Characteristic Selection ---
    const processData = {
        qMin: minCase.q, qMax: designCase.q,
        dpMin: designCase.p1 - designCase.p2, 
        dpMax: minCase.p1 - minCase.p2       
    };

    // Get the Controlled Variable
    const controlledVar = UserInputs.getProperty("getInputs.controlledVariable");

    // Get the Recommended Flow Characteristic
    const senaiRec = LTS_ValveSizing_RecommendFlowCharacteristic(processData, controlledVar);

    let preferredChar = senaiRec.type;
    let charReason = senaiRec.reason;
    
    // --- SCORING HELPER FUNCTION ---
    const calculateScore = (valve, opData) => {
        let score = 0;
        const { openOp, noise, rangeability, isPreferredChar } = opData;

        // 1. OPERATING OPENING (Max 30)
        const open = openOp;
        let dist = 0;
        if (open < 60) dist = 60 - open;
        else if (open > 80) dist = open - 80;
        
        if (dist === 0) score += 30; // Ideal range (60-80%)
        else score += Math.max(0, 30 - (dist * 0.6)); // Decay

        // 2. CHARACTERISTIC MATCH (Max 20)  <-- Increased by 5 points
        if (isPreferredChar) score += 20;
        else score += 5; 

        // 3. HYDRODYNAMIC NOISE (Max 15)
        if (noise === 0) score += 0; 
        else if (noise <= 75) score += 15;
        else if (noise >= 85) score += 0;
        else score += 15 - ((noise - 75) * 1.5);

        // 4. RANGEABILITY (Max 10)
        // Keep Max 10, but saturate at 25:1
        if (rangeability >= 25) score += 10;
        else if (rangeability <= 10) score += 0;
        else score += ((rangeability - 10) / 15) * 10; // 15 is the span (25-10)

        // 5. PRESSURE CLASS (Max 15)
        const classes = [150, 300, 600, 900, 1500, 2500, 4500];
        const reqIdx = classes.indexOf(minProcessClass);
        const valIdx = classes.indexOf(valve.pressureClass);
        
        if (valIdx !== -1 && reqIdx !== -1) {
            const diff = valIdx - reqIdx;
            if (diff === 0) score += 15;      // Exact match
            else if (diff === 1) score += 10; // One class higher
            else score += 5;                  // Significantly higher/oversized
        } else score += 15; // Fallback for legacy data

        // 6. VALVE SIZE vs PIPE SIZE (Max 10)      
        if (inputs.pipeDiameter > 0) {
            const dratio = valve.size / inputs.pipeDiameter;

            if (dratio <= 0.5) {
                // Valve is <= 50% of pipe size (Full Score)
                score += 10;
            } else if (dratio < 1.0) {
                // Progressive score between 0.5 and 1.0
                // Formula: (1.0 - ratio) * 20 scales it from 0 to 10 points
                score += (1.0 - dratio) * 20; 
            }
            // If dratio >= 1.0 (Valve >= Pipe), score remains 0
        }
        
        return Math.round(score);
    };

    let recommendations = [];

    // ============================================================
    // BRANCH A: CUSTOM VALVE ANALYSIS
    // ============================================================
    if (inputs.valveType === 99) {
        const customCandidate = {
            manufacturer: "User Defined", model: "Custom Analysis",
            size: inputs.customValveSize || inputs.pipeDiameter, 
            Fl: inputs.customFl, Fd: inputs.customFd,
            Kc: UserInputs.getSelectedValve().Kc || (0.65 * Math.pow(inputs.customFl, 2)),
            pressureClass: "N/A", characteristic: "User Defined"
        };

        // Temporarily override UserInputs for the sizing engine
        UserInputs.setValveProperty("Fl", customCandidate.Fl);
        UserInputs.setValveProperty("Fd", customCandidate.Fd);

        // Run sizing on Corrected Scenarios
        const resOp = LTS_ValveSizing_SizingOrchestrator(customCandidate.size, opCase.q, opCase.p1, opCase.p2);
        let resMin = null, resMax = null;
        try {
            resMin = LTS_ValveSizing_SizingOrchestrator(customCandidate.size, minCase.q, minCase.p1, minCase.p2);
            resMax = LTS_ValveSizing_SizingOrchestrator(customCandidate.size, designCase.q, designCase.p1, designCase.p2);
        } catch(e) {}

        const valveStyleDefaults = UserInputs.getSelectedValve();

        if (resOp) {
            const noiseVal = LTS_ValveSizing_CalculateFinalSPL({
                InletPressure: opCase.p1,
                OutletPressure: opCase.p2,
                Pv: Pv,
                Fl: customCandidate.Fl,
                Kc: customCandidate.Kc || valveStyleDefaults.Kc, 
                Cv: resOp.Cv
            });

            let processTurndown = (resMin && resMax && resMin.Cv > 0) ? (resMax.Cv / resMin.Cv).toFixed(2) : "N/A";
            
            recommendations.push({
                name: "Custom Valve Analysis",
                size: customCandidate.size,
                ratedCv: 0, characteristic: "N/A", pressureClass: "N/A",
                noise: noiseVal.toFixed(2),
                sigma: processSigmaData.sigma, cavSeverity: processSigmaData.severity, cavConsequence: processSigmaData.consequence, cavColor: processSigmaData.color,
                cvOp: resOp.Cv, cvMin: resMin ? resMin.Cv : 0, cvMax: resMax ? resMax.Cv : 0,
                openOp: 0, openMin: 0, openMax: 0, 
                fit: "Custom", score: 100, 
                isPreferredChar: false, turndownPass: true,
                valveRangeability: 0, processTurndown: processTurndown,
                regime: (resOp.choked ? "Choked" : "Non-Choked") + (resOp.Re > 10000 ? " Turbulent" : " Non-Tubulent")
            });
        }
    }
    // ============================================================
    // BRANCH B: STANDARD DATABASE SEARCH
    // ============================================================
    else {
        const originalFl = UserInputs.getProperty("getSelectedValve.Fl");
        const originalFd = UserInputs.getProperty("getSelectedValve.Fd");
        // Ensure valveDatabase is loaded globally
        const candidates = typeof valveDatabase !== 'undefined' ? valveDatabase.filter(valve => valve.type === inputs.valveType) : [];

        for (const originalCandidate of candidates) {
            let candidate = { ...originalCandidate };

            // FILTER: Physical Constraint (Valve cannot be larger than pipe)
            if (candidate.size > inputs.pipeDiameter) continue; 

            // FILTER: Check against baseCvMax (True Design Cv)
            if (candidate.ratedCv < (baseCvMax * 0.90)) continue; 

            // FILTER: Characteristic Constraint (Strict Eq% logic)
            const valveChar = candidate.characteristic || "Linear"; 
            if (preferredChar === "Equal Percentage") {
                if (valveChar !== "Equal Percentage") continue;
            } else {
                if (valveChar !== preferredChar && valveChar !== "Equal Percentage") continue;
            }

            // FILTER: Pressure Class Constraint
            let selectedClass = null;
            let availableClasses = candidate.pressureClasses || (candidate.pressureClass ? [candidate.pressureClass] : []);
            
            availableClasses.sort((a, b) => a - b);
            for (const cls of availableClasses) {
                const check = LTS_ValveSizing_CheckPressureClass(cls, pMaxProcess, tempF);
                if (check.pass) { selectedClass = cls; break; }
            }
            if (selectedClass === null) continue; 
            candidate.pressureClass = selectedClass;

            // Set valve-specific geometry for the sizing engine
            UserInputs.setValveProperty("Fl", candidate.Fl);
            UserInputs.setValveProperty("Fd", candidate.Fd);

            // Run sizing on Corrected Scenarios
            const resOp = LTS_ValveSizing_SizingOrchestrator(candidate.size, opCase.q, opCase.p1, opCase.p2);
            if (!resOp) continue;
            const resMin = LTS_ValveSizing_SizingOrchestrator(candidate.size, minCase.q, minCase.p1, minCase.p2);
            if (!resMin) continue;
            const resMax = LTS_ValveSizing_SizingOrchestrator(candidate.size, designCase.q, designCase.p1, designCase.p2);
            if (!resMax) continue;

            const valveStyleDefaults = UserInputs.getSelectedValve();

            // Calculate Noise
            const noiseVal = LTS_ValveSizing_CalculateFinalSPL({
                InletPressure: opCase.p1,
                OutletPressure: opCase.p2,
                Pv: Pv,
                Fl: candidate.Fl,
                Kc: candidate.Kc || valveStyleDefaults.Kc || (0.65 * Math.pow(candidate.Fl, 2)), 
                Cv: resOp.Cv,
            });

            // Calculate Openings
            const getOpen = (cv, val) => {
                if (val.cvCurve && val.cvCurve.length === 10) return LTS_ValveSizing_CalculateOpeningFromCv(cv, val.cvCurve);
                return (cv / val.ratedCv) * 100;
            };
            const openOp = getOpen(resOp.Cv, candidate);
            const openMin = getOpen(resMin.Cv, candidate);
            const openMax = getOpen(resMax.Cv, candidate);
            const isPreferredChar = (candidate.characteristic === preferredChar);

            // Classify Fit
            let fit = "Marginal";
            if (openMin >= 0 && openMax <= 100) {
                if (openOp >= 60 && openOp <= 80 && openMax <= 95 && openMin > 10 && isPreferredChar) fit = "Ideal";
                else if (openOp >= 20 && openOp <= 80 && openMax <= 95 && openMin > 10) fit = "Acceptable";
            } else fit = "Invalid"; 
            
            if (fit === "Invalid") continue;

            // Turndown Analysis
            const processTurndown = resMax.Cv / resMin.Cv;
            let cvAt90 = candidate.cvCurve ? candidate.cvCurve[8] : candidate.ratedCv * 0.9;
            let cvAt10 = candidate.cvCurve ? candidate.cvCurve[0] : candidate.ratedCv * 0.1;
            const valveRangeability = cvAt90 / cvAt10;
            const turndownPass = (valveRangeability >= processTurndown);

            // Calculate Total Score
            const opData = { openOp, openMin, openMax, noise: noiseVal, rangeability: valveRangeability, isPreferredChar };
            const score = calculateScore(candidate, opData);

            recommendations.push({
                name: `${candidate.manufacturer} ${candidate.model}`,
                size: candidate.size, ratedCv: candidate.ratedCv,
                characteristic: candidate.characteristic || "Linear",
                pressureClass: candidate.pressureClass,
                noise: noiseVal.toFixed(2),
                sigma: processSigmaData.sigma, cavSeverity: processSigmaData.severity, cavConsequence: processSigmaData.consequence, cavColor: processSigmaData.color,
                cvOp: resOp.Cv, cvMin: resMin.Cv, cvMax: resMax.Cv,
                openOp: openOp, openMin: openMin, openMax: openMax,
                fit: fit, score: Math.round(score),
                isPreferredChar: isPreferredChar, turndownPass: turndownPass,
                valveRangeability: valveRangeability, processTurndown: processTurndown.toFixed(2),
                regime: (resOp.choked ? "Choked" : "Non-Choked") + (resOp.Re > 10000 ? " Turbulent" : " Non-Turbulent")
            });
        }

        // Restore global defaults
        UserInputs.setValveProperty("Fl", originalFl);
        UserInputs.setValveProperty("Fd", originalFd);

        // Define weights for the Fit categories
        const fitWeights = { "Ideal": 3, "Acceptable": 2, "Marginal": 1, "Invalid": 0, "Custom": 4 };

        // Refined sorting logic
        recommendations.sort((a, b) => {
            // Turndown Capability (Safety First)
            if (a.turndownPass && !b.turndownPass) return -1;
            if (!a.turndownPass && b.turndownPass) return 1;

            // Qualitative Fit (Ideal vs Marginal)
            const weightA = fitWeights[a.fit] || 0;
            const weightB = fitWeights[b.fit] || 0;
            if (weightA !== weightB) return weightB - weightA;

            // Numerical Score (The tie-breaker)
            return b.score - a.score;
        });
    }

    return { 
        recommendations, 
        processInfo: { 
            preferredChar, charReason, baseCvOp, baseCvMin, baseCvMax, minProcessClass,
            sigma: processSigmaData.sigma, sigmaSeverity: processSigmaData.severity, 
            sigmaColor: processSigmaData.color, sigmaLabel: processSigmaData.consequence,
            sortWarning: sortWarning 
        } 
    };
}

//------------------------------------------------------------------------------------
// SECTION 9: UI INTEGRATION & EVENT LISTENERS
//------------------------------------------------------------------------------------
/**
 * Reads all input values from the HTML document (DOM) and updates the application state.
 * It also triggers necessary fluid property calculations (e.g., density, viscosity).
 */
function LTS_ValveSizing_ParameterUpdate() {
    // --- 1. Read Standard DOM Elements ---
    const inletPressureEl = document.getElementById("inletPressure");
    const outletPressureEl = document.getElementById("outletPressure");
    const temperatureEl = document.getElementById("temperature");
    const flowRateEl = document.getElementById("flowRate");
    const pipeDiameterEl = document.getElementById("pipeDiameter");
    const fluidTypeEl = document.getElementById("fluidType");
    const valveTypeEl = document.getElementById("valveType");
    const controlledVarEl = document.getElementById("controlledVariable");

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
    const controlledVariable = (controlledVarEl ? parseInt(controlledVarEl.value, 10) : 0);

    UserInputs.setInput("inletPressure", inletPressure);
    UserInputs.setInput("outletPressure", outletPressure);
    UserInputs.setInput("temperature", temperature);
    UserInputs.setInput("flowRate", flowRate);
    UserInputs.setInput("pipeDiameter", pipeDiameter);
    UserInputs.setInput("fluidType", fluidType);
    UserInputs.setInput("valveType", valveType);
    UserInputs.setInput("controlledVariable", controlledVariable);
    
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
        const waterVaporPressure = LTS_ValveSizing_EstimateWaterVaporPressure(temperature);
        UserInputs.setFluidProperty("Pv", waterVaporPressure);
        LTS_ValveSizing_CalculateSpecificGravity();
        const estViscosity = LTS_ValveSizing_EstimateWaterKinematicViscosityCSt(temperature);
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
    const resultEl = document.getElementById("result");
    
    // --- FIX #3: Clear the previous analysis immediately ---
    if (resultEl) {
        resultEl.innerHTML = `<div style="color:var(--muted); font-size:14px; padding:12px;">Calculating...</div>`;
    }

    // Show debug log container if debug is enabled
    if (typeof DEBUG !== 'undefined' && DEBUG) {
      const debugContainer = document.getElementById("debugContainer");
      if (debugContainer) debugContainer.style.display = "block";
    }

    // Use setTimeout to yield a fraction of a second so the browser can visually clear the UI
    setTimeout(() => {
        // --- FIX #2: Validate Custom Fluid Inputs ---
        const fluidTypeEl = document.getElementById("fluidType");
        const isCustomFluid = fluidTypeEl && parseInt(fluidTypeEl.value, 10) === UserInputs.Fluids.CUSTOM;
        
        if (isCustomFluid) {
            const pcVal = document.getElementById("Pc")?.value;
            const pvVal = document.getElementById("Pv")?.value;
            const rhoVal = document.getElementById("density")?.value;
            const viscVal = document.getElementById("viscosity")?.value;

            // If any of the strings are empty
            if (!pcVal || !pvVal || !rhoVal || !viscVal) {
                if (resultEl) {
                    resultEl.innerHTML = `
                        <div style="background:rgba(239, 68, 68, 0.1); border-left:3px solid var(--accent-red); padding:12px; border-radius:4px; font-size:14px; color:#fca5a5;">
                            <strong>Missing Custom Fluid Data:</strong><br>
                            <span style="font-size:13px; color:var(--muted);">Please ensure Critical Pressure (Pc), Vapor Pressure (Pv), Density, and Viscosity are all filled out.</span>
                        </div>`;
                }
                return; // Stop execution here
            }
        }

        // 1. Update Parameters & Run Calculation
        LTS_ValveSizing_ParameterUpdate(); 
        const { recommendations, processInfo } = LTS_ValveSizing_CalculateAndRecommend();
        
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

                ${processInfo.sortWarning ? `
                    <div style="background:rgba(234, 179, 8, 0.15); border-left:3px solid #eab308; padding:8px 10px; border-radius:4px; font-size:12px; color:#fef08a; line-height:1.4;">
                        <strong style="color:#facc15;">Note:</strong> ${processInfo.sortWarning}
                    </div>
                ` : ''}

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

                    // Footer (Rangeability (Turndown) + Regime)
                    const rangeText = valve.valveRangeability > 0 ? valve.valveRangeability.toFixed(0) + ":1" : "N/A";
                    footerHtml = `
                        <div style="margin-top:8px; font-size:11px; color:var(--muted); display:flex; justify-content:space-between; align-items:center;">
                            <span style="font-style:italic;" title="Ratio of Cv at 90% to Cv at 10% open">
                                Installed Rangeability: ${rangeText}
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
                                    ${valve.fit !== "Custom" ? `
                                    <span style="margin-left: 4px; font-size:10px; color:#f3f4f6; font-weight:700; background:rgba(255,255,255,0.05); padding:2px 6px; border-radius:4px; text-transform:uppercase; letter-spacing:0.5px;">
                                        Rated Cv: ${Number(valve.ratedCv).toFixed(1)}
                                    </span>` : ''}
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

    }, 20); // 20ms timeout is enough for the browser to render the "Calculating..." text
  });
}

//------------------------------------------------------------------------------------
// SECTION 10: AUXILIARY UI INTEGRATION & EVENT LISTENERS
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
    document.getElementById("flowRate").value = 535;
    document.getElementById("inletPressure").value = 85;
    document.getElementById("outletPressure").value = 30;
    document.getElementById("pipeDiameter").value = 4;
    document.getElementById("temperature").value = 77;
    
    // Reset Selectors (Back to defaults)
    document.getElementById("fluidType").value = 0; // Reset selection to Water
    document.getElementById("valveType").value = 0; // Reset selection to Globe
    document.getElementById("customValveStyle").value = 0; // Reset Custom Style
    document.getElementById("controlledVariable").value = 0; //Reset the Controlled Variable

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
    LTS_ValveSizing_ParameterUpdate();
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
        document.getElementById("flowRate").value = 10;
        document.getElementById("inletPressure").value = 100;
        document.getElementById("outletPressure").value = 90; // dp = 10 psi
        document.getElementById("pipeDiameter").value = 1;
        document.getElementById("temperature").value = 68; // 20C
    } else if (preset === 2) {
        document.getElementById("flowRate").value = 535;
        document.getElementById("inletPressure").value = 85;
        document.getElementById("outletPressure").value = 30; // dp = 1 psi
        document.getElementById("pipeDiameter").value = 4;
        document.getElementById("temperature").value = 77; // 20C
    } else {
        document.getElementById("flowRate").value = 390;
        document.getElementById("inletPressure").value = 20;
        document.getElementById("outletPressure").value = 10; // dp = 100 psi
        document.getElementById("pipeDiameter").value = 4;
        document.getElementById("temperature").value = 68; // 20C
    }
    // LTS_ValveSizing_ParameterUpdate(); // This is called by the click()
    document.getElementById("calculateBtn").click();
}

document.getElementById('test1').addEventListener('click', () => runPreset(1));
document.getElementById('test2').addEventListener('click', () => runPreset(2));
document.getElementById('test3').addEventListener('click', () => runPreset(3));



