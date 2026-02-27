// This variable is now globally accessible
const valveDatabase = [
    // =========================================================================
    // FISHER CAV4 (Cavitrol IV) - CL1500 - GLOBE - Flow Down
    // Source: Catalog 12, Page CAV4-1
    // Characteristic: Pure Linear
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "CAV4 Cavitrol™ IV Trim",
        type: 1, // GLOBECLOSE (Flow Down)
        size: 2,
        ratedCv: 1.6,
        Fl: 0.99, // Very high recovery (Anti-Cavitation)
        Fd: 1.00, // Standard Globe approximation
        pressureClasses: [1500], // Specifically listed for CL1500
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.1, 1.3, 1.5, 1.6]
    },
    // =========================================================================
    // FISHER CHP - CL2500 - FLOW DOWN (Standard Trims)
    // Source: Fisher Catalog 12, Page 79 (CHP-1)
    // Model Name: "CHP" (Standard)
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "CHP", 
        type: 1, // GLOBECLOSE (Flow Down)
        size: 8,
        ratedCv: 519,
        Fl: 0.86,
        Fd: 1.0, 
        pressureClasses: [2500],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [57.9, 134, 205, 262, 308, 363, 417, 458, 500, 519]
    },
    {
        manufacturer: "Fisher",
        model: "CHP",
        type: 1, // GLOBECLOSE
        size: 8,
        ratedCv: 514,
        Fl: 0.86,
        Fd: 1.0,
        pressureClasses: [2500],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [19.4, 44.1, 68.0, 110, 196, 307, 409, 449, 487, 514]
    },
    {
        manufacturer: "Fisher",
        model: "CHP",
        type: 1, // GLOBECLOSE
        size: 8,
        ratedCv: 381,
        Fl: 0.84,
        Fd: 1.0,
        pressureClasses: [2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [11.5, 27.9, 44.1, 60.6, 75.4, 110, 164, 236, 307, 381]
    },
    // =========================================================================
    // FISHER CHP - CL2500 - CAVITROL III TRIM - FLOW DOWN
    // Source: Fisher Catalog 12, Page 79 (CHP-1)
    // Model Name: "CHP Cavitrol III" (Special Trim)
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "CHP Cavitrol III",
        type: 1, // GLOBECLOSE
        size: 8,
        ratedCv: 169, 
        Fl: 0.98, // High recovery for anti-cavitation
        Fd: 1.0, 
        pressureClasses: [2500],
        characteristic: "Linear", 
        flowDirection: "Flow Down",
        cvCurve: [9.56, 28.4, 47.0, 65.6, 83.8, 102, 119, 136, 153, 169]
    },
    // =========================================================================
    // FISHER CHP - CL2500 - FLOW UP (Standard Trims)
    // Source: Fisher Catalog 12, Page 80 (CHP-2)
    // Model Name: "CHP" (Standard)
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "CHP",
        type: 0, // GLOBEOPEN (Flow Up)
        size: 8,
        ratedCv: 482,
        Fl: 0.84, 
        Fd: 0.46, // Standard Fd for Flow Up
        pressureClasses: [2500],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [58.7, 133, 206, 267, 313, 343, 382, 421, 456, 482]
    },
    {
        manufacturer: "Fisher",
        model: "CHP",
        type: 0, // GLOBEOPEN
        size: 8,
        ratedCv: 489,
        Fl: 0.87,
        Fd: 0.46,
        pressureClasses: [2500],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [20.3, 45.9, 72.5, 115, 200, 301, 408, 459, 487, 489]
    },
    {
        manufacturer: "Fisher",
        model: "CHP",
        type: 0, // GLOBEOPEN
        size: 8,
        ratedCv: 375, 
        Fl: 0.80,
        Fd: 0.46,
        pressureClasses: [2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [11.3, 28.5, 45.9, 63.5, 82.7, 115, 169, 237, 301, 375]
    },
    // =========================================================================
    // FISHER CV500 (V-Notch Ball) - FORWARD FLOW
    // Source: Catalog 12, Page 85 (CV500-1)
    // Characteristic: Modified Equal Percentage
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "CV500 (Forward)",
        type: 3, // BALL
        size: 3,
        ratedCv: 166,
        Fl: 0.69,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Forward",
        cvCurve: [4.74, 14.1, 34.6, 60.1, 84.0, 107, 133, 163, 166, 166]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Forward)",
        type: 3, 
        size: 4,
        ratedCv: 346,
        Fl: 0.62,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Forward",
        cvCurve: [11.1, 27.1, 61.7, 106, 149, 193, 252, 324, 346, 346]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Forward)",
        type: 3, 
        size: 6,
        ratedCv: 809,
        Fl: 0.57,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Forward",
        cvCurve: [15.7, 33.0, 86.1, 154, 229, 330, 497, 718, 809, 809]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Forward)",
        type: 3, 
        size: 8,
        ratedCv: 1440,
        Fl: 0.58,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Forward",
        cvCurve: [21.5, 82.4, 156, 259, 402, 592, 832, 1120, 1440, 1440]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Forward)",
        type: 3, 
        size: 10,
        ratedCv: 2360,
        Fl: 0.54,
        Fd: 1.00,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Forward",
        cvCurve: [41.4, 162, 301, 455, 699, 995, 1300, 1820, 2360, 2360]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Forward)",
        type: 3, 
        size: 12,
        ratedCv: 3050,
        Fl: 0.51,
        Fd: 1.00,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Forward",
        cvCurve: [60.4, 215, 443, 699, 1020, 1390, 1850, 2560, 3050, 3050]
    },
    // =========================================================================
    // FISHER CV500 (V-Notch Ball) - REVERSE FLOW
    // Source: Catalog 12, Page 86 (CV500-2)
    // Characteristic: Modified Equal Percentage
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "CV500 (Reverse)",
        type: 3, // BALL
        size: 3,
        ratedCv: 181,
        Fl: 0.53,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Reverse",
        // Note: Table shows 181 at both 80 and 90 degrees
        cvCurve: [3.25, 14.2, 34.2, 61.8, 94.5, 129, 160, 181, 181, 181]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Reverse)",
        type: 3, 
        size: 4,
        ratedCv: 300,
        Fl: 0.61,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Reverse",
        cvCurve: [7.20, 27.2, 64.8, 116, 172, 223, 263, 290, 300, 300]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Reverse)",
        type: 3, 
        size: 6,
        ratedCv: 808,
        Fl: 0.49,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Reverse",
        cvCurve: [5.20, 33.3, 88.5, 170, 268, 372, 476, 600, 808, 808]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Reverse)",
        type: 3, 
        size: 8,
        ratedCv: 1240,
        Fl: 0.58,
        Fd: 0.99,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Reverse",
        cvCurve: [8.68, 61.1, 156, 293, 463, 656, 856, 1050, 1240, 1240]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Reverse)",
        type: 3, 
        size: 10,
        ratedCv: 2140,
        Fl: 0.49,
        Fd: 1.00,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Reverse",
        cvCurve: [37.0, 137, 288, 505, 752, 1080, 1460, 1710, 2140, 2140]
    },
    {
        manufacturer: "Fisher",
        model: "CV500 (Reverse)",
        type: 3, 
        size: 12,
        ratedCv: 3080,
        Fl: 0.50,
        Fd: 1.00,
        pressureClasses: [150, 300, 600],
        characteristic: "Modified Equal Percentage",
        flowDirection: "Reverse",
        cvCurve: [39.0, 192, 411, 703, 1090, 1560, 2040, 2490, 3080, 3080]
    },
    // =========================================================================
    // FISHER DESIGN D - MICRO-FORM (Flow Up)
    // Source: Catalog 12, Page 87 (D-1)
    // Fits: NPS 3/4, 1, and 2 (Listed here as Size 1)
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 1/4\")",
        type: 0, // GLOBEOPEN
        size: 1, 
        ratedCv: 1.66,
        Fl: 0.87,
        Fd: 0.68, 
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.070, 0.115, 0.164, 0.224, 0.315, 0.450, 0.641, 0.921, 1.28, 1.66]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 3/8\")",
        type: 0, 
        size: 1,
        ratedCv: 4.03,
        Fl: 0.84,
        Fd: 0.61,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.155, 0.260, 0.407, 0.596, 0.858, 1.21, 1.65, 2.22, 3.00, 4.03]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 1/2\")",
        type: 0, 
        size: 1,
        ratedCv: 6.51,
        Fl: 0.84,
        Fd: 0.56,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.273, 0.436, 0.631, 0.911, 1.30, 1.84, 2.57, 3.65, 5.08, 6.51]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 3/4\")",
        type: 0, 
        size: 1,
        ratedCv: 12.3,
        Fl: 0.92,
        Fd: 0.49,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.483, 0.775, 1.25, 1.97, 2.89, 4.13, 5.87, 8.16, 10.9, 12.3]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 1/4\")",
        type: 0, // GLOBEOPEN
        size: 2, 
        ratedCv: 1.66,
        Fl: 0.87,
        Fd: 0.68, 
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.070, 0.115, 0.164, 0.224, 0.315, 0.450, 0.641, 0.921, 1.28, 1.66]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 3/8\")",
        type: 0, 
        size: 2,
        ratedCv: 4.03,
        Fl: 0.84,
        Fd: 0.61,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.155, 0.260, 0.407, 0.596, 0.858, 1.21, 1.65, 2.22, 3.00, 4.03]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 1/2\")",
        type: 0, 
        size: 2,
        ratedCv: 6.82,
        Fl: 0.81,
        Fd: 0.56,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.348, 0.505, 0.709, 0.998, 1.38, 1.92, 2.69, 3.82, 5.25, 6.82]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 3/4\")",
        type: 0, 
        size: 2,
        ratedCv: 14.1,
        Fl: 0.81,
        Fd: 0.49,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.613, 0.952, 1.44, 2.06, 2.92, 4.13, 5.87, 8.16, 11.1, 14.1]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 1\")",
        type: 0, 
        size: 2,
        ratedCv: 23.7,
        Fl: 0.82,
        Fd: 0.46,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.2, 1.68, 2.44, 3.53, 5.05, 7.28, 10.5, 14, 18.4, 23.7]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Form 1-1/4\")",
        type: 0, 
        size: 2,
        ratedCv: 34.5,
        Fl: 0.85,
        Fd: 0.44,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.32, 1.76, 2.50, 3.66, 5.42, 8.25, 12.7, 20.6, 29, 34.5]
    },
    // =========================================================================
    // FISHER DESIGN D - MICRO-FLUTE (Flow Up)
    // Source: Catalog 12, Page 88 (D-2)
    // Fits: NPS 1 and 2 (Listed here as Size 1)
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "D (Micro-Flute 1-Flute)",
        type: 0, // GLOBEOPEN
        size: 1,
        ratedCv: 0.354,
        Fl: 0.87,
        Fd: 0.46,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.0385, 0.0455, 0.0560, 0.0719, 0.0942, 0.124, 0.162, 0.212, 0.278, 0.354]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Flute 3-Flutes)",
        type: 0, 
        size: 1,
        ratedCv: 1.07,
        Fl: 0.90,
        Fd: 0.46,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.0562, 0.0725, 0.101, 0.146, 0.216, 0.312, 0.433, 0.588, 0.802, 1.07]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Flute 1-Flute)",
        type: 0, // GLOBEOPEN
        size: 2,
        ratedCv: 0.354,
        Fl: 0.87,
        Fd: 0.46,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.0385, 0.0455, 0.0560, 0.0719, 0.0942, 0.124, 0.162, 0.212, 0.278, 0.354]
    },
    {
        manufacturer: "Fisher",
        model: "D (Micro-Flute 3-Flutes)",
        type: 0, 
        size: 2,
        ratedCv: 1.07,
        Fl: 0.90,
        Fd: 0.46,
        pressureClasses: [150, 300, 600, 900, 1500, 2500],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.0562, 0.0725, 0.101, 0.146, 0.216, 0.312, 0.433, 0.588, 0.802, 1.07]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/8\")",
        type: 0, // GLOBEOPEN
        size: 1, // NPS 1 Body
        ratedCv: 3.70,
        Fl: 0.94,
        Fd: 0.46, 
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.266, 0.480, 0.698, 0.984, 1.31, 1.69, 2.10, 2.63, 3.23, 3.70]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/4\")",
        type: 0, 
        size: 1, 
        ratedCv: 10.4,
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.708, 1.34, 1.97, 2.60, 3.58, 4.95, 6.44, 7.95, 9.30, 10.4]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 1\")",
        type: 0, 
        size: 1, 
        ratedCv: 13.6,
        Fl: 0.86,
        Fd: 0.46,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.38, 2.70, 3.98, 5.30, 6.80, 8.00, 9.60, 11.2, 12.5, 13.6]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/8\")",
        type: 1, // GLOBECLOSE
        size: 1, 
        ratedCv: 5.00,
        Fl: 0.59, 
        Fd: 1.0, 
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [0.600, 1.02, 1.40, 1.90, 2.37, 2.74, 2.90, 3.78, 4.58, 5.00]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/4\")",
        type: 1, 
        size: 1, 
        ratedCv: 9.07,
        Fl: 0.91,
        Fd: 1.0,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [1.31, 2.30, 3.20, 2.90, 3.50, 5.10, 6.58, 7.60, 8.50, 9.07]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 1\")",
        type: 1, 
        size: 1, 
        ratedCv: 11.5,
        Fl: 0.90,
        Fd: 1.0,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [1.50, 2.90, 4.85, 6.58, 7.25, 8.00, 8.70, 9.84, 10.8, 11.5]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/8\")",
        type: 0, // GLOBEOPEN
        size: 2,
        ratedCv: 3.74,
        Fl: 0.95,
        Fd: 0.46, 
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.27, 0.48, 0.69, 0.99, 1.32, 1.7, 2.13, 2.66, 3.29, 3.74]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/4\")",
        type: 0, 
        size: 2, 
        ratedCv: 11.5, // Higher than NPS 1 equivalent (10.4)
        Fl: 0.92,
        Fd: 0.46,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.653, 1.28, 1.9, 2.55, 3.5, 4.86, 6.58, 8.4, 10.1, 11.5]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 1\")",
        type: 0, 
        size: 2, 
        ratedCv: 16.8, // Significantly higher than NPS 1 equivalent (13.6)
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.71, 3.03, 4.54, 6, 7.35, 9.1, 11.1, 13.2, 15.3, 16.8]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/8\")",
        type: 1, // GLOBECLOSE
        size: 2, 
        ratedCv: 4.7, 
        Fl: 0.62, 
        Fd: 1.0, 
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [0.48, 0.9, 1.3, 1.75, 2.2, 2.5, 2.85, 3.5, 4.4, 4.7]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 3/4\")",
        type: 1, 
        size: 2, 
        ratedCv: 10.4,
        Fl: 0.91,
        Fd: 1.0,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [1.08, 2, 3, 3.2, 3.50, 5, 6.75, 8.26, 9.4, 10.4]
    },
    {
        manufacturer: "Fisher",
        model: "D3 (Micro-Form 1\")",
        type: 1, 
        size: 2, 
        ratedCv: 13.8,
        Fl: 0.86,
        Fd: 1.0,
        pressureClasses: [600, 900],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [1.55, 3.4, 5.35, 6.85, 7.9, 8.7, 10.3, 11.9, 12.9, 13.8]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, // GLOBECLOSE (Flow Down)
        size: 1,
        ratedCv: 22.1,
        Fl: 0.81,
        Fd: 0.36, 
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [4.86, 9.39, 13.4, 16.9, 18.9, 20.3, 21.1, 21.8, 21.9, 22.1]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, // GLOBECLOSE (Flow Down)
        size: 1.25,
        ratedCv: 22.1,
        Fl: 0.81,
        Fd: 0.36, 
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [4.86, 9.39, 13.4, 16.9, 18.9, 20.3, 21.1, 21.8, 21.9, 22.1]
    },
    {
        manufacturer: "Fisher",
        model: "ED (Port 1-7/8\")",
        type: 1, 
        size: 1.5,
        ratedCv: 44.0,
        Fl: 0.79,
        Fd: 0.36,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [7.79, 14.4, 20.5, 26.8, 32.0, 36.6, 39.4, 41.3, 42.7, 44.0]
    },
    {
        manufacturer: "Fisher",
        model: "ED (Port 1-5/16\")",
        type: 1, 
        size: 1.5,
        ratedCv: 29.9,
        Fl: 0.88,
        Fd: 0.36,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [5.05, 9.99, 14.7, 20, 24, 25.7, 26.2, 27.4, 28.6, 29.9]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 2,
        ratedCv: 77.6,
        Fl: 0.77,
        Fd: 0.36,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [13.4, 26.8, 39.9, 51.3, 62.9, 70.6, 73.7, 75.6, 76.8, 77.6]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 2.5,
        ratedCv: 109,
        Fl: 0.81,
        Fd: 0.36,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [20.9, 39.6, 58.8, 74.2, 84.9, 97.0, 103, 106, 108, 109]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 3,
        ratedCv: 161,
        Fl: 0.77,
        Fd: 0.36,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [27.2, 52.2, 77.9, 99.5, 124, 140, 149, 154, 158, 161]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 4,
        ratedCv: 251,
        Fl: 0.79,
        Fd: 0.30,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [37.7, 75.0, 125, 163, 193, 220, 238, 247, 251, 251]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 6,
        ratedCv: 460, // Standard 7" Port
        Fl: 0.82,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [73.6, 150, 232, 306, 353, 389, 416, 441, 451, 460]
    },
    {
        manufacturer: "Fisher",
        model: "ED (Port 4-3/8\")", // Restricted Trim
        type: 1, 
        size: 6,
        ratedCv: 358, 
        Fl: 0.87,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [52.3, 101, 150, 199, 247, 284, 310, 329, 345, 358]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 8,
        ratedCv: 744, // 2" Travel
        Fl: 0.87,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [80.3, 188, 290, 389, 480, 554, 615, 658, 705, 744]
    },
    {
        manufacturer: "Fisher",
        model: "ED", // 3" Travel
        type: 1, 
        size: 8,
        ratedCv: 863, 
        Fl: 0.85,
        Fd: 0.27,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Down",
        cvCurve: [135, 291, 434, 551, 639, 706, 759, 807, 841, 863]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 1,
        ratedCv: 20.6,
        Fl: 0.84,
        Fd: 0.34, 
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [3.21, 5.50, 8.18, 10.9, 13.2, 15.0, 16.9, 18.6, 19.9, 20.6]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 1.5,
        ratedCv: 39.2,
        Fl: 0.82,
        Fd: 0.34,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [4.23, 7.84, 11.8, 15.8, 20.4, 25.3, 30.3, 34.7, 37.2, 39.2]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 2,
        ratedCv: 72.9,
        Fl: 0.77,
        Fd: 0.33,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [7.87, 16.0, 24.9, 33.4, 42.1, 51.8, 62.0, 68.1, 70.6, 72.9]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 2.5,
        ratedCv: 108,
        Fl: 0.81,
        Fd: 0.27,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [9.34, 21.6, 35.5, 49.5, 62.7, 74.1, 83.6, 93.5, 102, 108]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 3,
        ratedCv: 148,
        Fl: 0.82,
        Fd: 0.30,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [14.5, 32.9, 52.1, 70.4, 88.5, 105, 118, 133, 142, 148]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 4,
        ratedCv: 236,
        Fl: 0.82,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [23.3, 50.3, 78.1, 105, 127, 152, 181, 203, 223, 236]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 6,
        ratedCv: 433,
        Fl: 0.84,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [46.3, 107, 171, 228, 279, 327, 367, 402, 420, 433]
    },
    {
        manufacturer: "Fisher",
        model: "ED (Port 4-3/8\")", // Restricted Trim
        type: 1, 
        size: 6,
        ratedCv: 322,
        Fl: 0.88,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [16.7, 38.6, 65.4, 93.7, 123, 156, 194, 244, 290, 322]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 8,
        ratedCv: 688, // 2" Travel
        Fl: 0.87,
        Fd: 0.31,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [60.2, 129, 206, 285, 363, 444, 526, 581, 640, 688]
    },
    {
        manufacturer: "Fisher",
        model: "ED", // 3" Travel
        type: 1, 
        size: 8,
        ratedCv: 846, 
        Fl: 0.87,
        Fd: 0.31,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [91.4, 207, 325, 440, 550, 639, 711, 760, 795, 846]
    },
{
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 1,
        ratedCv: 17.2,
        Fl: 0.88,
        Fd: 0.38, 
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [0.783, 1.54, 2.20, 2.89, 4.21, 5.76, 7.83, 10.9, 14.1, 17.2]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 1.5,
        ratedCv: 35.8,
        Fl: 0.84,
        Fd: 0.38,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [1.52, 2.63, 3.87, 5.41, 7.45, 11.2, 17.4, 24.5, 30.8, 35.8]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 2,
        ratedCv: 59.7,
        Fl: 0.85,
        Fd: 0.31,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [1.66, 2.93, 4.66, 6.98, 10.8, 16.5, 25.4, 37.3, 50.7, 59.7]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 2.5,
        ratedCv: 99.4,
        Fl: 0.84,
        Fd: 0.27,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [3.43, 7.13, 10.8, 15.1, 22.4, 33.7, 49.2, 71.1, 89.5, 99.4]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 3,
        ratedCv: 136,
        Fl: 0.82,
        Fd: 0.32,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [4.32, 7.53, 10.9, 17.1, 27.2, 43.5, 66.0, 97.0, 120, 136]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 4,
        ratedCv: 224,
        Fl: 0.82,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [5.85, 11.6, 18.3, 30.2, 49.7, 79.7, 125, 171, 205, 224]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 6,
        ratedCv: 394,
        Fl: 0.85,
        Fd: 0.26,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [12.9, 25.8, 43.3, 67.4, 104, 162, 239, 316, 368, 394]
    },
    {
        manufacturer: "Fisher",
        model: "ED", // Restricted Trim
        type: 1, 
        size: 6,
        ratedCv: 274,
        Fl: 0.88,
        Fd: 0.26,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [5.40, 10.1, 15.8, 26.7, 45.2, 71.2, 111, 169, 232, 274]
    },
    {
        manufacturer: "Fisher",
        model: "ED",
        type: 1, 
        size: 8,
        ratedCv: 567, // 2" Travel
        Fl: 0.85,
        Fd: 0.26,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [18.5, 38.0, 58.4, 86.7, 130, 189, 268, 371, 476, 567]
    },
    {
        manufacturer: "Fisher",
        model: "ED", // 3" Travel
        type: 1, 
        size: 8,
        ratedCv: 818,
        Fl: 0.86,
        Fd: 0.26,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Down",
        cvCurve: [27.0, 58.1, 105, 188, 307, 478, 605, 695, 761, 818]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, // GLOBECLOSE
        size: 12,
        ratedCv: 1570,
        Fl: 0.88,
        Fd: 0.28, // Standard approximation (not listed in table for Size 12)
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        // Curve: Min, 10%, 20% ... 100% [cite: 4]
        cvCurve: [40, 206, 415, 630, 852, 1079, 1295, 1465, 1557, 1570]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 12,
        ratedCv: 1560,
        Fl: 0.88,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        // Curve: Min, 10% ... 100% [cite: 4]
        cvCurve: [40, 106, 222, 368, 541, 727, 922, 1109, 1283, 1437, 1560]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 14,
        ratedCv: 1860,
        Fl: 0.88,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [40, 206, 415, 629, 838, 1110, 1420, 1630, 1760, 1850, 1860]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 14,
        ratedCv: 1601, // Note: Lower Cv at higher travel due to cage window/plug interaction
        Fl: 0.88,
        Fd: 0.28,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [40, 73, 169, 314, 495, 698, 909, 1115, 1305, 1468, 1601]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 16,
        ratedCv: 2772,
        Fl: 0.85,
        Fd: 0.30, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [76, 202, 456, 709, 963, 1232, 1534, 1900, 2234, 2526, 2772]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 16,
        ratedCv: 3312,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [124, 298, 649, 1000, 1376, 1846, 2324, 2701, 2988, 3188, 3312]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 18,
        ratedCv: 2810,
        Fl: 0.85,
        Fd: 0.21, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [82, 218, 492, 765, 1038, 1310, 1622, 1989, 2311, 2583, 2810]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 18,
        ratedCv: 3351,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [133, 322, 701, 1078, 1460, 1935, 2396, 2745, 3010, 3208, 3351]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 16,
        ratedCv: 3253,
        Fl: 0.85,
        Fd: 0.21, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [80, 213, 480, 747, 1000, 1283, 1627, 2065, 2490, 2888, 3253]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 18,
        ratedCv: 3253,
        Fl: 0.85,
        Fd: 0.21, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [80, 213, 480, 747, 1000, 1283, 1627, 2065, 2490, 2888, 3253]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 20,
        ratedCv: 3253,
        Fl: 0.85,
        Fd: 0.21, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [80, 213, 480, 747, 1000, 1283, 1627, 2065, 2490, 2888, 3253]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 24,
        ratedCv: 3253,
        Fl: 0.85,
        Fd: 0.21, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [80, 213, 480, 747, 1000, 1283, 1627, 2065, 2490, 2888, 3253]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 16,
        ratedCv: 3721,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [130, 314, 683, 1079, 1472, 1984, 2514, 2944, 3286, 3543, 3721]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 18,
        ratedCv: 3721,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [130, 314, 683, 1079, 1472, 1984, 2514, 2944, 3286, 3543, 3721]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 20,
        ratedCv: 3721,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [130, 314, 683, 1079, 1472, 1984, 2514, 2944, 3286, 3543, 3721]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 16,
        ratedCv: 4265,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [130, 314, 683, 1041, 1445, 1998, 2609, 3145, 3606, 3982, 4265]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 18,
        ratedCv: 4265,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [130, 314, 683, 1041, 1445, 1998, 2609, 3145, 3606, 3982, 4265]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 20,
        ratedCv: 4265,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [130, 314, 683, 1041, 1445, 1998, 2609, 3145, 3606, 3982, 4265]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 24,
        ratedCv: 4265,
        Fl: 0.85,
        Fd: 0.31, // [cite: 4]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [130, 314, 683, 1041, 1445, 1998, 2609, 3145, 3606, 3982, 4265]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 20,
        ratedCv: 5196,
        Fl: 0.85,
        Fd: 0.30, // [cite: 11]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        // Note: Curve has 5196 at both 90% and 100% [cite: 11]
        cvCurve: [270, 647, 1377, 2089, 2720, 3274, 3747, 4309, 4787, 5196, 5196]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 20,
        ratedCv: 5841,
        Fl: 0.85,
        Fd: 0.30, // [cite: 11]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [271, 649, 1394, 2123, 2827, 3455, 4010, 4589, 5157, 5626, 5841]
    },
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 24,
        ratedCv: 5841,
        Fl: 0.85,
        Fd: 0.30, // [cite: 11]
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [271, 649, 1394, 2123, 2827, 3455, 4010, 4589, 5157, 5626, 5841]
    },    
    {
        manufacturer: "Fisher",
        model: "Large ED",
        type: 1, 
        size: 30,
        ratedCv: 9530,
        Fl: 0.99,
        Fd: 0.3,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Down",
        cvCurve: [100, 906, 2000, 3080, 4230, 5290, 6690, 7710, 8450, 9260, 9530]
    },
    // =========================================================================
    // FISHER ES
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "ES 1\"",
        type: 0, // GLOBEOPEN
        size: 1,
        ratedCv: 20.1,
        Fl: 0.89,
        Fd: 0.46, // Standard Flow-Up Default
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [2.27, 4.12, 6.23, 8.54, 11.0, 13.4, 15.8, 17.8, 19.3, 20.1]
    },
    {
        manufacturer: "Fisher",
        model: "ES 1.1/4\"",
        type: 0, // GLOBEOPEN
        size: 1.25,
        ratedCv: 20.1,
        Fl: 0.89,
        Fd: 0.46, // Standard Flow-Up Default
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [2.27, 4.12, 6.23, 8.54, 11.0, 13.4, 15.8, 17.8, 19.3, 20.1]
    },
    {
        manufacturer: "Fisher",
        model: "ES 1.1/2\"",
        type: 0, 
        size: 1.5,
        ratedCv: 34.9,
        Fl: 0.92,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [3.56, 7.01, 11.1, 15.1, 19.0, 22.9, 26.7, 30.0, 33.1, 34.9]
    },
    {
        manufacturer: "Fisher",
        model: "ES 1.1/2\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 1.5,
        ratedCv: 26.9,
        Fl: 0.95,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [2.42, 4.30, 6.40, 8.77, 11.5, 14.6, 17.8, 21.1, 24.3, 26.9]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2\"",
        type: 0, 
        size: 2,
        ratedCv: 65.3,
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [8.49, 17.1, 25.9, 35.3, 44.4, 52.9, 59.2, 62.0, 63.9, 65.3]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 2,
        ratedCv: 30.9,
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [2.22, 4.11, 6.06, 8.25, 11.0, 14.3, 18.0, 21.8, 26.0, 30.9]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2.1/2\"",
        type: 0, 
        size: 2.5,
        ratedCv: 86.5,
        Fl: 0.93,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [10.4, 22.2, 34.9, 47.1, 58.2, 66.6, 73.7, 79.3, 84.4, 86.5]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2.1/2\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 2.5,
        ratedCv: 48.6,
        Fl: 0.93,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [3.50, 6.85, 10.8, 14.8, 18.9, 23.3, 28.2, 34.1, 41.1, 48.6]
    },
    {
        manufacturer: "Fisher",
        model: "ES 3\"",
        type: 0, 
        size: 3,
        ratedCv: 135,
        Fl: 0.89,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [15.3, 34.3, 52.8, 71.4, 87.8, 101, 112, 121, 129, 135]
    },
    {
        manufacturer: "Fisher",
        model: "ES 3\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 3,
        ratedCv: 88.8,
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [6.39, 13.0, 20.7, 29.1, 38.2, 47.9, 58.0, 68.4, 79.3, 88.8]
    },
    {
        manufacturer: "Fisher",
        model: "ES 4\"",
        type: 0, 
        size: 4,
        ratedCv: 212,
        Fl: 0.89,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [23.7, 46.4, 72.9, 98.2, 122, 145, 165, 183, 199, 212]
    },
    {
        manufacturer: "Fisher",
        model: "ES 4\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 4,
        ratedCv: 139,
        Fl: 0.93,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [10.6, 22.5, 35.0, 47.5, 60.2, 73.1, 88.0, 103, 120, 139]
    },
    {
        manufacturer: "Fisher",
        model: "ES 6\"",
        type: 0, 
        size: 6,
        ratedCv: 417,
        Fl: 0.81,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [55.0, 118, 180, 235, 280, 312, 341, 368, 390, 417]
    },
    {
        manufacturer: "Fisher",
        model: "ES 6\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 6,
        ratedCv: 271,
        Fl: 0.89,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [15.7, 35.8, 60.2, 86.2, 115, 146, 179, 215, 247, 271]
    },
    {
        manufacturer: "Fisher",
        model: "ES 8\"",
        type: 0, 
        size: 8,
        ratedCv: 701,
        Fl: 0.84,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [66.6, 147, 221, 292, 375, 450, 522, 592, 652, 701]
    },
    {
        manufacturer: "Fisher",
        model: "ES 8\"",
        type: 0, 
        size: 8,
        ratedCv: 836,
        Fl: 0.85,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [100, 213, 330, 451, 553, 648, 719, 773, 809, 836]
    },
    // =========================================================================
    // FISHER ES - EQUAL PERCENTAGE
    // =========================================================================
    {
        manufacturer: "Fisher",
        model: "ES 1\"",
        type: 0, // GLOBEOPEN
        size: 1,
        ratedCv: 17.4,
        Fl: 0.95,
        Fd: 0.46, 
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.783, 1.29, 1.86, 2.71, 4.18, 6.44, 9.54, 13.1, 15.7, 17.4]
    },
    {
        manufacturer: "Fisher",
        model: "ES 1.1/4\"",
        type: 0, // GLOBEOPEN
        size: 1.25,
        ratedCv: 17.4,
        Fl: 0.95,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.783, 1.29, 1.86, 2.71, 4.18, 6.44, 9.54, 13.1, 15.7, 17.4]
    },
    {
        manufacturer: "Fisher",
        model: "ES 1.1/2\"",
        type: 0, 
        size: 1.5,
        ratedCv: 33.4,
        Fl: 0.94,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.54, 2.52, 3.57, 4.94, 7.41, 11.6, 17.2, 23.5, 28.7, 33.4]
    },
    {
        manufacturer: "Fisher",
        model: "ES 1.1/2\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 1.5,
        ratedCv: 21.0,
        Fl: 0.96,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.882, 1.35, 1.89, 2.52, 3.68, 5.52, 8.13, 12.0, 16.6, 21.0]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2\"",
        type: 0, 
        size: 2,
        ratedCv: 56.2,
        Fl: 0.92,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.74, 3.15, 4.72, 6.91, 10.6, 16.3, 25.0, 36.7, 47.8, 56.2]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 2,
        ratedCv: 20.8,
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.849, 1.34, 1.83, 2.39, 3.43, 5.12, 7.49, 11.2, 15.8, 20.8]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2.1/2\"",
        type: 0, 
        size: 2.5,
        ratedCv: 82.7,
        Fl: 0.93,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [4.05, 7.19, 10.6, 14.5, 21.2, 31.6, 45.5, 64.2, 77.7, 82.7]
    },
    {
        manufacturer: "Fisher",
        model: "ES 2.1/2\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 2.5,
        ratedCv: 40.3,
        Fl: 0.95,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.43, 2.37, 3.34, 4.76, 7.25, 11.3, 17.3, 24.2, 31.8, 40.3]
    },
    {
        manufacturer: "Fisher",
        model: "ES 3\"",
        type: 0, 
        size: 3,
        ratedCv: 121,
        Fl: 0.89,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [4.05, 6.84, 10.0, 15.0, 23.8, 37.8, 59.0, 87.1, 110, 121]
    },
    {
        manufacturer: "Fisher",
        model: "ES 3\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 3,
        ratedCv: 67.5,
        Fl: 0.94,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [2.74, 3.44, 4.86, 6.95, 10.6, 16.5, 25.0, 37.7, 52.7, 67.5]
    },
    {
        manufacturer: "Fisher",
        model: "ES 4\"",
        type: 0, 
        size: 4,
        ratedCv: 203,
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [6.56, 11.4, 17.3, 27.0, 42.2, 66.4, 103, 146, 184, 203]
    },
    {
        manufacturer: "Fisher",
        model: "ES 4\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 4,
        ratedCv: 121,
        Fl: 0.94,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [3.96, 7.14, 10.6, 14.5, 21.1, 31.7, 48.0, 69.7, 95.6, 121]
    },
    {
        manufacturer: "Fisher",
        model: "ES 6\"",
        type: 0, 
        size: 6,
        ratedCv: 357,
        Fl: 0.86,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [13.2, 24.6, 41.1, 62.5, 97.1, 155, 223, 286, 326, 357]
    },
    {
        manufacturer: "Fisher",
        model: "ES 6\" (Restricted Trim)", // Restricted Trim
        type: 0, 
        size: 6,
        ratedCv: 233,
        Fl: 0.91,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [4.96, 9.02, 14.0, 24.2, 39.4, 60.8, 94.6, 144, 199, 233]
    },
    {
        manufacturer: "Fisher",
        model: "ES 8\"",
        type: 0, 
        size: 8,
        ratedCv: 570,
        Fl: 0.85,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [18.8, 33.6, 53.6, 79.8, 114, 168, 242, 345, 467, 570]
    },
    {
        manufacturer: "Fisher",
        model: "ES 8\"",
        type: 0, 
        size: 8,
        ratedCv: 808,
        Fl: 0.85,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [25.9, 53.3, 97.8, 178, 299, 461, 618, 727, 768, 808]
    },
    // =============================================================================
    // FISHER EZ - QUICK OPENING - FLOW UP
    // Source: Catalog 12, Page EZ-1
    // =============================================================================
    {
        manufacturer: "Fisher",
        model: "EZ 1/2\"",
        type: 0, // GLOBEOPEN
        size: 0.5,
        ratedCv: 4.44,
        Fl: 0.83,
        Fd: 0.46, // Not explicitly listed for 0.5", assuming default
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [1.76, 3.29, 4.29, 4.44, 4.44, 4.44, 4.44, 4.44, 4.44, 4.44]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3/4\"",
        type: 0, 
        size: 0.75,
        ratedCv: 9.72,
        Fl: 0.88,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [3.85, 7.19, 9.40, 9.72, 9.72, 9.72, 9.72, 9.72, 9.72, 9.72]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 1\"",
        type: 0, 
        size: 1,
        ratedCv: 16.9,
        Fl: 0.94,
        Fd: 0.50, // Listed as 0.50 at max travel
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [4.39, 10.3, 14.0, 15.5, 16.2, 16.6, 16.8, 16.8, 16.9, 16.9]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 1.1/2\"",
        type: 0, 
        size: 1.5,
        ratedCv: 34.2,
        Fl: 0.96,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [5.64, 11.9, 20.6, 27.4, 30.5, 32.4, 33.4, 33.7, 34.1, 34.2]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 2\" (Restricted Trim)",
        type: 0, 
        size: 2,
        ratedCv: 19.4,
        Fl: 0.90,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [4.17, 8.94, 14.6, 17.4, 18.3, 18.8, 18.9, 19.0, 19.1, 19.4]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 2\"",
        type: 0, 
        size: 2,
        ratedCv: 58.6,
        Fl: 0.94,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [13.0, 30.1, 44.3, 52.4, 56.4, 57.8, 58.4, 58.5, 58.6, 58.6]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\" (Restricted Trim)",
        type: 0, 
        size: 3,
        ratedCv: 17.9,
        Fl: 0.86,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [4.35, 9.79, 14.9, 16.6, 17.3, 17.5, 17.5, 17.6, 17.7, 17.9]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\"",
        type: 0, 
        size: 3,
        ratedCv: 129,
        Fl: 0.91,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [30.8, 53.8, 65.1, 92.4, 110, 118, 123, 126, 128, 129]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 4\" (Restricted Trim)",
        type: 0, 
        size: 4,
        ratedCv: 88.4,
        Fl: 0.95,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [9.99, 27.6, 44.9, 61.0, 71.9, 78.4, 83.1, 86.2, 87.5, 88.4]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 4\"",
        type: 0, 
        size: 4,
        ratedCv: 223,
        Fl: 0.88,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [50.8, 68.2, 116, 159, 185, 201, 212, 219, 222, 223]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 4\" (Restricted Trim)", // 2nd listing for 4" body with 2" trim
        type: 0, 
        size: 4,
        ratedCv: 86.7,
        Fl: 0.85,
        Fd: 0.50,
        pressureClasses: [150, 300, 600],
        characteristic: "Quick Opening",
        flowDirection: "Flow Up",
        cvCurve: [13.5, 32.3, 52.2, 66.2, 74.4, 81.1, 85.0, 85.8, 86.3, 86.7]
    },
{
        manufacturer: "Fisher",
        model: "EZ 1\"",
        type: 0, 
        size: 1,
        ratedCv: 13.6,
        Fl: 0.96,
        Fd: 0.46, 
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [2.21, 3.87, 5.29, 6.56, 8.20, 9.82, 11.1, 12.1, 13.0, 13.6]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 1.1/2\"",
        type: 0, 
        size: 1.5,
        ratedCv: 31.9,
        Fl: 0.96,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [3.99, 7.53, 11.1, 14.8, 18.7, 22.5, 25.8, 29.2, 31.2, 31.9]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 2\" (Restricted Trim)",
        type: 0, 
        size: 2,
        ratedCv: 16.7,
        Fl: 0.96,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [1.96, 3.42, 4.94, 6.11, 7.80, 9.30, 10.9, 13.0, 15.1, 16.7]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 2\"",
        type: 0, 
        size: 2,
        ratedCv: 52.4,
        Fl: 0.95,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [6.08, 11.9, 18.0, 24.1, 30.1, 36.4, 42.8, 49.9, 52.0, 52.4]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\" (Restricted Trim)",
        type: 0, 
        size: 3,
        ratedCv: 15.7,
        Fl: 0.94,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [1.88, 3.41, 4.95, 6.49, 8.06, 9.67, 11.23, 12.79, 14.35, 15.7]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\"",
        type: 0, 
        size: 3,
        ratedCv: 110.4,
        Fl: 0.92,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [15.4, 29.6, 43.4, 58.3, 71.8, 83.9, 93.8, 103, 108, 110.4]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\" (Restricted Trim)",
        type: 0, 
        size: 3,
        ratedCv: 80.4,
        Fl: 0.94,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [6.59, 13.3, 20.7, 28.1, 36.0, 44.0, 55.6, 67.5, 76.2, 80.4]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 4\" (Restricted Trim)",
        type: 0, 
        size: 4,
        ratedCv: 86.8,
        Fl: 0.90,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [6.16, 12.8, 20.0, 27.8, 36.1, 45.1, 58.8, 67.5, 78.8, 86.8]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 4\"",
        type: 0, 
        size: 4,
        ratedCv: 209,
        Fl: 0.89,
        Fd: 0.46,
        pressureClasses: [150, 300, 600],
        characteristic: "Linear",
        flowDirection: "Flow Up",
        cvCurve: [21.3, 39.7, 57.5, 75.8, 100, 129, 157, 180, 199, 209]
    },
{
        manufacturer: "Fisher",
        model: "EZ 1\"",
        type: 0, 
        size: 1,
        ratedCv: 13.2,
        Fl: 0.96,
        Fd: 0.50, // Listed: 0.09...0.50
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.79, 1.25, 1.80, 2.53, 3.63, 5.28, 7.59, 10.7, 12.7, 13.2]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 1.1/2\"",
        type: 0, 
        size: 1.5,
        ratedCv: 28.1,
        Fl: 0.97,
        Fd: 0.40, // Listed: 0.077...0.40
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.795, 1.23, 1.91, 2.95, 4.30, 6.46, 9.84, 16.4, 22.2, 28.1]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 2\" (Restricted Trim)",
        type: 0, 
        size: 2,
        ratedCv: 17.3,
        Fl: 0.98,
        Fd: 0.50, // Listed: 0.069...0.50
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [0.770, 1.23, 1.78, 2.58, 3.67, 5.54, 8.30, 12.0, 15.1, 17.3]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 2\"",
        type: 0, 
        size: 2,
        ratedCv: 53.8,
        Fl: 0.95,
        Fd: 0.46, // Listed: 0.062...0.46
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.65, 2.61, 4.30, 6.62, 11.1, 20.7, 32.8, 44.7, 50.0, 53.8]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\" (Restricted Trim)",
        type: 0, 
        size: 3,
        ratedCv: 13.8, // Note: Drops from 15.9 (90%) to 13.8 (100%)
        Fl: 0.92,
        Fd: 0.46, 
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.02, 1.50, 2.05, 2.78, 3.90, 5.57, 8.16, 11.8, 14.5, 13.8]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\"",
        type: 0, 
        size: 3,
        ratedCv: 114,
        Fl: 0.92,
        Fd: 0.46, 
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [3.11, 5.77, 9.12, 13.7, 21.7, 36.0, 60.4, 86.4, 104, 114]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 3\" (Restricted Trim)",
        type: 0, 
        size: 3,
        ratedCv: 71.6,
        Fl: 0.92,
        Fd: 0.44, // Listed: 0.052...0.44
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [2.11, 3.11, 4.58, 6.76, 10.7, 20.7, 34.3, 48.3, 61.5, 71.6]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 4\" (Restricted Trim)",
        type: 0, 
        size: 4,
        ratedCv: 72.7,
        Fl: 0.92,
        Fd: 0.44, 
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [1.96, 3.05, 4.43, 6.98, 11.9, 22.3, 36.7, 50.9, 61.8, 72.7]
    },
    {
        manufacturer: "Fisher",
        model: "EZ 4\"",
        type: 0, 
        size: 4,
        ratedCv: 190,
        Fl: 0.90,
        Fd: 0.46, // Standard Default (Not explicitly listed in this snippet)
        pressureClasses: [150, 300, 600],
        characteristic: "Equal Percentage",
        flowDirection: "Flow Up",
        cvCurve: [4.90, 8.19, 13.5, 20.1, 31.2, 52.6, 96.7, 140, 170, 190]
    }
];
