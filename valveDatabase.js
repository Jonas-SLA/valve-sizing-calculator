// This variable is now globally accessible
const valveDatabase = [
    {
        manufacturer: "Fisher",
        model: "ET Globe",
        type: 0, // 0 = GLOBEOPEN
        size: 4,
        ratedCv: 150,
        Fl: 0.90,
        Fd: 0.46,
        characteristic: "Equal Percentage",
        pressureClasses: [150, 300, 600, 900],
        cvCurve: [8, 14, 23, 35, 52, 72, 95, 120, 138, 150]
    },
    {
        manufacturer: "Masoneilan",
        model: "21000 Series",
        type: 0, // 0 = GLOBEOPEN
        size: 4,
        ratedCv: 135,
        Fl: 0.89,
        Fd: 0.46,
        pressureClasses: [150, 300, 600, 900],
        cvCurve: [8, 14, 23, 35, 52, 72, 95, 110, 120, 135]
    },
    { // A smaller valve (reduced port)
        manufacturer: "Masoneilan",
        model: "21000 Series",
        type: 0, // 0 = GLOBEOPEN
        size: 3,
        ratedCv: 90,
        Fl: 0.92,
        Fd: 0.47,
        pressureClasses: [300, 600, 900],
        cvCurve: [8, 14, 23, 35, 52, 62, 75, 80, 85, 90]
    },
    {
        manufacturer: "Generic",
        model: "V-Port Ball",
        type: 3, // 3 = BALL
        size: 4,
        ratedCv: 250,
        Fl: 0.70,
        Fd: 0.95,
        pressureClasses: [900, 1500, 2500, 4500],
        cvCurve: [8, 14, 23, 35, 52, 72, 95, 140, 188, 250]
    }
    // ... add all your other valves here
];
