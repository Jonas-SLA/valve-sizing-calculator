// This variable is now globally accessible
const valveDatabase = [
    {
        manufacturer: "Fisher",
        model: "ET Globe",
        type: 0, // 0 = GLOBEOPEN
        size: 4,
        ratedCv: 150,
        Fl: 0.90,
        Fd: 0.46
    },
    {
        manufacturer: "Masoneilan",
        model: "21000 Series",
        type: 0, // 0 = GLOBEOPEN
        size: 4,
        ratedCv: 135,
        Fl: 0.89,
        Fd: 0.46
    },
    { // A smaller valve (reduced port)
        manufacturer: "Masoneilan",
        model: "21000 Series",
        type: 0, // 0 = GLOBEOPEN
        size: 3,
        ratedCv: 90,
        Fl: 0.92,
        Fd: 0.47
    },
    {
        manufacturer: "Generic",
        model: "V-Port Ball",
        type: 3, // 3 = BALL
        size: 4,
        ratedCv: 250,
        Fl: 0.70,
        Fd: 0.95
    }
    // ... add all your other valves here
];