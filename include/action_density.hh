#ifndef ACTION_DENSITY_HH
#define ACTION_DENSITY_HH

// Wilson gauge action density in 4D from the normalized average plaquette:
// s = S / V = 6 * beta * (1 - <P>).
inline double avg_action_density(double beta, double avg_plaquette) {
    return 6.0 * beta * (1.0 - avg_plaquette);
}

#endif // ACTION_DENSITY_HH
