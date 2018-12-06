#define phs_rambo void *

void * declare_phs_rambo (const int n_particles, const double total_momentum[4], const double masses[]);

void free_phs_rambo (phs_rambo phs);

void write_phs_rambo (phs_rambo phs);

void generate_phs_rambo (phs_rambo phs, const int n_r_in, const double r_in[]);

void invert_phs_rambo (phs_rambo phs, double r_out[]);

double get_weight_phs_rambo (phs_rambo phs);

void get_p_phs_rambo (phs_rambo phs, const int i, double p[4]);
