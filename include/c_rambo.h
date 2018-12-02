void * declare_phs_rambo (const int n_particles, const double total_momentum[4], const double masses[]);

void free_phs_rambo (void * phs);

void write_phs_rambo (void * phs);

void generate_phs_rambo (void * phs, const int n_r_in, const double r_in[]);

double get_weight_phs_rambo (void * phs);
