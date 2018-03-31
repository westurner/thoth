gsl_rng *rng();

void norm_prob(prob *a);
void norm_jprob(j_prob *a);

void print_prob(prob *a);
void print_jprob(j_prob *a);

double entropy(prob *a);
double entropy_ww(prob *a);

double jsd(prob *a, prob *b, double alpha);

double mi(j_prob *a);
double mi_ww(j_prob *pij, double beta);

void entropy_bs_wrapper(prob *a, long n_samp, gsl_rng *r, double* out_results, int quantities);
void jsd_bs_wrapper(prob *a, prob *b, double alpha, long n_samp, gsl_rng *r, double* out_results, int quantities);
void mi_bs_wrapper(j_prob *a, long n_samp, gsl_rng *r, double* out_results, int quantities);

void do_entropy(int n_states, int n_fact, char *filename);

void do_mi(int n_states, int m_states, int n_fact, char *filename);

void do_jsd(int n_states, double alpha, int n_fact, char *filename);

double entropy_nsb(prob *p);
double mi_nsb(j_prob *p);

prob* import_prob(char* file_name);
prob* prob_from_array(double* in_array, int n);
j_prob* jprob_from_array(double* j_in_array, int n1, int n2);

