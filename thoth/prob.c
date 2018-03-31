#include "prob.h"
// extern void rancom_(int *in, int *ik, int *ip, int *iseed);

double cumulative_distribution(double beta, void *p) {
	struct nsb_params *params = (struct nsb_params *)p;
    
	return (-gsl_sf_psi_n(0, 1.0+beta) + gsl_sf_psi_n(0, 1.0+beta*(double)params->n))/log(params->n) - (double)params->choice;
}

void sample_prob(prob *a, double k, gsl_rng *r) {
	double *alpha, *alpha_breaker, *breaker, *p_long, choice;
	int i, j, start, finish, iter=0, status;
	int true_bins, dist_bins;
	int *ip, iseed, long_i;
	gsl_root_fsolver *s;
	gsl_function F;

	double root = 0, x_lo = 1e-5, x_hi = 1e3;
	
	if (k > 0) {
		
		alpha=(double *)malloc(a->n*sizeof(double));
		for(i=0;i<true_bins;i++) {
			alpha[i]=k;
		}
		gsl_ran_dirichlet(r, true_bins, alpha, a->p);
		free(alpha);

		if (a->pre_proc != NULL) {
			gsl_ran_discrete_free(a->pre_proc);
			a->pre_proc=NULL;
		}
				
	} else { // uniform requested
		choice=gsl_rng_uniform(r);
		struct nsb_params params={(double)a->n, choice};
		F.function = &cumulative_distribution;
		F.params = &params;
		
		if (cumulative_distribution(x_hi, &params) < 0) {
			sample_prob(a, x_hi, r);
			return;
		}
		if (cumulative_distribution(x_lo, &params) > 0) {
			sample_prob(a, x_lo, r);
			return;
		}

		s = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
		gsl_root_fsolver_set(s, &F, x_lo, x_hi);
		
		iter=0;
		do {
           iter++;
           	status = gsl_root_fsolver_iterate(s);
           	root = gsl_root_fsolver_root(s);
           	x_lo = gsl_root_fsolver_x_lower(s);
           	x_hi = gsl_root_fsolver_x_upper(s);
			// printf("%lf (%lf %lf) -- value: %lf\n", root, x_lo, x_hi, cumulative_distribution(root, &params));
           	status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-5);
		} while (status == GSL_CONTINUE && iter < 1000);

		gsl_root_fsolver_free(s);
		sample_prob(a, root, r);
	}
	
	a->norm=1.0;
}

// void sample_prob(prob *a, double k, gsl_rng *r) {
// 	double *alpha, *alpha_breaker, *breaker, *p_long, choice;
// 	int i, j, start, finish, iter=0, status;
// 	int true_bins, dist_bins;
// 	int *ip, iseed, long_i;
// 	gsl_root_fsolver *s;
// 	gsl_function F;
// 
// 	double root = 0, x_lo = 1e-5, x_hi = 1e3;
// 	
// 	if (k > 0) {
// 		// choose a random number of bins between a->n and a->n/2
// 		// true_bins=(int)gsl_rng_uniform_int(r, (unsigned long int)(a->n/2+1))+a->n; //
// 		true_bins=a->n;
// 		
// 		alpha=(double *)malloc(true_bins*sizeof(double));
// 		// p_long=(double *)malloc(true_bins*sizeof(double));
// 		for(i=0;i<true_bins;i++) {
// 			alpha[i]=k;
// 		}
// 		// gsl_ran_dirichlet(r, true_bins, alpha, p_long);
// 		gsl_ran_dirichlet(r, true_bins, alpha, a->p);
// 		free(alpha);
// 		
// 		// printf("Beta is %lf\n", k);
// 		
// 		// printf("ORIGINAL (%i)\n", true_bins);
// 		// for(i=0;i<true_bins;i++) {
// 		// 	printf("%lf ", p_long[i]);
// 		// }
// 		// printf("\n");
// 		
// 		// if (true_bins > a->n+1) {
// 		// 	dist_bins=true_bins-a->n;
// 		// 	
// 		// 	iseed=gsl_rng_get(r);
// 		// 	ip=(int *)malloc(a->n*sizeof(int));
// 		// 	rancom_(&dist_bins, &a->n, &ip[0], &iseed); // sort left over bins do the a->n you have
// 		// 
// 		// 	i=0;
// 		// 	long_i=0;
// 		// 
// 		// 	while(long_i < true_bins) {
// 		// 		// printf("INDEX %i\n", i);
// 		// 		a->p[i]=0;
// 		// 		// printf("Add %i to %i (%i units)\n", long_i, long_i+ip[i], ip[i]+1);
// 		// 		for(j=0;j<ip[i]+1;j++) {
// 		// 			a->p[i] += p_long[long_i+j];
// 		// 		}
// 		// 		long_i += ip[i]+1;
// 		// 		i += 1;
// 		// 	}
// 		// 	free(ip);
// 		// 	
// 		// } else {
// 		// 	for(i=0;i<a->n;i++) {
// 		// 		a->p[i]=p_long[i];
// 		// 	}
// 		// 	if (true_bins == a->n+1) {
// 		// 		a->p[a->n-1]=p_long[a->n];
// 		// 	}
// 		// }
// 		// free(p_long);
// 
// 		// printf("FINAL: ");
// 		// for(i=0;i<a->n;i++) {
// 		// 	printf("%lf ", a->p[i]);
// 		// }
// 		// printf("\n");
// 		
// 		// you just changed the sampling distribution -- so go and delete the pre-process
// 		if (a->pre_proc != NULL) {
// 			gsl_ran_discrete_free(a->pre_proc);
// 			a->pre_proc=NULL;
// 		}
// 				
// 	} else { // uniform requested
// 		choice=gsl_rng_uniform(r);
// 		struct nsb_params params={(double)a->n, choice};
// 		F.function = &cumulative_distribution;
// 		F.params = &params;
// 		
// 		if (cumulative_distribution(x_hi, &params) < 0) {
// 			sample_prob(a, x_hi, r);
// 			return;
// 		}
// 		if (cumulative_distribution(x_lo, &params) > 0) {
// 			sample_prob(a, x_lo, r);
// 			return;
// 		}
// 
// 		s = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
// 		gsl_root_fsolver_set(s, &F, x_lo, x_hi);
// 		
// 		iter=0;
// 		do {
//            iter++;
//            	status = gsl_root_fsolver_iterate(s);
//            	root = gsl_root_fsolver_root(s);
//            	x_lo = gsl_root_fsolver_x_lower(s);
//            	x_hi = gsl_root_fsolver_x_upper(s);
// 			// printf("%lf (%lf %lf) -- value: %lf\n", root, x_lo, x_hi, cumulative_distribution(root, &params));
//            	status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-5);
// 		} while (status == GSL_CONTINUE && iter < 1000);
// 
// 		gsl_root_fsolver_free(s);
// 		sample_prob(a, root, r);
// 	}
// 	
// 	a->norm=1.0;
// }
// 

void sample_jprob(j_prob *a, double k, gsl_rng *r) { // this one is easy -- call the above, and then reform
	prob *a_temp;
	int i;
	
	a_temp=new_prob(a->n*a->m);
	sample_prob(a_temp, k, r);

	for(i=0;i<a->n*a->m;i++) {
		a->p[i/a->m][i % a->m]=a_temp->p[i];
	}
	if (a->pre_proc != NULL) {
		gsl_ran_discrete_free(a->pre_proc);
		a->pre_proc=NULL;
	}
	
	delete_prob(a_temp);
	a->norm=1.0;
}

prob *mixture(prob *a, prob *b, double alpha) {
	prob *m;
	int i;
	
	m=new_prob(a->n);
	(m->norm)=1;
	
	norm_prob(a);
	norm_prob(b);
	
	for(i=0;i<a->n;i++) {
		m->p[i]=alpha*a->p[i]/a->norm + (1.0-alpha)*b->p[i]/b->norm;
	}
	m->norm=1.0;
	
	return m;
}

double fast_mixture(prob *a, prob *b, double alpha) { // takes the entropy of the mixture of a and b, a subcomputation for JSD
	int i;
	double ans, temp_prob;
		
	norm_prob(a);
	norm_prob(b);
	
	ans=0.0;
	
	for(i=0;i<a->n;i++) {
		temp_prob=alpha*a->p[i]/a->norm + (1.0-alpha)*b->p[i]/b->norm;
		if (temp_prob > 0) {
			ans += temp_prob*log(temp_prob);
		}
	}
	
	return -ans/log(2);
}

double entropy(prob *a) {
	int i;
	double ent;
	
	norm_prob(a);
	
	ent=0;
	for(i=0;i<(a->n);i++) {
		if ((a->p)[i] > 0) {
			ent += (a->p)[i]*log(a->p[i]);
		}
	}
	
	ent = (ent/a->norm - log(a->norm))/log(2.0);
	
	return -ent;
}

double jsd(prob *a, prob *b, double alpha) {
	
	return fast_mixture(a, b, alpha) - alpha*entropy(a) - (1.0-alpha)*entropy(b);
}

double mi(j_prob *pij) {
	int i, j;
	prob *pi, *pj;
	double total, total_mi;
	
	norm_jprob(pij);
	
	pi=new_prob(pij->n);
	pj=new_prob(pij->m);
	
	for(i=0;i<pij->n;i++) {
		total=0;
		for(j=0;j<pij->m;j++) {
			total += pij->p[i][j];
		}
		pi->p[i]=total;
	}
	
	for(j=0;j<pij->m;j++) {
		total=0;
		for(i=0;i<pij->n;i++) {
			total += pij->p[i][j];
		}
		pj->p[j]=total;
	}
	
	total_mi=0;
	for(i=0;i<pij->n;i++) {
		for(j=0;j<pij->m;j++) {
			if ((pij->p[i][j] > 0) && (pi->p[i]*pj->p[j] > 0)) {
				total_mi += pij->p[i][j]*log(pij->p[i][j]*pij->norm/(pi->p[i]*pj->p[j]));
			}
		}
	}
	total_mi = total_mi/pij->norm;
	
	delete_prob(pi);
	delete_prob(pj);
	
	return total_mi/log(2);
}

double *entropy_bs(prob *a, long n_samp, gsl_rng *r) {
	long i;
	double *ent_list, *ent_out;
	double ent_mean, ent_naieve;
	double lo_1sig, hi_1sig, lo_2sig, hi_2sig;
	prob *samp;
	
	norm_prob(a);
	// for safety, clear the pre-processed draw variable
	if (a->pre_proc != NULL) {
		gsl_ran_discrete_free(a->pre_proc);
		a->pre_proc=NULL;
	}
	
	samp=new_prob(a->n);
	ent_list=(double *)malloc(sizeof(double)*n_samp);
	
	ent_mean=0.0;
	for(i=0;i<n_samp;i++) {
		draw(a, samp, a->norm, r);
		ent_list[i]=entropy(samp);
		ent_mean += ent_list[i];
	}
	ent_mean = ent_mean/(double)n_samp;
	
	delete_prob(samp);
	
	qsort(ent_list, n_samp, sizeof(double), d_compare);
	
	lo_1sig=(double)n_samp*(0.136+0.021+0.001);
	lo_2sig=(double)n_samp*(0.021+0.001);
	hi_1sig=(double)n_samp*(1.0-(0.136+0.021+0.001));
	hi_2sig=(double)n_samp*(1.0-(0.021+0.001));
    
	ent_naieve=entropy(a);
	
	ent_out=(double *)malloc(5*sizeof(double));
	ent_out[0]=2.0*ent_naieve-ent_mean;
	ent_out[1]=2.0*ent_naieve-ent_list[(long)hi_1sig];
	ent_out[2]=2.0*ent_naieve-ent_list[(long)lo_1sig];
	ent_out[3]=2.0*ent_naieve-ent_list[(long)hi_2sig];
	ent_out[4]=2.0*ent_naieve-ent_list[(long)lo_2sig];
	
	free(ent_list);
	
	return ent_out;
}

void entropy_bs_wrapper(prob *a, long n_samp, gsl_rng *r, double* out_results, int quantities) {
	double* ent_out;
	int i;	
	ent_out = (double*)malloc(5*sizeof(double));
	ent_out = entropy_bs(a, n_samp, r);
	for (i = 0; i < 5; i++) {
		out_results[i] = ent_out[i];
	}
}


double *jsd_bs(prob *a, prob *b, double alpha, long n_samp, gsl_rng *r) {
	long i;
	double *jsd_list, *jsd_out;
	double jsd_mean, jsd_naieve;
	double lo_1sig, hi_1sig, lo_2sig, hi_2sig;
	prob *a_samp, *b_samp;
	
	norm_prob(a);
	norm_prob(b);
	
	a_samp=new_prob(a->n);
	b_samp=new_prob(a->n);
	
	jsd_list=(double *)malloc(sizeof(double)*n_samp);
	jsd_mean=0;
	for(i=0;i<n_samp;i++) {
		draw(a, a_samp, a->norm, r);
		draw(b, b_samp, b->norm, r);
		
		jsd_list[i]=jsd(a_samp, b_samp, alpha);
		jsd_mean += jsd_list[i];
	}
	jsd_mean = jsd_mean/(double)n_samp;
	
	qsort(jsd_list, n_samp, sizeof(double), d_compare);
	
	lo_1sig=(double)n_samp*(0.136+0.021+0.001);
	lo_2sig=(double)n_samp*(0.021+0.001);
	hi_1sig=(double)n_samp*(1.0-(0.136+0.021+0.001));
	hi_2sig=(double)n_samp*(1.0-(0.021+0.001));
    
	jsd_naieve=jsd(a, b, alpha);
	
	jsd_out=(double *)malloc(5*sizeof(double));
	jsd_out[0]=2.0*jsd_naieve-jsd_mean;
	jsd_out[1]=2.0*jsd_naieve-jsd_list[(long)hi_1sig];
	jsd_out[2]=2.0*jsd_naieve-jsd_list[(long)lo_1sig];
	jsd_out[3]=2.0*jsd_naieve-jsd_list[(long)hi_2sig];
	jsd_out[4]=2.0*jsd_naieve-jsd_list[(long)lo_2sig];
	
	free(jsd_list);
	
	return jsd_out;
}

void jsd_bs_wrapper(prob *a, prob *b, double alpha, long n_samp, gsl_rng *r, double* out_results, int quantities) {
   double* jsd_out;
	int i;	
	jsd_out = (double*)malloc(5*sizeof(double));
	jsd_out = jsd_bs(a, b, alpha, n_samp, r);
	for (i = 0; i < 5; i++) {
		out_results[i] = jsd_out[i];
	}
}

double *mi_bs(j_prob *a, long n_samp, gsl_rng *r) {
	long i;
	double *mi_list, *mi_out;
	double mi_mean, mi_naieve;
	double lo_1sig, hi_1sig, lo_2sig, hi_2sig;
	j_prob *samp;

	norm_jprob(a);
	// for safety, clear the pre-processed draw variable
	if (a->pre_proc != NULL) {
		gsl_ran_discrete_free(a->pre_proc);
		a->pre_proc=NULL;
	}
	
	samp=new_jprob(a->n, a->m);
	mi_list=(double *)malloc(sizeof(double)*n_samp);
	mi_mean=0;
	for(i=0;i<n_samp;i++) {
		draw_joint(a, samp, a->norm, r);
		mi_list[i]=mi(samp);
		mi_mean += mi_list[i];
	}
	delete_jprob(samp);
	
	mi_mean = mi_mean/(double)n_samp;

	qsort(mi_list, n_samp, sizeof(double), d_compare);

	lo_1sig=(double)n_samp*(0.136+0.021+0.001);
	lo_2sig=(double)n_samp*(0.021+0.001);
	hi_1sig=(double)n_samp*(1.0-(0.136+0.021+0.001));
	hi_2sig=(double)n_samp*(1.0-(0.021+0.001));

	mi_naieve=mi(a);

	mi_out=(double *)malloc(5*sizeof(double));
	mi_out[0]=2.0*mi_naieve-mi_mean;
	mi_out[1]=2.0*mi_naieve-mi_list[(long)hi_1sig];
	mi_out[2]=2.0*mi_naieve-mi_list[(long)lo_1sig];
	mi_out[3]=2.0*mi_naieve-mi_list[(long)hi_2sig];
	mi_out[4]=2.0*mi_naieve-mi_list[(long)lo_2sig];

	free(mi_list);
	
	return mi_out;
}

void mi_bs_wrapper(j_prob *a, long n_samp, gsl_rng *r, double* out_results, int quantities) {
   double* mi_out;
	int i;	
	mi_out = (double*)malloc(5*sizeof(double));
	mi_out = mi_bs(a, n_samp, r);
	for (i = 0; i < 5; i++) {
		out_results[i] = mi_out[i];
	}
}

double entropy_ww(prob *a) {
	double ent;
	int i;
	
	norm_prob(a);
	
	ent=0;
	for(i=0;i<a->n;i++) {
		ent += ((a->p[i] + 1.0)/(a->norm + a->n))*(gsl_sf_psi_n(0, a->p[i]+2) - gsl_sf_psi_n(0, a->norm+a->n+1));
	}
	
	return -ent/log(2.0);
}

double mi_ww(j_prob *pij, double beta) {
	double total_mi, total, nu;
	prob *pi, *pj;
	int i, j;

	norm_jprob(pij);
	
	pi=new_prob(pij->n);
	pj=new_prob(pij->m);
	
	// NU CONVENTION
	for(i=0;i<pij->n;i++) {
		total=0;
		for(j=0;j<pij->m;j++) {
			total += pij->p[i][j]+beta;
		}
		pi->p[i]=total;
	}

	// NU CONVENTION	
	for(j=0;j<pij->m;j++) {
		total=0;
		for(i=0;i<pij->n;i++) {
			total += pij->p[i][j]+beta;
		}
		pj->p[j]=total;
	}

	norm_prob(pi);
	nu=pi->norm;
	
	// have to add on the betas here...
	total_mi=0;
	for(i=0;i<pij->n;i++) {
		for(j=0;j<pij->m;j++) {
			total_mi += (pij->p[i][j]+beta)*(gsl_sf_psi_n(0, pij->p[i][j]+beta+1) - gsl_sf_psi_n(0, nu+1));
		}
	}
	// but these are already accounted for...
	for(i=0;i<pij->n;i++) {
		total_mi -= pi->p[i]*(gsl_sf_psi_n(0, pi->p[i]+1) - gsl_sf_psi_n(0, nu+1));
	}
	for(j=0;j<pij->m;j++) {
		total_mi -= pj->p[j]*(gsl_sf_psi_n(0, pj->p[j]+1) - gsl_sf_psi_n(0, nu+1));
	}
	
	total_mi = total_mi/nu;
	
	delete_prob(pi);
	delete_prob(pj);

	return total_mi/log(2.0);
}

