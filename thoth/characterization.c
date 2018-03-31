#include "prob.h"

void do_entropy(int n_states, int n_fact, char *filename) {
	prob *a, *a_samp;
	prob *a_cg_1, *a_cg_2;
	FILE *f_out;
	int n_samples, i, j;
	double *bs, *bs_a, *bs_b, beta;
	int one_sig_ok, two_sig_ok;
	gsl_rng *r_global;
	
	f_out=fopen(filename, "w");
	
	//
	r_global=rng();
	//
	
	// allocate the true and sample space
	a=new_prob(n_states);
	a_samp=new_prob(n_states);
	
	// if you have three states, allocate the coarse-grained case
	if (n_states == 3) {
		a_cg_1=new_prob(2);
		a_cg_2=new_prob(2);
		n_samples=100000;
	} else {
		n_samples=100000;
	}
	
	for(i=0;i<10000;i++) {
		beta=-1;
		// draw a new true probability (uniform)
		sample_prob(a, beta, r_global);
		// sample from it
		draw(a, a_samp, n_states*n_fact, r_global);
		while(!(entropy(a_samp) > 0)) {
			sample_prob(a, beta, r_global);
			draw(a, a_samp, n_states*n_fact, r_global);			
		};
		
		// run the sampler
		bs=entropy_bs(a_samp, n_samples, r_global);

		// if we are checking coarse-graining as well...
		if (n_states == 3) {
			
			for(j=0;j<n_states-1;j++) {
				a_cg_1->p[j]=a_samp->p[j];
			}
			a_cg_1->p[n_states-2] += a_samp->p[n_states-1];
			a_cg_2->p[0]=a_samp->p[n_states-2];
			a_cg_2->p[1]=a_samp->p[n_states-1];
			
			// we need to renormalize...
			a_cg_1->norm=-1;
			a_cg_2->norm=-1;
			
			bs_a=entropy_bs(a_cg_1, n_samples, r_global);
			
			// if there are non-zero hits...
			if ((a_cg_2->p[0]+a_cg_2->p[1]) > 0) {
				bs_b=entropy_bs(a_cg_2, n_samples, r_global);
			} else {
				bs_b=(double *)malloc(sizeof(double));
				bs_b[0]=0.0;
			}	
		}
		
		one_sig_ok=0;
		two_sig_ok=0;
		if ((entropy(a) >= (bs[1]-1e-4)) && (entropy(a) <= (bs[2]+1e-4))) {
			one_sig_ok=1;
		}
		if ((entropy(a) >= (bs[3]-1e-4)) && (entropy(a) <= (bs[4]+1e-4))) {
			two_sig_ok=1;
		}						

		
		fprintf(f_out, "%12.15lf %12.15lf %12.15lf %12.15lf %12.15lf %12.15lf %12.15lf %i %i", entropy(a), entropy_ww(a_samp), entropy_nsb(a_samp), entropy(a_samp), bs[0], bs[1], bs[2], one_sig_ok, two_sig_ok);
		
		if (n_states == 3) {
			fprintf(f_out, " %12.15lf ", bs[0]-(bs_a[0]+bs_b[0]*(double)(a_samp->p[n_states-2]+a_samp->p[n_states-1])/(1e-64+(double)a_samp->norm)));
			fprintf(f_out, " %12.15lf ", entropy_ww(a_samp) - (entropy_ww(a_cg_1) + entropy_ww(a_cg_2)*(double)(a_samp->p[n_states-2]+a_samp->p[n_states-1])/(1e-64+(double)a_samp->norm)));
			fprintf(f_out, " %12.15lf\n", entropy_nsb(a_samp) - (entropy_nsb(a_cg_1) + entropy_nsb(a_cg_2)*(double)(a_samp->p[n_states-2]+a_samp->p[n_states-1])/(1e-64+(double)a_samp->norm)));
			
			free(bs_a);
			free(bs_b);
		} else {
			fprintf(f_out, "\n");
		}
		
		free(bs);
	}
	
	delete_prob(a);
	delete_prob(a_samp);
	if (n_states == 3) {
		delete_prob(a_cg_1);
		delete_prob(a_cg_2);
	}
	
	gsl_rng_free(r_global);
	fclose(f_out);
}

void do_mi(int n_states, int m_states, int n_fact, char *filename) {
	j_prob *pij, *pij_samp, *pij_cg_1, *pij_cg_2;
	int one_sig_ok, two_sig_ok, n_samples, i;
	gsl_rng *r_global;
	double *bs, *bs_a, *bs_b, beta;
	FILE *f_out;

	f_out=fopen(filename, "w");
	
	//
	r_global=rng();
	//

	// allocate the samples space
	pij=new_jprob(n_states, m_states);
	pij_samp=new_jprob(pij->n, pij->m);
	
	// if you have two/three state case, allocate the coarse-grained
	if ((n_states == 3) && (m_states == 2)) {
		pij_cg_1=new_jprob(2, 2);
		pij_cg_2=new_jprob(2, 2);
		n_samples=100000;
	} else {
		n_samples=100000;
	}
		
	for(i=0;i<10000;i++) {
		beta=-1;
		// draw a new true probability (uniform)
		sample_jprob(pij, beta, r_global);
		// sample from it
		draw_joint(pij, pij_samp, n_states*m_states*n_fact, r_global);
		while(!(mi(pij_samp) > 0)) {
			sample_jprob(pij, beta, r_global);
			draw_joint(pij, pij_samp, n_states*m_states*n_fact, r_global);
		};
		
		// run the sampler
		bs=mi_bs(pij_samp, n_samples, r_global);

		// check the error bars
		one_sig_ok=0;
		two_sig_ok=0;
		if ((mi(pij) > bs[1]) && (mi(pij) < bs[2])) {
			one_sig_ok=1;
		}
		if ((mi(pij) > bs[3]) && (mi(pij) < bs[4])) {
			two_sig_ok=1;
		}						
		
		fprintf(f_out, "%12.15lf %12.15lf %12.15lf %12.15lf %12.15lf %12.15lf %12.15lf %i %i", mi(pij), mi_ww(pij_samp, 1.0), mi_nsb(pij_samp), mi(pij_samp), bs[0], bs[1], bs[2], one_sig_ok, two_sig_ok);

		if ((n_states == 3) && (m_states == 2)) {
			
			pij_cg_1->p[0][0]=pij_samp->p[0][0];
			pij_cg_1->p[0][1]=pij_samp->p[0][1];

			pij_cg_1->p[1][0]=pij_samp->p[1][0]+pij_samp->p[2][0];
			pij_cg_1->p[1][1]=pij_samp->p[1][1]+pij_samp->p[2][1];
			
			pij_cg_2->p[0][0]=pij_samp->p[1][0];
			pij_cg_2->p[1][0]=pij_samp->p[2][0];
			
			pij_cg_2->p[0][1]=pij_samp->p[1][1];
			pij_cg_2->p[1][1]=pij_samp->p[2][1];
						
			// we need to renormalize...
			pij_cg_1->norm=-1;
			pij_cg_2->norm=-1;
			norm_jprob(pij_cg_1);
			norm_jprob(pij_cg_2);
			
			bs_a=mi_bs(pij_cg_1, n_samples, r_global);
			
			// if there are non-zero hits...
			if (pij_cg_2->norm > 0) {
				bs_b=mi_bs(pij_cg_2, n_samples, r_global);
			} else {
				bs_b=(double *)malloc(sizeof(double));
				bs_b[0]=0.0;
			}	
			
			fprintf(f_out, " %12.15lf ", bs[0]-(bs_a[0]+bs_b[0]*pij_cg_2->norm/(1e-64+pij_samp->norm)));
			fprintf(f_out, "%12.15lf ", mi_ww(pij_samp, 1.0/(n_states*m_states)) - (mi_ww(pij_cg_1, 1.0/(n_states*m_states)) + mi_ww(pij_cg_2, 1.0/(n_states*m_states))*pij_cg_2->norm/(1e-64+pij_samp->norm)));
			fprintf(f_out, "%12.15lf\n", mi_nsb(pij_samp) - (mi_nsb(pij_cg_1) + mi_nsb(pij_cg_2)*pij_cg_2->norm/(1e-64+pij_samp->norm)));
			
			free(bs_a);
			free(bs_b);
			
		} else {
			fprintf(f_out, "\n");
		}
	
		free(bs);
	}
	
	
	delete_jprob(pij);
	delete_jprob(pij_samp);
	if ((n_states == 3) && (m_states == 2)) {
		delete_jprob(pij_cg_1);
		delete_jprob(pij_cg_2);
	}
	
	gsl_rng_free(r_global);	
	fclose(f_out);
	
}

void do_jsd(int n_states, double alpha, int n_fact, char *filename) {
	gsl_rng *r_global;
	prob *a, *a_samp, *b, *c, *c_samp;
	double *bs;
	
	a=new_prob(n_states);
	sample_prob(a, -1, r_global);

	b=new_prob(n_states);
	sample_prob(b, -1, r_global);

	c=mixture(a, b, alpha);
	delete_prob(b);
	
	a_samp=new_prob(a->n);
	draw(a, a_samp, n_states*n_fact, r_global);
	c_samp=new_prob(a->n);
	draw(c, c_samp, n_states*n_fact, r_global);

	bs=jsd_bs(a_samp, c_samp, alpha, 10000, r_global);
	
}
