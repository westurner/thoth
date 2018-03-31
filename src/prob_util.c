#include "prob.h"

int d_compare(const void *elem1, const void *elem2) {
	double f = *((double *)elem1);
	double s = *((double *)elem2);

	return (f > s) - (f < s);
}

void norm_prob(prob *a) { // for speed, we just save the overall normalizations
	int i;
	
	if (a->norm < 0) {
		a->norm=0;
		for(i=0;i<a->n;i++) {
			a->norm += a->p[i];
		}
	}
}

void norm_jprob(j_prob *a) { // for speed, we just save the overall normalizations
	int i,j;
	
	if (a->norm < 0) {
		a->norm=0;
		for(i=0;i<a->n;i++) {
			for(j=0;j<a->m;j++) {
				a->norm += a->p[i][j];
			}
		}
	}
}
 
prob *new_prob(int n) {
	prob *m;
	
	m=(prob *)malloc(sizeof(prob));
	(m->n)=n;
	(m->norm)=-1;
	(m->p)=(double *)malloc(sizeof(double)*n);
	(m->pre_proc)=NULL;
	
	return m;
}
	
void delete_prob(prob *a) {
	free(a->p);
	if (a->pre_proc != NULL) {
		gsl_ran_discrete_free(a->pre_proc);
	}
	free(a);
}

j_prob *new_jprob(int n, int m) {
	j_prob *q;
	int i;
	
	q=(j_prob *)malloc(sizeof(prob));
	(q->n)=n;
	(q->m)=m;
	(q->norm)=-1;
	
	(q->p)=(double **)malloc(sizeof(double *)*n);
	for(i=0;i<n;i++) {
		(q->p)[i]=(double *)malloc(sizeof(double)*m);
	}
	
	(q->pre_proc)=NULL;
	
	return q;
}

void delete_jprob(j_prob *a) {
	int i;
	
	for(i=0;i<a->n;i++) {
		free(a->p[i]);
	}
	free(a->p);
	
	if (a->pre_proc != NULL) {
		gsl_ran_discrete_free(a->pre_proc);
	}
	free(a);
}

void draw(prob *a, prob *samp, int k, gsl_rng *r) { // draw 
	int i;
	
	if (a->pre_proc == NULL) {
		a->pre_proc=gsl_ran_discrete_preproc(a->n, a->p);
	}
	
	for(i=0;i<a->n;i++) {
		samp->p[i]=0;
	}
	
	for(i=0;i<k;i++) {
		samp->p[gsl_ran_discrete(r, a->pre_proc)] += 1;
	}
	samp->norm=(double)k;
		
}

void draw_joint(j_prob *a, j_prob *samp, int k, gsl_rng *r) {
	int i, j, choice;
	double *flatten;
		
	if (a->pre_proc == NULL) {
		flatten=(double *)malloc(sizeof(double)*a->m*a->n);
		for(i=0;i<a->m*a->n;i++) {
			flatten[i]=a->p[i/a->m][i % a->m];
		}
		a->pre_proc=gsl_ran_discrete_preproc(a->n*a->m, flatten);
		free(flatten);
	}
	
	for(i=0;i<a->n;i++) {
		for(j=0;j<a->m;j++) {
			samp->p[i][j]=0;
		}
	}

	
	for(i=0;i<k;i++) {
		choice=gsl_ran_discrete(r, a->pre_proc);
		samp->p[choice/a->m][choice % a->m]++;
	}
	samp->norm=(double)k;
	
}

void print_prob(prob *a) {
	int i;
	
	printf("%i states; %lf norm\n", a->n, a->norm);
	for(i=0;i<a->n;i++) {
		printf("%lf ", a->p[i]);
	}
	printf("\n");
}

void print_jprob(j_prob *a) {
	int i,j;
	
	printf("%i x %i states; %lf norm\n", a->n, a->m, a->norm);
	for(i=0;i<a->n;i++) {
		for(j=0;j<a->m;j++) {
			printf("%le ", a->p[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

gsl_rng *rng() {
	unsigned long s;
	FILE *fn;
	gsl_rng *r_global;
	
	fn = fopen("/dev/urandom", "rb"); 		
	if (fread(&s, sizeof(s), 1, fn) != 1) 
		exit(-1); /* Failed! */
	r_global=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r_global, s);
	fclose(fn);
	
	return r_global;
}

prob* import_prob(char* file_name) {
        int i, ik, k;
        //double ik;
        FILE* fin;
        fin = fopen(file_name, "r");
        if (!fin) {
                printf("Error reading input file: %c \n", file_name);
                exit(-1);
        }
        fscanf(fin, "%i", &k);
        printf("There are %d state(s) \n", k);
        prob* pi_samp = new_prob(k);
        for (i = 0; i < k; i++) {
                fscanf(fin, "%i", &ik);
                printf("Reading %dth value: %f \n", i+1, (double)ik);
                pi_samp->p[i] = (double)ik;
        }
        return pi_samp;
}

prob* prob_from_array(double* in_array, int n) {
	prob* a_prob = new_prob(n);
	int i;
	for (i = 0; i < n; i++) {
		a_prob->p[i] = (double)in_array[i];
	}
	return a_prob;
}

void convert(double* in, int n1, int n2, double*** out) {
    double** data;
    int i;

    data = malloc(sizeof(*data) * n1);
    for(i = 0; i < n1; ++i) {
            data[i] = in + i * n2;
    }       

    *out = data;
}

j_prob* jprob_from_array(double* j_in_array, int n1, int n2) {
   double** tmp_p;
   j_prob* a_jprob = new_jprob(n1, n2);
   convert(j_in_array, n1, n2, &tmp_p);

   int i,j;
   for (i = 0; i < n1; i++) {
      for (j = 0; j < n2; j++) {
         a_jprob->p[i][j] = tmp_p[i][j];
      }
   }
   //printf("Done reading probabilities \n");
   //print_jprob(a_jprob);
   return a_jprob;
}


void return_array(double* array, int n) {
       	int i;
	for (i = 0; i < n; i++) {
		array[i] = array[i];
	} 
}

void print_bs(double* bs, prob* a_prob) {
	printf("Entropy: %lf (%lf, %lf) (NSB: %lf)\n", bs[0], bs[1], bs[2], entropy_nsb(a_prob)); 
}
