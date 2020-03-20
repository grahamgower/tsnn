#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glob.h>
#include <errno.h>
#include <sys/types.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include <tskit.h>
#include "kann.h"

int
ts_num_samples(const char *fn, tsk_table_collection_t *tables)
{
	tsk_treeseq_t ts;
	int ret;
	int tret;

	if ((tret = tsk_table_collection_load(tables, fn, TSK_NO_INIT)) != 0) {
		fprintf(stderr, "tsk_table_collection_load: %s: %s\n",
				fn, tsk_strerror(tret));
		ret = -1;
		goto err0;
	}

	if ((tret = tsk_treeseq_init(&ts, tables, 0)) != 0) {
		fprintf(stderr, "tsk_treeseq_init: %s: %s\n",
				fn, tsk_strerror(tret));
		ret = -1;
		goto err0;
	}

	ret = tsk_treeseq_get_num_samples(&ts);

	tsk_treeseq_free(&ts);
err0:
	return ret;
}

int
ts_load(const char *fn, tsk_table_collection_t *tables,
		float *m,
		uint n_rows, // squash all sites to this size
		uint n_cols) // number of samples
{
	tsk_treeseq_t ts;
	tsk_vargen_t vg;
	tsk_variant_t *variant;
	int ret;
	int tret;
	int i, j;
	int n_samples;
	double sequence_length;

	if ((tret = tsk_table_collection_load(tables, fn, TSK_NO_INIT)) != 0) {
		fprintf(stderr, "tsk_table_collection_load: %s: %s\n",
				fn, tsk_strerror(tret));
		ret = 1;
		goto err0;
	}

	if ((tret = tsk_treeseq_init(&ts, tables, 0)) != 0) {
		fprintf(stderr, "tsk_treeseq_init: %s: %s\n",
				fn, tsk_strerror(tret));
		ret = -1;
		goto err0;
	}

	if ((tret = tsk_vargen_init(&vg, &ts, NULL, 0, 0)) != 0) {
		fprintf(stderr, "tsk_vargen_init: %s: %s\n",
				fn, tsk_strerror(tret));
		ret = -1;
		goto err1;
	}

	n_samples = tsk_treeseq_get_num_samples(&ts);
	sequence_length = tsk_treeseq_get_sequence_length(&ts);

	if (n_cols != n_samples) {
		fprintf(stderr, "%s: has %d samples, not %d\n",
				fn, n_samples, n_cols);
		ret = 1;
		goto err2;
	}

	for (j=0; (tret = tsk_vargen_next(&vg, &variant)) == 1; j++) {
		uint8_t *v8 = variant->genotypes.u8;
		double pos = variant->site->position;
		uint row = n_rows * pos / sequence_length;
		for (i=0; i<n_samples; i++) {
			uint idx = row * n_samples + i;
			assert(idx < n_rows * n_samples);
			m[idx] += v8[i];
		}
	}

	ret = 0;
err2:
	tsk_vargen_free(&vg);
err1:
	tsk_treeseq_free(&ts);
err0:
	return ret;
}

int
_glob_err(const char *epath, int _errno)
{
	fprintf(stderr, "%s: %s\n", epath, strerror(_errno));
	return 0;
}

int
tsdir_load(const char *dir, float **_genotypes, int *_n_files,
		uint n_rows,
		uint *_n_cols)
{
	tsk_table_collection_t tables;
	char pat[strlen(dir)+3];
	glob_t pglob = {0,};
	float *genotypes;
	int ret;
	int tret;
	int i;
	int n_files;
	int n_samples;

	*_genotypes = NULL;
	*_n_files = 0;
	*_n_cols = 0;

	sprintf(pat, "%s/*", dir);

	switch (glob(pat, GLOB_NOSORT, _glob_err, &pglob)) {
		case GLOB_NOSPACE:
			errno = ENOMEM;
			perror("glob");
			ret = 1;
			goto err0;
		case GLOB_ABORTED:
			fprintf(stderr,
				"glob: %s: aborting due to earlier errors.\n",
				dir);
			ret = 1;
			goto err0;
		case GLOB_NOMATCH:
			fprintf(stderr,
				"glob: %s: no files found.\n",
				dir);
			ret = 1;
			goto err0;
	}

	if ((tret = tsk_table_collection_init(&tables, 0)) != 0) {
		fprintf(stderr, "tsk_table_collection_init: %s\n",
				tsk_strerror(tret));
		ret = 1;
		goto err0;
	}

	n_samples = ts_num_samples(pglob.gl_pathv[0], &tables);
	if (n_samples < 0) {
		ret = 1;
		goto err1;
	}

	genotypes = calloc(pglob.gl_pathc * n_rows * n_samples, sizeof(float));
	if (genotypes == NULL) {
		perror("calloc");
		ret = 1;
		goto err1;
	}

	n_files = 0;
	for (i=0; i<pglob.gl_pathc; i++) {
		float *m = &genotypes[n_files * (n_rows * n_samples)];
		switch (ts_load(pglob.gl_pathv[i], &tables, m,
					n_rows, n_samples)) {
			case -1:
				ret = 1;
				goto err2;
			case 1:
				continue;
		}

		n_files++;
	}

	if (n_files == 0) {
		fprintf(stderr, "%s: no valid tskit files found.\n", dir);
		ret = 1;
		goto err2;
	}

	if (n_files < pglob.gl_pathc) {
		// TODO: is it worth realloc'ing here?
	}

	*_genotypes = genotypes;
	*_n_files = n_files;
	*_n_cols = n_samples;
	ret = 0;
err2:
	if (ret) {
		for (i=0; i<n_files; i++)
			free(genotypes+i);
		free(genotypes);
	}
err1:
	tsk_table_collection_free(&tables);
err0:
	globfree(&pglob);
	return ret;
}

// Mean centre and divide by standard deviation.
void
normalise(float **x, int n, int dim)
{
	int i, j;
	float mean, std_inv, sum = 0;
	for (i=0; i<n; i++) {
		float *xx = x[i];
		for (j=0; j<dim; j++)
			sum += xx[j];
	}
	mean = sum / (n * dim);
	sum = 0;
	for (i=0; i<n; i++) {
		float *xx = x[i];
		float tmp;
		for (j=0; j<dim; j++) {
			tmp = xx[j] - mean;
			sum += tmp*tmp;
			xx[j] = tmp;
		}
	}
	std_inv = 1.0 / sqrt(sum / (n*dim - 1));
	for (i=0; i<n; i++) {
		float *xx = x[i];
		for (j=0; j<dim; j++) {
			xx[j] *= std_inv;
		}
	}
	printf("transforming data: mean=%f, stddev=%f\n", mean, 1/std_inv);
}

int
prepare_data(float *gts_a, float *gts_b, int n_files_a, int n_files_b, int in_dim,
		float **_labels, float ***_x,  float ***_y)
{
	float **x;
	float **y;
	float *labels;
	int n = n_files_a + n_files_b;
	int i;

	x = malloc(n * sizeof(float*));
	if (x == NULL) {
		perror("malloc");
		return 1;
	}

	for (i=0; i<n_files_a; i++)
		x[i] = &gts_a[i * in_dim];
	for (i=0; i<n_files_b; i++)
		x[n_files_a + i] = &gts_b[i * in_dim];

	normalise(x, n, in_dim);

	labels = malloc(n * sizeof(float));
	if (labels == NULL) {
		perror("malloc");
		return 1;
	}

	for (i=0; i<n_files_a; i++)
		labels[i] = 0.0;
	for (i=0; i<n_files_b; i++)
		labels[n_files_a + i] = 1.0;

	y = malloc(n * sizeof(float*));
	if (y == NULL) {
		perror("malloc");
		return 1;
	}

	for (i=0; i<n; i++)
		y[i] = &labels[i];

	*_x = x;
	*_y = y;
	*_labels = labels;
	return 0;
}

void
kad_print_dot(FILE *fp, int n, kad_node_t **v)
{
	int i, j;
	for (i = 0; i < n; ++i) v[i]->tmp = i;
	fprintf(fp, "digraph {\n");
	for (i = n - 1; i >= 0; --i) {
		kad_node_t *p = v[i];
		if (p->op > 0) fprintf(fp, "\t%d [label=\"%s\"]\n", i, kad_op_name[p->op]);
		for (j = 0; j < p->n_child; ++j)
			fprintf(fp, "\t%d -> %d\n", p->child[j]->tmp, i);
		if (p->pre) fprintf(fp, "\t%d -> %d [style=dotted,weight=0,constraint=false]\n", i, p->pre->tmp);
	}
	fprintf(fp, "}\n");
	for (i = 0; i < n; ++i) v[i]->tmp = 0;
}

kann_t *
model_gen(int n_rows, int n_cols, int n_conv, int n_filt, int k_size, float dropout)
{
	kad_node_t *t;
	int i;

	/* The input is a 4D array with the four dimensions being:
	 * mini-batch size, number of channels, height and width.
	 */
	t = kad_feed(4, 1, 1, n_rows, n_cols), t->ext_flag += KANN_F_IN;

	for (i=0; i<n_conv; i++) {
		// `k_size`x1 kernel; 1x1 stride; 0x0 padding
		t = kad_relu(kann_layer_conv2d(t, n_filt, k_size, 1, 1, 1, 0, 0));
		if (dropout > 0)
			t = kann_layer_dropout(t, dropout);
	}

	// combine haplotypes
	t = kad_reduce_mean(t, 3);

	for (i=0; i<n_conv; i++) {
		t = kad_relu(kann_layer_conv1d(t, n_filt, k_size, 1, 0));
		if (dropout > 0)
			t = kann_layer_dropout(t, dropout);
	}

	// combine loci
	t = kad_reduce_mean(t, 2);
	
	t = kad_relu(kann_layer_dense(t, 16));
	return kann_new(kann_layer_cost(t, 1, KANN_C_CEB), 0);
}

int
main(int argc, char **argv)
{
	const int verbose = 1;
	float *gts_a, *gts_b;
	float *labels;
	int n_files_a, n_files_b;
	kann_t *ann;
	float **x, **y;
	const uint n_rows = 32; // n_sites dimension
	uint n_cols, n_cols2; // n_samples dimension

	if (argc != 3) {
		fprintf(stderr, "usage: %s tsdir1 tsdir2\n", argv[0]);
		return 1;
	}

	{
		// Verify that the model doesn't depend on the input image size.
		kann_t *nn1, *nn2;
		nn1 = model_gen(100, 100, 3, 16, 4, 0.3);
		nn2 = model_gen(200, 200, 3, 16, 4, 0.3);
		assert(kann_size_var(nn1) == kann_size_var(nn2));
		assert(kann_size_const(nn1) == kann_size_const(nn2));
		kann_delete(nn1);
		kann_delete(nn2);
	}

	const char *dir_a = argv[1];
	const char *dir_b = argv[2];

	if (tsdir_load(dir_a, &gts_a, &n_files_a, n_rows, &n_cols) != 0)
		return 2;
	if (tsdir_load(dir_b, &gts_b, &n_files_b, n_rows, &n_cols2) != 0)
		return 3;

	if (n_cols != n_cols2) {
		fprintf(stderr, "%s and %s having differing numbers of samples\n",
				dir_a, dir_b);
		return 4;
	}

	if (prepare_data(gts_a, gts_b, n_files_a, n_files_b, n_rows*n_cols,
				&labels, &x, &y) != 0)
		return 5;

	kad_trap_fe();
	kann_srand(31415);
	ann = model_gen(n_rows, n_cols, 2, 12, 3, 0.3);

	const int n = n_files_a + n_files_b;
	if (verbose) {
		uint64_t cells = n * n_rows * n_cols;
		printf("data: %.2f MiB across %d ts files\n",
				cells * sizeof(float) / 1024.0 / 1024.0, n);
		printf("model: %d variables\n", kann_size_var(ann));
		printf("dim(in) = %d, dim(out) = %d\n",
				kann_dim_in(ann), kann_dim_out(ann));

		kad_print_graph(stdout, ann->n, ann->v);
		//kad_print_dot(stdout, ann->n, ann->v);
	}

	const float lr = 0.01;
	const int mini_size = 50;
	const int max_epoch = 10;
	const int max_drop_streak = 10;
	const float frac_val = 0.1;
	const int n_threads = 1;

	if (n_threads > 1) kann_mt(ann, n_threads, mini_size);
	kann_train_fnn1(ann, lr, mini_size, max_epoch, max_drop_streak,
			frac_val, n, x, y);

//	kann_save("blah", ann);
	kann_delete(ann);

	free(x);
	free(y);
	free(gts_a);
	free(gts_b);
	free(labels);

	return 0;
}
