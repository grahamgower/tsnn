#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <glob.h>
#include <errno.h>
#include <sys/types.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include <tskit.h>
#include "kann.h"

static int verbose = 0;

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
	uint i, j;
	uint n_samples;
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
#define PATTERN_SUFFIX "/*.trees"
	char pat[strlen(dir) + strlen(PATTERN_SUFFIX) + 1];
	glob_t pglob = {0,};
	float *genotypes;
	int ret;
	int tret;
	uint i;
	uint n_files;
	int n_samples;

	*_genotypes = NULL;
	*_n_files = 0;
	*_n_cols = 0;

	sprintf(pat, "%s%s", dir, PATTERN_SUFFIX);

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
			case 1:
				ret = 1;
				goto err2;
		}

		n_files++;
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
	if (verbose)
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

/*
 * Copy/paste/modify of kann_train_fnn1().
 *  - Print training info as tab-separated columns.
 *  - Accept vector of learning rates, and use lr[i % lr_size] in epoch i.
 */
int kann_train_tsnn1(kann_t *ann,
		float *lr, int lr_size,
		int mini_size, int max_epoch, int max_drop_streak, float frac_val,
		int n, float **_x, float **_y)
{
	int i, j, *shuf, n_train, n_val, n_in, n_out, n_var, n_const, drop_streak = 0, min_set = 0;
	float **x, **y, *x1, *y1, *r, min_val_cost = FLT_MAX, *min_x, *min_c;

	n_in = kann_dim_in(ann);
	n_out = kann_dim_out(ann);
	if (n_in < 0 || n_out < 0) return -1;
	n_var = kann_size_var(ann);
	n_const = kann_size_const(ann);
	r = (float*)calloc(n_var, sizeof(float));
	shuf = (int*)malloc(n * sizeof(int));
	x = (float**)malloc(n * sizeof(float*));
	y = (float**)malloc(n * sizeof(float*));
	kann_shuffle(n, shuf);
	for (j = 0; j < n; ++j)
		x[j] = _x[shuf[j]], y[j] = _y[shuf[j]];
	n_val = (int)(n * frac_val);
	n_train = n - n_val;
	min_x = (float*)malloc(n_var * sizeof(float));
	min_c = (float*)malloc(n_const * sizeof(float));

	x1 = (float*)malloc(n_in  * mini_size * sizeof(float));
	y1 = (float*)malloc(n_out * mini_size * sizeof(float));
	kann_feed_bind(ann, KANN_F_IN,    0, &x1);
	kann_feed_bind(ann, KANN_F_TRUTH, 0, &y1);

	printf("epoch\ttraining_cost\ttraining_error\tvalidation_cost\tvalidation_error\n");

	for (i = 0; i < max_epoch; ++i) {
		int n_proc = 0, n_train_err = 0, n_val_err = 0, n_train_base = 0, n_val_base = 0;
		double train_cost = 0.0, val_cost = 0.0;
		kann_shuffle(n_train, shuf);
		kann_switch(ann, 1);
		while (n_proc < n_train) {
			int b, c, ms = n_train - n_proc < mini_size? n_train - n_proc : mini_size;
			for (b = 0; b < ms; ++b) {
				memcpy(&x1[b*n_in],  x[shuf[n_proc+b]], n_in  * sizeof(float));
				memcpy(&y1[b*n_out], y[shuf[n_proc+b]], n_out * sizeof(float));
			}
			kann_set_batch_size(ann, ms);
			train_cost += kann_cost(ann, 0, 1) * ms;
			c = kann_class_error(ann, &b);
			n_train_err += c, n_train_base += b;
			kann_RMSprop(n_var, lr[i % lr_size], 0, 0.9f, ann->g, ann->x, r);
			n_proc += ms;
		}
		train_cost /= n_train;
		kann_switch(ann, 0);
		n_proc = 0;
		while (n_proc < n_val) {
			int b, c, ms = n_val - n_proc < mini_size? n_val - n_proc : mini_size;
			for (b = 0; b < ms; ++b) {
				memcpy(&x1[b*n_in],  x[n_train+n_proc+b], n_in  * sizeof(float));
				memcpy(&y1[b*n_out], y[n_train+n_proc+b], n_out * sizeof(float));
			}
			kann_set_batch_size(ann, ms);
			val_cost += kann_cost(ann, 0, 0) * ms;
			c = kann_class_error(ann, &b);
			n_val_err += c, n_val_base += b;
			n_proc += ms;
		}
		if (n_val > 0) val_cost /= n_val;
		if (n_train_base && n_val_base && n_val > 0) {
			printf("%d\t%g\t%g\t%g\t%g\n", i+1,
					train_cost, ((float)n_train_err) / n_train,
					val_cost, ((float)n_val_err) / n_val);
		}
		if (i >= max_drop_streak && n_val > 0) {
			if (val_cost < min_val_cost) {
				min_set = 1;
				memcpy(min_x, ann->x, n_var * sizeof(float));
				memcpy(min_c, ann->c, n_const * sizeof(float));
				drop_streak = 0;
				min_val_cost = (float)val_cost;
			} else if (++drop_streak >= max_drop_streak)
				break;
		}
	}
	if (min_set) {
		memcpy(ann->x, min_x, n_var * sizeof(float));
		memcpy(ann->c, min_c, n_const * sizeof(float));
	}

	free(min_c); free(min_x); free(y1); free(x1); free(y); free(x); free(shuf); free(r);
	return i;
}

kann_t *
model_gen(int n_rows, int n_cols, int n_conv2d, int n_conv1d, int n_filt,
		int k_size, int stride, int fc_size, float dropout)
{
	kad_node_t *t;
	int i;

	/* The input is a 4D array with the four dimensions being:
	 * mini-batch size, number of channels, height and width.
	 */
	t = kad_feed(4, 1, 1, n_rows, n_cols), t->ext_flag += KANN_F_IN;

	for (i=0; i<n_conv2d; i++) {
		// `k_size`x1 kernel; 1x1 stride; 0x0 padding
		t = kad_relu(kann_layer_conv2d(t, n_filt, k_size, 1, stride, 1, 0, 0));
		if (dropout > 0)
			t = kann_layer_dropout(t, dropout);
	}

	// combine haplotypes
	t = kad_reduce_mean(t, 3);

	for (i=0; i<n_conv1d; i++) {
		t = kad_relu(kann_layer_conv1d(t, n_filt, k_size, stride, 0));
		if (dropout > 0)
			t = kann_layer_dropout(t, dropout);
	}

	// combine loci
	t = kad_relu(kann_layer_dense(t, fc_size));

	if (dropout > 0)
		t = kann_layer_dropout(t, dropout);

	return kann_new(kann_layer_cost(t, 1, KANN_C_CEB), 0);
}

void
usage(const char *argv0)
{
	fprintf(stderr, "usage: %s [...] tsdir1 tsdir2\n", argv0);
	fprintf(stderr, "  -c INT     Number of conv2d layers [2].\n");
	fprintf(stderr, "  -C INT     Number of conv1d layers [1].\n");
	fprintf(stderr, "  -f INT     Number of filters per conv layer [16].\n");
	fprintf(stderr, "  -k INT     Kernel size of each conv layer [8].\n");
	fprintf(stderr, "  -S INT     Stride of each conv layer [2]\n");
	fprintf(stderr, "  -F INT     Number of nodes in fully connected layer [4].\n");
	fprintf(stderr, "  -d FLOAT   Dropout after conv and fc layers [0.2].\n");
	//fprintf(stderr, "  -i FILE    Input model file [NULL]\n");
	fprintf(stderr, "  -o FILE    Output model file [NULL]\n");
	//fprintf(stderr, "  -l FLOAT   Learning rate [0.01]\n");
	fprintf(stderr, "  -m INT     Max number of training epochs [100]\n");
	fprintf(stderr, "  -r INT     Row length (haplotype squashing) [128]\n");
	fprintf(stderr, "  -s INT     Seed for random number generator [31415]\n");
	fprintf(stderr, "  -t INT     Number of threads to use [1]\n");
	fprintf(stderr, "  -v         Increase verbosity\n");
	exit(1);
}

int
main(int argc, char **argv)
{
	kann_t *ann;
	float *gts_a, *gts_b;
	float *labels;
	float **x, **y;
	char *fn_in = NULL, *fn_out = NULL;
	char *dir_a, *dir_b;
	int n_files_a, n_files_b;
	int opt;
	int seed = 31415;
	uint n_rows = 128; // n_sites dimension
	uint n_cols, n_cols2; // n_samples dimension

	// model params
	int n_conv2d = 2;
	int n_conv1d = 1;
	int n_filt = 16;
	int k_size = 8;
	int stride = 2;
	int fc_size = 4;
	float dropout = 0.2f;

	// training params
	int n_threads = 1;
	int mini_size = 64;
	int max_epoch = 100;
	int max_drop_streak = 10;
	float lr[] = {0.01f, 0.005f, 0.001f};  // sawtooth learning rate schedule
	float frac_val = 0.1f;

	while ((opt = getopt(argc, argv, "c:C:f:F:k:d:i:o:l:m:r:s:S:t:v")) != -1) {
		switch(opt) {
			case 'c':
				n_conv2d = atoi(optarg);
				break;
			case 'C':
				n_conv1d = atoi(optarg);
				break;
			case 'f':
				n_filt = atoi(optarg);
				break;
			case 'F':
				fc_size = atoi(optarg);
				break;
			case 'k':
				k_size = atoi(optarg);
				break;
			case 'd':
				dropout = atof(optarg);
				break;
			// TODO: do something useful with input model
			/*case 'i':
				fn_in = optarg;
				break;*/
			case 'o':
				fn_out = optarg;
				break;
			/*case 'l':
				lr = atof(optarg);
				break;*/
			case 'm':
				max_epoch = atoi(optarg);
				break;
			case 'r':
				n_rows = atoi(optarg);
				break;
			case 's':
				seed = atoi(optarg);
				break;
			case 'S':
				stride = atoi(optarg);
				break;
			case 't':
				n_threads = atoi(optarg);
				break;
			case 'v':
				verbose++;
				break;
			default:
				usage(argv[0]);
				// NORETURN
		}
	}

	if (argc - optind != 2 || (fn_in && argc - optind != 1)) {
		usage(argv[0]);
	}

	dir_a = argv[optind];
	dir_b = argv[optind+1];

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
	kann_srand(seed);
	ann = model_gen(n_rows, n_cols, n_conv2d, n_conv1d, n_filt, k_size,
			stride, fc_size, dropout);

	if (verbose) {
		uint64_t cells = (n_files_a + n_files_b) * n_rows * n_cols;
		fprintf(stderr, "data: %.2f MiB across %d ts files\n",
				cells * sizeof(float) / 1024.0 / 1024.0,
				n_files_a + n_files_b);
		fprintf(stderr, "model: %d variables\n", kann_size_var(ann));
		fprintf(stderr, "dim(in) = %d, dim(out) = %d\n",
				kann_dim_in(ann), kann_dim_out(ann));
		kad_print_graph(stderr, ann->n, ann->v);
	}

	if (n_threads > 1)
		kann_mt(ann, n_threads, mini_size);

	kann_train_tsnn1(ann, lr, sizeof(lr)/sizeof(lr[0]),
			mini_size, max_epoch, max_drop_streak, frac_val,
			n_files_a + n_files_b, x, y);

	if (fn_out)
		kann_save(fn_out, ann);
	kann_delete(ann);

	free(x);
	free(y);
	free(gts_a);
	free(gts_b);
	free(labels);

	return 0;
}
