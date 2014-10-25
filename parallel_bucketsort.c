#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <mpi.h>
#include <prand.h>

/* 
 * Sequential bucketsort for randomly generated uint64_tegers.
 *
 */

#define VERBOSE 1

const uint64_t DEBUG_LEVEL = DEBUG;
uint64_t count_array_grows = 0;

int compareTo(const void *, const void *);

const uint64_t TRUE = 1;
const uint64_t Fnum_listLSE = 0;

static uint64_t min_range,max_range,box_share,n; 
static int	pid,nproc;
static uint64_t *num_list;
static uint64_t **bucket; // bucket[i] = array holding the ith bucket
static uint64_t *capacity; // capacity[i] = capacity of ith bucket
static uint64_t *size; // size[i] = next free location in ith bucket

void checkIfSorted(uint64_t *array, uint64_t n){
	uint64_t i;
	uint64_t sorted;

	sorted = TRUE;
	for (i=0; i<n-1; i++){
		if (array[i] > array[i+1]) {
				sorted = Fnum_listLSE;
				break;
		}
	}

	if (sorted) {
		if (DEBUG_LEVEL >= 1){
			fprintf(stderr, "array is sorted\n");
		}
	} else {
			fprintf(stderr, "Error: array is not sorted!\n");
	}
}



/*
 * Generate n numbers using the given seed
 */
void generateInput(uint64_t *num_list, uint64_t n, uint64_t seed){
	uint64_t i;

	srandom(seed);
	for (i=0; i<n; i++) {
		num_list[i] = random();
	}
}


/*
 * Pruint64_t the array, one element per line
 */
void printnum_listrray(uint64_t *num_list, uint64_t n){
	uint64_t i;

	for (i=0; i<n; i++){
		printf(" %16"PRIu64" \n", num_list[i]);
	}
}


/*
 *
 * Insert a given value uint64_to the specified bucket
 */

void insertInBucket(uint64_t value, uint64_t bucketIndex){
	uint64_t *tmp;

	if (size[bucketIndex] == capacity[bucketIndex]){
		//grow the bucket array
		tmp = malloc(sizeof(uint64_t)*(2*capacity[bucketIndex]));
		memcpy(tmp, bucket[bucketIndex], capacity[bucketIndex]*sizeof(uint64_t));
		free(bucket[bucketIndex]);
		bucket[bucketIndex] = tmp;
		capacity[bucketIndex] = 2 * capacity[bucketIndex];
		count_array_grows++;
		if (DEBUG_LEVEL >= 1) {
			fprintf(stderr, "Growing bucket %"PRIu64" from %"PRIu64" to %"PRIu64" elements\n",
					bucketIndex, capacity[bucketIndex]/2, capacity[bucketIndex]);
		}
	}

	bucket[bucketIndex][size[bucketIndex]] = value;
	size[bucketIndex]++;
}


/*
 * compareTo function for using qsort
 * returns  -ve if *x < *y, 0 if *x == *y, +ve if *x > *y
 */
int compareTo(const void *x, const void *y){
	return ((*(uint64_t *)x) - (*(uint64_t *)y));
}


/*
 * Sort indiviual bucket using quick sort from std C library
 */

void sortEachBucket(uint64_t numBuckets){
	uint64_t i;

	for (i=0; i<numBuckets; i++)
	{
		qsort(bucket[i], size[i], sizeof(uint64_t), compareTo);
	}
	if (DEBUG_LEVEL >= 2)
		for (i=0; i<numBuckets; i++) {
			fprintf(stderr, "bucket %"PRIu64" has %"PRIu64" elements\n", i, size[i]);
		}
}


/* 
 * Combine all buckets back uint64_to the original array to finish the sorting
 *
 */
void combineBuckets(uint64_t *num_list, uint64_t n, uint64_t numBuckets){
	uint64_t i;

	uint64_t start = 0;
	for (i=0; i<numBuckets; i++) {
		memcpy(num_list+start, bucket[i], sizeof(uint64_t)*size[i]);
		start = start + size[i];
		free(bucket[i]);
	}
	free(bucket);
}

/*
 * Use bucketsort to sort n uniformly distributed numbers in the range [0..2^31-1].
 * Input: uint64_t *num_list: array of uint64_ts num_list[0..n-1]
 *        uint64_t n: number of elements in the input array
 *        uint64_t numBuckets: number of buckets to use
 *
 */

void parallelBucketsort(uint64_t *num_list, uint64_t n, uint64_t numBuckets){
	uint64_t share;
	uint64_t i;
	uint64_t bucketRange;
	uint64_t bucketIndex;

	share = n / numBuckets;
	share = share + (share * 11)/100; // 11% extra for overflow

	capacity = (uint64_t *) malloc(sizeof(uint64_t)*numBuckets);
	size = (uint64_t *) malloc(sizeof(uint64_t)*numBuckets);
	bucket = (uint64_t **) malloc(sizeof(uint64_t *)* numBuckets);

	for (i=0; i<numBuckets; i++) {
		bucket[i] = (uint64_t *) malloc(sizeof(uint64_t)*share);
		capacity[i] = share;
		size[i] = 0;
	}

	bucketRange = RAND_MAX/numBuckets;
	for (i=0; i<n; i++){
		bucketIndex = num_list[i]/bucketRange;
		if (bucketIndex > numBuckets - 1)
				bucketIndex = numBuckets - 1;
		insertInBucket(num_list[i], bucketIndex);
	}

	sortEachBucket(numBuckets);
	combineBuckets(num_list, n, numBuckets);
	free(capacity);
	free(size);
}


void pruint64_t_usage(char * program){
	fprintf(stderr, "Usage %s <n, must be > 1> <#buckets, must between 1 and n> <random seed>\n", program);
}

void gen_ranges(){
	if(0 == pid){
		min_range = 0;
	}else{
		min_range = pid*box_share;
	}
	if(pid == (nproc-1)){
		max_range = n;
	}else{
		max_range = (pid+1)*box_share;
	}
#if VERBOSE == 1
	printf("pid %d has range %"PRIu64"-%"PRIu64"\n",pid,min_range,max_range);
#endif
}

int main(int argc, char **argv){
	uint64_t numBuckets;
	uint64_t seed;
	uint64_t start_time,total_time;
        MPI_Status status;

	if (argc != 4) {
		pruint64_t_usage(argv[0]);
		exit(1);
	}

	n = atoi(argv[1]);
	numBuckets = atoi(argv[2]);
	seed = atoi(argv[3]);

	MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&pid);
        start_time = MPI_Wtime();

	box_share = n/nproc;
	srandom(seed);
	unrankRand(pid*box_share);
	if(pid == (nproc - 1)){
		box_share += numBuckets%nproc;
	}
	num_list = malloc(sizeof(uint64_t)*box_share);
	gen_ranges();
	uint64_t x;
	for(x=0;x<box_share;x++){
		uint64_t temp = random();
			
		num_list[x] = temp;
	}
	if ((numBuckets < 1) || (n < 1) || (n < numBuckets)) {
		pruint64_t_usage(argv[0]);
		exit(1);
	}
			
	num_list = (uint64_t *) malloc(sizeof(uint64_t) * n);

	generateInput(num_list, n, seed);
 	if (DEBUG_LEVEL >= 3){ 
		printnum_listrray(num_list,n);
	}

	parallelBucketsort(num_list, n, numBuckets);
	checkIfSorted(num_list, n);
	if (DEBUG_LEVEL >= 1){
		printf("Number of array grows is %"PRIu64"\n", count_array_grows);
	}
 	if (DEBUG_LEVEL >= 3) {
 		printnum_listrray(num_list,n);
	}
	total_time = start_time - MPI_Wtime();
	printf("bucketsort: n = %"PRIu64"  m = %"PRIu64" buckets seed = %"PRIu64" time = %"PRIu64" seconds\n",
		   n, numBuckets, seed,	total_time);

	free(num_list);

	MPI_Finalize();
	exit(0);
}

/* vim: set ts=4: */
