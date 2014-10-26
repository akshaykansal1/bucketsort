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
const uint64_t FALSE = 0;

static uint64_t min_range,max_range,box_share,n; 
static int	pid,nproc;
static uint64_t *num_list;
static uint64_t **bucket; // bucket[i] = array holding the ith bucket
static uint64_t *capacity; // capacity[i] = capacity of ith bucket
static uint64_t *size; // size[i] = next free location in ith bucket

typedef struct send_list{
	int pid;
	uint64_t *num_list;
}send_list;

void checkIfSorted(uint64_t *array, uint64_t n){
	uint64_t i;
	uint64_t sorted;

	sorted = TRUE;
	for (i=0; i<n-1; i++){
		if (array[i] > array[i+1]) {
				sorted = FALSE;
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
void generateInput(){
	uint64_t i,k=0;
	uint64_t pnum = box_share/nproc;
	for (i=0; i<pnum; i++) {
		num_list[i] = random()%max_range;
		if(num_list[i]<min_range){
			num_list[i] += min_range;
		}
#if VERBOSE > 1
		printf("generated %"PRIu64" for range %"PRIu64"-%"PRIu64"\n",num_list[i],min_range,max_range);
#endif
	}
	send_list *to_send = malloc(sizeof(send_list)*(nproc-1));
	uint64_t tpnum = box_share + box_share%nproc;
	to_send[0].num_list = malloc(sizeof(uint64_t)*pnum);
	if(-1==(pid-1)){
		to_send[0].pid = nproc-1;
	}else{	
		to_send[0].pid = (pid-1)%nproc;
	}
	for(i=0;i<tpnum;i++){
		to_send[0].num_list[i] = random();
	}
	//uint64_t **to_send = malloc(sizeof(uint64_t*)*(nproc-1));
	for(i=1;i<(nproc-1);i++){
		to_send[i].num_list = malloc(sizeof(uint64_t)*pnum);
		to_send[i].pid = (pid+i)%nproc;
		printf("\t%d: pid to send to is %d\n",pid,to_send[i].pid);
		for(k=0;k<pnum;k++){
			to_send[i].num_list[k] = random();
		}
	}
	//need to send the data
	for(i=0;i<(nproc-1);i++){ 
		/*MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			int dest, int sendtag, void *recvbuf, int recvcount,
			MPI_Datatype recvtype, int source, int recvtag,
			MPI_Comm comm, MPI_Status *status);*/
	}
	free(to_send);
}


/*
 * Pruint64_t the array, one element per line
 */
void print_array(uint64_t *num_list, uint64_t n){
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


void print_usage(char * program){
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
		print_usage(argv[0]);
		exit(1);
	}

	n = atoi(argv[1]);
	numBuckets = atoi(argv[2]);
	seed = atoi(argv[3]);

	if ((numBuckets < 1) || (n < 1) || (n < numBuckets)) {
		print_usage(argv[0]);
		exit(1);
	}
			
	MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&pid);
        start_time = MPI_Wtime();

	srandom(seed);
	unrankRand(pid*box_share);
	
	box_share = n/nproc;
	if(pid == (nproc - 1)){
		box_share += n%nproc;
	}
#if VERBOSE ==1
	printf("pid %d has share %"PRIu64"\n",pid,box_share);
#endif
	num_list = malloc(sizeof(uint64_t)*box_share);
	gen_ranges();
	
	generateInput();
 	
	if (DEBUG_LEVEL >= 3){ 
		print_array(num_list,n);
	}

	parallelBucketsort(num_list, n, numBuckets);
	checkIfSorted(num_list, n);
	if (DEBUG_LEVEL >= 1){
		printf("Number of array grows is %"PRIu64"\n", count_array_grows);
	}
 	if (DEBUG_LEVEL >= 3) {
 		print_array(num_list,n);
	}
	total_time = start_time - MPI_Wtime();
	printf("bucketsort: n = %"PRIu64"  m = %"PRIu64" buckets seed = %"PRIu64" time = %"PRIu64" seconds\n",
		   n, numBuckets, seed,	total_time);

	free(num_list);

	MPI_Finalize();
	exit(0);
}

/* vim: set ts=4: */
