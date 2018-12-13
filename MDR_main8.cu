#include <iostream>
#include <math.h>
//#include "MDR.h"
//#include "MDR_kernel.cu"
//#include "MDR.cu"

#define ORDER 3
/*

#define BSx 2
#define BSy 1
#define BSz 1

#define GSx 2
#define GSy 1
#define GSz 1
*/
#define IDX(i,j,ld) (((i)*(ld))+(j))
#define imin(a,b) (a<b?a:b)

const int NUMCOMBS = 3200000;
const int BSx = 256;
//const int GSx = imin( 32, (NUMCOMBS+BSx-1) / BSx );
const int GSx = (NUMCOMBS+BSx-1) / BSx ;
#define NSNPS 20000 //rows of input file
#define NIND 6000  //columns of input file
#define THR 1

#define CV 3

#define TESTCOMB 999

#define mat_SNP_size (NIND * NSNPS * sizeof(int))
#define v_pheno_size (NIND * sizeof(int))
#define output_size (NIND * sizeof(int))
#define combinations_size (NUMCOMBS * ORDER * sizeof(int))
#define indices_size (NIND * sizeof(int))


char phenoFile[] = "pheno7";
char genoFile[] = "geno7";
char output[] = "output7";
char combFile[] = "combinations7";


struct controlscases {
int controls;
int cases;
};

/*
__device__ int dev_rand() {
	return rand();
}
*/

//#include "MDR.h"
//#include "MDR.cu"

__constant__ int dev_v_pheno_values[NIND];
__constant__ int dev_cv_indices[NIND];

__global__ void MDR( int* dev_SNP_values, int* dev_output_values, int* dev_combinations ) {
    
    
	//printf(" %d + %d * %d :", threadIdx.x, blockIdx.x, blockDim.x);
	//__shared__ float cache[BS][threadsPerBlock];
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	int thread_combination[ORDER]; //a combination (thread level)
	//retrieve the combination indices
	for (int i=0; i< ORDER; i++) {
		thread_combination[i] = *(dev_combinations + tid * ORDER + i);
	}
	//printf("thread with tid %d is assigned combination: <%d, %d, %d>\n", tid, thread_combination[0], 
	//thread_combination[1], thread_combination[2]); 
	
	//retrieve the genotype of each snp in the combination, from SNPvalues, for ALL individuals
	int thread_geno[NIND * ORDER];
	for (int i=0; i< ORDER; i++) {
		for (int j=0; j< NIND; j++) {
			thread_geno[i * NIND + j] = *(dev_SNP_values + NIND * thread_combination[i] + j);
		}
	}
	
	if (tid == TESTCOMB){
			for (int i=0; i< 20; i++) {
				printf("thread %d %d-th-snp's geno: %d\n", tid, i, thread_geno[i]);
			}
		}
	
	
	//CV loop
	for (int cv=0; cv<CV; cv++){
		

		if (tid == TESTCOMB){
			printf("\n******************************\nthread %d, iteration of CV %d/%d:\n\n", tid, cv, CV);
			printf("train interval:(0;%d) U (%d;%d)\n", int((cv/float(CV))*NIND), int(((cv+1)/float(CV))*NIND), NIND);
			printf("test interval: [%d;%d]\n", int((cv/float(CV))*NIND), int(((cv+1)/float(CV))*NIND));
		}
	
		struct controlscases thread_table[3][3][3];
	
		//replace this initialization?
		for (int i=0; i< 3; i++) {
			for (int j=0; j< 3; j++) {
				for (int k=0; k< 3; k++) {
					thread_table[i][j][k].controls = 0;
					thread_table[i][j][k].cases = 0;
					}}}
	
		//populate the 3^ORDER-tot-entries table
		int f,s,t,ind;
		for (int n=0; n< NIND; n++) { //first NIND_TRAIN of NIND are for train
			 if ((n >= int((cv/float(CV))*NIND)) && (n <= ((cv+1)/float(CV))*NIND )) //reserved for test
			 	continue;
			 ind = *(dev_cv_indices + n);
			 f = thread_geno[0 * NIND + ind]; //1st snp geno
			 s = thread_geno[1 * NIND + ind]; //2nd snp geno
			 t = thread_geno[2 * NIND + ind]; //3rd snp geno
			 if (int(*(dev_v_pheno_values + ind))) //get the pheno
			 	thread_table[f][s][t].cases += 1;
			 else
			 	thread_table[f][s][t].controls += 1;
		}
	
	
		//only a print
	
		if (tid == TESTCOMB){
		printf("\n***************\nthread %d\n:", tid);
		for (int i=0; i< 3; i++) {
			for (int j=0; j< 3; j++) {
				for (int k=0; k< 3; k++) {
					printf("thread_table[%d][%d][%d].controls=%d\n",i,j,k,thread_table[i][j][k].controls);
					printf("thread_table[%d][%d][%d].cases=%d\n",i,j,k,thread_table[i][j][k].cases);
					printf("\n");
				
				}
			}
		}
		}
	
	
	
		//moving two a two-dim variable
		int high_cases = 0;
		int high_controls = 0;
		int low_cases = 0;
		int low_controls = 0;
		int high_genos[3*3*3][3]; //content: 001,021,112,..,9xx,..: 3**ORDER strings, each 3 chars
		int c = 0;
		for (int i=0; i< 3; i++) {
			for (int j=0; j< 3; j++) {
				for (int k=0; k< 3; k++) {
					if ((thread_table[i][j][k].cases)/(thread_table[i][j][k].controls + 0.01) >= THR){
						high_cases += thread_table[i][j][k].cases;
						high_controls += thread_table[i][j][k].controls;
						if (tid == TESTCOMB){
							printf("tid %d (comb. <%d, %d, %d>),"
							" geno %d%d%d is HIGH\n",
							tid, thread_combination[0], thread_combination[1],
							thread_combination[2], i, j ,k);
						}
					
						high_genos[c][0] = i;
						high_genos[c][1] = j;
						high_genos[c][2] = k;
						c+=1;
					}
					else{
						//here in LOW also the case 0 controls 0 cases
						low_cases += thread_table[i][j][k].cases;
						low_controls += thread_table[i][j][k].controls;
						if (tid == TESTCOMB){
							printf("tid %d (comb. <%d, %d, %d>),"
							" geno %d%d%d is LOW\n",
							tid, thread_combination[0], thread_combination[1],
							thread_combination[2], i, j ,k);
						}
					
					}
				}
			}
		}
		high_genos[c][0] = 9; //end sequence, since high_genos only reports the high ones
	
		//printf("******************\n");
		float train_error = float(high_controls + low_cases)/float(high_cases + high_controls + low_cases + low_controls);
	
	
		if (tid == TESTCOMB){
		printf("snp comb. <%d, %d, %d> (tid %d) TRAIN error %1.5f = (%d+%d)/(%d+%d+%d+%d)\n", 
				thread_combination[0], thread_combination[1], thread_combination[2], tid,
				train_error, high_controls, low_cases, high_cases, high_controls, low_cases, low_controls);
		}

	
		//*****************
		//TESTING
		//*****************
		int high_cases_test = 0;
		int high_controls_test = 0;
		int low_cases_test = 0;
		int low_controls_test = 0;
		for (int n=0; n< NIND; n++) {
			 if ((n < int(cv/float(CV)*NIND)) || (n > int((cv+1)/float(CV)*NIND)) )//reserved for training
			 	continue;
			 ind = *(dev_cv_indices + n);
			 f = thread_geno[0 * NIND + ind]; //1st snp geno
			 s = thread_geno[1 * NIND + ind]; //2nd snp geno
			 t = thread_geno[2 * NIND + ind]; //3rd snp geno
			 int ph = int(*(dev_v_pheno_values + ind));
			 //check if fst is in high or low
			 for (int i=0; i< (3*3*3); i++) {
			 	 if (high_genos[i][0] == 9)
			 	 	break;
				 if  (high_genos[i][0] == f && high_genos[i][1] == s && high_genos[i][2] == t){
				 	if (ph)
				 		high_cases_test += 1;
				 	else
				 		high_controls_test += 1;
				 	break;
				 	
				 }
			}
		 	 //not in high_cases, it's low
		 	 if (ph)
			 	low_cases_test += 1;
			 else
			 	low_controls_test += 1;
			 
			
	
		}
		
		//printf("******************\n");
		float test_error = float(high_controls_test + low_cases_test)/float(high_cases_test + high_controls_test + low_cases_test + low_controls_test);
	
	
		if (tid == TESTCOMB){
		printf("snp comb. <%d, %d, %d> (tid %d) TEST error %1.5f = (%d+%d)/(%d+%d+%d+%d)\n", 
				thread_combination[0], thread_combination[1], thread_combination[2], tid,
				test_error, high_controls_test, low_cases_test, high_cases_test, high_controls_test, low_cases_test, low_controls_test);
		}
	

	}
}


static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "ERROR HANDLED: %s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
  	}
  }
  
void readintData(char *dataFile, unsigned int rows, unsigned int cols, int * data){
  FILE *fp;
  int *dp = data;
  int i;

  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    //exit(1);
  } 
  
  for (i=0; i<rows*cols; ++i){
	  fscanf(fp, "%d", dp);
	  dp++;
  } 
  fclose(fp);
  return;
}


void readCombinations(char *dataFile, int rows, int cols, int * data){
  FILE *fp;
  int *dp = data;
  int i;

  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    //exit(1);
  } 
  
  for (i=0; i<rows*cols; ++i){
	  fscanf(fp, "%d", dp);
	  dp++;
  } 
  fclose(fp);
  return;
}
  	
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


#define HANDLE_NULL( a ) {if (a == NULL) { \
                            printf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
                            exit( EXIT_FAILURE );}}

void print_cudaGetDeviceProperties(){
	cudaDeviceProp prop;
	int count;
	HANDLE_ERROR( cudaGetDeviceCount( &count ) );
	
	if (count == 0) {
        	printf("error in print_cudaGetDeviceProperties: no devices supporting CUDA.\n");
        	return;
    	}
	
	for (int i=0; i< count; i++) {
		HANDLE_ERROR( cudaGetDeviceProperties( &prop, i ) );
		printf( "   --- General Information for device %d ---\n", i );
		printf( "Name:  %s\n", prop.name );
		printf( "Compute capability:  %d.%d\n", prop.major, prop.minor );
		printf( "Clock rate:  %d\n", prop.clockRate );
		printf( "Device copy overlap:  " );
		if (prop.deviceOverlap)
		    printf( "Enabled\n" );
		else
		    printf( "Disabled\n");
		printf( "Kernel execution timeout :  " );
		if (prop.kernelExecTimeoutEnabled)
		    printf( "Enabled\n" );
		else
		    printf( "Disabled\n" );

		printf( "   --- Memory Information for device %d ---\n", i );
		printf( "Total global mem:  %ld\n", prop.totalGlobalMem );
		printf( "Total constant Mem:  %ld\n", prop.totalConstMem );
		printf( "Max mem pitch:  %ld\n", prop.memPitch );
		printf( "Texture Alignment:  %ld\n", prop.textureAlignment );

		printf( "   --- MP Information for device %d ---\n", i );
		printf( "Multiprocessor count:  %d\n",
		            prop.multiProcessorCount );
		printf( "Shared mem per mp:  %ld\n", prop.sharedMemPerBlock );
		printf( "Registers per mp:  %d\n", prop.regsPerBlock );
		printf( "Threads in warp:  %d\n", prop.warpSize );
		printf( "Max threads per block:  %d\n",
		            prop.maxThreadsPerBlock );
		printf( "Max thread dimensions:  (%d, %d, %d)\n",
		            prop.maxThreadsDim[0], prop.maxThreadsDim[1],
		            prop.maxThreadsDim[2] );
		printf( "Max grid dimensions:  (%d, %d, %d)\n",
		            prop.maxGridSize[0], prop.maxGridSize[1],
		            prop.maxGridSize[2] );
		printf( "\n" );
	}
	return;
}

/************************/
//MAIN
/************************/

int main(void){
	 
  	//print_cudaGetDeviceProperties(); 
  	
  	printf("*****************\n");
	printf("Multifactor Dimensionality Reduction\n");
	printf("*****************\n");	
	
	int* combinations = (int*)malloc(NUMCOMBS * ORDER * sizeof(int));
	int* dev_combinations;
	//int* dev_cv_indices;
  	int* cv_indices = (int*)malloc(NIND * sizeof(int));
  	int* dev_mat_SNP;
	int* dev_output;

	//generate all indices for CV loop
	for(int i=0;i<NIND;++i){
        	*(cv_indices + i) = i;
    		}
    		
		//permute r with Fisher-Yates shuffling algorithm
	for (int i = NIND; i >= 0; --i){
		//generate a random number [0, n-1]
		int j = rand() % (i+1);

		//swap the last element with element at random index
		int temp = *(cv_indices + i);
		*(cv_indices + i) = *(cv_indices + j);
		*(cv_indices + j) = temp;
	}

	
	//Allocate host memory 
	int* mat_SNP = (int*)malloc(mat_SNP_size); 
	int* v_pheno = (int*)malloc(v_pheno_size);
	int* output = (int*)malloc(output_size); //TODO
  	
  	//Read the matrix in host data
	readintData(genoFile, NSNPS, NIND, mat_SNP);
	printf("geno file read..\n");
	readintData(phenoFile, NIND, 1, v_pheno);
	printf("pheno file read..\n");
	
	//Read combinations
	readCombinations(combFile, NUMCOMBS, ORDER, combinations);
	printf("combinations file read..\n");

  	
  	//Allocate device memory
  	//constant memory?
	cudaMalloc((void**)&dev_mat_SNP, mat_SNP_size);
	//cudaMalloc((void**)&dev_v_pheno.values, dev_v_pheno.mem_size);
	
	cudaMalloc((void**)&dev_output, output_size);
  	cudaMalloc((void**)&dev_combinations, combinations_size  );
  	cudaMalloc((void**)&dev_cv_indices, indices_size);
  	
  	// Copy host memory to device
  	//HANDLE_ERROR( cudaMemcpy(dev_v_pheno.values, v_pheno.values, dev_v_pheno.mem_size, cudaMemcpyHostToDevice));
  	cudaMemcpyToSymbol( dev_v_pheno_values,  v_pheno,  v_pheno_size );
  	cudaMemcpyToSymbol( dev_cv_indices,  cv_indices,  indices_size );
	HANDLE_ERROR( cudaMemcpy(dev_mat_SNP, mat_SNP, mat_SNP_size, cudaMemcpyHostToDevice));
  	HANDLE_ERROR( cudaMemcpy(dev_combinations, combinations, combinations_size, cudaMemcpyHostToDevice));
  	//HANDLE_ERROR( cudaMemcpy(dev_cv_indices, cv_indices, indices_size, cudaMemcpyHostToDevice));
  	fprintf(stderr,"matrices copied  to GPU \n");
  	
  	
  	
  	//cudaHostAlloc((void**)&output.values,output.mem_size,cudaHostAllocDefault);
  	
  	// kernel call
	dim3 dimBlock(BSx);//,BSy,BSz);
	dim3 dimGrid(GSx);//,GSy,GSz);
	
	printf("\ncalling the kernel with this configuration:\n");
	printf(" interaction order: %d\n NSNPS: %d\n NIND: %d\n # cross validations: %d\n THRESHOLD: %d\n #BLOCKS: %d\n BLOCK SIZE: %d\n",ORDER, NSNPS, NIND, CV, THR, GSx, BSx);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float elapsedTime;
	cudaEventRecord(start, 0);

	MDR<<< dimGrid, dimBlock >>>(dev_mat_SNP, dev_output, dev_combinations);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime,start,stop);
	printf("kernel computation terminated, Time required (ms): %4.5f\n", elapsedTime);
	
	HANDLE_ERROR( cudaEventDestroy( start ) );
	HANDLE_ERROR( cudaEventDestroy( stop ) );


  	//free
  	cudaFree(dev_mat_SNP);
	//cudaFree(dev_v_pheno.values);
	cudaFree(dev_output);
	cudaFree(dev_combinations);
	//cudaFree(dev_cv_indices);

	free(mat_SNP);
	free(v_pheno);
	free(output);
	free(combinations);
	free(cv_indices);

  	
  	
  	
  	
  	
  	
  	
 	return 0;
}
