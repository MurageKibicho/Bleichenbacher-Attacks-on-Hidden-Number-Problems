#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpzi.h>
#include <flint/fmpq.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod.h>
#define MPFR_PRECISION 1024
//clear && gcc UnderstandingBias.c -o m.o -lm -lgmp -lmpfr -lflint && ./m.o 

/**
 * Compute |B_n(V)| - the magnitude of the sampled bias.
 * Using real arithmetic: |sum(cosθ + i·sinθ)| = sqrt((sum cosθ)² + (sum sinθ)²)
 */
void ComputeBiasMagnitude(int length, fmpz_t *values,fmpz_t modulo)
{
	mpfr_t tempFloat,sum_cos, sum_sin, angle, n_mpfr, two_pi;
	mpfr_init2(tempFloat, MPFR_PRECISION);mpfr_init2(sum_cos, MPFR_PRECISION);mpfr_init2(sum_sin, MPFR_PRECISION);mpfr_init2(angle, MPFR_PRECISION);mpfr_init2(n_mpfr, MPFR_PRECISION);mpfr_init2(two_pi, MPFR_PRECISION);
	int stringBase = 2;mpfr_set_zero(sum_cos, 1);mpfr_set_zero(sum_sin, 1);
	//2π
	mpfr_const_pi(two_pi, MPFR_RNDN);mpfr_mul_ui(two_pi, two_pi, 2, MPFR_RNDN);
	//Set modulo to mpfr
	char *tempString_modulo = fmpz_get_str(NULL, stringBase, modulo);mpfr_set_str(n_mpfr, tempString_modulo, stringBase, MPFR_PRECISION);free(tempString_modulo);
	for(int i = 0; i < length; i++)
	{
		//Cast fmpz_t to mpfr
		char *tempString = fmpz_get_str(NULL, stringBase, values[i]);mpfr_set_str(tempFloat, tempString, stringBase, MPFR_PRECISION);free(tempString);
		// angle = 2π * (values[i]) / n
		mpfr_mul(angle, two_pi, tempFloat, MPFR_RNDN);mpfr_div(angle, angle, n_mpfr, MPFR_RNDN);
		// sum_cos += cos(angle)
		mpfr_cos(tempFloat, angle, MPFR_RNDN);mpfr_add(sum_cos, sum_cos, tempFloat, MPFR_RNDN);
		// sum_sin += sin(angle)
		mpfr_sin(tempFloat, angle, MPFR_RNDN);mpfr_add(sum_sin, sum_sin, tempFloat, MPFR_RNDN);
	}
	//Divide by length
	mpfr_div_ui(sum_cos, sum_cos, length, MPFR_RNDN);
	mpfr_div_ui(sum_sin, sum_sin, length, MPFR_RNDN);
	// result = sqrt(sum_cos² + sum_sin²)
	mpfr_sqr(tempFloat, sum_cos, MPFR_RNDN);
	mpfr_sqr(angle, sum_sin, MPFR_RNDN); //reusing angle as temp2
	mpfr_add(tempFloat, tempFloat, angle, MPFR_RNDN);
	mpfr_sqrt(tempFloat, tempFloat, MPFR_RNDN);
	printf("Bias: ");mpfr_printf("%.10Rf\n", tempFloat);
	mpfr_clear(sum_cos);mpfr_clear(sum_sin);mpfr_clear(angle);mpfr_clear(n_mpfr);mpfr_clear(two_pi);mpfr_clear(tempFloat);
}

void ComputeBiasMagnitude_Frequency(int length, fmpz_t *values,fmpz_t modulo,fmpz_t frequency)
{
	mpfr_t tempFloat,sum_cos, sum_sin, angle, n_mpfr, two_pi,frequency_mpfr;
	mpfr_init2(frequency_mpfr, MPFR_PRECISION);mpfr_init2(tempFloat, MPFR_PRECISION);mpfr_init2(sum_cos, MPFR_PRECISION);mpfr_init2(sum_sin, MPFR_PRECISION);mpfr_init2(angle, MPFR_PRECISION);mpfr_init2(n_mpfr, MPFR_PRECISION);mpfr_init2(two_pi, MPFR_PRECISION);
	int stringBase = 2;mpfr_set_zero(sum_cos, 1);mpfr_set_zero(sum_sin, 1);
	//2π
	mpfr_const_pi(two_pi, MPFR_RNDN);mpfr_mul_ui(two_pi, two_pi, 2, MPFR_RNDN);
	//Set modulo to mpfr
	char *tempString_modulo = fmpz_get_str(NULL, stringBase, modulo);mpfr_set_str(n_mpfr, tempString_modulo, stringBase, MPFR_PRECISION);free(tempString_modulo);
	//Set frequency to mpfr
	char *tempString_frequency = fmpz_get_str(NULL, stringBase, frequency);mpfr_set_str(frequency_mpfr, tempString_frequency, stringBase, MPFR_PRECISION);free(tempString_frequency);
	
	for(int i = 0; i < length; i++)
	{
		//Cast fmpz_t to mpfr
		char *tempString = fmpz_get_str(NULL, stringBase, values[i]);mpfr_set_str(tempFloat, tempString, stringBase, MPFR_PRECISION);free(tempString);
		// angle = 2π * frequency* (values[i]) / n
		mpfr_mul(angle, two_pi, tempFloat, MPFR_RNDN);mpfr_mul(angle, angle, frequency_mpfr, MPFR_RNDN);mpfr_div(angle, angle, n_mpfr, MPFR_RNDN);
		// sum_cos += cos(angle)
		mpfr_cos(tempFloat, angle, MPFR_RNDN);mpfr_add(sum_cos, sum_cos, tempFloat, MPFR_RNDN);
		// sum_sin += sin(angle)
		mpfr_sin(tempFloat, angle, MPFR_RNDN);mpfr_add(sum_sin, sum_sin, tempFloat, MPFR_RNDN);
	}
	//Divide by length
	mpfr_div_ui(sum_cos, sum_cos, length, MPFR_RNDN);
	mpfr_div_ui(sum_sin, sum_sin, length, MPFR_RNDN);
	// result = sqrt(sum_cos² + sum_sin²)
	mpfr_sqr(tempFloat, sum_cos, MPFR_RNDN);
	mpfr_sqr(angle, sum_sin, MPFR_RNDN); //reusing angle as temp2
	mpfr_add(tempFloat, tempFloat, angle, MPFR_RNDN);
	mpfr_sqrt(tempFloat, tempFloat, MPFR_RNDN);
	printf("Bias: ");mpfr_printf("%.10Rf\n", tempFloat);
	mpfr_clear(frequency_mpfr);mpfr_clear(sum_cos);mpfr_clear(sum_sin);mpfr_clear(angle);mpfr_clear(n_mpfr);mpfr_clear(two_pi);mpfr_clear(tempFloat);
}
void SampleBiasListTestMSB(fmpz_t modulo,flint_rand_t rngState)
{
	int numberOfSamples = 10000;
	int totalFixedMSB = 8;
	int moduloBitLength = fmpz_sizeinbase(modulo,2);
	//Allocate memory
	fmpz_t maxVal;fmpz_init(maxVal);
	fmpz_t *biasedSamples = malloc(numberOfSamples * sizeof(fmpz_t));for(int currentSample = 0; currentSample < numberOfSamples; currentSample++){fmpz_init(biasedSamples[currentSample]);}
	printf("Prime:");fmpz_print(modulo);printf(" (%d) bits\n",moduloBitLength);
	for(int numberOfFixedMSB = 0; numberOfFixedMSB <= totalFixedMSB; numberOfFixedMSB++)
	{
		fmpz_fdiv_q_2exp(maxVal, modulo, numberOfFixedMSB);
		for(int currentSample = 0; currentSample < numberOfSamples; currentSample++)
		{
			//Generate biased nonce with MSBs equal to 0
			fmpz_randm(biasedSamples[currentSample], rngState, maxVal);
			//printf("(%4d):",currentSample);fmpz_print(biasedSamples[currentSample]);printf("\n");
		}
		printf("(%3d): ",numberOfFixedMSB);ComputeBiasMagnitude(numberOfSamples, biasedSamples, modulo);	
	}
	//Free memory
	for(int currentSample = 0; currentSample < numberOfSamples; currentSample++){fmpz_clear(biasedSamples[currentSample]);}free(biasedSamples);
	fmpz_clear(maxVal);
}

void SampleBiasListTestLSB(fmpz_t modulo,flint_rand_t rngState)
{
	int numberOfSamples = 10000;
	int totalFixedMSB = 8;
	int moduloBitLength = fmpz_sizeinbase(modulo,2);
	//Allocate memory
	fmpz_t maxVal,optimalFrequency;fmpz_init(maxVal);fmpz_init(optimalFrequency);
	fmpz_t *biasedSamples = malloc(numberOfSamples * sizeof(fmpz_t));for(int currentSample = 0; currentSample < numberOfSamples; currentSample++){fmpz_init(biasedSamples[currentSample]);}
	printf("Prime:");fmpz_print(modulo);printf(" (%d) bits\n",moduloBitLength);
	for(int numberOfFixedMSB = 0; numberOfFixedMSB <= totalFixedMSB; numberOfFixedMSB++)
	{
		fmpz_fdiv_q_2exp(maxVal, modulo, numberOfFixedMSB);
		for(int currentSample = 0; currentSample < numberOfSamples; currentSample++)
		{
			//Generate biased nonce with MSBs equal to 0
			fmpz_randm(biasedSamples[currentSample], rngState, maxVal);
			//Shift To make lsb fixed
			fmpz_mul_2exp(biasedSamples[currentSample], biasedSamples[currentSample], numberOfFixedMSB);
			//printf("(%4d):",currentSample);fmpz_print(biasedSamples[currentSample]);printf("\n");
		}
		printf("(%3d): ",numberOfFixedMSB);ComputeBiasMagnitude_Frequency(numberOfSamples, biasedSamples, modulo,maxVal);		
	}
	//Free memory
	for(int currentSample = 0; currentSample < numberOfSamples; currentSample++){fmpz_clear(biasedSamples[currentSample]);}free(biasedSamples);
	fmpz_clear(maxVal);fmpz_clear(optimalFrequency);
}

void GenerateTable1()
{
	fmpz_t modulo;fmpz_init(modulo);
	flint_rand_t rngState;flint_rand_init(rngState);
	
	int targetBitCount = 13;
	fmpz_randprime(modulo, rngState, targetBitCount, 1);
	
	SampleBiasListTestLSB(modulo, rngState);
	
	flint_rand_clear(rngState);
	fmpz_clear(modulo);
}

void GenerateDatasetLSB(int equationsToGenerate, int numberOfFixedMSB, flint_rand_t rngState, fmpz_t secretKey, fmpz_t modulo,fmpz_t *multipliers,fmpz_t *observations)
{
	fmpz_t maxVal,temp0,temp1;fmpz_init(maxVal);fmpz_init(temp0);fmpz_init(temp1);
	for(int i = 0; i < equationsToGenerate; i++)
	{
		fmpz_fdiv_q_2exp(maxVal, modulo, numberOfFixedMSB);
		//Generate biased nonce with MSBs equal to 0
		fmpz_randm(temp0, rngState, maxVal);
		//Shift To make lsb fixed
		fmpz_mul_2exp(temp0, temp0, numberOfFixedMSB);
			
		//Generate random multiplier
		fmpz_randm(multipliers[i], rngState, maxVal);
		//Multiply secret and multiplier
		fmpz_mul(temp1, multipliers[i], secretKey);
		//Set Observation
		fmpz_sub(observations[i], temp0, temp1);
		fmpz_mod(observations[i],observations[i],modulo);
				
	}
	fmpz_clear(maxVal);fmpz_clear(temp0);fmpz_clear(temp1);
}

void GenerateDataset(int equationsToGenerate, int numberOfFixedMSB, flint_rand_t rngState, fmpz_t secretKey, fmpz_t modulo,fmpz_t *multipliers,fmpz_t *observations)
{
	fmpz_t maxVal,temp0,temp1;fmpz_init(maxVal);fmpz_init(temp0);fmpz_init(temp1);
	for(int i = 0; i < equationsToGenerate; i++)
	{
		fmpz_fdiv_q_2exp(maxVal, modulo, numberOfFixedMSB);
		//Generate biased nonce with MSBs equal to 0
		fmpz_randm(temp0, rngState, maxVal);
		//Generate random multiplier
		fmpz_randm(multipliers[i], rngState, maxVal);
		//Multiply secret and multiplier
		fmpz_mul(temp1, multipliers[i], secretKey);
		//Set Observation
		fmpz_sub(observations[i], temp0, temp1);
		fmpz_mod(observations[i],observations[i],modulo);
				
	}
	fmpz_clear(maxVal);fmpz_clear(temp0);fmpz_clear(temp1);
}

void PrintKnownDataset(int equationsToGenerate, fmpz_t secretKey, fmpz_t modulo,fmpz_t *multipliers,fmpz_t *observations)
{
	for(int i = 0; i < equationsToGenerate; i++)
	{
		printf("%3d: k = (",i);fmpz_print(observations[i]);printf("+");fmpz_print(multipliers[i]);printf("*");fmpz_print(secretKey);printf(") mod(");fmpz_print(modulo);printf(")\n");
	}
}

void BruteforceBiasMagnitude(fmpz_t modulo,flint_rand_t rngState)
{
	int equationsToGenerate = 30;
	int numberOfFixedMSB = 1;
	fmpz_t secretKey;fmpz_init(secretKey);fmpz_t *multipliers = malloc(equationsToGenerate * sizeof(fmpz_t));fmpz_t *observations = malloc(equationsToGenerate * sizeof(fmpz_t));fmpz_t *kGuess = malloc(equationsToGenerate * sizeof(fmpz_t));for(int i = 0; i < equationsToGenerate; i++){fmpz_init(multipliers[i]);fmpz_init(observations[i]);fmpz_init(kGuess[i]);}
	//Generate secret and dataset
	fmpz_randm(secretKey, rngState, modulo);
	GenerateDataset(equationsToGenerate, numberOfFixedMSB, rngState, secretKey, modulo, multipliers,observations);
	//PrintKnownDataset(equationsToGenerate, secretKey, modulo, multipliers, observations);
	int moduloLength = fmpz_sizeinbase(modulo,2);
	assert(moduloLength < 32);
	//Get modulo as integer
	int moduloInt = fmpz_get_ui(modulo);
	printf("Target:");fmpz_print(secretKey);printf("\n");
	for(int secretKeyGuess = 0; secretKeyGuess < moduloInt - 1; secretKeyGuess++)
	{
		//Generate k's for this secret key
		for(int i = 0; i < equationsToGenerate; i++)
		{
			fmpz_mul_ui(kGuess[i], multipliers[i],secretKeyGuess);
			fmpz_add(kGuess[i],kGuess[i],observations[i]);fmpz_mod(kGuess[i],kGuess[i],modulo);
		}
		//Calculate bias magnitude for the generated k's
		printf("(%3d): ",secretKeyGuess);ComputeBiasMagnitude(equationsToGenerate, kGuess, modulo);	
	}
	
	for(int i = 0; i < equationsToGenerate; i++){fmpz_clear(multipliers[i]);fmpz_clear(observations[i]);fmpz_clear(kGuess[i]);}
	free(multipliers);free(observations);free(kGuess);		
	fmpz_clear(secretKey);
}
void BruteforceBiasMagnitudeLSB(fmpz_t modulo,flint_rand_t rngState)
{
	int equationsToGenerate = 30;
	int numberOfFixedMSB = 1;
	fmpz_t secretKey,optimalFrequency;fmpz_init(secretKey);;fmpz_init(optimalFrequency);fmpz_t *multipliers = malloc(equationsToGenerate * sizeof(fmpz_t));fmpz_t *observations = malloc(equationsToGenerate * sizeof(fmpz_t));fmpz_t *kGuess = malloc(equationsToGenerate * sizeof(fmpz_t));for(int i = 0; i < equationsToGenerate; i++){fmpz_init(multipliers[i]);fmpz_init(observations[i]);fmpz_init(kGuess[i]);}
	//Generate secret and dataset
	fmpz_fdiv_q_2exp(optimalFrequency, modulo, numberOfFixedMSB);
	fmpz_randm(secretKey, rngState, modulo);
	GenerateDatasetLSB(equationsToGenerate, numberOfFixedMSB, rngState, secretKey, modulo, multipliers,observations);
	//PrintKnownDataset(equationsToGenerate, secretKey, modulo, multipliers, observations);
	int moduloLength = fmpz_sizeinbase(modulo,2);
	assert(moduloLength < 32);
	//Get modulo as integer
	int moduloInt = fmpz_get_ui(modulo);
	printf("Target:");fmpz_print(secretKey);printf("\n");
	for(int secretKeyGuess = 0; secretKeyGuess < moduloInt - 1; secretKeyGuess++)
	{
		//Generate k's for this secret key
		for(int i = 0; i < equationsToGenerate; i++)
		{
			fmpz_mul_ui(kGuess[i], multipliers[i],secretKeyGuess);
			fmpz_add(kGuess[i],kGuess[i],observations[i]);fmpz_mod(kGuess[i],kGuess[i],modulo);
		}
		//Calculate bias magnitude for the generated k's
		printf("|(%3d): ",secretKeyGuess);ComputeBiasMagnitude_Frequency(equationsToGenerate, kGuess, modulo,optimalFrequency);	
	}
	
	for(int i = 0; i < equationsToGenerate; i++){fmpz_clear(multipliers[i]);fmpz_clear(observations[i]);fmpz_clear(kGuess[i]);}
	free(multipliers);free(observations);free(kGuess);		
	fmpz_clear(secretKey);fmpz_clear(optimalFrequency);
}
void BruteforceTest()
{
	fmpz_t modulo;fmpz_init(modulo);
	flint_rand_t rngState;flint_rand_init(rngState);
	
	int targetBitCount = 13;
	fmpz_randprime(modulo, rngState, targetBitCount, 1);
	//fmpz_mul_ui(modulo, modulo, 2);
	BruteforceBiasMagnitude(modulo, rngState);
	//BruteforceBiasMagnitudeLSB(modulo, rngState);
	
	flint_rand_clear(rngState);
	fmpz_clear(modulo);
}

int main()
{
	//GenerateTable1();
	BruteforceTest();
	flint_cleanup();
	return 0;
}
