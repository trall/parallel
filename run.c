#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include "mpiRandomWalk.h"


int main (int argc, char **argv) {
	
	InputParams* inputParams = getInputParams (argc, argv);
	if (!inputParams) {
		perror ("Invalid input data!\n");
		return 1;
	}
	struct timeval start, end;
	gettimeofday(&start, NULL);
	Environment* env = processInit (&argc, &argv, inputParams);
	int count = 0;
	while (havePointsToMove (env->points, env->maxSteps)) {
		movePoints (env->maxSteps, env->points, env->square, 0.25, inputParams->p);
		pointsExchange (env);
	}

	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	resultsOutput (env->rank, env->points, inputParams, delta);
	processFinalize (inputParams, env);
	return 0;
}


	InputParams* getInputParams (int argc, char** argv) {
	if (argc != 10) {
		printf ("You must pass 9 arguments, given %d\n", argc - 1);
		return NULL;
	}

	InputParams* params = (InputParams*) calloc (1, sizeof (InputParams));

	params->l = atoi (argv[1]);
	params->a = atoi (argv[2]);
	params->b = atoi (argv[3]);
	params->n = atoi (argv[4]);
	params->N = atoi (argv[5]);
	if (params->l <= 0 || params->a <= 0 || params->b <= 0 || params->n <= 0 || params->N <= 0){
		printf("Incorrect parameters\n");
		return NULL;
	}

	params->p[LEFT] = strtod (argv[6], NULL);
	params->p[UP] = strtod (argv[8], NULL) + params->p[LEFT];
	params->p[RIGHT] = strtod (argv[7], NULL) + params->p[UP];
	params->p[DOWN] = strtod (argv[9], NULL) + params->p[RIGHT];
	params->argc = argc;
	params->argv = argv;
	return params;
}