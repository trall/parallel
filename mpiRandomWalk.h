#pragma once 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

typedef char bool;
enum BOOL {FALSE, TRUE};

enum DIRECTIONS {LEFT, UP, RIGHT, DOWN};
enum MPI_TAGS {TAG_NUMBER, TAG_VECTOR};
enum SWAP_ACTION_ORDER {SEND_FIRST, RECEIVE_FIRST};

typedef struct InputParams {
	int argc;
	char** argv;
	int l; 
	int a; 
	int b; 
	int n; 
	int N; 

	double p[4]; 
} InputParams;

typedef struct Point {
	int id;
	int x;
	int y;
	int lifetime;
} Point;


typedef struct PointsVector {
	Point* content;
	int length;
	int size;
} PointsVector;

typedef struct Square {
	
	int length;
	int bounds[4];
	int fieldSizeX;
	int fieldSizeY;
	int neighbours[4];
	int posX;
	int posY;
} Square;

typedef struct Environment {
	int rank;
	MPI_Datatype* MyMPIpoint;
	int maxSteps;
	Square* square;
	PointsVector* points;
	PointsVector* bufferStay;
	PointsVector* buffersSend[4];
	PointsVector* buffersReceive[4];
	void* mpiBuffer;
} Environment;

Environment* processInit (int* argc, char*** argv, InputParams* input);
void pointsExchange (Environment* env);
void returnPointToTheField (Point* point, Square* square);
int getNeighbourRank (int thisRank, int offsetX, int offsetY, int squaresInRow, int squaresInColumn);
PointsVector* pointsVectorInit (int size);
void movePoints (int maxSteps, PointsVector* points, Square* square, double extra, double* probabilities);
int getPointTransitionDirection (Square* square, Point* point);
bool havePointsToMove (PointsVector* points, int maxSteps);
bool inBounds (Point* point, Square* square, int addition);
void pointsSwap (Environment* env, int direction, int actionOrder);
void fillBufferSend (Square* square, PointsVector* points, PointsVector* bufferSend, int direction);
void fillBufferStay (Square* square, PointsVector* points, PointsVector* bufferStay);
void pointsVectorResize (PointsVector* vector, int newSize);
void processFinalize (InputParams* inputParams, Environment* env);
void mergeBuffers (PointsVector* points, PointsVector* bufferStay, PointsVector* buffersReceive[]);
InputParams* getInputParams (int argc, char** argv);
void resultsOutput (int rank, PointsVector* points, InputParams* inputParams, double time);