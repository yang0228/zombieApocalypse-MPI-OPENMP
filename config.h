/**
 * Just a config file to keep everything clean
 */
#define THREADS				64
#define STEPS               3650
#define WIDTH 	 			1024
#define HEIGHT 	 			768
#define INIT    			0.5
#define MOVE_PROB    		0.1
#define MOVE_YOUNG      	0.08
#define MOVE_ELDERLY    	0.05

// POPULATION
#define PROB_EMPTY  		0.8 // 0.2 density
#define PROB_BIRTH          5.5/1000.0/365.0
#define PROB_DEATH          5.5/1000.0/365.0
#define STARTING_POPULATION (int)floor(((VERTICAL-2)*HORIZONTAL)*POPULATION_DENSITY)
#define INITIAL_ZOMBIES     200

// DEMOGRAPHICS
#define PROB_MALE           0.5
#define PROB_FEMALE         1-MALE
#define PROB_YOUNG          0.19
#define PROB_ADULT          0.67
#define PROB_ELDERLY        0.14

// INFECTION PARAMETERS
#define PROB_INFECTION          0.6
#define PROB_INFECTION_YOUNG    PROB_INFECTION + 0.25
#define PROB_INFECTION_ELDERLY	PROB_INFECTION + 0.35

// EPIDEMIC
#define HUMAN 			0
#define INFECTED 		1
#define ZOMBIE 			2
#define DEAD 			3
#define ZOMBIE_LIFE 	28

// CONSTANTS
#define MALE 	0
#define FEMALE 	1
#define YOUNG   0
#define ADULT   1
#define ELDERLY 2

// MPI
#define NORTH 0
#define SOUTH 1