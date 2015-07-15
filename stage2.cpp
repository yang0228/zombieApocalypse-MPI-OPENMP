#include <cstdlib>
#include <mpi.h>
#include <omp.h>
#include <iostream>

#include "config.h"
#include "locks.h"
#include "Person.h"
#include "MersenneTwister.h"

using namespace std;

const int nThreads = THREADS;
MTRand *twisters[nThreads];

Person spawnHuman(bool birth = false){
	
	Person p = Person();

	double genderRand = twisters[omp_get_thread_num()]->rand();

	// gender
	if (genderRand <= PROB_MALE){
		p.setGender(MALE);
	} else {
		p.setGender(FEMALE);
	}

	// age
	if (birth){
		p.setAge(YOUNG);
	} else {
		double ageRand = twisters[omp_get_thread_num()]->rand();
		if (ageRand <= PROB_YOUNG){
			p.setAge(YOUNG);
		} else if (ageRand > 1-PROB_ELDERLY){
			p.setAge(ELDERLY);
		} else {
			p.setAge(ADULT);
		}
	}

	p.setStage(HUMAN);
	p.setEmpty(false);

	return p;
}

void initialiseZombies(Person **Mesh){
	int zombies = 0;
	while (zombies < INITIAL_ZOMBIES){
		int randX = twisters[omp_get_thread_num()]->randInt(WIDTH)+1;
		int randY = twisters[omp_get_thread_num()]->randInt(HEIGHT)+1;
		Mesh[randY][randX] = spawnHuman();
		Mesh[randY][randX].setStage(ZOMBIE);
		zombies += 1;
	}
}

void setOwnership(Person *row, int owner){
	for (int i = 0; i < WIDTH+2; i++){
		row[i].setOwner(owner);
	}
}

/**
 * Fills is the ghost cells for the model
 * @param Mesh   [description]
 * @param region [description]
 */
void setGhostCells(Person **Mesh, int region){

	// if north, open south and duplicate top
	if (region == NORTH){
		Mesh[0] = Mesh[1];
	// if south, open north and duplicate bottom
	} else {
		Mesh[HEIGHT+1] = Mesh[HEIGHT];
	}

	// copy east/west
	for (int i = 0; i < HEIGHT+2; i++){
		Mesh[i][0] = Mesh[i][1];
		Mesh[i][WIDTH+1] = Mesh[i][WIDTH];
	}

}

/**
 * Clears all of the ghost cells from the meas
 * @param Mesh [description]
 */
void clearGhostCells(Person **Mesh){

	for (int i = 0; i < HEIGHT+2; i++){
		Mesh[i][0].setEmpty(true);
		Mesh[i][WIDTH+1].setEmpty(true);
	}	

	for (int i = 0; i < WIDTH+2; i++){
		Mesh[0][i].setEmpty(true);
		Mesh[HEIGHT+1][i].setEmpty(true);
	}

}

void printStatistics(Person **mesh){

	int male = 0;
	int female = 0;
	
	int young = 0;
	int adult = 0;
	int elderly = 0;

	int human = 0;
	int infected = 0;
	int zombie = 0;

	#pragma omp parallel for default(none) collapse(2) shared(mesh) reduction(+:male,female,young,adult,elderly,human,infected,zombie) num_threads(nThreads)
	for (int i = 1; i < HEIGHT + 1; i++){
		for (int j = 1; j < WIDTH + 1; j++){

			if (mesh[i][j].isEmpty()){continue;}

			if (mesh[i][j].getStage() == HUMAN){
				if (mesh[i][j].getGender() == MALE){
					male += 1;
				} else {
					female += 1;
				}
			}

			if (mesh[i][j].getAge() == YOUNG){
				young += 1;
			} else if (mesh[i][j].getAge() == ADULT){
				adult += 1;
			} else {
				elderly += 1;
			}

			if (mesh[i][j].getStage() == HUMAN){
				human += 1;
			} else if (mesh[i][j].getStage() == INFECTED){
				infected += 1;
			} else {
				zombie += 1;
			}

		}
	}

	printf("%d,%d,%d,%d,%d,%d,%d,%d\n",male,female,young,adult,elderly,human,infected,zombie);
}

int main(int argc, char** argv){

	printf("male,female,young,adult,elderly,human,infected,zombie\n");

	// create the twister array
	for (int i = 0; i < nThreads; i++){
		twisters[i] = new MTRand(i);
	}	
	
	// initialise the locks
	bool *locks = new bool[HEIGHT+2];
	for (int i = 0; i < HEIGHT+2; i++){
		locks[i] = false;
	}

	// create the meshes
	Person **MeshA = new Person *[HEIGHT+2];
	Person **MeshB = new Person *[HEIGHT+2];
	for (int i = 0; i < HEIGHT+2; i++){
		MeshA[i] = new Person[WIDTH+2];
		MeshB[i] = new Person[WIDTH+2];			
		for (int j = 0; j < WIDTH+2; j++){
			double populate = twisters[0]->rand();
			if (populate > PROB_EMPTY){
				MeshA[i][j] = spawnHuman();
			} else {
				MeshA[i][j].setEmpty(true);
			}
			MeshB[i][j].setEmpty(true);
		}
	}

	// initialise the zombies
	initialiseZombies(MeshA);

	// start the clock
	double startTime = omp_get_wtime();

	int region;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&region);
	MPI_Status status;

	// create the MPI datatypes
	MPI_Datatype personDataType;
	MPI_Type_contiguous(sizeof(Person),MPI_BYTE,&personDataType);
	MPI_Type_commit(&personDataType);
	MPI_Datatype personVectorDataType;
	MPI_Type_vector(WIDTH+2,1,1,personDataType,&personVectorDataType);
	MPI_Type_commit(&personVectorDataType);

	// start the simulation
	for (int s = 0; s < STEPS; s++){

		// MOVEMENT

		// NORTH
		if (region == NORTH){

			Person rowSwap[WIDTH+2];
			MPI_Send(&MeshA[HEIGHT],1,personVectorDataType,SOUTH,0,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			setOwnership(rowSwap,SOUTH);
			MeshA[HEIGHT+1] = rowSwap;
			setGhostCells(MeshA,NORTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if(MeshA[i][j].isEmpty()){
						continue;
					}

					MeshA[i][j].setOwner(NORTH);
					MeshA[i][j].setEmpty(true);

					double move = twisters[omp_get_thread_num()]->rand();
					double moveProb = MeshA[i][j].getMoveAdjustment(MOVE_PROB);

					// up
					if (move < 1.0*moveProb && MeshA[i-1][j].isEmpty() && MeshB[i-1][j].isEmpty()){
						MeshB[i-1][j] = MeshA[i][j];
					}

					// down
					else if (move < 2.0*moveProb && MeshA[i+1][j].isEmpty() && MeshB[i+1][j].isEmpty()){
						MeshB[i+1][j] = MeshA[i][j];
					}

					// left
					else if (move < 3.0*moveProb && MeshA[i][j-1].isEmpty() && MeshB[i][j-1].isEmpty()){
						MeshB[i][j-1] = MeshA[i][j];
					}

					// right
					else if (move < 4.0*moveProb && MeshA[i][j+1].isEmpty() && MeshB[i][j+1].isEmpty()){
						MeshB[i][j+1] = MeshA[i][j];
					}
					
					// stay still
					else {
						MeshB[i][j] = MeshA[i][j];
						
					}

					MeshB[i][j].setEmpty(false);

				}
				unlock(i,locks);
			}

			// transfer ownership
			MPI_Send(&MeshB[HEIGHT],1,personVectorDataType,SOUTH,0+10,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,1+10,MPI_COMM_WORLD,&status);
			MeshB[HEIGHT] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);
		}

		// SOUTH
		else {

			Person rowSwap[WIDTH+2];
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshA[1],1,personVectorDataType,NORTH,1,MPI_COMM_WORLD);
			setOwnership(rowSwap,NORTH);
			MeshA[0] = rowSwap;
			setGhostCells(MeshA,SOUTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if(MeshA[i][j].isEmpty()){
						continue;
					}

					// set ownership to ourselves
					MeshA[i][j].setOwner(SOUTH);
					MeshA[i][j].setEmpty(true);

					double move = twisters[omp_get_thread_num()]->rand();
					double moveProb = MeshA[i][j].getMoveAdjustment(MOVE_PROB);

					// up
					if (move < 1.0*moveProb && MeshA[i-1][j].isEmpty() && MeshB[i-1][j].isEmpty()){
						MeshB[i-1][j] = MeshA[i][j];
					}

					// down
					else if (move < 2.0*moveProb && MeshA[i+1][j].isEmpty() && MeshB[i+1][j].isEmpty()){
						MeshB[i+1][j] = MeshA[i][j];
					}

					// left
					else if (move < 3.0*moveProb && MeshA[i][j-1].isEmpty() && MeshB[i][j-1].isEmpty()){
						MeshB[i][j-1] = MeshA[i][j];
					}

					// right
					else if (move < 4.0*moveProb && MeshA[i][j+1].isEmpty() && MeshB[i][j+1].isEmpty()){
						MeshB[i][j+1] = MeshA[i][j];
					}
					
					// stay still
					else {
						MeshB[i][j] = MeshA[i][j];
						
					}

					MeshB[i][j].setEmpty(false);

				}

				unlock(i,locks);
			}

			// transfer ownership
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshB[0],1,personVectorDataType,NORTH,1+10,MPI_COMM_WORLD);
			MeshB[1] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);
		}

		// wait up
		MPI_Barrier(MPI_COMM_WORLD);

		// ZOMBIFICATION
		if (region == NORTH){

			Person rowSwap[WIDTH+2];
			MPI_Send(&MeshA[HEIGHT],1,personVectorDataType,SOUTH,0,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			setOwnership(rowSwap,SOUTH);
			MeshA[HEIGHT+1] = rowSwap;
			setGhostCells(MeshA,NORTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if (MeshA[i][j].isEmpty()){
						continue;
					}

					MeshA[i][j].setEmpty(true);

					// turn the infected
					if (MeshA[i][j].getStage() == INFECTED){
						MeshB[i][j].setStage(ZOMBIE);
						MeshB[i][j].setEmpty(false);
						continue;
					}

					if (MeshA[i][j].getStage() == ZOMBIE){
						MeshA[i][j].decreaseZombieLife();
						MeshB[i][j] = MeshA[i][j];
						MeshB[i][j].setEmpty(false);
						continue;
					}

					if (MeshA[i-1][j].isZombie() || MeshB[i-1][j].isZombie() ||
						MeshA[i+1][j].isZombie() || MeshB[i+1][j].isZombie() ||
						MeshA[i][j-1].isZombie() || MeshB[i][j-1].isZombie() ||
						MeshA[i][j+1].isZombie() || MeshB[i][j+1].isZombie()){

						double infectRand = twisters[omp_get_thread_num()]->rand();
						double infectProb = MeshA[i][j].getInfectionProbability();

						if (infectRand <= infectProb){
							MeshA[i][j].setStage(INFECTED);
						}

					}

					MeshB[i][j] = MeshA[i][j];
					MeshB[i][j].setEmpty(false);
				}

				unlock(i,locks);

			}

			// transfer ownership
			MPI_Send(&MeshB[HEIGHT],1,personVectorDataType,SOUTH,0+10,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,1+10,MPI_COMM_WORLD,&status);
			MeshB[HEIGHT] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);

		} else {

			Person rowSwap[WIDTH+2];
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshA[1],1,personVectorDataType,NORTH,1,MPI_COMM_WORLD);
			setOwnership(rowSwap,NORTH);
			MeshA[0] = rowSwap;
			setGhostCells(MeshA,SOUTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if (MeshA[i][j].isEmpty()){
						continue;
					}

					MeshA[i][j].setEmpty(true);

					// turn the infected
					if (MeshA[i][j].getStage() == INFECTED){
						MeshB[i][j].setStage(ZOMBIE);
						MeshB[i][j].setEmpty(false);
						continue;
					}

					if (MeshA[i][j].getStage() == ZOMBIE){
						MeshA[i][j].decreaseZombieLife();
						MeshB[i][j] = MeshA[i][j];
						MeshB[i][j].setEmpty(false);
						continue;
					}

					if (MeshA[i-1][j].isZombie() || MeshB[i-1][j].isZombie() ||
						MeshA[i+1][j].isZombie() || MeshB[i+1][j].isZombie() ||
						MeshA[i][j-1].isZombie() || MeshB[i][j-1].isZombie() ||
						MeshA[i][j+1].isZombie() || MeshB[i][j+1].isZombie()){

						double infectRand = twisters[omp_get_thread_num()]->rand();
						double infectProb = MeshA[i][j].getInfectionProbability();

						if (infectRand <= infectProb){
							MeshA[i][j].setStage(INFECTED);
						}

					}

					MeshB[i][j] = MeshA[i][j];
					MeshB[i][j].setEmpty(false);
				}

				unlock(i,locks);

			}

			// transfer ownership
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshB[0],1,personVectorDataType,NORTH,1+10,MPI_COMM_WORLD);
			MeshB[1] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);

		}

		MPI_Barrier(MPI_COMM_WORLD);

		// BIRTHS
		if (region == NORTH){

			Person rowSwap[WIDTH+2];
			MPI_Send(&MeshA[HEIGHT],1,personVectorDataType,SOUTH,0,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			setOwnership(rowSwap,SOUTH);
			MeshA[HEIGHT+1] = rowSwap;
			setGhostCells(MeshA,NORTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if (MeshA[i][j].isEmpty()){
						continue;
					}

					MeshA[i][j].setEmpty(true);

					if (MeshA[i][j].canBreed() == false){
						MeshB[i][j] = MeshA[i][j];
						MeshB[i][j].setEmpty(false);
						continue;
					}

					// look in the local network
					if (MeshA[i-1][j].canBreed() || MeshB[i-1][j].canBreed() ||
						MeshA[i+1][j].canBreed() || MeshB[i+1][j].canBreed() ||
						MeshA[i][j-1].canBreed() || MeshB[i][j-1].canBreed() ||
						MeshA[i][j+1].canBreed() || MeshB[i][j+1].canBreed()){

						double breedRand = twisters[omp_get_thread_num()]->rand();
						if (breedRand <= PROB_BIRTH){

							// look for a space
							if (MeshA[i-1][j].isEmpty() && MeshB[i-1][j].isEmpty()){
								MeshA[i-1][j] = spawnHuman(true);
							} else if (MeshA[i+1][j].isEmpty() && MeshB[i+1][j].isEmpty()){
								MeshA[i+1][j] = spawnHuman(true);
							} else if (MeshA[i][j-1].isEmpty() && MeshB[i][j-1].isEmpty()){
								MeshA[i][j-1] = spawnHuman(true);
							}  else if (MeshA[i][j+1].isEmpty() && MeshB[i][j+1].isEmpty()){
								MeshA[i][j+1] = spawnHuman(true);
							}

							// nowehere to go, no babies
						}
					}

					MeshB[i][j] = MeshA[i][j];
					MeshB[i][j].setEmpty(false);
				}

				unlock(i,locks);

			}

			// transfer ownership
			MPI_Send(&MeshB[HEIGHT],1,personVectorDataType,SOUTH,0+10,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,1+10,MPI_COMM_WORLD,&status);
			MeshB[HEIGHT] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);

		} else {

			Person rowSwap[WIDTH+2];
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshA[1],1,personVectorDataType,NORTH,1,MPI_COMM_WORLD);
			setOwnership(rowSwap,NORTH);
			MeshA[0] = rowSwap;
			setGhostCells(MeshA,SOUTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if (MeshA[i][j].isEmpty()){
						continue;
					}

					MeshA[i][j].setEmpty(true);

					if (MeshA[i][j].canBreed() == false){
						MeshB[i][j] = MeshA[i][j];
						MeshB[i][j].setEmpty(false);
						continue;
					}

					// look in the local network
					if (MeshA[i-1][j].canBreed() || MeshB[i-1][j].canBreed() ||
						MeshA[i+1][j].canBreed() || MeshB[i+1][j].canBreed() ||
						MeshA[i][j-1].canBreed() || MeshB[i][j-1].canBreed() ||
						MeshA[i][j+1].canBreed() || MeshB[i][j+1].canBreed()){

						double breedRand = twisters[omp_get_thread_num()]->rand();
						if (breedRand <= PROB_BIRTH){

							// look for a space
							if (MeshA[i-1][j].isEmpty() && MeshB[i-1][j].isEmpty()){
								MeshA[i-1][j] = spawnHuman(true);
							} else if (MeshA[i+1][j].isEmpty() && MeshB[i+1][j].isEmpty()){
								MeshA[i+1][j] = spawnHuman(true);
							} else if (MeshA[i][j-1].isEmpty() && MeshB[i][j-1].isEmpty()){
								MeshA[i][j-1] = spawnHuman(true);
							}  else if (MeshA[i][j+1].isEmpty() && MeshB[i][j+1].isEmpty()){
								MeshA[i][j+1] = spawnHuman(true);
							}

							// nowehere to go, no babies
						}
					}

					MeshB[i][j] = MeshA[i][j];
					MeshB[i][j].setEmpty(false);
				}

				unlock(i,locks);

			}

			// transfer ownership
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshB[0],1,personVectorDataType,NORTH,1+10,MPI_COMM_WORLD);
			MeshB[1] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);

		}

		MPI_Barrier(MPI_COMM_WORLD);

		// DEATHS
		if (region == NORTH){

			Person rowSwap[WIDTH+2];
			MPI_Send(&MeshA[HEIGHT],1,personVectorDataType,SOUTH,0,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			setOwnership(rowSwap,SOUTH);
			MeshA[HEIGHT+1] = rowSwap;
			setGhostCells(MeshA,NORTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if (MeshA[i][j].isEmpty()){
						continue;
					}

					MeshA[i][j].setEmpty(true);

					// check if they can, and will be killed
					if (MeshA[i][j].isHuman()){
						double deathRand = twisters[omp_get_thread_num()]->rand();
						if (deathRand <= PROB_DEATH){
							MeshB[i][j] = MeshA[i][j];
							MeshB[i][j].setStage(DEAD);
							continue;
						}
					}

					MeshB[i][j] = MeshA[i][j];
					MeshB[i][j].setEmpty(false);
				}

				unlock(i,locks);

			}

			// transfer ownership
			MPI_Send(&MeshB[HEIGHT],1,personVectorDataType,SOUTH,0+10,MPI_COMM_WORLD);
			MPI_Recv(&rowSwap,1,personVectorDataType,SOUTH,1+10,MPI_COMM_WORLD,&status);
			MeshB[HEIGHT] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);

		} else {

			Person rowSwap[WIDTH+2];
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshA[1],1,personVectorDataType,NORTH,1,MPI_COMM_WORLD);
			setOwnership(rowSwap,NORTH);
			MeshA[0] = rowSwap;
			setGhostCells(MeshA,SOUTH);

			#pragma omp parallel for default(none) shared(s,MeshA,MeshB,locks,twisters) num_threads(nThreads)
			for (int i = 1; i < HEIGHT; i++){

				lock(i,locks);

				for (int j = 1; j < WIDTH; j++){

					if (MeshA[i][j].isEmpty()){
						continue;
					}

					MeshA[i][j].setEmpty(true);

					// check if they can, and will be killed
					if (MeshA[i][j].isHuman()){
						double deathRand = twisters[omp_get_thread_num()]->rand();
						if (deathRand <= PROB_DEATH){
							MeshB[i][j] = MeshA[i][j];
							MeshB[i][j].setStage(DEAD);
							continue;
						}
					}

					MeshB[i][j] = MeshA[i][j];
					MeshB[i][j].setEmpty(false);
				}

				unlock(i,locks);

			}

			// transfer ownership
			MPI_Recv(&rowSwap,1,personVectorDataType,NORTH,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Send(&MeshB[0],1,personVectorDataType,NORTH,1+10,MPI_COMM_WORLD);
			MeshB[1] = rowSwap;
			swap(MeshA,MeshB);
			clearGhostCells(MeshA);

		}

		// wait up everyone
		MPI_Barrier(MPI_COMM_WORLD);

		// gather statistics
		if (region == NORTH){
			printStatistics(MeshA);
		}
	}

	// wait up everyone
	MPI_Barrier(MPI_COMM_WORLD);

	// destroy MPI
	MPI_Finalize();

	// stop the clock
	double endTime = omp_get_wtime();
	double totalTime  = endTime - startTime;

	// done
	printf("Total runtime @ %d cores = %f\n", THREADS, totalTime);
}