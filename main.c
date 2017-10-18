/*
 * Solving the traveling salesman problem using a genetic algorithm.
 * Author: Theo Morales <theo.morales.fr@gmail.com>
 * 2017
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <zconf.h>


#define POPULATION_SIZE 100
#define GENERATION_SIZE 10
#define NB_GENES 9
#define GENE_SIZE 3
#define MUTATION_RATE 0.015


typedef struct {
    int x;
    int y;
    int binaryIndex[3];
    char *name;
} city;


const city cities[8] = {
    { .x = 4, .y = 2, .binaryIndex = { 0, 0, 0 }, .name = "Lille" },
    { .x = 4, .y = 6, .binaryIndex = { 0, 0, 1 }, .name = "Paris" },
    { .x = 7, .y = 7, .binaryIndex = { 0, 1, 0 }, .name = "Reims" },
    { .x = 6, .y = 14, .binaryIndex = { 0, 1, 1 }, .name =  "Lyon" },
    { .x = 10, .y = 22, .binaryIndex = { 1, 0, 0 }, .name = "Marseille" },
    { .x = 2, .y = 9, .binaryIndex = { 1, 0, 1 }, .name = "Nantes" },
    { .x = 1, .y = 11, .binaryIndex = { 1, 1, 0 }, .name = "La Rochelle" },
    { .x = 3, .y = 12, .binaryIndex = { 1, 1, 1 }, .name = "Bordeaux" }
};

/*
 * Describes the order to visit cities:
 * A gene is a 3-bit array, a chromosome contains 9 genes for the 8 cities + the first one
 * */
typedef struct {
    int genes[NB_GENES][GENE_SIZE];
    float fitness;
} chromosome;


float maxFitness = 0;

int bin_to_dec(int bin[GENE_SIZE]) {

    int dec = 0;

    for (int i = 0; i < GENE_SIZE; i ++)
        dec += bin[i] * pow(2, (GENE_SIZE - 1 - i));

    return dec;
}


/* Generate bits to fill a gene */
void make_gene(int *gene) {

    for (int i = 0; i < GENE_SIZE; i++)
        gene[i] = rand() % 2;
}


/*
 * Map a gene from its binary index to the corresponding city
 * and return its coordinates
 */
void get_gene_coord(int gene[], int coord[]) {

    int geneValue = bin_to_dec(gene);

    coord[0] = cities[geneValue].x;
    coord[1] = cities[geneValue].y;
}


int validate_chromosome(chromosome c) {

    int isThereTwice = 0;

    /* Check if a gene is present twice */
    for (int i = 0; i < NB_GENES; i++) {
        for (int j = 1; j < NB_GENES - 1; j++) { //Ignore the last one, which SHOULD be equal to the first one !
            if (i != j && c.genes[i] == c.genes[j]) {
                isThereTwice = 1;
                break;
            }
        }
    }

    return !isThereTwice;
}

/* Objective function: total length of the trip */
float score(chromosome c) {

    float score = 0;

    for (int i =  1; i < NB_GENES; i++) {
        int geneCoord[2] = {0}, prevGeneCoord[2] = {0};

        get_gene_coord(c.genes[i], geneCoord);
        get_gene_coord(c.genes[i - 1], prevGeneCoord);

        score += sqrt(
                pow(fabs(geneCoord[0] - prevGeneCoord[0]), 2)
                + pow(fabs(geneCoord[1] - prevGeneCoord[1]), 2)
        );
    }

    if (!validate_chromosome(c))
        score = 0.1;

    return score;
}

/* Mean objective function value over population */
float mean(chromosome *generation, int size) {

    /* TODO: Optimize by caching the mean for the current generation */
    float mean = 0;

    for (int i = 0; i < size; i++)
        mean += score(generation[i]);

    mean /= size;

    return mean;
}

/* Objective score function: 1 / score */
float fitness(chromosome c) {

    float fitness = 1 / score(c);

    /* Store the max fitness of current generation to adapt our random pick */
    if (fitness > maxFitness)
        maxFitness = fitness;

    return fitness;
}


chromosome fittest(chromosome generation[], int size) {

    chromosome fittest;
    fittest.fitness = 0;

    for (int i = 0; i < size; i++)
        if (generation[i].fitness > fittest.fitness)
            fittest = generation[i];

    return fittest;
}



float rand_a_b(float a, float b) {

    return (rand() / (float) RAND_MAX) * (b - a) + a;
}


void select_chromosomes(chromosome generation[], int toSize, chromosome population[], int fromSize) {

    srand((unsigned)time(NULL));
    int picked = 0;

    for (int i = 0; i < toSize; i++) {
        picked = 0;

        while (!picked) {
            for (int j = i; j < fromSize; j++) {
                float wheel = rand_a_b(0, 1);

                if (wheel <= population[j].fitness) {
                    generation[i] = population[j];
                    picked = 1;
                }
            }
        }
    }
}


/*
 * Crossover method based on subsets of parents.
 * The offspring will be a valid chromosome.
 */
chromosome crossover(chromosome *candidates, int size) {

    chromosome dad, mom, kid;
    kid.fitness = 0;

    int firstPick = (int) random() % size;
    int secondPick = firstPick;

    dad = candidates[firstPick];

    while (secondPick == firstPick) secondPick = (int) random() % size;

    mom = candidates[secondPick];

    int a = (int) rand_a_b(1, NB_GENES - 1);
    int b = (int) rand_a_b(a, NB_GENES - 1);

    while (b < a) b = (int) rand_a_b(a, NB_GENES - 1);

    int dadSubset[b - a + 1];
    int k = 0;

    for (int i = a; i <= b; i++) {
        for (int j = 0; j < GENE_SIZE; j++)
            kid.genes[i][j] = dad.genes[i][j]; //That's for dad

        dadSubset[k++] = bin_to_dec(dad.genes[i]); //Keep track of the used cities
    }

    int skip = 0;
    k = 0;

    for (int i = 0; i < NB_GENES; i++) {
        if (i >= a && i <= b) continue;

        do {
            skip = 0;

            for (int j = 0; j < (b - a + 1); j++) {
                if (bin_to_dec(mom.genes[k]) == dadSubset[j]) {
                    skip = 1;
                    break;
                }
            }

            if (skip) k++;
        } while (skip);

        for (int j = 0; j < GENE_SIZE; j++)
            kid.genes[i][j] = mom.genes[k][j]; //That's for mom

        k++;
    }

    return kid;
}





/* Alternate mutation method: swap two points at random indexes
 * This will ensure a valid mutated chromosome */
chromosome mutate(chromosome c) {

    srand((unsigned)time(NULL));

    /* Spin the wheel to mutate */
    if (rand_a_b(0, 1) > MUTATION_RATE) {
        c.fitness = fitness(c);
        return c;
    }

    int a = (int) rand_a_b(1, NB_GENES - 1); //Only swap between the start and end points
    int b = a;
    int temp[GENE_SIZE] = {0};

    while (b == a) b = (int) rand_a_b(1, NB_GENES - 1); //Make sure the two indexes are different

    /* Copy gene a into temp, b into a, and a into b*/
    for (int i = 0; i < GENE_SIZE; i++) {
        temp[i] = c.genes[a][i];
        c.genes[a][i] = c.genes[b][i];
        c.genes[b][i] = temp[i];
    }

    c.fitness = fitness(c);

    return c;
}



int main() {

    chromosome population[POPULATION_SIZE];
    chromosome generation[GENERATION_SIZE];

    printf("* Generating %d candidates for base population...\n", POPULATION_SIZE);

    for (int p = 0; p < POPULATION_SIZE; p++) {
        chromosome c;
        int invalid = 1;

        /* Start from Lille, return to Lille */
        /* TODO: Make this scalable (if we add new cities) */
        c.genes[0][0] = c.genes[NB_GENES - 1][0] = 0;
        c.genes[0][1] = c.genes[NB_GENES - 1][1] = 0;
        c.genes[0][2] = c.genes[NB_GENES - 1][2] = 0;

        /* Reset C */
        for (int i = 1; i < NB_GENES - 1; i++)
            for (int j = 0; j < GENE_SIZE; j++)
                c.genes[i][j] = 0;


        for (int i = 1; i < NB_GENES - 1; i++) {

            invalid = 1;

            while (invalid) {
                invalid = 0;
                make_gene(c.genes[i]);

                /* Check validity */
                for (int j = 0; j < i; j++) {
                    if (j != i
                        && (c.genes[j][0] == c.genes[i][0])
                        && (c.genes[j][1] == c.genes[i][1])
                        && (c.genes[j][2] == c.genes[i][2])) {
                        invalid = 1;
                        break;
                    }
                }
            }
        }

        c.fitness = fitness(c);

        population[p] = c;
    }

    /* Selecting original candidates */
    select_chromosomes(generation, GENERATION_SIZE, population, POPULATION_SIZE);

    float firstSelectionScore = mean(generation, GENERATION_SIZE);
    float scoreOfFittest = score(fittest(generation, GENERATION_SIZE));

    printf("* First selection average score: %.2f\n", firstSelectionScore);
    printf("* Generation's fittest chromosome's score: %.2f\n", scoreOfFittest);

    int iteration = 0;

    while(iteration < 200) {

        system("clear");
        printf("* Generating %d candidates for base population...\n", POPULATION_SIZE);
        printf("* First selection average score: %.2f\n", firstSelectionScore);
        printf("* Generation's fittest chromosome's score: %.2f\n", scoreOfFittest);
        printf("* Iteration: %d\n", iteration++);

        /* Select GENERATION_SIZE individuals from population, based on their fitness or randomly */
        printf("* Selecting %d individuals from generation %d...\n", GENERATION_SIZE, iteration);

        chromosome nextGeneration[GENERATION_SIZE] = {0}; //The offsprings of the (intermediate) generation

        maxFitness = 0;

        /* Crossover and mutate from the selection */
        for (int i = 0; i < GENERATION_SIZE; i++)
            nextGeneration[i] = mutate(crossover(generation, GENERATION_SIZE));

        chromosome selection[GENERATION_SIZE] = {0}; //Selection of nextGeneration (offsprings) + base generation

        /*
         * Select new generation, based on the offsprings and the parent generation,
         * and the fitness of their chromosomes
         */
        select_chromosomes(selection, GENERATION_SIZE, nextGeneration, GENERATION_SIZE);
        //select_chromosomes(&selection[GENERATION_SIZE / 2], GENERATION_SIZE / 2, generation, GENERATION_SIZE);

        /* Replace generation by the selection -> cross-breed of old generation + next generation */
        for (int i = 0; i < GENERATION_SIZE; i++)
            generation[i] = selection[i];

        /* Update final mean score to give feedback */

        printf("* Generation average score: %.2f\n", mean(generation, GENERATION_SIZE));
        printf("* Generation's fittest chromosome's score: %.2f\n", score(fittest(generation, GENERATION_SIZE)));

        usleep(8000);
    }

    return 0;
}