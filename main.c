/*
 * Solving the traveling salesman problem using a genetic algorithm.
 * Author: Theo Morales <theo.morales.fr@gmail.com>
 * 2017
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>


#define POPULATION_SIZE 100
#define GENERATION_SIZE 15
#define ITERATIONS 50
#define NB_GENES 8
#define GENE_SIZE 3
#define CROSSOVER_MODE "SP" // SP:  Single Point, 2P: Two Point, U: Uniform
#define MUTATION_RATE 0.1

#define GENE_X(gene) gene_coord(gene)[0]
#define GENE_Y(gene) gene_coord(gene)[1]

typedef struct {
    int x;
    int y;
    int binaryIndex[3];
    char *name;
} city;


const city cities[8] = {
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
 * A gene is a 3-bit array, a chromosome contains 8 genes for the 7 cities + the first one
 * */
typedef struct {
    int genes[8][GENE_SIZE];
    float fitness;
} chromosome;



int bin_to_dec(int bin[GENE_SIZE]) { //TODO: Fix bug - received bin is not an array

    int dec = 0;

    for (int i = 0; i < GENE_SIZE; i ++)
        dec += bin[i] * pow(2, (GENE_SIZE - 1 - i));

    return dec;
}


/* Generate bits to fill a gene */
void make_gene(int *gene) {

    for (int i = 0; i < GENE_SIZE; i++) {
        //Currently generates both valid and invalid genes (same city 2+ times)
        gene[i] = rand() % 2;
    }
}


/*
 * Map a gene from its binary index to the corresponding city
 * and return its coordinates
 */
void get_gene_coord(int gene[], int coord[]) { //TODO: Fix bug - received bin is not an array

    int geneValue = bin_to_dec(gene);

    coord[0] = cities[geneValue].x;
    coord[1] = cities[geneValue].y;
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

/* Objective score function: score / mean */
float fitness(chromosome c, chromosome *generation, int generationSize) {

    return score(c) / mean(generation, generationSize);
}

/* Select one chromosome from the given generation */
chromosome select_chromosome(chromosome *generation, int size) {

    chromosome c = { 0 };

    for (int i = 0; i < size; i++) {
        float wheel = (rand()%100) / 100;

        /* Need to make sure that the fitness is under the form: 0.XX */
        if (wheel <= generation[i].fitness) {
            c = generation[i];
            break;
        }
    }

    return c;
}

/*
 * Crossover of two randomly selected chromosomes from the given population
 * Returns the offspring
 */
chromosome make_random_offspring(chromosome *candidates, int size) {

    chromosome dad, mom, kid;
    kid.fitness = 0;

    int firstPick = random() % size;
    int secondPick = random() % size;

    dad = candidates[firstPick];

    while (secondPick == firstPick)
        secondPick = random() % size;

    mom = candidates[secondPick];


    if (strcmp(CROSSOVER_MODE, "SP") == 0) {
        /* Single point crossover */
        int crossover_point = 1 + rand() % 7;

        for (int i = 0; i < crossover_point; i++)
            for (int j = 0; j < GENE_SIZE; j++)
                kid.genes[i][j] = dad.genes[i][j];

        for (int i = crossover_point; i < 8; i++)
            for (int j = 0; j < GENE_SIZE; j++)
                kid.genes[i][j] = mom.genes[i][j];

    } else if (strcmp(CROSSOVER_MODE, "DP") == 0) {
        /* Double point crossover */
        //TODO

    } else if (strcmp(CROSSOVER_MODE, "U") == 0) {
        /* Uniform crossover */
        //TODO
    }
}


chromosome mutate(chromosome c) {
    //TODO: Implement the mutation
}


int main() {

    chromosome population[POPULATION_SIZE];
    chromosome generation[GENERATION_SIZE];

    printf("* Generating %d candidates for base population...\n", POPULATION_SIZE);

    for (int i = 0; i < POPULATION_SIZE; i++) {
        chromosome c;
        c.fitness = 0;

        for (int i = 0; i < 8; i++)
            make_gene(c.genes[i]);

        population[i] = c;
    }

    /* Calculate fitness for each candidate in the population */
    for (int i = 0; i < POPULATION_SIZE; i++)
        population[i].fitness = fitness(population[i], population, POPULATION_SIZE);

    /* Selecting original candidates */
    for (int i = 0; i < GENERATION_SIZE; i++)
        generation[i] = select_chromosome(population, POPULATION_SIZE);


    int iteration = ITERATIONS;

    while(iteration--) {

        system("clear");
        printf("* Generating %d candidates for base population...\n", POPULATION_SIZE);
        printf("* Iteration: %d", iteration);

        /* Select GENERATION_SIZE individuals from population, based on their fitness or randomly */
        printf("* Selecting %d individuals from generation %d...", GENERATION_SIZE, ITERATIONS - iteration);

        chromosome nextGeneration[GENERATION_SIZE] = {0}; //The offsprings of the (intermediate) generation

        /* Crossover and mutate from the selection */
        for (int i = 0; i < GENERATION_SIZE; i++)
            nextGeneration[i] = make_random_offspring(generation, GENERATION_SIZE);

        /*
         * Update the fitness of the offsprings
         * Needs to be done after the crossover is complete because the fitness
         * is calculated depending on the mean of the whole generation
         */
        for (int i = 0; i < GENERATION_SIZE; i++)
            nextGeneration[i].fitness = fitness(nextGeneration[i], nextGeneration, GENERATION_SIZE);

        chromosome selection[GENERATION_SIZE] = {0}; //Selection of nextGeneration (offsprings) + base generation

        /*
         * Select new generation, based on the offsprings and the parent generation,
         * and the fitness of their chromosomes
         */
        for (int i = 0; i < GENERATION_SIZE / 2; i++)
            selection[i] = select_chromosome(generation, GENERATION_SIZE);


        for (int i = GENERATION_SIZE / 2; i < GENERATION_SIZE; i++)
            selection[i] = select_chromosome(generation, GENERATION_SIZE);
    }

    return 0;
}