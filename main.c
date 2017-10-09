/*
 * Solving the traveling salesman problem using a genetic algorithm.
 * Author: Theo Morales <theo.morales.fr@gmail.com>
 * 2017
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define POPULATION_SIZE 100
#define GENERATION_SIZE 15
#define ITERATIONS 50
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

/* Describes the order to visit cities:
 * A gene is a 3-bit array, a chromosome contains 8 genes for the 7 cities + the first one
 * */
typedef struct {
    int genes[8][GENE_SIZE];
    float fitness;
} chromosome;



int bin_to_dec(int *bin) {

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


/* Map a gene from its binary index to the corresponding city
 * and return its coordinates
 */
int * gene_coord(int gene[GENE_SIZE]) {

    int coord[2] = { 0, 0 };

    coord[0] = cities[bin_to_dec(gene)].x;
    coord[1] = cities[bin_to_dec(gene)].y;

    return coord;
}

/* Objective function: total length of the trip */
float score(chromosome c) {

    float score = 0;

    for (int i =  1; i < sizeof(c.genes) - 1; i++) {
        score += sqrt(
                pow(GENE_X(c.genes[i]) - GENE_X(c.genes[i - 1]), 2)
                + pow(GENE_Y(c.genes[i]) - GENE_Y(c.genes[i - 1]), 2)
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


    int iteration = ITERATIONS;

    while(iteration--) {

        system("clear");
        printf("* Generating %d candidates for base population...\n", POPULATION_SIZE);
        printf("* Iteration: %d", iteration);

        /* Select GENERATION_SIZE individuals from population, based on their fitness or randomly */
        printf("* Selecting %d individuals from generation %d...", GENERATION_SIZE, ITERATIONS - iteration);

        /* TODO: Clean up this mess and refactor this spaghetti code */
        if (iteration == ITERATIONS - 1) { //Select from population if first iteration

            for (int i = 0; i < GENERATION_SIZE; i++) {
                for (int j = 0; j < POPULATION_SIZE; j++) {
                    float wheel = (rand()%100) / 100;

                    /* Need to make sure that the fitness is under the form: 0.XX
                     * and that the right comparison sign is used
                     */
                    if ((population[j].fitness <= 0.5 && wheel <= population[j].fitness)
                        || (population[j].fitness > 0.5 && wheel > population[j].fitness)){

                        generation[i] = population[j];
                    }
                }
            }

        } else { //Select from current generation

            chromosome nextGeneration[GENERATION_SIZE] = {0};

            for (int i = 0; i < GENERATION_SIZE; i++) {
                float wheel = (rand()%100) / 100;

                /* Need to make sure that the fitness is under the form: 0.XX
                     * and that the right comparison sign is used
                     */
                if ((generation[i].fitness <= 0.5 && wheel <= generation[i].fitness)
                    || (generation[i].fitness > 0.5 && wheel > generation[i].fitness)){

                    nextGeneration[i] = generation[i];
                }
            }

            /* Kill and replace previous generation */
            for (int i = 0; i < GENERATION_SIZE; i++)
                generation[i] = nextGeneration[i];
        }


        /* Crossover from the selection */

        /* Mutation from the offsprings */

        /* Replace generation by mutated offsprings */
    }

    return 0;
}