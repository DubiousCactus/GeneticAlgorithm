/* Solving the traveling salesman problem using a genetic algorithm.
 * Author: Theo Morales <theo.morales.fr@gmail.com>
 * 2017
*/

#include <stdio.h>
#include <math.h>


#define POPULATION_SIZE 100
#define GENERATION_SIZE 15
#define ITERATIONS 50

#define GENE_X(gene) gene_coord(gene)[0]
#define GENE_Y(gene) gene_coord(gene)[1]

typedef struct {
    int x;
    int y;
    int binaryIndex[3];
    char *name;
} city;


const city cities[8] = {
    { .x = 3, .y = 1, .binaryIndex = { 0, 0, 0}, .name = "Lille" },
    { .x = 4, .y = 6, .binaryIndex = { 0, 0, 1 }, .name = "Paris" },
    { .x = 7, .y = 7, .binaryIndex = { 0, 1, 0 }, .name = "Reims" },
    { .x = 6, .y = 14, .binaryIndex = { 0, 1, 1 }, .name =  "Lyon" },
    { .x = 10, .y = 22, .binaryIndex = { 1, 0, 0 }, .name = "Marseille" },
    { .x = 2, .y = 9, .binaryIndex = { 1, 0, 1 }, .name = "Nantes" },
    { .x = 1, .y = 11, .binaryIndex = { 1, 1, 0 }, .name = "La Rochelle" },
    { .x = 3, .y = 12, .binaryIndex = { 1, 1, 1 }, .name = "Bordeaux" }
};

/* Describes the order to visit cities:
 * A gene is a 3-bit array, a chromosome contains 8 genes for the 8 cities
 * */
typedef struct {
    int genes[8][3];
} chromosome;



int bin_to_dec(int *bin, int length) {

    int dec = 0;

    for (int i = 0; i < length; i ++)
        dec += bin[i] * pow(2, (length - 1 - i));

    return dec;
}


/* Map a gene from its binary index to the corresponding city
 * and return its coordinates
 */
int * gene_coord(int gene[3]) {

    int coord[2] = { 0, 0 };

    coord[0] = cities[bin_to_dec(gene, 3)].x;
    coord[1] = cities[bin_to_dec(gene, 3)].y;

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

    /* While ITERATIONS--

        /* Generate random population */

        /* Select GENERATION_SIZE individuals from population, based on their fitness or randomly */

        /* Crossover from the selection */

        /* Mutation from the offsprings */

        /* Replace generation by mutated offsprings */

    return 0;
}