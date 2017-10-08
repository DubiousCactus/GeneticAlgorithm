/* Solving the traveling salesman problem using a genetic algorithm.
 * Author: Theo Morales <theo.morales.fr@gmail.com>
 * 2017
*/

#include <stdio.h>
#include <math.h>


#define POPULATION_SIZE 100
#define GENERATION_SIZE 15
#define ITERATIONS 50

typedef struct {
    int x;
    int y;
    int binaryIndex[3];
} city;

typedef enum cities {
    PARIS = (city) { .x = 4, .y = 6, .binaryIndex = { 0, 0, 1 } },
    REIMS = (city) { .x = 7, .y = 7, .binaryIndex = { 0, 1, 0 } },
    LYON = (city) { .x = 6, .y = 14, .binaryIndex = { 0, 1, 1 } },
    MARSEILLE = (city) { .x = 10, .y = 22, .binaryIndex = { 1, 0, 0 } },
    NANTES = (city) { .x = 2, .y = 9, .binaryIndex = { 1, 0, 1 } },
    LA_ROCHELLE = (city) { .x = 1, .y = 11, .binaryIndex = { 1, 1, 0 } },
    BORDEAUX = (city) { .x = 3, .y = 12, .binaryIndex = { 1, 1, 1 } }
} cities;

/* Describes the order to visit cities */
typedef struct {
    cities genes[7];
} chromosome;

/* Objective function: total length of the trip */
float score(chromosome c) {

    float score = 0;

    for (int i =  1; i < sizeof(c.genes) - 1; i++) {
        score += sqrt(
                pow(c.genes[i].x - c.genes[i - 1].x, 2)
                + pow(c.genes[i].y - c.genes[i - 1].y, 2)
        );
    }

    return score;
}

/* Mean objective function value over population */
float mean(chromosome *generation, int size) {

    float mean = 0;

    for (int i = 0; i < size; i++)
        mean += score(generation[i]);

    mean /= size;

    return mean;
}

/* Objective score function: score / mean */
float fitness() {

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