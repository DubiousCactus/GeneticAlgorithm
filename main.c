#include <stdio.h>

#define POPULATION_SIZE 100
#define GENERATION_SIZE 15
#define ITERATIONS 50

typedef struct {
    int x;
    int y;
    int binaryIndex[3];
} city;

typedef enum {
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


/* Mean objective function value over population */
float mean() {

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