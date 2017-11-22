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
#include <ncurses.h>


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


city cities[8];

/*
 * Describes the order to visit cities:
 * A gene is a 3-bit array, a chromosome contains 9 genes for the 8 cities + the first one
 * */
typedef struct {
    int genes[NB_GENES][GENE_SIZE];
    float fitness;
} chromosome;


float maxFitness = 0;
WINDOW *visualization_window, *details_window;

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

WINDOW *create_newwin(int height, int width, int starty, int startx)
{	WINDOW *local_win;

	local_win = newwin(height, width, starty, startx);
	box(local_win, 0 , 0);		/* 0, 0 gives default characters
					 * for the vertical and horizontal
					 * lines			*/
    refresh();
	wrefresh(local_win);		/* Show that box 		*/

	return local_win;
}

void destroy_win(WINDOW *local_win)
{
	/* box(local_win, ' ', ' '); : This won't produce the desired
	 * result of erasing the window. It will leave it's four corners
	 * and so an ugly remnant of window.
	 */
	wborder(local_win, ' ', ' ', ' ',' ',' ',' ',' ',' ');
	/* The parameters taken are
	 * 1. win: the window on which to operate
	 * 2. ls: character to be used for the left side of the window
	 * 3. rs: character to be used for the right side of the window
	 * 4. ts: character to be used for the top side of the window
	 * 5. bs: character to be used for the bottom side of the window
	 * 6. tl: character to be used for the top left corner of the window
	 * 7. tr: character to be used for the top right corner of the window
	 * 8. bl: character to be used for the bottom left corner of the window
	 * 9. br: character to be used for the bottom right corner of the window
	 */
	wrefresh(local_win);
	delwin(local_win);
}


/* Draw the TSP graph based on the given chromosome */
void visualize(chromosome journey)
{

    int prev_coord[2] = {-1};
    int prev_gene[GENE_SIZE];
    for (int j = 0; j < NB_GENES; j++) {
        int coord[2] = {0};
        int gene[GENE_SIZE];
        int not_drawing = 1;
        
        for (int k = 0; k < GENE_SIZE; k++)
           gene[k] = journey.genes[j][k];

        get_gene_coord(gene, coord);

        /* Draw the point (city) */
        wattron(visualization_window, COLOR_PAIR(1));
        mvwprintw(visualization_window, coord[1], coord[0], cities[bin_to_dec(gene)].name);
        mvwprintw(visualization_window, coord[1] + 1, coord[0] + (sizeof cities[bin_to_dec(gene)].name / sizeof *cities[bin_to_dec(gene)].name) / 2 - 1, "%d", j);
        wattroff(visualization_window, COLOR_PAIR(1));

        if (prev_coord[0] != -1 && prev_coord[1] != -1) {
            int from_y, to_y, from_x, to_x, reverse_x = 0;
            float slope, intercept, dx, dy;
            
            from_y = prev_coord[1];
            to_y = coord[1];
            from_x = prev_coord[0];
            to_x = coord[0];

            if (from_x > to_x) {
                from_x = coord[0];
                to_x = prev_coord[0];
                from_y = coord[1];
                to_y = prev_coord[1];
                reverse_x = 1;
            }

            dx = to_x - from_x;
            dy = to_y - from_y;
            slope = dy / dx;
            intercept = from_y - slope * from_x;

            int prev_y = -1, y = 0;

            wattron(visualization_window, COLOR_PAIR(2));
            for (int x = from_x; x < to_x; x++) {
                y = slope * x + intercept;
                if ((mvwinch(visualization_window, y, x) & A_CHARTEXT) == ' ' && y != prev_y) {
                    mvwprintw(visualization_window, y, x, "o");
                    prev_y = y;
                }
            }
            wattroff(visualization_window, COLOR_PAIR(2));
        }
        
        for (int k = 0; k < GENE_SIZE; k++)
           prev_gene[k] = gene[k];

        prev_coord[0] = coord[0];
        prev_coord[1] = coord[1];
    }
    
}



int main() {

    /* Init ncurses */
    initscr();
    noecho();
    cbreak();

    if(has_colors() == FALSE)
    {
        endwin();
        printf("Your terminal does not support color\n");
        exit(1);
    }

    start_color();
    init_pair(1, COLOR_BLUE, COLOR_YELLOW);
    init_pair(2, COLOR_GREEN, COLOR_BLACK);

    int yMax, xMax;
    getmaxyx(stdscr, yMax, xMax);

    memcpy(cities, (city[]) {
        { .x = xMax * 0.47, .y = yMax * 0.03, .binaryIndex = { 0, 0, 0 }, .name = "Lille" },
        { .x = xMax * 0.51, .y = yMax * 0.15, .binaryIndex = { 0, 0, 1 }, .name = "Paris" },
        { .x = xMax * 0.8, .y = yMax * 0.06, .binaryIndex = { 0, 1, 0 }, .name = "Reims" },
        { .x = xMax * 0.7, .y = yMax * 0.65, .binaryIndex = { 0, 1, 1 }, .name =  "Lyon" },
        { .x = xMax * 0.9, .y = yMax * 0.8, .binaryIndex = { 1, 0, 0 }, .name = "Marseille" },
        { .x = xMax * 0.15, .y = yMax * 0.3, .binaryIndex = { 1, 0, 1 }, .name = "Nantes" },
        { .x = xMax * 0.05, .y = yMax * 0.42, .binaryIndex = { 1, 1, 0 }, .name = "La Rochelle" },
        { .x = xMax * 0.2, .y = yMax * 0.59, .binaryIndex = { 1, 1, 1 }, .name = "Bordeaux" }
    }, sizeof cities);


    visualization_window = create_newwin(yMax - 7, xMax - 2, 0, 1);
    details_window = create_newwin(7, xMax - 2, yMax - 7, 1);

    chromosome population[POPULATION_SIZE];
    chromosome generation[GENERATION_SIZE];

    /*printf("* Generating %d candidates for base population...\n", POPULATION_SIZE);*/

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

    /*printf("* First selection average score: %.2f\n", firstSelectionScore);
    printf("* Generation's fittest chromosome's score: %.2f\n", scoreOfFittest);*/

    int iteration = 0;

    while(iteration < 200) {

        /*system("clear");
        printf("* Generating %d candidates for base population...\n", POPULATION_SIZE);
        printf("* First selection average score: %.2f\n", firstSelectionScore);
        printf("* Generation's fittest chromosome's score: %.2f\n", scoreOfFittest);
        printf("* Iteration: %d\n", iteration++);

        [> Select GENERATION_SIZE individuals from population, based on their fitness or randomly <]
        printf("* Selecting %d individuals from generation %d...\n", GENERATION_SIZE, iteration);*/

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

        mvwprintw(details_window, 1, 1, "* Generation %d -> chromosome 0", iteration);
        visualize(selection[0]); //Visualize each chromosome of the generation
        /* Replace generation by the selection -> cross-breed of old generation + next generation */
        /*for (int i = 0; i < GENERATION_SIZE; i++) {
            generation[i] = selection[i];
            visualize(generation[i]); //Visualize each chromosome of the generation
            usleep(50000);
        }*/

        /* Update final mean score to give feedback */

        /*printf("* Generation average score: %.2f\n", mean(generation, GENERATION_SIZE));
        printf("* Generation's fittest chromosome's score: %.2f\n", score(fittest(generation, GENERATION_SIZE)));*/

        mvwprintw(details_window, 2, 1, "* Generation %d -> average score: %.2f", iteration++, mean(generation, GENERATION_SIZE));
        mvwprintw(details_window, 3, 1, "               -> fittest chromosome's score: %.2f", score(fittest(generation, GENERATION_SIZE)));

        wrefresh(details_window);
        wrefresh(visualization_window);
        refresh();

        sleep(1);
        wclear(visualization_window);
        usleep(1000);
    }

    endwin();

    return 0;
}
