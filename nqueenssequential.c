// Genetic algorithm to solve the n-queens problem for unspecified n 
#include "nqueens.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

int* crossover(int* strand1, int* strand2, int length) {
    int cutoff = (int) (rand() % length);
    
    int* new_strand = malloc(sizeof(int) * length);
    int* s1; 
    int* s2;

    // Determine (pseudo) randomly which strand to use for the first half
    int s = (int) (rand() % 2);

    if (s == 0) {
        s1 = strand1;
        s2 = strand2;
    } else {
        s1 = strand2;
        s2 = strand1;
    }

    for (int i = 0; i < cutoff; i++) {
        new_strand[i] = s1[i];
    }

    for (int i = cutoff; i < length; i++) {
        new_strand[i] = s2[i];
    }

    return new_strand;
}

/*
*   @effects - randomly mutates a given strand of DNA by changing a character
*/
void random_mut(int* strand, int length) {

    // Position of random mutation
    int mut_position = (int) (rand() % length);
    // Value of random mutation
    int mutation = (int) (rand() % length);
    strand[mut_position] = mutation;
}


/*
*   @effects - generates a number representing the fitness of a strand of DNA - in our case, it is the number of attacks
*/
int fitness(int* strand, int board_dim) {
    
    // Six possible attack directions for each queen
    // Strand specifies position within column - we need to create a 2D array, since it will be easier to work with
    int num_attacks = 0;

    // 0's will mean empty, 1 occupied
    int** board = malloc(sizeof(int*) * board_dim);

    // Construct the board
    for (int i = 0; i < board_dim; i++) {
        int* row = calloc(board_dim, sizeof(int));
        int pos = strand[i];
        row[pos] = 1;
        board[i] = row;
    }

    // Check and count all possible attacks
    int queenSeen = 0;

    // Horizontal attacks
    for (int i = 0; i < board_dim; i++) {
        for (int j = 0; j < board_dim; j++) {
            if (board[i][j] != 0 && queenSeen == 1) {
                num_attacks++;
            } else if (board[i][j] != 0) {
                queenSeen = 1;
            }
        }
        queenSeen = 0;
    }

    // Vertical attacks
    for (int i = 0; i < board_dim; i++) {
        for (int j = 0; j < board_dim; j++) {
            if (board[j][i] != 0 && queenSeen == 1) {
                num_attacks++;
            } else if (board[j][i] != 0) {
                queenSeen = 1;
            }
        }
        queenSeen = 0;
    }

    // Diagonal attacks 
    // TODO: These are kinda more difficult, I don't wanna do them right now 

    // Free resources
    for (int i = 0; i < board_dim; i++) {
        free(board[i]);
    }
    free(board);

    return num_attacks;
}

/*
*   @effects - generates a random assignment of n values to a strand
*/
int* random_assignment(int n) {
    
    // Allocate space for the strand of DNA
    int* strand = calloc(n, sizeof(int));

    for (int i = 0; i < n; i++) {
        int r = (rand() % n);
        strand[i] = r;
    }
    printf("Strand is ");
    for (int i = 0; i < n; i++) {
        printf("%i", strand[i]);
    }
    printf("\n");

    return strand;
}

// Custom comparator for use with bultin qsort()
int comparator(const void* lhs, const void* rhs) {
    return ((*(fit_to_index*)lhs).fitness - (*(fit_to_index*)rhs).fitness);
}

int** get_next_gen(int** strands, int num_strands, int length) {
    fit_to_index fitness_scores[num_strands];

    int** new_gen = malloc(sizeof(int*) * num_strands);

    // Calculate all fitness scores for the current generation
    for (int i = 0; i < num_strands; i++) {
        fitness_scores[i].index = i;
        fitness_scores[i].fitness = fitness(strands[i], length);
    }


    // Keep best 1/4 of the parents, use them to create new children
    // Start by sorting 
    qsort(fitness_scores, num_strands, sizeof(fit_to_index), &comparator);

    int cutoff_point = (int) floor(num_strands/4);

    // Add first 1/4 most fit 
    for (int i = 0; i < cutoff_point; i++) {
        new_gen[i] = strands[fitness_scores[i].index];
    }

    // Generate some new strands using the fittest 1/4 of gene pool
    for (int i = cutoff_point; i < num_strands; i++) {
        // Generate two random numbers to select parents from the 'fit' pool
        int parent1 = (int) (rand() % cutoff_point);
        int parent2 = (int) (rand() % cutoff_point);

        int* strand1 = new_gen[parent1];
        int* strand2 = new_gen[parent2];

        // Crossover and random mutation
        int* new_strand = crossover(strand1, strand2, length);
        random_mut(new_strand, length);

        new_gen[i] = new_strand;
    }

    return new_gen;

}

// Used at the end to return the single fittest strand from the current gene pool
int* get_fittest_strand(int** strands, int num_strands, int dimension) {
    int * fitnesses = malloc(sizeof(int) * num_strands);
    for (int i = 0; i < num_strands; i++) {
        fitnesses[i] = fitness(strands[i], dimension);
    }

    int index = -1;
    int fittest = INT_MAX;
    for (int i = 0; i < num_strands; i++) {
        if (fitnesses[i] < fittest) {
            index = i;
            fittest = fitnesses[i];
        }
    }

    if (index != -1) {
        return strands[index];
    }

    return NULL;
}


// Used for debugging to print out current genepool
void printstrands(int** strands, int num_strands, int length) {

    for (int i = 0; i < num_strands; i++) {
        printf("Strand %i: ", i);
        for (int j = 0; j < length; j++) {
            printf("%i", strands[i][j]);
        }
        printf("\n");
    }

}

int main(int argc, char* argv[]) {
    
    char* ptr;

    if (argc != 4) {
        printf("Usage: ./genalgoseq.out <dimension> <num_strands> <num_generations>");
        exit(EXIT_FAILURE);
        return -1;
    }

    srand(time(0));
    int n = (int) strtol(argv[1], &ptr, 10);
    int num_strands = (int) strtol(argv[2], &ptr, 10);
    int num_generations = (int) strtol(argv[3], &ptr, 10);

    // Allocate space for strands of DNA
    printf("Num_strands is %i\n", num_strands);
    printf("n is %i\n", n);
    printf("num_generations is %i\n", num_generations);
    int** strands = malloc(sizeof(int*) * num_strands);

    // Generate n strands of DNA
    for (int i = 0; i < num_strands; i++) {
        strands[i] = random_assignment(n);
    }

    printstrands(strands, num_strands, n);

    // TODO: Loop until an 'optimal' solution is found?
    // What defines an optimal solution for this problem? 
    // Should this loop until a SIGINT is sent?
    while(num_generations != 0) {
        strands = get_next_gen(strands, num_strands, n);
        printstrands(strands, num_strands, n);
        num_generations--;
    }

    int* fittest_strand = get_fittest_strand(strands, num_strands, n);

    num_generations = (int) strtol(argv[3], &ptr, 10);

    if (fittest_strand != NULL) {
        printf("The fittest strand is ");
        for (int i = 0; i < n; i++) {
            printf("%i", fittest_strand[i]);
        }
        int fit = fitness(fittest_strand, n);
        printf(" with fitness score %i", fit);
        printf(" for an %i x %i board with %i strands and %i generations.\n", n, n, num_strands, num_generations);
        
        exit(EXIT_SUCCESS);
    } else {
        printf("Exiting with a failure...\n");
        exit(EXIT_FAILURE);
        return -1;
    }

    return 0;

}