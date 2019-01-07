// Function prototypes
int* crossover(int*, int*, int);
void random_mut(int*, int);
int fitness(int*, int);
int* random_assignment(int);
int** get_next_gen(int**, int, int);

typedef struct {
    int index;
    int fitness;
} fit_to_index;