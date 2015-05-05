#include <stdlib.h>
#include "utils.h"

float* create1Darray(int dimention)
{
    int i;
    float* vect = malloc(dimention * sizeof(float));
    for ( i = 0; i < dimention; i++)
    {
        //vect[i] = rand();
        vect[i] = i;
    }
    return vect;
}
