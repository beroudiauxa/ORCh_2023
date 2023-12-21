#include "tools.h"

//---tools---

void give_random_number_from_0_to_1 (double &random_from_0_to_1) 
{
   int random = rand()%10000;
   random_from_0_to_1 = double(random)/10000.0;
}

int compare_numbers (const void *x, const void *y) 
{
   double xx = *(double*)x;
   double yy = *(double*)y;
   if (xx < yy) return -1;
   if (xx > yy) return 1;
   return 0;
}

void print_to_screen_with_newline (string a, int loglevel, int debuglevel)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      if (debuglevel > loglevel)
         cout << a << endl;
   }
}

void print_to_screen (string a, int loglevel, int debuglevel)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      if (debuglevel > loglevel)
         cout << a;
   }
}

void print_to_screen_with_newline (double a, int loglevel, int debuglevel)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      if (debuglevel > loglevel)
         cout << a << endl;
   }
}

void print_to_screen (double a, int loglevel, int debuglevel)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      if (debuglevel > loglevel)
         cout << a;
   }
}

void print_to_screen_with_newline (int a, int loglevel, int debuglevel)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      if (debuglevel > loglevel)
         cout << a << endl;
   }
}

void print_to_screen (int a, int loglevel, int debuglevel)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      if (debuglevel > loglevel)
         cout << a;
   }
}
