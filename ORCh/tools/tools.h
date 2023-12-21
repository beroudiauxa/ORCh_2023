#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

using namespace std;

//---------------------
   void give_random_number_from_0_to_1(double &random_0_to_1);
   int compare_numbers (const void *x, const void *y);

   void print_to_screen_with_newline (string a, int loglevel, int debuglevel);
   void print_to_screen (string a, int loglevel, int debuglevel);
   void print_to_screen_with_newline (double a, int loglevel, int debuglevel);
   void print_to_screen (double a, int loglevel, int debuglevel);
   void print_to_screen_with_newline (int a, int loglevel, int debuglevel);
   void print_to_screen (int a, int loglevel, int debuglevel);

#endif

