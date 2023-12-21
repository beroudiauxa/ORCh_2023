#ifndef Write_QSS_FORTRAN_H
#define Write_QSS_FORTRAN_H

#include "read.h"

using namespace std;


//--------------------
class Write_QSS_FORTRAN
{
   public:
   //constructeur
   Write_QSS_FORTRAN();

   virtual void Write_QSS_file_in_FORTRAN(string mech, string write_mech, vector<bool> &QSS_Species) const;

   //destructeur
   virtual ~Write_QSS_FORTRAN();


   private:




};

#endif
















