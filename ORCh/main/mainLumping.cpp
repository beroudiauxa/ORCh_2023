#include "mainLumping.h"

void Lumping(int debuglevel, string configuration, string mech, string trajectory_ref, string mech_desc)
{




   if (configuration == "MultipleInlet")
   {
      cout << "Lumping on MultipleInlet trajectory :" << endl << endl;
      cout << "Make sure the concentration files for lumping are up to date! " << endl << endl;
  //    getchar();
      
      lumping *l = new lumping();
      l->lumping_mech(mech, trajectory_ref, mech_desc);
   }

   if (configuration == "AutoIgnition")
   {
      cout << "Lumping on AutoIgnition trajectory " << endl;
      cout << "Make sure the concentration files for lumping are up to date! " << endl;
  //    getchar();

      lumping *l = new lumping();
      l->lumping_mech(mech, trajectory_ref, mech_desc);
   }





}







