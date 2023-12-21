#ifndef QSSscenario_H
#define QSSscenario_H

#include <iostream>
#include <vector>

using namespace std;


//--------------------
class QSSscenario
{
   public:
   //constructeur
   QSSscenario(vector<string> QSSspecies); 

   //destructeur
   virtual ~QSSscenario();

   public:
   vector<string> m_QSSspecies;

};


#endif

