#include <fstream>
#include <sstream>
#include <cstdlib>
#include "read.h"

//---Read---

Read::Read() //Constructeur
{}

void Read::Read_species(string mech, vector<Species_ORCh*> &listSpecies) const
{


   string mech_string;
   Input_file_into_string(mech, mech_string);

   int nSpecies = 0;

//---Species Array contains the list of species involved---
   string species_array;
   find_XML_key(mech_string, species_array, "speciesArray");

   stringstream s(species_array);

   while (!s.eof())
   {
      string Name;
      s >> Name;

      if (Name != "")
      {
         listSpecies.push_back(new Species_ORCh(Name, "", vector<double> (), 0, 0, vector<double> (), 0, 0, 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, false));
         nSpecies += 1;
      }
   }


//---Species Data contains the characteristics of the species involved---
   string species_data;
   find_XML_key(mech_string, species_data, "speciesData");

   for (int n=0; n<nSpecies; n++)
   {
      string search_init = "<species name=\"";
      search_init.append(listSpecies[n]->m_Name).append("\">");
      size_t init_left = species_data.find(search_init);

      if (init_left != string::npos)
      {
         string search_end = "</species>";
         size_t end_right = species_data.find(search_end, init_left);
         listSpecies[n]->m_Description = species_data.substr( init_left, 
                                                            (end_right+search_end.length()-(init_left)) );
         string NasaArray_1;
         string NasaArray_2;
    
         string Nasa_lowT_highT_1;
         string Nasa_lowT_highT_2;

         string getSpeciesDescription = listSpecies[n]->m_Description;

         //get species LJ_welldepth
         string LJ_welldepth_description;
         find_XML_key(getSpeciesDescription, LJ_welldepth_description, "LJ_welldepth");
         listSpecies[n]->m_LJ_welldepth = atof(LJ_welldepth_description.c_str());

         //get species LJ_diameter
         string LJ_diameter_description;
         find_XML_key(getSpeciesDescription, LJ_diameter_description, "LJ_diameter");
         listSpecies[n]->m_LJ_diameter = atof(LJ_diameter_description.c_str());

         //get species dipoleMoment
         string dipoleMoment_description;
         find_XML_key(getSpeciesDescription, dipoleMoment_description, "dipoleMoment");
         listSpecies[n]->m_dipoleMoment = atof(dipoleMoment_description.c_str());

         //get species polarizability
         string polarizability_description;
         find_XML_key(getSpeciesDescription, polarizability_description, "polarizability");
         listSpecies[n]->m_polarizability = atof(polarizability_description.c_str());

         //get species rotRelax
         string rotRelax_description;
         find_XML_key(getSpeciesDescription, rotRelax_description, "rotRelax");
         listSpecies[n]->m_rot_relax = atof(rotRelax_description.c_str());

         //get species geometry
         string geometry_description;
         find_XML_key(getSpeciesDescription, geometry_description, "string");
         listSpecies[n]->m_geometry = geometry_description;





         //get species number of atoms
         string AtomDescription;
         vector<int> AtomCoeffs;
         vector<string> AtomList;

         find_XML_key(getSpeciesDescription, AtomDescription, "atomArray");

         Atom_into_Species(AtomDescription, AtomList, AtomCoeffs);

         for (int i=0; i<AtomList.size(); i++)
         {
            if (AtomList[i] == "C")
               listSpecies[n]->m_C = AtomCoeffs[i];
            if (AtomList[i] == "H")
               listSpecies[n]->m_H = AtomCoeffs[i];
            if (AtomList[i] == "O")
               listSpecies[n]->m_O = AtomCoeffs[i];
            if (AtomList[i] == "N")
               listSpecies[n]->m_N = AtomCoeffs[i];
            if (AtomList[i] == "Ar")
               listSpecies[n]->m_Ar = AtomCoeffs[i];
         }




         find_XML_key(getSpeciesDescription, NasaArray_1, "NASA", Nasa_lowT_highT_1);
         find_XML_key(getSpeciesDescription, NasaArray_2, "NASA", Nasa_lowT_highT_2);

         string Tmin1;
         string Tmin2;
         string Tmax1;
         string Tmax2;

         translate_type(Nasa_lowT_highT_1, "Tmin", Tmin1);
         translate_type(Nasa_lowT_highT_2, "Tmin", Tmin2);

         translate_type(Nasa_lowT_highT_1, "Tmax", Tmax1);
         translate_type(Nasa_lowT_highT_2, "Tmax", Tmax2);

         vector<double> NasaCoeffs1;
         vector<double> NasaCoeffs2;
         
         string NasaCoefficients1;
         string NasaCoefficients2;

         find_XML_key(NasaArray_1, NasaCoefficients1, "floatArray");
         find_XML_key(NasaArray_2, NasaCoefficients2, "floatArray");

         //get the Nasa coefficients 1
         int nCoeff = 0;
 
         stringstream s1(NasaCoefficients1);
         while (!s1.eof())
         {
            string getCoeff;
            s1 >> getCoeff;
 
            if (getCoeff != "")
            {
               NasaCoeffs1.push_back(atof(getCoeff.c_str()));
               nCoeff += 1;
            }
         }
 
         //get the Nasa coefficients 2
         nCoeff = 0;
 
         stringstream s2(NasaCoefficients2);
         while (!s2.eof())
         {
            string getCoeff;
            s2 >> getCoeff;
 
            if (getCoeff != "")
            {
               NasaCoeffs2.push_back(atof(getCoeff.c_str()));
               nCoeff += 1;
            }
         }

         double Tmin_1 = atof(Tmin1.c_str());
         double Tmin_2 = atof(Tmin2.c_str());
 
         double Tmax_1 = atof(Tmax1.c_str());
         double Tmax_2 = atof(Tmax2.c_str());
 
         if (Tmin_1 > Tmin_2)
         {

            //2 is the Low temperature region
            listSpecies[n]->m_NASA_lowT_min = Tmin_2;
            listSpecies[n]->m_NASA_lowT_max = Tmax_2;
            for (unsigned int c=0; c<NasaCoeffs2.size(); c++)
               listSpecies[n]->m_NASACoeffs_lowT.push_back(NasaCoeffs2[c]);
 
            //1 is the High temperature region
            listSpecies[n]->m_NASA_highT_min = Tmin_1;
            listSpecies[n]->m_NASA_highT_max = Tmax_1;
            for (unsigned int c=0; c<NasaCoeffs1.size(); c++)
               listSpecies[n]->m_NASACoeffs_highT.push_back(NasaCoeffs1[c]);
         }
         else
         {
            //1 is the Low temperature region
            listSpecies[n]->m_NASA_lowT_min = Tmin_1;
            listSpecies[n]->m_NASA_lowT_max = Tmax_1;
            for (unsigned int c=0; c<NasaCoeffs1.size(); c++)
               listSpecies[n]->m_NASACoeffs_lowT.push_back(NasaCoeffs1[c]);
 
            //2 is the High temperature region
            listSpecies[n]->m_NASA_highT_min = Tmin_2;
            listSpecies[n]->m_NASA_highT_max = Tmax_2;
            for (unsigned int c=0; c<NasaCoeffs2.size(); c++)
               listSpecies[n]->m_NASACoeffs_highT.push_back(NasaCoeffs2[c]);
         }
      }
   }

}









void Read::Read_reactions(string mech, vector<Reaction_ORCh*> &listReactions) const
{

  string mech_string;
  Input_file_into_string(mech, mech_string);

  string Allreactions_data;
  find_XML_key(mech_string, Allreactions_data, "reactionData");

  string reaction_description;
  string reaction_type_duplicate_reversible_id;

  bool foundReaction = true;

  while (foundReaction)
  {
     find_XML_key(Allreactions_data, reaction_description, "reaction", 
                  reaction_type_duplicate_reversible_id, foundReaction);

     if (foundReaction)
     {
         string type;
         string duplicate_str;
         string reversible_str;
         string ID;

         translate_type(reaction_type_duplicate_reversible_id, "type", type);
         translate_type(reaction_type_duplicate_reversible_id, "duplicate", duplicate_str);
         translate_type(reaction_type_duplicate_reversible_id, "reversible", reversible_str);
         translate_type(reaction_type_duplicate_reversible_id, "id", ID);

         bool reversible = false;
         if (reversible_str == "yes")
            reversible = true;

         bool duplicate = false;
         if (duplicate_str == "yes")
            duplicate = true;

         string equation;
         string Arrhenius_description;
         string reactant_list;
         string product_list;

         find_XML_key(reaction_description, equation, "equation");
         find_XML_key(reaction_description, reactant_list, "reactants");
         find_XML_key(reaction_description, product_list, "products");
         find_XML_key(reaction_description, Arrhenius_description, "Arrhenius");

         double A;
         double b;
         double E;
         string Etype;

         vector<double> ReactantsCoeffs;
         vector<double> ProductsCoeffs;

         vector<string> ReactantsSpeciesList;
         vector<string> ProductsSpeciesList;

         get_Arrhenius_coefficients(Arrhenius_description, A, b, E, Etype);

         Species_into_Reactants_Products(reactant_list, ReactantsSpeciesList, ReactantsCoeffs);
         Species_into_Reactants_Products(product_list, ProductsSpeciesList, ProductsCoeffs);

         if (type == "")
         {
            listReactions.push_back(new Simple(A, b, E, ID, 
                        equation, Etype, reversible, duplicate,
                        ReactantsSpeciesList, ProductsSpeciesList,
                        ReactantsCoeffs, ProductsCoeffs)); 
         }


         if (type == "threeBody" || type == "falloff")
         {
            string efficiencies;
            vector<double> EfficienciesCoeffs;
            vector<string> EfficienciesSpeciesList;

            find_XML_key(reaction_description, efficiencies, "efficiencies");
            Species_into_Efficiencies(efficiencies, EfficienciesSpeciesList, EfficienciesCoeffs);

            if (type == "threeBody")
            {
               listReactions.push_back(new ThreeBody(A, b, E, ID, 
                           equation, Etype, reversible, duplicate,
                           ReactantsSpeciesList, ProductsSpeciesList,
                           ReactantsCoeffs, ProductsCoeffs, 
                           EfficienciesCoeffs, EfficienciesSpeciesList));
            }

            if (type == "falloff")
            {

               double A_lowP;
               double b_lowP;
               double E_lowP;
               string Etype_lowP;

               string Arrhenius_falloff_description;
               find_XML_key(reaction_description, Arrhenius_falloff_description, 
                            "Arrhenius");
               get_Arrhenius_coefficients(Arrhenius_falloff_description, 
                                          A_lowP, b_lowP, E_lowP, Etype_lowP);

               string falloff_type;
               string Troe_Lindemann;
               string falloff_description;

               find_XML_key(reaction_description, falloff_description, "falloff", falloff_type);
               translate_type(falloff_type, "type", Troe_Lindemann); 

               if (Troe_Lindemann == "Troe")
               {

                  vector<double> TroeCoefficients;
                  getTroeCoefficients(falloff_description, TroeCoefficients);

                  listReactions.push_back(new Troe(A, b, E, ID, 
                              equation, Etype, reversible, duplicate,
                              ReactantsSpeciesList, ProductsSpeciesList,
                              ReactantsCoeffs, ProductsCoeffs, 
                              EfficienciesCoeffs, EfficienciesSpeciesList,
                              A_lowP, b_lowP, E_lowP, Etype_lowP, TroeCoefficients));
               }
               else if (Troe_Lindemann == "Lindemann")
               {
                  listReactions.push_back(new Lindemann(A, b, E, ID, 
                              equation, Etype, reversible, duplicate,
                              ReactantsSpeciesList, ProductsSpeciesList,
                              ReactantsCoeffs, ProductsCoeffs, 
                              EfficienciesCoeffs, EfficienciesSpeciesList,
                              A_lowP, b_lowP, E_lowP, Etype_lowP));
               }
               else
               {
                  cout << "The type of the falloff is neither Troe nor Lindemann" << endl;
               }

            } //if type == falloff
         } //if type == threeBody || type == falloff
      }
   } //end while foundReaction



}


















void Read::Input_file_into_string(string mech, string& mech_string) const
{
   struct stat buf;
   if (stat(mech.c_str(), &buf) != -1) //Check if the mech exists
   {
      ifstream xml_file;
      xml_file.open(mech.c_str(), ios::in);

      string line_string = "";

      while (!xml_file.eof())
      {
         getline(xml_file,line_string);
         mech_string.append(line_string);
         mech_string.append("\n");
      }
   }
   else
   {
      cout << "ERROR, the mech isn't present" << endl;
   }
}

void Read::find_XML_key(string& initial_string, string& substracted_part, string keyword) const
{
     string search_init = "<";
     search_init.append(keyword);

     size_t init_left = initial_string.find(search_init);

     //If we found a string corresponding to the key
     if (init_left != string::npos)
     {
        size_t init_right = initial_string.find(">",init_left);

        string search_end = "</";
        search_end.append(keyword);
        size_t end_left = initial_string.find(search_end);
        size_t end_right = initial_string.find(">",end_left);

        substracted_part = initial_string.substr( init_right+1, (end_left-(init_right+1)) );
        initial_string.erase( init_left, end_right+1-(init_left) );
     }
}

void Read::find_XML_key(string& initial_string, string& substracted_part, 
                        string keyword, string& type) const
{
     string search_init = "<";
     search_init.append(keyword);

     size_t init_left = initial_string.find(search_init);

     //If we found a string corresponding to the key
     if (init_left != string::npos)
     {
        size_t init_right = initial_string.find(">",init_left);

        string search_end = "</";
        search_end.append(keyword);
        size_t end_left = initial_string.find(search_end);
        size_t end_right = initial_string.find(">",end_left);

        type = initial_string.substr( init_left+search_init.length(), 
                                            (init_right-(init_left+search_init.length())) );

        substracted_part = initial_string.substr( init_right+1, (end_left-(init_right+1)) );
        initial_string.erase( init_left, end_right+1-(init_left) );
     }
}

void Read::find_XML_key(string& initial_string, string& substracted_part, 
                        string keyword, string& type, bool &found) const
{
     string search_init = "<";
     search_init.append(keyword);

     size_t init_left = initial_string.find(search_init);

     //If we found a string corresponding to the key
     if (init_left != string::npos)
     {
        found = true;
      
        size_t init_right = initial_string.find(">",init_left);

        string search_end = "</";
        search_end.append(keyword);
        size_t end_left = initial_string.find(search_end);
        size_t end_right = initial_string.find(">",end_left);

        type = initial_string.substr( init_left+search_init.length(), 
                                            (init_right-(init_left+search_init.length())) );

        substracted_part = initial_string.substr( init_right+1, (end_left-(init_right+1)) );
        initial_string.erase( init_left, end_right+1-(init_left) );
     }
     else
     {
        found = false;
     }
}

void Read::translate_type(string initial_string, string type, string &associated_text) const
{
    string search_init = "=\"";
    type.append(search_init);
    size_t left = initial_string.find(type);

    if (left != string::npos)
    {
       size_t right = initial_string.find("\"",left+type.length());
       associated_text = initial_string.substr(left+type.length(), right-(left+type.length()));
    }
}

void Read::Species_into_Reactants_Products(string initial_string, vector<string>& Species_List, 
                                           vector<double>& Coeffs) const
{
    string search_species = ":";
    size_t species_place = 0;
    while (species_place != initial_string.length())
    {
       size_t species_place = initial_string.find(search_species);
       size_t blank = initial_string.find(" ");

       if (blank != 0)
       {
          string name = initial_string.substr(0, species_place);
          Species_List.push_back(name);
          string a  = initial_string.substr(species_place+1, blank-species_place);
          Coeffs.push_back(atof(a.c_str()));
          initial_string.erase(0, blank);
       }
       else
       {
          initial_string.erase(0, 1);
       }
    }
}

void Read::Atom_into_Species(string initial_string, vector<string>& Atom_List,
                                           vector<int>& Coeffs) const
{
    string search_atom = ":";
    size_t atom_place = 0;
    while (atom_place != initial_string.length())
    {
       size_t atom_place = initial_string.find(search_atom);
       size_t blank = initial_string.find(" ");

       if (blank != 0)
       {
          string name = initial_string.substr(0, atom_place);
          Atom_List.push_back(name);
          string a  = initial_string.substr(atom_place+1, blank-atom_place);
          Coeffs.push_back(atof(a.c_str()));
          initial_string.erase(0, blank);
       }
       else
       {
          initial_string.erase(0, 1);
       }
    }
}



void Read::Species_into_Efficiencies(string initial_string, vector<string>& Species_List, 
                                     vector<double>& Coeffs) const
{
    string search_species = ":";
    size_t species_place = 0;
    while (species_place != initial_string.length())
    {
       size_t species_place = initial_string.find(search_species);
       size_t blank = initial_string.find(" ");

       if (blank != 0)
       {
          string name = initial_string.substr(0, species_place);
          Species_List.push_back(name);
          string a  = initial_string.substr(species_place+1, blank-species_place);
          Coeffs.push_back(atof(a.c_str()));
          initial_string.erase(0, blank);
       }
       else
       {
          initial_string.erase(0, 1);
       }
    }
}


void Read::get_Arrhenius_coefficients(string Arrhenius_description, 
                                      double &A, double &b, double &E, string &Etype) const
{
   string get_A;
   string get_b;
   string get_E;

   find_XML_key(Arrhenius_description, get_A, "A");
   find_XML_key(Arrhenius_description, get_b, "b");
   find_XML_key(Arrhenius_description, get_E, "E", Etype);
    
   A = atof(get_A.c_str());
   b = atof(get_b.c_str());
   E = atof(get_E.c_str());
}


void Read::getTroeCoefficients(string initial_string, vector<double> &Coefficients) const
{
   size_t coeff_place = 0;
   while (coeff_place != initial_string.length())
   {
      size_t blank = initial_string.find(" ");

      if (blank != 0)
      {
         string a = initial_string.substr(0,blank);
         Coefficients.push_back(atof(a.c_str()));
         initial_string.erase(0, blank);
      }
      else
      {
         initial_string.erase(0, 1);
      }
   }
}


Read::~Read() //Destructeur
{}



