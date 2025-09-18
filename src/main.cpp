#include"class.h"




int main()
{
    srand (static_cast <unsigned> (time(0))); // Initialize the randomness

    string path ="output/data/";

    Electric_Field E(0.,1.33,800); // Champ Electrique avec : Angle Polaris, indice du milieu , longueur d'onde
    Setup S(100,90,90);  // Paramtres de la simulation : Nb de Frames, Angles de collection.

    int radius = 10; // Rayon de la sphre de dipoles
    
    Population Pop; //Declaration de la population
    int n=12; //Puis Paramtrage :
    for(double i=0;i<=n;i++)
        Pop.Add_Dipole(0,0,0,"a",radius*cos(i*2*M_PI/n),radius*sin(i*2*M_PI/n),0);

    Pop.Orient_Radial_out();

    S.Balise_Save_config(path+to_string(radius)+"_"+to_string(n));
    S.PolarPattern(E,Pop,5,path+to_string(radius)+"_"+to_string(n));

    return 0;
}
/*
 int main()
{
    srand (static_cast <unsigned> (time(0))); // Initialize the randomness

    int radius =10;

    Population P(0,0,0,0,0,0);
    P.ReadConfig_just_position("Thomson-Spheres/32-Thomson.txt");
    P.Orient_Radial_out();
    do{
        P.Place_Radial_out(radius);

        Setup S(10000,90,90);
        Electric_Field Ew(0.,1.33,800);



        S.Balise_Save_config("C:/Users/zacharie.behel/Desktop/HRS-f/Thomson-Spheres/Polar/32Dipoles_varR/"+to_string(radius));

        S.PolarPattern(Ew,P,10,"C:/Users/zacharie.behel/Desktop/HRS-f/Thomson-Spheres/Polar/32Dipoles_varR/"+to_string(radius));
        radius+=5;
    }while(radius<=1000);

    return 0;
}*/
/*
 int main()
{
    srand (static_cast <unsigned> (time(0))); // Initialize the randomness

    double fact=0.55;
    do{
        Population P(0,0,0,0,0,0);
        P.ReadConfig_just_position("Thomson-Spheres/32-Thomson.txt");

        P.Orient_Radial_out();
        P.Place_Radial_out(50);

        for(unsigned int i=0;i<P.Get_Nb_Dip();i++)
            P.Get_Dip(i)->Move_to_(P.Get_Dip(i)->getPosition().x(),P.Get_Dip(i)->getPosition().y(),fact*P.Get_Dip(i)->getPosition().z());

        P.Orient_Radial_out();

        Electric_Field Ew(0.,1.33,800);


        Setup S(10000,90,90);
        S.Select_EME(0);
        S.Select_MEE(0);
        S.Select_QEE(0);
        S.Balise_Save_config("Ellipse_Cigare/"+to_string(fact));

        S.PolarPattern(Ew,P,10,"Ellipse_Cigare/"+to_string(fact));
        fact+=0.05;
    }while(fact<=2);
    return 0;
}*/

/**********     MAIN FOR HELIX   int main()
{
    srand (static_cast <unsigned> (time(0))); // Initialize the randomness


    bool enantiomere = 0;
    for(unsigned int rayon = 165;rayon<200;rayon+=5)
    {
        Population Dipole_Mill(0,0,0,-25,-25,-50);
        Dipole_Mill.Build_Helix(100,rayon,50,2,enantiomere);


        Eliptic_Electric_Field Ew(0,0.,1.33,800);


        Setup S(1000,90,90);
        S.Balise_Save_config("");

        S.PolarPattern(Ew,Dipole_Mill,10,("Helix/test/R="+to_string(rayon)+"_e"+to_string(enantiomere)+"_"));
        cout<<"Done : "<<rayon<<"nm enantiomere "<<enantiomere<<"."<<endl;
        if(enantiomere == 0)
        {
            enantiomere = 1;
            rayon -= 5;
        }
        else enantiomere = 0;
//        Dipole_Mill.~Population();
    }
    return 0;
}
*/
