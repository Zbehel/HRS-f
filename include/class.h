#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <fstream>
#include <chrono>
#include <functional>
#include <random>

#ifndef CLASS_H_INCLUDED
#define CLASS_H_INCLUDED

// #define M_PI 3.14159265358979323846264338327915288

using namespace std;

complex<double> const i_C(0.,1.);

double enoise(double a, double b);
vector<complex<double>> operator*(complex<double> const & a,vector<complex<double>>const & b);
complex<double> D(double w);

/****************************************************************************************/
/*******************     Classes constituant les �l�ments de base            ************/
/****************************************************************************************/
class Position
{
protected:
    double x_,y_,z_;
public:
    Position(double x=0,double y=0,double z=0)
    :x_(x),y_(y),z_(z)
    {}

    double x()const;
    double y()const;
    double z()const;
    double norme()const;

    void Deplacer_de_(double x,double y,double z);
    void Deplacer_en_(double x,double y,double z);

    void operator+=(vector<double> V);
    //Position operator+(Position p,vector<double>const& V);

    ~Position(){}
};

class Euler_Angles
{
protected:
    double Psi_,Theta_,Phi_;

public:
    Euler_Angles(double Psi=0,double Theta=0,double Phi=0)
    :Psi_(Psi),Theta_(Theta),Phi_(Phi)
    {}

    double Psi()const;
    double Theta()const;
    double Phi()const;

    void Set_Angles(double psi,double theta,double phi);
    vector<vector<double>>const Matrice_Passage() const;
    ~Euler_Angles(){}
};

class T_Hyperpolarizabilite
{
/****************************************************************************************/
protected:
    vector<vector<vector<complex<double>>>> Beta;

public:
    T_Hyperpolarizabilite( string name="" )
    {
        Beta.resize(3);
        for(unsigned int i=0;i<3;i++)
        {
            Beta[i].resize(3);
            for(unsigned int j=0;j<3;j++)
                Beta[i][j].resize(3);
        }
        if(name == "LN")
        {
            Beta[1][1][1] = 2.4;
            Beta[2][0][0] = -4.52;
            Beta[2][2][2] = 31.5;
        }
        else
        {
            Beta[2][2][2] = 1;
            //Beta[2][0][0] = 1;
            //Beta[2][1][1] = 1;
        }
    }

    T_Hyperpolarizabilite(vector<vector<vector<complex<double>>>> const& B)
    {
        Beta.resize(3);
        for(unsigned int i=0;i<3;i++)
        {
            Beta[i].resize(3);
            for(unsigned int j=0;j<3;j++)
            {
                Beta[i][j].resize(3);
                for(unsigned int k=0;k<3;k++)
                    {Beta[i][j][k] = B[i][j][k];
                    //cout<<i<<" "<<j<<" "<<k<<endl;
                    }
            }
        }
    }

    T_Hyperpolarizabilite& operator=(T_Hyperpolarizabilite B);
    T_Hyperpolarizabilite operator=(vector<vector<vector<complex<double>>>> const& B);
    complex<double> Get_Element(double i,double j, double k)const;
    void Change_element(int i,int j,int k,complex<double>new_val);
    ~T_Hyperpolarizabilite(){}

};
/*******************     Fin des Classes constituant les �l�ments de base    ************/
/****************************************************************************************/

vector<complex<double>>const TM_DotProduct(T_Hyperpolarizabilite const TensA,vector<vector<complex<double>>>const & MatB);
vector<double>const MVProduct(vector<vector<double>>const& M,vector<double> V);
vector<vector<double>>const MM_Product(vector<vector<double>>const& M1,vector<vector<double>>const& M2);
vector<vector<complex<double>>>const VV_Dyadic_Prod(vector<complex<double>>const &VecA,vector<complex<double>>const &VecB);
T_Hyperpolarizabilite const FrameTransformation(T_Hyperpolarizabilite const T1,vector<vector<double>>const &M);
vector<complex<double>>const VV_Cross_Product(vector<complex<double>>const& VecA,vector<complex<double>>const& VecB);
complex<double>const VV_Scal_Prod(vector<complex<double>> const& A, Position const& B);
double norm(vector<double>const& A);


class Electric_Field
{
protected :
vector<complex<double>> Composantes_;

complex<double> k;

public:
    Electric_Field()
    {
        Composantes_.resize(3);
    }

    Electric_Field(double gamma, complex<double> n,double lambda) //Constructor Fundamental Field for now must propagates along Z...
    {
        Composantes_.resize(3);
        Composantes_[0] = cos(gamma*M_PI/180.);
        Composantes_[1] = sin(gamma*M_PI/180.);
        Composantes_[2] = 0.;

        Set_k(n,lambda);
    }

    Electric_Field(vector<complex<double>> E)
    {
        Composantes_.resize(3);
        if(E.size()==3)
        {
            Composantes_[0] = E[0];
            Composantes_[1] = E[1];
            Composantes_[2] = E[2];
        }

    }

    Electric_Field(Electric_Field const* E)
    {
        Composantes_.resize(3);

        Composantes_[0] = E->Comp()[0];
        Composantes_[1] = E->Comp()[1];
        Composantes_[2] = E->Comp()[2];
    }

    void Rotate_Field(double gamma);


    void Set_k(complex<double> n,double lambda);
    complex<double> Get_k() const;
    vector<complex<double>> const& Comp()const;
    complex<double> Amplitude(Position const& P);

    void operator+=(Electric_Field const& E);
    void operator-=(Electric_Field const& E);

    void operator*=(complex<double> const& a);

    ~Electric_Field(){}

};
 //Electric_Field operator+(Electric_Field const& E1,Electric_Field const& E2);
 //Electric_Field operator*(Electric_Field const& E1,complex<double> const& a)const;

class Eliptic_Electric_Field : public Electric_Field
{
protected :
    double Alpha_;
public :
    Eliptic_Electric_Field(double Alpha,double gamma, complex<double> n,double lambda) //Constructor Fundamental Field for now must propagates along Z...
    :Electric_Field(gamma,n,lambda),Alpha_(Alpha*M_PI/180.)
    {}
    void Rotate_Field(double gamma,double Alpha);
    ~Eliptic_Electric_Field(){}
};



 class Dipole
 {
 protected:
    Euler_Angles orientation_ ;
    string name_;
    Position coord_;
    T_Hyperpolarizabilite Beta_;
    T_Hyperpolarizabilite Beta_eme;
    T_Hyperpolarizabilite Beta_mee;
    vector<T_Hyperpolarizabilite> Beta_qee;

 public:
     Dipole(double Psi=0,double Theta=0,double Phi=0,string nom="",double x=0,double y=0,double z=0)
     :orientation_(Psi,Theta,Phi),name_(nom),coord_(x,y,z),Beta_(nom),Beta_eme(nom),Beta_mee(nom)
     {}

    const Position& getPosition()const;
    const Euler_Angles& getOrientation()const;
    const string& getNom()const;

    void SetPosition(double x,double y, double z);
    void SetOrientation(double psi,double theta,double phi);
    void Randomize_Orientation();
    void Translate_(double x,double y,double z);
    void Move_to_(double x,double y,double z);

    void Randomize_Position(double bornes);

    void Set_Beta(vector<vector<vector<complex<double>>>> const& H);

    void Build_Beta_MEE(double r,double pitch, bool enan);
    void Build_Beta_QEE(double r,double pitch, bool enan);
    void Build_Beta_EEE_EME(double r,double pitch, bool enan);

    T_Hyperpolarizabilite const Beta_Dip()const;
    T_Hyperpolarizabilite const Beta_EME()const;
    T_Hyperpolarizabilite const Beta_MEE()const;
    T_Hyperpolarizabilite const Beta_QEE(unsigned int i)const;
    void Save_Beta(string file);
    ~Dipole(){}
 };

class Population
{
protected:
    vector<Dipole*> liste_;
    vector<Dipole*> conv_liste_;
    Euler_Angles Orient_pop;
    Position Syst_coord_;


public:
    Population(double psi=0.,double theta=0.,double phi=0.,double x=0.,double y=0.,double z=0.)
    :Orient_pop(psi,theta,phi),Syst_coord_(x,y,z)
    {}

    void Add_Dipole(double Psi=0,double Theta=0,double Phi=0,string nom="",double x=0,double y=0,double z=0){liste_.push_back(new Dipole(Psi,Theta,Phi,nom,x,y,z));}
    void Add_Dipole(Position const& A, Euler_Angles const& B, string nom="");
    unsigned int Get_Nb_Dip()const;
    Euler_Angles Get_Orientation()const;
    Position Get_Position()const{return Syst_coord_;}
    T_Hyperpolarizabilite const Get_Beta_Lab (unsigned int dip_no)const;
    T_Hyperpolarizabilite const Get_Beta_MEE_Lab (unsigned int dip_no)const;
    T_Hyperpolarizabilite const Get_Beta_EME_Lab (unsigned int dip_no)const;
    T_Hyperpolarizabilite const Get_Beta_QEE_Lab (unsigned int dip_no,unsigned int n)const;

    void Change_Orientation(double psi,double theta,double phi);

    Dipole* Get_Dip(unsigned int i)const;
    Dipole* Get_Dip_in_Lab(unsigned int i)const;

    void Redraw_Coord();
    void Place_Element_in_Lab_Frame();

    void Randomize_Orientation();
    void Randomize_Population(double bornes);

    void Translate_Pop(double x,double y, double z);
    void Move_Pop(double x,double y, double z);
    void Place_Radial_out(double distance);

    void Orient_Radial_in();
    void Orient_Radial_out();

    void ReadConfig_just_position(string file);
    void ReadConfig_with_Angles(string file);

    void Build_Helix(unsigned int resolution,double radius,double pitch_height, int N_loop, bool Chir);
    void Build_Hemisphere(string file);
    void Build_Cylinder(double radius, double length,double d_phi,double d_l);

    void SaveConfig(string file);
    /*~Population()
    {
        for(unsigned int i=0;i<liste_.size();i++)
        {
            delete liste_[i];
        }
        for(unsigned int i=0;i<conv_liste_.size();i++)
        {
            delete conv_liste_[i];
        }
    }*/

    // Custom copy constructor
    Population(const Population& other)
    : liste_(other.liste_), Orient_pop(other.Orient_pop), Syst_coord_(other.Syst_coord_)
    {
        // Do not copy conv_liste_ pointers to avoid double-free or sharing issues.
        // conv_liste_ starts empty.
    }

    void operator+=(Population P);
    Population& operator=(Population P);

    ~Population()
    {
        // We do not own liste_ elements (shared), so we don't delete them.
        // We own conv_liste_ elements, so we delete them.
        for(unsigned int i=0;i<conv_liste_.size();i++)
        {
            delete conv_liste_[i];
        }
        conv_liste_.clear();
    }
 };


class Setup
{
 protected:
    unsigned int Frame_;
    double Ang_Col_u;   // Angle autour de l'axe X' ( aussi theta ).
    double Ang_Col_v;   // Angle autour de l'axe Z ( aussi phi )

    vector<double> SommeV;
    vector<double> SommeH;
    vector<double> Gamma;

    bool Save_config;
    bool eme,MEE,QEE;

 public:
    Setup(int nF=1,double u=90.,double v=90.)
    :Frame_(nF),Ang_Col_u(u),Ang_Col_v(v)
    {
        Save_config = 0;
        eme = false;
        MEE = false;
        QEE = false;
    }

    unsigned int get_nFrame()const;

    void Mov_u(double a);
    void Mov_v(double a);

    void Balise_Save_config(string name);
    bool Check_save_config();

    Electric_Field*const Get_Field_Amplitude(const Population& P,Electric_Field& Ew)const;
    Electric_Field*const Full_p_developement(const Population& P,Electric_Field& Ew)const;
    Electric_Field*const Get_M_Contribution(const Population& P,Electric_Field& Ew)const;
    Electric_Field*const Get_RetardationContrib(const Population& P,Electric_Field& Ew)const;

    void PolarPattern(Eliptic_Electric_Field& Ew,Population& Pop,double dGamma=10,string path="");
    void PolarPattern(Electric_Field& Ew,Population& Pop,double dGamma=10,string path="");
    void PolarPattern(Electric_Field& Ew,vector<Population> Pop,double dGamma,string path);
    
    // Enhanced optimized versions
    void PolarPattern_Optimized(Electric_Field& Ew, Population& Pop, double dGamma, string path);
    void ecrire_with_metadata(const vector<double>& t1, const vector<double>& t2, 
                             const vector<double>& t3, string name, 
                             int frames, double dGamma, long duration_ms) const;

    void Select_eme(bool val);
    void Select_MEE(bool val);
    void Select_QEE(bool val);

    bool Treat_eme() const;
    bool Treat_MEE() const;
    bool Treat_QEE() const;

    //void AngularPattern

    void ecrire(vector<double> t1, vector<double> t2, vector<double> t3, string name) const;

    // Refactored core simulation logic
    template <typename PopT, typename FieldT>
    void RunParallelSimulation(
        PopT& Pop, 
        FieldT& Ew, 
        double dGamma, 
        double start_angle, 
        const string& path,
        std::function<void(PopT&)> update_pop,
        std::function<void(PopT&, FieldT&, std::vector<Electric_Field>&, int, double)> calc_field
    );

    ~Setup(){}
 };


/******************************************************************************************************/
/****************                    FONCTIONS ANNEXES              ***********************************/
/******************************************************************************************************/

void Save_Config(vector<Population>const& Pop,string path);



#endif // CLASS_H_INCLUDED
