#include"class.h"
#include<utility>



using namespace std;

/******     M�thodes classe Position     *****/
    double Position::x()const{return x_;}
    double Position::y()const{return y_;}
    double Position::z()const{return z_;}
    double Position::norme()const{return sqrt(x_*x_+y_*y_+z_*z_);}

    void Position::Deplacer_de_(double x,double y,double z)
    {
        x_+=x;
        y_+=y;
        z_+=z;
    }

    void Position::Deplacer_en_(double x,double y,double z)
    {
        x_=x;
        y_=y;
        z_=z;
    }

    void Position::operator+=(vector<double>const V)
{
    if(V.size() == 3)
    {
        x_+=V[0];
        y_+=V[1];
        z_+=V[2];
    }
}
Position operator+(Position p,vector<double>const& V)
{
    p+=V;
    return p;
}
/*********************************************/
vector<double> operator+(vector<double>const& A,vector<double>const& B)
{
    vector<double> V;
    V.push_back(A[0]+B[0]);
    V.push_back(A[1]+B[1]);
    V.push_back(A[2]+B[2]);
    return V;
}

vector<double> operator-(vector<double>const& A,vector<double>const& B)
{
    vector<double> V;
    V.push_back(A[0]-B[0]);
    V.push_back(A[1]-B[1]);
    V.push_back(A[2]-B[2]);
    return V;
}

double norm(vector<double>const& A)
{
    return sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
}
/******     M�thodes classe Euler_Angles     *****/

    double Euler_Angles::Psi()const{return Psi_;}
    double Euler_Angles::Theta()const{return Theta_;}
    double Euler_Angles::Phi()const{return Phi_;}

    void Euler_Angles::Set_Angles(double psi,double theta,double phi)
    {
        Psi_ = psi;
        Theta_ = theta;
        Phi_ = phi;
    }

    vector<vector<double>>const Euler_Angles::Matrice_Passage() const
    {
        vector<vector<double>> M,M_Psi,M_Theta,M_Phi;
        M.resize(3);M_Psi.resize(3);M_Theta.resize(3); M_Phi.resize(3);
        for(unsigned int i=0;i<3;i++)
        {
            M[i].resize(3);
            M_Psi[i].resize(3);
            M_Theta[i].resize(3);
            M_Phi[i].resize(3);
        }
        //cout<<Psi()<<"    "<<Theta()<<"   "<<Phi()<<endl;
        double cpsi = cos(Psi()*M_PI/180.),  spsi = sin(Psi()*M_PI/180.);
        double ctheta = cos(Theta()*M_PI/180.),  stheta = sin(Theta()*M_PI/180.);
        double cphi = cos(Phi()*M_PI/180.),  sphi = sin(Phi()*M_PI/180.);


        //////
        M_Psi[0][0] =  cpsi ;
        M_Psi[0][1] =  -spsi ;
        M_Psi[0][2] =  0 ;

        M_Psi[1][0] =  spsi ;
        M_Psi[1][1] =  cpsi ;
        M_Psi[1][2] =   0;

        M_Psi[2][0] = 0;
        M_Psi[2][1] =  0;
        M_Psi[2][2] =  1;

        //////
        M_Theta[0][0] =  ctheta;
        M_Theta[0][1] =   0;
        M_Theta[0][2] =   stheta;

        M_Theta[1][0] =   0;
        M_Theta[1][1] =   1;
        M_Theta[1][2] =   0;

        M_Theta[2][0] = -stheta;
        M_Theta[2][1] =  0;
        M_Theta[2][2] =  ctheta;

        ///////
        M_Phi[0][0] =   cphi;
        M_Phi[0][1] =   -sphi;
        M_Phi[0][2] =   0;

        M_Phi[1][0] =   sphi;
        M_Phi[1][1] =   cphi;
        M_Phi[1][2] =   0;

        M_Phi[2][0] = 0;
        M_Phi[2][1] = 0;
        M_Phi[2][2] = 1;

        M = MM_Product(M_Psi,MM_Product(M_Theta,M_Phi));
        return M;
    }
/*************************************************/


/******     M�thodes classe T_Hyperpoarizabilite     *****/

complex<double> T_Hyperpolarizabilite::Get_Element(double i,double j, double k)const{return Beta[i][j][k];}

/*********************************************************/



/******     M�thodes classe Electric_Field     *****/
    void Electric_Field::Rotate_Field(double gamma)
    {
        gamma *= M_PI/180.;
        if(Composantes_.size()==3)
        {
            Composantes_[0] = cos(gamma);
            Composantes_[1] = sin(gamma);
            Composantes_[2] = 0.;
        }
        else cout<<"Fund Field Dim ERROR"<<endl;
    }

    void Eliptic_Electric_Field::Rotate_Field(double gamma,double Alpha)
    {
        gamma *= M_PI/180.;
        Alpha *= M_PI/180.;

        if(Composantes_.size()==3)
        {
            Composantes_[0] = cos(gamma)*((cos(Alpha))*(cos(Alpha)) + i_C*(sin(Alpha))*sin(Alpha)) + (1.-i_C)*sin(gamma)*sin(Alpha)*cos(Alpha);
            Composantes_[1] = sin(gamma)*(+i_C*(cos(Alpha))*(cos(Alpha)) + (sin(Alpha))*sin(Alpha)) + (1.-i_C)*cos(gamma)*sin(Alpha)*cos(Alpha);
            Composantes_[2] = 0.;
        }
        else cout<<"Fund Field Dim ERROR"<<endl;
    }

    void Electric_Field::Set_k(complex<double> n,double lambda){k = 2*M_PI*n/lambda;}
    complex<double> Electric_Field::Get_k() const{return k;}
    vector<complex<double>> const& Electric_Field::Comp()const{return Composantes_;}
    complex<double> Electric_Field::Amplitude(Position const& P)
    {
        complex<double> A(0.,0.);


        A = exp(i_C *Get_k()* P.z());
        return A;
    }

/******     M�thodes classe Dipole     *****/
    const Position& Dipole::getPosition()const{return coord_;}
    const Euler_Angles& Dipole::getOrientation()const{return orientation_;}
    const string& Dipole::getNom()const{return name_;}

    void Dipole::SetPosition(double x,double y, double z)
    {
        coord_.Deplacer_en_(x,y,z);
    }

    void Dipole::SetOrientation(double psi,double theta,double phi)
    {
        orientation_.Set_Angles(psi,theta,phi);
    }

    void Dipole::Randomize_Orientation()
    {
        orientation_.Set_Angles(enoise(0.,360.),180.*acos(enoise(-1.,1.))/M_PI,enoise(0.,360.));
    }

    void Dipole::Translate_(double x,double y,double z)
    {
        coord_.Deplacer_de_(x,y,z);
    }

    void Dipole::Move_to_(double x,double y,double z)
    {
        coord_.Deplacer_en_(x,y,z);
    }

    void Dipole::Randomize_Position(double bornes)
    {
        coord_.Deplacer_en_(enoise(-bornes,bornes),enoise(-bornes,bornes),enoise(-bornes,bornes));
    }


    T_Hyperpolarizabilite& T_Hyperpolarizabilite::operator=(T_Hyperpolarizabilite B)
    {
        swap(*this,B);
        return *this;
    }

    T_Hyperpolarizabilite T_Hyperpolarizabilite::operator=(vector<vector<vector<complex<double>>>> const& B)
    {
        Beta.resize(3);
        for(unsigned int i=0;i<3;i++)
        {
            Beta[i].resize(3);
            for(unsigned int j=0;j<3;j++)
            {
                Beta[i][j].resize(3);
                for(unsigned int k=0;k<3;k++)
                    Beta[i][j][k] = B[i][j][k];
            }
        }
        return Beta;
    }
    void T_Hyperpolarizabilite::Change_element(int i,int j,int k,complex<double>new_val)
    {
        Beta[i][j][k] = new_val;
    }


    void Dipole::Set_Beta(vector<vector<vector<complex<double>>>> const& H)
    {
        Beta_= T_Hyperpolarizabilite(H);
    }
    void Dipole::Build_Beta_EEE_EME(double r,double pitch, bool enan)
    {
         /************Init Tensor size************/
        vector<vector<vector<complex<double>>>> N,M;
        N.resize(3);
        M.resize(3);
        for(unsigned int i=0;i<3;i++)
        {
            N[i].resize(3);
            M[i].resize(3);
            for(unsigned int j=0;j<3;j++)
            {
                N[i][j].resize(3);
                M[i][j].resize(3);
                for(unsigned int k=0;k<3;k++)
                    {
                        N[i][j][k] = 0;
                        M[i][j][k] = 0;
                    }
            }
        }
        /*****************************************/
        /****** List of variables needed    ******/
        double w = 2*M_PI/800.;
        double e = 1;//-1.602e-19; //C
        double m = 1;//9.109e-31;  //kg
        double c = 1;//  nm/s !!

        double L = sqrt(r*r+pitch*pitch/(4*M_PI*M_PI));
        double a;
        if(enan) a = 1;
        else a = -1;
        double b = a*pitch*pitch/(16*pow(M_PI,4)*pow(L,3));

        N[0][1][1] = +(pow(e,3)*pow(r,3))/(2.*pow(D(w),2)*pow(L,4)*pow(m,2));
        N[0][1][2] = +(pitch*pow(e,3)*pow(r,2))/(2*M_PI*pow(D(w),2)*pow(L,4)*pow(m,2));
        N[0][2][2] = +(pow(pitch,2)*pow(e,3)*r)/(8.*pow(M_PI,2)*pow(D(w),2)*pow(L,4)*pow(m,2));

        N[1][1][1] = +(b*pow(e,3)*pow(r,3))/(D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        N[1][0][1] = +(pow(e,3)*pow(r,3))/(D(2*w)*D(w)*pow(L,4)*pow(m,2));
        N[1][1][2] = +(pitch*b*pow(e,3)*pow(r,2))/(M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        N[1][0][2] = +(pitch*pow(e,3)*pow(r,2))/(2*M_PI*D(2*w)*D(w)*pow(L,4)*pow(m,2));
        N[1][2][2] = +(pow(pitch,2)*b*pow(e,3)*r)/(4.*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));

        N[2][1][1] = +(pitch*b*pow(e,3)*pow(r,2))/(2*M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        N[2][0][1] = +(pitch*pow(e,3)*pow(r,2))/(2*M_PI*D(2*w)*D(w)*pow(L,4)*pow(m,2));
        N[2][1][2] = +(pow(pitch,2)*b*pow(e,3)*r)/(2.*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        N[2][0][2] = +(pow(pitch,2)*pow(e,3)*r)/(4.*pow(M_PI,2)*D(2*w)*D(w)*pow(L,4)*pow(m,2));
        N[2][2][2] = +(pow(pitch,3)*b*pow(e,3))/(8.*pow(M_PI,3)*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));

        Beta_= N;


        M[0][1][2] = +(i_C*pow(e,3)*pow(r,4)*w)/(2.*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[0][2][2] = +(i_C*pitch*pow(e,3)*pow(r,3)*w)/(4*M_PI*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[0][1][1] = -(i_C*pitch*pow(e,3)*pow(r,3)*w)/(4*M_PI*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[0][2][1] = -(i_C*pow(pitch,2)*pow(e,3)*pow(r,2)*w)/(8*pow(M_PI,2)*pow(D(w),2)*pow(L,4)*c*pow(m,2));

        M[1][1][2] = +(i_C*b*pow(e,3)*pow(r,4)*w)/(D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[1][0][2] = +(i_C*pow(e,3)*pow(r,4)*w)/(2.*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        M[1][2][2] = +(i_C*pitch*b*pow(e,3)*pow(r,3)*w)/(2*M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[1][1][1] = -(i_C*pitch*b*pow(e,3)*pow(r,3)*w)/(2*M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[1][1][0] = -(i_C*pitch*pow(e,3)*pow(r,3)*w)/(4*M_PI*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        M[1][0][1] = -(i_C*pitch*pow(e,3)*pow(r,3)*w)/(4*M_PI*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        M[1][2][1] = -(i_C*pow(pitch,2)*b*pow(e,3)*pow(r,2)*w)/(4*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[1][2][0] = -(i_C*pow(pitch,2)*pow(e,3)*pow(r,2)*w)/(8*pow(M_PI,2)*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));

        M[2][1][2] = +(i_C*pitch*b*pow(e,3)*pow(r,3)*w)/(2*M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[2][0][2] = +(i_C*pitch*pow(e,3)*pow(r,3)*w)/(4*M_PI*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        M[2][2][2] = +(i_C*pow(pitch,2)*b*pow(e,3)*pow(r,2)*w)/(4*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[2][1][1] = -(i_C*pow(pitch,2)*b*pow(e,3)*pow(r,2)*w)/(4*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[2][1][0] = -(i_C*pow(pitch,2)*pow(e,3)*pow(r,2)*w)/(8*pow(M_PI,2)*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        M[2][0][1] = -(i_C*pow(pitch,2)*pow(e,3)*pow(r,2)*w)/(8*pow(M_PI,2)*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        M[2][2][1] = -(i_C*pow(pitch,3)*b*pow(e,3)*r*w)/(8*pow(M_PI,3)*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        M[2][2][0] = -(i_C*pow(pitch,3)*pow(e,3)*r*w)/(16*pow(M_PI,3)*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));


        Beta_eme = M;

    }

    void Dipole::Build_Beta_MEE(double r,double pitch, bool enan)
    {
        /************Init Tensor size************/
        vector<vector<vector<complex<double>>>> N;
        N.resize(3);
        for(unsigned int i=0;i<3;i++)
        {
            N[i].resize(3);
            for(unsigned int j=0;j<3;j++)
            {
                N[i][j].resize(3);
                for(unsigned int k=0;k<3;k++)
                    N[i][j][k] = 0;
            }
        }
        /*****************************************/
        /****** List of variables needed    ******/
        double w = 2*M_PI/800.;
        double e = 1; //C
        double m = 1;  //kg
        double c = 1;//  nm/s !!

        double L = sqrt(r*r+pitch*pitch/(4*M_PI*M_PI));
        double a;
        if(enan) a = 1;
        else a = -1;
        double b = a*pitch*pitch/(16*pow(M_PI,4)*pow(L,3));

        N[1][1][1] = +(pitch*b*pow(e,3)*i_C*pow(r,3)*w)/(2*M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        N[1][0][1] = +(pitch*pow(e,3)*i_C*pow(r,3)*w)/(2*M_PI*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        N[1][1][2] = +(pow(pitch,2)*b*pow(e,3)*i_C*pow(r,2)*w)/(2.*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        N[1][0][2] = +(pow(pitch,2)*pow(e,3)*i_C*pow(r,2)*w)/(4.*pow(M_PI,2)*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        N[1][2][2] = +(pow(pitch,3)*b*pow(e,3)*i_C*r*w)/(8.*pow(M_PI,3)*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));

        N[2][1][1] = -(b*pow(e,3)*i_C*pow(r,4)*w)/(D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        N[2][0][1] = -(pow(e,3)*i_C*pow(r,4)*w)/(D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        N[2][1][2] = -(pitch*b*pow(e,3)*i_C*pow(r,3)*w)/(M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        N[2][0][2] = -(pitch*pow(e,3)*i_C*pow(r,3)*w)/(2*M_PI*D(2*w)*D(w)*pow(L,4)*c*pow(m,2));
        N[2][2][2] = -(pow(pitch,2)*b*pow(e,3)*i_C*pow(r,2)*w)/(4.*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*c*pow(m,2));
        Beta_mee= N;
    }


    void Dipole::Build_Beta_QEE(double r,double pitch, bool enan)
    {
        Beta_qee.resize(3);
        /************Init Tensor size************/
        vector<vector<vector<vector<complex<double>>>>> Q;
        Q.resize(3);
        for(unsigned int i=0;i<3;i++)
        {
            Q[i].resize(3);
            for(unsigned int j=0;j<3;j++)
            {
                Q[i][j].resize(3);
                for(unsigned int k=0;k<3;k++)
                    {
                        Q[i][j][k].resize(3);
                        for(unsigned int l=0;l<3;l++)
                            Q[i][j][k][l] = 0;
                    }
            }
        }
        /*****************************************/
        /****** List of variables needed    ******/
        double w = 2*M_PI/800.;
        double e = 1; //C
        double m = 1;  //kg
        double c = 1;//  nm/s !!

        double L = sqrt(r*r+pitch*pitch/(4*M_PI*M_PI));
        double a;
        if(enan) a = 1;
        else a = -1;
        double b = a*pitch*pitch/(16*pow(M_PI,4)*pow(L,3));

        Q[1][0][1][1] = +(b*pow(e,3)*pow(r,4))/(D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[1][0][0][1] = +(pow(e,3)*pow(r,4))/(D(2*w)*D(w)*pow(L,4)*pow(m,2));
        Q[1][0][1][2] = +(pitch*b*pow(e,3)*pow(r,3))/(M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[1][0][0][2] = +(pitch*pow(e,3)*pow(r,3))/(2*M_PI*D(2*w)*D(w)*pow(L,4)*pow(m,2));
        Q[1][0][2][2] = +(pow(pitch,2)*b*pow(e,3)*pow(r,2))/(4*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
/*
        Q[1][1][][] = (Ey*e^2*pow(r,3))/(D(w)*L^2*m);
        Q[1][1][][] = +(Ez*pitch*e^2*pow(r,2))/(2*M_PI*D(w)*L^2*m);

        Q[2][1][][] = -(pitch*e)/(2*M_PI);
*/

        Q[0][1][1][1] = +b*pow(e,3)*pow(r,4)/(D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[0][1][0][1] = +(pow(e,3)*pow(r,4))/(D(2*w)*D(w)*pow(L,4)*pow(m,2));
        Q[0][1][1][2] = +(pitch*b*pow(e,3)*pow(r,3))/(M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[0][1][0][2] = +(pitch*pow(e,3)*pow(r,3))/(2*M_PI*D(2*w)*D(w)*pow(L,4)*pow(m,2));
        Q[0][1][2][2] = +(pow(pitch,2)*b*pow(e,3)*pow(r,2))/(4*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));


        Q[0][2][1][1] = +(pitch*b*pow(e,3)*pow(r,3))/(2*M_PI*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[0][2][0][1] = +(pitch*pow(e,3)*pow(r,3))/(2*M_PI*D(2*w)*D(w)*pow(L,4)*pow(m,2));
        Q[0][2][1][2] = +(pow(pitch,2)*b*pow(e,3)*pow(r,2))/(2*pow(M_PI,2)*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[0][2][0][2] = +(pow(pitch,2)*pow(e,3)*pow(r,2))/(4*pow(M_PI,2)*D(2*w)*D(w)*pow(L,4)*pow(m,2));
        Q[0][2][2][2] = +(pow(pitch,3)*b*pow(e,3)*r)/(8*pow(M_PI,3)*D(2*w)*pow(D(w),2)*pow(L,4)*pow(m,2));

        Q[1][2][1][1] = -(pitch*pow(e,3)*pow(r,3))/(2*M_PI*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[1][2][1][2] = -(pow(pitch,2)*pow(e,3)*pow(r,2))/(2*pow(M_PI,2)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[1][2][2][2] = -(pow(pitch,3)*pow(e,3)*r)/(8*pow(M_PI,3)*pow(D(w),2)*pow(L,4)*pow(m,2));

        Q[2][2][1][1] = -(pow(pitch,2)*pow(e,3)*pow(r,2))/(4*pow(M_PI,2)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[2][2][1][2] = -(pow(pitch,3)*pow(e,3)*r)/(4*pow(M_PI,3)*pow(D(w),2)*pow(L,4)*pow(m,2));
        Q[2][2][2][2] = -(pow(pitch,4)*pow(e,3))/(16*pow(M_PI,4)*pow(D(w),2)*pow(L,4)*pow(m,2));

        for(unsigned int i=0;i<3;i++)
        {
            for(unsigned int j=0;j<3;j++)
            {
                for(unsigned int k=0;k<3;k++)
                    {
                        for(unsigned int l=0;l<3;l++)
                            Beta_qee[i].Change_element(j,k,l,Q[j][i][k][l]);
                    }
            }
        }
    }

    T_Hyperpolarizabilite const Dipole::Beta_Dip()const{return FrameTransformation(Beta_,getOrientation().Matrice_Passage());}
    T_Hyperpolarizabilite const Dipole::Beta_MEE()const{return FrameTransformation(Beta_mee,getOrientation().Matrice_Passage());}
    T_Hyperpolarizabilite const Dipole::Beta_EME()const{return FrameTransformation(Beta_eme,getOrientation().Matrice_Passage());}
    T_Hyperpolarizabilite const Dipole::Beta_QEE(unsigned int i)const{return FrameTransformation(Beta_qee[i],getOrientation().Matrice_Passage());}

    complex<double> D(double w)
    {
        return pow((1.e16/3e17),2)- w*w - 2.*i_C*w;
    }
/******     M�thodes classe Population     *****/

    void Population::Add_Dipole(Position const& A, Euler_Angles const& B, string nom)
    {
        Add_Dipole(B.Psi(),B.Theta(),B.Phi(),nom,A.x(),A.y(),A.z());
    }

    unsigned int Population::Get_Nb_Dip()const{return liste_.size();}
    Euler_Angles Population::Get_Orientation()const{ return Orient_pop;}

    T_Hyperpolarizabilite const Population::Get_Beta_Lab (unsigned int dip_no)const{return FrameTransformation(liste_[dip_no]->Beta_Dip(),Orient_pop.Matrice_Passage());}
    T_Hyperpolarizabilite const Population::Get_Beta_MEE_Lab (unsigned int dip_no)const{return FrameTransformation(liste_[dip_no]->Beta_MEE(),Orient_pop.Matrice_Passage());}
    T_Hyperpolarizabilite const Population::Get_Beta_EME_Lab (unsigned int dip_no)const{return FrameTransformation(liste_[dip_no]->Beta_EME(),Orient_pop.Matrice_Passage());}
    T_Hyperpolarizabilite const Population::Get_Beta_QEE_Lab (unsigned int dip_no,unsigned int n)const{return FrameTransformation(liste_[dip_no]->Beta_QEE(n),Orient_pop.Matrice_Passage());}

    void Population::Change_Orientation(double psi,double theta,double phi){Orient_pop.Set_Angles(psi,theta,phi);}
    Dipole* Population::Get_Dip(unsigned int i)const{return liste_[i];}
    Dipole* Population::Get_Dip_in_Lab(unsigned int i)const{return conv_liste_[i];}

    void Population::Randomize_Orientation()
    {
        //Orient_pop.Set_Angles(enoise(0.,360.),enoise(0.,180.),enoise(0.,360.));
        Orient_pop.Set_Angles(enoise(0.,360.),180.*acos(enoise(-1.,1.))/M_PI,enoise(0.,360.));

    }

    void Population::Randomize_Population(double bornes)
    {
        for(unsigned int i=0;i<liste_.size();i++)
        {
            liste_[i]->Randomize_Orientation();
            liste_[i]->Randomize_Position(bornes);
        }
    }

    void Population::Translate_Pop(double x,double y, double z)
    {
        for(unsigned int i=0;i<liste_.size();i++)
        {
            liste_[i]->Translate_(x,y,z);
        }
    }


    void Population::Build_Hemisphere(string file)
    {
        ifstream fichier(file.c_str(),ios::in);
        unsigned int i_size;
        fichier>>i_size;
        double x,y,z;
        string name;
        for(unsigned int i = 0;i<i_size;i++)
        {
            fichier>>name>>x>>y>>z;
            if(z>=0)    Add_Dipole(0,0,0,name,x,y,z);

        }
    }

    void Population::Build_Cylinder(double radius, double length,double d_phi,double d_l)
    {
        for(double i=0;i<=length;i+=d_l)
        {
            for(double j=0;j<360;j+=d_phi)
            {
                Add_Dipole(j,90,0,"Cyl",radius*cos(j*M_PI/180.),radius*sin(j*M_PI/180.),i);
            }
        }
    }

    void Population::Move_Pop(double x,double y, double z)
    {
        Syst_coord_.Deplacer_en_(x,y,z);
    }

    void Population::Redraw_Coord()
    {
        while(!conv_liste_.empty())
        {
            delete conv_liste_.back();
            conv_liste_.pop_back();
        }
        vector<double> temoin,temoin2; temoin.resize(3);temoin2.resize(3);
        vector<double> Pos_Lab; Pos_Lab.resize(3);
        for(unsigned int i=0;i<liste_.size();i++)
        {
            Pos_Lab[0] = Syst_coord_.x() + liste_[i]->getPosition().x();
            Pos_Lab[1] = Syst_coord_.y() + liste_[i]->getPosition().y();
            Pos_Lab[2] = Syst_coord_.z() + liste_[i]->getPosition().z();

            temoin[0] =  Syst_coord_.x() + liste_[i]->getPosition().x();
            temoin[1] =  Syst_coord_.y() + liste_[i]->getPosition().y();
            temoin[2] =  Syst_coord_.z() + liste_[i]->getPosition().z();

            temoin2 = MVProduct(liste_[i]->getOrientation().Matrice_Passage(),vector<double> {0.,0.,1.});
            temoin2 = MVProduct(Orient_pop.Matrice_Passage(),temoin2 + Pos_Lab);

            temoin = MVProduct(Orient_pop.Matrice_Passage(),Pos_Lab);

            double psi =  (180./M_PI)*atan2((temoin2[1]-temoin[1]),(temoin2[0]-temoin[0]));
            double theta = (180./M_PI)*acos((temoin2[2]-temoin[2])/norm(temoin2-temoin));
            conv_liste_.push_back(new Dipole(psi,theta,0,liste_[i]->getNom(),temoin[0],temoin[1],temoin[2]));
            //cout<<psi<<"\t"<<theta<<endl;
        }
    }

    void Population::Place_Element_in_Lab_Frame()
    {
        while(!conv_liste_.empty())
        {
            delete conv_liste_.back();
            conv_liste_.pop_back();
        }
        vector<double> temoin,temoin2; temoin.resize(3);temoin2.resize(3);
        vector<double> Pos_Lab; Pos_Lab.resize(3);
        for(unsigned int i=0;i<liste_.size();i++)
        {
            Pos_Lab[0] = Syst_coord_.x() + liste_[i]->getPosition().x();
            Pos_Lab[1] = Syst_coord_.y() + liste_[i]->getPosition().y();
            Pos_Lab[2] = Syst_coord_.z() + liste_[i]->getPosition().z();

            temoin[0] =  Syst_coord_.x() + liste_[i]->getPosition().x();
            temoin[1] =  Syst_coord_.y() + liste_[i]->getPosition().y();
            temoin[2] =  Syst_coord_.z() + liste_[i]->getPosition().z();

            temoin2 = MVProduct(liste_[i]->getOrientation().Matrice_Passage(),vector<double> {0.,0.,1.});
            temoin2 = MVProduct(Orient_pop.Matrice_Passage(),temoin2 + Pos_Lab);

            temoin = MVProduct(Orient_pop.Matrice_Passage(),Pos_Lab);

            //double psi =  (180./M_PI)*atan2((temoin2[1]-temoin[1]),(temoin2[0]-temoin[0]));
            //double theta = (180./M_PI)*acos((temoin2[2]-temoin[2])/norm(temoin2-temoin));
            conv_liste_.push_back(new Dipole(liste_[i]->getOrientation().Psi(),liste_[i]->getOrientation().Theta(),0,liste_[i]->getNom(),temoin[0],temoin[1],temoin[2]));
            //cout<<psi<<"\t"<<theta<<endl;
        }
    }

    void Population::Place_Radial_out(double distance)
    {
        for(unsigned int i=0;i<liste_.size();i++)
        {
            double cpsi = cos(liste_[i]->getOrientation().Psi()*M_PI/180.),  spsi = sin(liste_[i]->getOrientation().Psi()*M_PI/180.);
            double ctheta = cos(liste_[i]->getOrientation().Theta()*M_PI/180.),  stheta = sin(liste_[i]->getOrientation().Theta()*M_PI/180.);

            liste_[i]->SetPosition(distance*cpsi*stheta,distance*spsi*stheta,distance*ctheta);
        }

    }

    void Population::Orient_Radial_out()
    {
         for(unsigned int i=0;i<liste_.size();i++)
        {
            double x = liste_[i]->getPosition().x();    double y = liste_[i]->getPosition().y();    double z = liste_[i]->getPosition().z();

            double theta = 180*acos(z/sqrt(x*x+y*y+z*z))/M_PI;
            double psi = 180*atan2(y,x)/M_PI;
            liste_[i]->SetOrientation(psi,theta,0.);

        }
    }

    void Population::Orient_Radial_in()
    {
         for(unsigned int i=0;i<liste_.size();i++)
        {
            double x = liste_[i]->getPosition().x();    double y = liste_[i]->getPosition().y();    double z = liste_[i]->getPosition().z();

            double theta = 180.+ 180*acos(z/sqrt(x*x+y*y+z*z))/M_PI;
            double psi = 180*atan2(y,x)/M_PI;
            liste_[i]->SetOrientation(psi,theta,0.);

        }
    }

    void Population::ReadConfig_just_position(string file)
    {
        ifstream fichier(file.c_str(),ios::in);
        unsigned int i_size;
        fichier>>i_size;
        double x,y,z;
        string name;
        for(unsigned int i = 0;i<i_size;i++)
        {
            fichier>>name>>x>>y>>z;
            Add_Dipole(0,0,0,name,x,y,z);
        }
    }

    void Population::ReadConfig_with_Angles(string file)
    {
        ifstream fichier(file.c_str(),ios::in);
        unsigned int i_size;
        fichier>>i_size;
        double psi,theta,phi,x,y,z;
        string name="Helix";
        for(unsigned int i = 0;i<i_size;i++)
        {
            fichier>>x>>y>>z>>psi>>theta>>phi;

            Add_Dipole(psi,theta,phi,name,x,y,z);
        }
    }
    void Population::SaveConfig(string file)
    {
        vector<double> temoin; temoin.resize(3);
        vector<double> temoin2; temoin2.resize(3);
        ofstream fichier(file.c_str(), ios::out | ios::app);
        fichier<<2*liste_.size()<<endl<<"Disposition and Orientation"<<endl;
        
        if (conv_liste_.size() != liste_.size()) {
            cout << "Warning: conv_liste_ size mismatch in SaveConfig. Skipping." << endl;
            return;
        }

        for(unsigned int i=0;i<liste_.size();i++)
        {
            if(liste_[i]->getNom() =="") fichier<<"Dip\t";
            else fichier<<liste_[i]->getNom()<<"\t";

            temoin[0] = conv_liste_[i]->getPosition().x();
            temoin[1] = conv_liste_[i]->getPosition().y();
            temoin[2] = conv_liste_[i]->getPosition().z();
            //temoin = MVProduct(Orient_pop.Matrice_Passage(),temoin);

            fichier<<temoin[0]<<"\t"<<temoin[1]<<"\t"<<temoin[2]<<endl;
            //fichier<<conv_liste_[i]->getPosition().x()<<"\t"<<conv_liste_[i]->getPosition().y()<<"\t"<<conv_liste_[i]->getPosition().z()<<endl;
            temoin2 = MVProduct(Orient_pop.Matrice_Passage(), MVProduct(conv_liste_[i]->getOrientation().Matrice_Passage(),vector<double> {0.,0.,1.}));
            fichier<<"Top\t"<<temoin2[0]+temoin[0]<<"\t"<<temoin2[1]+temoin[1]<<"\t"<<temoin2[2]+temoin[2]<<endl;
        }
    }




    void Population::Build_Helix(unsigned int resolution,double radius,double pitch_height, int N_loop, bool Chir)
    {
        for(unsigned int i=0;i<resolution;i++)
        {
            double dip_x = (radius * cos(N_loop*2*M_PI*i/resolution)); //Position du dipole(N)
            double dip_y = (pow(-1,Chir) * radius * sin(N_loop*2*M_PI*i/resolution));
            double dip_z = ((-pitch_height*N_loop/2.)+(i*N_loop*pitch_height/resolution));


            double dip_x1 = (radius * cos(N_loop*2*M_PI*(i+1)/resolution));         //Position du dipole (N+1)  - Necessaire pour orienter les dipoles selon le cas 3.
            double dip_y1 = (pow(-1,Chir) * radius * sin(N_loop*2*M_PI*(i+1)/resolution));
            double dip_z1 = ((-pitch_height*N_loop/2.)+((i+1)*N_loop*pitch_height/resolution));

            double VecDip[3] = {dip_x1-dip_x,dip_y1-dip_y,dip_z1-dip_z};

            double psi = (180./M_PI) * acos((VecDip[0]/sqrt(VecDip[0]*VecDip[0]+VecDip[1]*VecDip[1])));
            double theta = (180./M_PI)*acos( VecDip[2] / sqrt(VecDip[0]*VecDip[0]+VecDip[1]*VecDip[1]+VecDip[2]*VecDip[2]) );
            if(VecDip[1]<=0.)   psi*=-1;
            Add_Dipole(psi,theta,0,"Helix",dip_x1,dip_y1,dip_z1);
        }
    }


/*********************************************************/
/******* Mathematical Functions **************************/
/*********************************************************/
T_Hyperpolarizabilite const FrameTransformation(T_Hyperpolarizabilite const T1,vector<vector<double>>const &M)
{
/*
    Calculates the new tensor in the new system of coordinates
    according to the matrix T for the transformation.
    TensA is the tensor in the OLD
    and TensB in the NEW FRAME.
*/
    vector<vector<vector<complex<double>>>> TR;
    TR.resize(3);
    for(int i=0;i<3;i++)
    {
        TR[i].resize(3);
        for(int j=0;j<3;j++)
        {
            TR[i][j].resize(3);
        }
    }

    int i, j, k, a, b, c;
    complex<double> somme;
    i = 0;
    do{
        j = 0;
        do{
            k = 0;
            do{
                somme = 0;
                a = 0;
                do{
                    b = 0;
                    do{
                        c = 0;
                        do{
                            somme += M[i][a] * M[j][b] * M[k][c] * T1.Get_Element(a,b,c);
                            c ++;
                        }while (c <= 2);
                        b ++;
                    }while (b <= 2);
                    a ++;
                }while ( a <= 2);

                TR[i][j][k] = somme;

                k ++;
            }while (k <= 2);
            j ++;
        }while (j <= 2);
        i ++;
    }while (i <= 2);

    return (T_Hyperpolarizabilite(TR));
}

vector<double>const MVProduct(vector<vector<double>>const& M,vector<double> V)
{
    vector<double> R;
    R.resize(3);
    for(int i =0;i<3;i++)
    {
        for(int j=0;j<3;j++)
            R[i]+=M[i][j]*V[j];
    }
    return R;
}

vector<vector<double>>const MM_Product(vector<vector<double>>const& M1,vector<vector<double>>const& M2)
{
    vector<vector<double>> M_Res;
    M_Res.resize(3);
    for(unsigned i=0;i<3;i++)M_Res[i].resize(3);

    unsigned int v;
    for(unsigned int i=0;i<M1.size();i++)
    {
        for(unsigned int j=0;j<M1[i].size();j++)
        {
            M_Res[i][j]=0;
            v=0;
            do
            {
                M_Res[i][j]+=M1[i][v]*M2[v][j];
                v++;
            }while(v<M2.size());
        }
    }
    return M_Res;
}

vector<complex<double>>const TM_DotProduct(T_Hyperpolarizabilite const TensA,vector<vector<complex<double>>>const & MatB)
{//Calculates the dot product of tensor and matrix

    vector<complex<double>> VecC;
    VecC.resize(3);
    complex<double>    somme;
    for(int i=0;i<=2;i++)
    {
        somme=0;
        for(int j=0;j<=2;j++)
        {
            for(int k=0;k<=2;k++)
            {
                somme += TensA.Get_Element(i,j,k)*MatB[j][k];
            }
        }
        VecC[i] = somme;
    }

    return VecC;
}

vector<vector<complex<double>>>const VV_Dyadic_Prod(vector<complex<double>>const &VecA,vector<complex<double>>const &VecB)
{//Calculates the dyadic product between two vectors

    vector<vector<complex<double>>> M_Dyad;
    M_Dyad.resize(3);
    for(unsigned int i=0;i<=2;i++)
    {
        M_Dyad[i].resize(3);
        for(unsigned int j=0;j<=2;j++)
             M_Dyad[i][j] = VecA[i] * VecB[j];
    }
    return M_Dyad;
}

vector<complex<double>>const VV_Cross_Product(vector<complex<double>>const& VecA,vector<complex<double>>const& VecB)
{
    vector<complex<double>> VecC;
    VecC.resize(3);

    VecC[0] = VecA[1] * VecB[2] - VecA[2] * VecB[1];
    VecC[1] = VecA[2] * VecB[0] - VecA[0] * VecB[2];
    VecC[2] = VecA[0] * VecB[1] - VecA[1] * VecB[0];

    return VecC;
}


complex<double>const VV_Scal_Prod(vector<complex<double>> const& A, Position const& B)
{
    return A[0]*B.x()+A[1]*B.y()+A[2]*B.z();
}






    void Save_Config(vector<Population>const& Pop,string path)
{
    int Total_n_dips = 0;
    for(unsigned int i=0;i<Pop.size();i++)
    {
        Total_n_dips += Pop[i].Get_Nb_Dip();
    }

    vector<double> temoin; temoin.resize(3);
    vector<double> temoin2; temoin2.resize(3);
    ofstream fichier(path.c_str(), ios::out | ios::app);
    fichier<<2*Total_n_dips<<endl<<"Disposition and Orientation"<<endl;

    for(unsigned int i=0;i<Pop.size();i++)
    {
        for(unsigned int j=0;j<Pop[i].Get_Nb_Dip();j++)
        {
            if(Pop[i].Get_Dip(j)->getNom() =="") fichier<<"Dip\t";
            else fichier<<Pop[i].Get_Dip(j)->getNom()<<"\t";

            temoin[0] = Pop[i].Get_Dip_in_Lab(j)->getPosition().x();
            temoin[1] = Pop[i].Get_Dip_in_Lab(j)->getPosition().y();
            temoin[2] = Pop[i].Get_Dip_in_Lab(j)->getPosition().z();
            temoin = MVProduct(Pop[i].Get_Orientation().Matrice_Passage(),temoin);
            fichier<<temoin[0]<<"\t"<<temoin[1]<<"\t"<<temoin[2]<<endl;
            //fichier<<conv_Pop[i].Get_Dip(j)->getPosition().x()<<"\t"<<conv_Pop[i].Get_Dip(j)->getPosition().y()<<"\t"<<conv_Pop[i].Get_Dip(j)->getPosition().z()<<endl;
            temoin2 = MVProduct(Pop[i].Get_Orientation().Matrice_Passage(), MVProduct(Pop[i].Get_Dip(j)->getOrientation().Matrice_Passage(),vector<double> {0.,0.,1.}));
            fichier<<"Top\t"<<temoin2[0]+temoin[0]<<"\t"<<temoin2[1]+temoin[1]<<"\t"<<temoin2[2]+temoin[2]<<endl;
        }
    }


}



/****************************************************************/
/******* Surcharge d'op�rateur pour Math. Op�rations ************/
/****************************************************************/
vector<complex<double>> operator*(complex<double> const & a,vector<complex<double>>const & b)
{
    vector<complex<double>> res;
    res.resize(b.size());
    for(unsigned int i=0;i<b.size();i++)
        res[i] = a*b[i];

    return res;
}

    void Electric_Field::operator+=(Electric_Field const& E)
    {
        Composantes_[0] += E.Comp()[0];
        Composantes_[1] += E.Comp()[1];
        Composantes_[2] += E.Comp()[2];
    }
    void Electric_Field::operator-=(Electric_Field const& E)
    {
        Composantes_[0] -= E.Comp()[0];
        Composantes_[1] -= E.Comp()[1];
        Composantes_[2] -= E.Comp()[2];
    }
     Electric_Field operator-(Electric_Field E1,Electric_Field const& E2)
    {
        E1-=E2;
        return (E1);
    }
         Electric_Field operator+(Electric_Field E1,Electric_Field const& E2)
    {
        E1+=E2;
        return (E1);
    }
    void Electric_Field::operator*=(complex<double> const& a)
    {
        Composantes_[0] *= a;
        Composantes_[1] *= a;
        Composantes_[2] *= a;
    }
    Electric_Field operator*(Electric_Field E1,complex<double> const& a)
    {
        E1*=a;
        return (E1);
    }

/*********************************************************/
/******* Surcharge d'op�rateur pour AFFICHAGE ************/
/*********************************************************/
ostream& operator<<(ostream& sortie,Electric_Field const& E)
{
    sortie<<"Ex = "<<E.Comp()[0]<<"\tEy = "<<E.Comp()[1]<<"\tEz = "<<E.Comp()[2];
    return sortie;
}

ostream& operator<<(ostream& sortie,Euler_Angles const& EA)
{
    sortie<<"("<<EA.Psi()<<", "<<EA.Theta()<<", "<<EA.Phi()<<")";
    return sortie;
}


ostream& operator<<(ostream& sortie,Position const& Pos)
{
    sortie<<"("<<Pos.x()<<", "<<Pos.y()<<", "<<Pos.z()<<")";
    return sortie;
}

ostream& operator<<(ostream& sortie,Dipole const& D)
{
    sortie<<"Dip : "<<D.getNom()<<"\tAngles : "<<D.getOrientation()<<"\tat : "<<D.getPosition();
    return sortie;
}

ostream& operator<<(ostream& sortie, T_Hyperpolarizabilite const& B)
{
    for(unsigned int i=0.;i<3;i++)
    {
        for(unsigned int j=0.;j<3;j++)
        {
            for(unsigned int k=0.;k<3;k++)
                sortie<<"Beta["<<i<<"]["<<j<<"]["<<k<<"]= "<<B.Get_Element(i,j,k)<<",\t";
            sortie<<endl;
        }
        sortie<<endl;
    }
    return sortie;
}


/*********************************************************/
/*******            Surcharge d'op�rateur     ************/
/*********************************************************/



void Population::operator+=(Population P)
{
    //Change_Orientation(P.Get_Orientation().Psi(),P.Get_Orientation().Theta(),P.Get_Orientation().Phi());
    //Move_Pop(P.Get_Position().x(),P.Get_Position().y(),P.Get_Position().z());
    P.Redraw_Coord();
    for(unsigned int i=0;i<P.Get_Nb_Dip();i++)
    {
        Add_Dipole(P.Get_Dip_in_Lab(i)->getPosition(),P.Get_Dip_in_Lab(i)->getOrientation(),P.Get_Dip_in_Lab(i)->getNom());
    }
}

Population operator+(Population a,Population const&b)
{
    a+=b;
    return a;
}

Population& Population::operator=(Population P)
{
    swap(*this,P);
    return (*this);
}

// Include optimized implementations
#ifdef HAS_OPENMP
#include <omp.h>
#endif
#include <chrono>

/******     M�thodes classe Setup     *****/
    unsigned int Setup::get_nFrame()const{return Frame_;}

    void Setup::Mov_u(double a){Ang_Col_u = a;}
    void Setup::Mov_v(double a){Ang_Col_v = a;}

    void Setup::Balise_Save_config(string name)
    {
        Save_config = true;
        name+="_Dip_Config.xyz";
        ofstream fichier(name.c_str(), ios::out |ios::trunc);
        fichier.close();
    }
    bool Setup::Check_save_config(){return Save_config;}

    void Setup::Select_eme(bool val){eme=val;}
    void Setup::Select_MEE(bool val){MEE=val;}
    void Setup::Select_QEE(bool val){QEE=val;}

    bool Setup::Treat_eme() const {return eme;}
    bool Setup::Treat_MEE() const {return MEE;}
    bool Setup::Treat_QEE() const {return QEE;}


    Electric_Field*const Setup::Get_Field_Amplitude(const Population& P,Electric_Field& Ew)const
    {
        vector<complex<double>> P2w,wN;
        vector<complex<double>> wD;
        P2w.resize(3);
        Electric_Field E2w;
        Electric_Field Ew_Cop(0.,1.33,800);
        E2w.Set_k(1.33,400.);
        wN.resize(3);
        wD.resize(3);

        wD[0] = 0.;wD[1] = 0.;wD[2] = 1.;

        wN[0] = cos(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[1] = sin(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[2] = cos(Ang_Col_u*M_PI / 180.);
        for(unsigned int i=0;i<P.Get_Nb_Dip();i++)
        {
            //cout<<*P.Get_Dip_in_Lab(i)<<endl;
            Ew_Cop = Ew*exp(i_C*Ew.Get_k()*VV_Scal_Prod(wD,P.Get_Dip_in_Lab(i)->getPosition()));
            P2w = TM_DotProduct(P.Get_Beta_Lab(i),VV_Dyadic_Prod(Ew_Cop.Comp(),Ew_Cop.Comp()));
            E2w += exp(-i_C*E2w.Get_k()*VV_Scal_Prod(wN,P.Get_Dip_in_Lab(i)->getPosition())) * VV_Cross_Product(wN,VV_Cross_Product(P2w,wN));
        }
        P2w.clear();
        wN.clear();
        wD.clear();
        return (new Electric_Field(E2w));
    }
    Electric_Field*const Setup::Full_p_developement(const Population& P,Electric_Field& Ew)const
    {
        vector<complex<double>> t;
        t.resize(3);
        t[0]=0;t[1]=0;t[2]=1;
        vector<complex<double>> B;B.resize(3);


        vector<complex<double>> P2w,wN;
        vector<complex<double>> wD;
        P2w.resize(3);
        Electric_Field E2w;
        Electric_Field Ew_Cop(0.,Ew.Get_k(),800);
        E2w.Set_k(1.33,400.);
        wN.resize(3);
        wD.resize(3);

        wD[0] = 0.;wD[1] = 0.;wD[2] = 1.;

        wN[0] = cos(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[1] = sin(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[2] = cos(Ang_Col_u*M_PI / 180.);
        for(unsigned int i=0;i<P.Get_Nb_Dip();i++)
        {
            Ew_Cop = Ew*exp(i_C*Ew.Get_k()*VV_Scal_Prod(wD,P.Get_Dip_in_Lab(i)->getPosition()));
            B = VV_Cross_Product(t,Ew.Comp());
            P2w = TM_DotProduct(P.Get_Beta_Lab(i),VV_Dyadic_Prod(Ew_Cop.Comp(),Ew_Cop.Comp()));
            E2w += exp(i_C*E2w.Get_k()*VV_Scal_Prod(wN,P.Get_Dip_in_Lab(i)->getPosition())) * VV_Cross_Product(wN,VV_Cross_Product(P2w,wN));
            P2w = TM_DotProduct(P.Get_Beta_EME_Lab(i),VV_Dyadic_Prod(Ew_Cop.Comp(),B));
            E2w += 2.*exp(i_C*E2w.Get_k()*VV_Scal_Prod(wN,P.Get_Dip_in_Lab(i)->getPosition())) * VV_Cross_Product(wN,VV_Cross_Product(P2w,wN));

        }
        P2w.clear();
        wN.clear();
        wD.clear();
        return (new Electric_Field(E2w));
    }

    Electric_Field*const Setup::Get_M_Contribution(const Population& P,Electric_Field& Ew)const
    {
        vector<complex<double>> M2w,wN;
        vector<complex<double>> wD;
        M2w.resize(3);
        Electric_Field E2w;
        Electric_Field Ew_Cop(0.,Ew.Get_k(),800);
        E2w.Set_k(1.33,400.);
        wN.resize(3);
        wD.resize(3);

        wD[0] = 0.;wD[1] = 0.;wD[2] = 1.;

        wN[0] = cos(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[1] = sin(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[2] = cos(Ang_Col_u*M_PI / 180.);

        for(unsigned int i=0;i<P.Get_Nb_Dip();i++)
        {
            Ew_Cop = Ew*exp(i_C*Ew.Get_k()*VV_Scal_Prod(wD,P.Get_Dip_in_Lab(i)->getPosition()));

            M2w = TM_DotProduct(P.Get_Beta_MEE_Lab(i),VV_Dyadic_Prod(Ew_Cop.Comp(),Ew_Cop.Comp()));
            E2w += exp(i_C*E2w.Get_k()*VV_Scal_Prod(wN,P.Get_Dip_in_Lab(i)->getPosition())) * VV_Cross_Product(wN,M2w);
        }
        M2w.clear();
        wN.clear();
        wD.clear();
        return (new Electric_Field(E2w));
    }

    Electric_Field*const Setup::Get_RetardationContrib(const Population& P,Electric_Field& Ew)const
    {
        vector<complex<double>> P2w,wN;
        vector<complex<double>> wD;
        P2w.resize(3);
        Electric_Field E2w;
        Electric_Field Ew_Cop(0.,Ew.Get_k(),800);
        E2w.Set_k(1.33,400.);
        wN.resize(3);
        wD.resize(3);

        wD[0] = 0.;wD[1] = 0.;wD[2] = 1.;

        wN[0] = cos(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[1] = sin(Ang_Col_v*M_PI / 180.) * sin(Ang_Col_u*M_PI / 180.);
        wN[2] = cos(Ang_Col_u*M_PI / 180.);
        for(unsigned int i=0;i<P.Get_Nb_Dip();i++)
        {
            Ew_Cop = Ew*exp(i_C*Ew.Get_k()*VV_Scal_Prod(wD,P.Get_Dip_in_Lab(i)->getPosition()));

            P2w = TM_DotProduct(P.Get_Beta_QEE_Lab(i,1),VV_Dyadic_Prod(Ew_Cop.Comp(),Ew_Cop.Comp()));
            E2w += (i_C*E2w.Get_k()/2.)* /*(1/P.Get_Dip_in_Lab(i)->getPosition().norme())*/exp(i_C*E2w.Get_k()*VV_Scal_Prod(wN,P.Get_Dip_in_Lab(i)->getPosition())) * VV_Cross_Product(wN,VV_Cross_Product(P2w,wN));
        }
        P2w.clear();
        wN.clear();
        wD.clear();
        return (new Electric_Field(E2w));
    }

// Thread-safe enoise implementation
double enoise(double a, double b)
{
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(a, b);
    return distribution(generator);
}


// Template implementation for RunParallelSimulation
template <typename PopT, typename FieldT>
void Setup::RunParallelSimulation(
    PopT& Pop, 
    FieldT& Ew, 
    double dGamma, 
    double start_angle, 
    const string& path,
    std::function<void(PopT&)> update_pop,
    std::function<void(PopT&, FieldT&, std::vector<Electric_Field>&, int, double)> calc_field
) {
    int nGamma = int(1 + (360. - start_angle) / dGamma);
    SommeV.assign(nGamma, 0.0);
    SommeH.assign(nGamma, 0.0);

    // Precompute Gamma values
    Gamma.resize(nGamma);
    for(int i=0; i<nGamma; ++i) Gamma[i] = start_angle + i*dGamma;

    int completed_frames = 0;

    #pragma omp parallel
    {
        // Thread-local copy of Population
        PopT Pop_local = Pop;
        
        // Thread-local accumulation
        vector<double> SommeV_local(nGamma, 0.0);
        vector<double> SommeH_local(nGamma, 0.0);
        
        // Reusable Electric_Field objects
        vector<Electric_Field> E2w_local(nGamma); 

        #pragma omp for schedule(dynamic)
        for(unsigned int f=0; f<Frame_; f++)
        {
            // Update population (randomize, move, etc.)
            update_pop(Pop_local);

            for(int i = 0; i < nGamma; ++i)
            {
                double resol = Gamma[i];
                
                // Reset field
                E2w_local[i] = Electric_Field();
                
                // Calculate field contribution
                calc_field(Pop_local, Ew, E2w_local, i, resol);

                SommeV_local[i] += norm(E2w_local[i].Comp()[0]);
                SommeH_local[i] += norm((-E2w_local[i].Comp()[1]*cos(Ang_Col_u*M_PI/180.)) + (E2w_local[i].Comp()[2]*sin(Ang_Col_u*M_PI/180.)));
            }
            
            // Atomic progress update
            int local_completed;
            #pragma omp atomic capture
            local_completed = ++completed_frames;
            
            if (local_completed % 100 == 0 || local_completed == Frame_) {
                #pragma omp critical
                {
                    cout<<"Avancement : "<<setprecision (3)<<setw(3)<<100.*local_completed/Frame_<<"%\t\t\r"<<flush;
                }
            }
        }
        
        // Accumulate results
        #pragma omp critical
        {
            for(int i=0; i<nGamma; ++i) {
                SommeV[i] += SommeV_local[i];
                SommeH[i] += SommeH_local[i];
            }
        }
    }
    
}


void Setup::PolarPattern(Electric_Field& Ew, Population& Pop, double dGamma, string path)
{
    RunParallelSimulation<Population, Electric_Field>(
        Pop, Ew, dGamma, 0.0, path,
        // Update Function
        [](Population& p) {
            p.Randomize_Orientation();
            p.Place_Element_in_Lab_Frame();
        },
        // Calculate Function
        [this](Population& p, Electric_Field& ew, vector<Electric_Field>& e2w, int idx, double angle) {
            Electric_Field Ew_rot = ew;
            Ew_rot.Rotate_Field(angle);

            if(Treat_eme()) {
                Electric_Field* ptr = Full_p_developement(p, Ew_rot);
                e2w[idx] += *ptr;
                delete ptr;
            } else {
                Electric_Field* ptr = Get_Field_Amplitude(p, Ew_rot);
                e2w[idx] += *ptr;
                delete ptr;
            }
            
            if(Treat_MEE()) {
                Electric_Field* ptr = Get_M_Contribution(p, Ew_rot);
                e2w[idx] -= *ptr;
                delete ptr;
            }
            if(Treat_QEE()) {
                Electric_Field* ptr = Get_RetardationContrib(p, Ew_rot);
                e2w[idx] -= *ptr;
                delete ptr;
            }
        }
    );
    
    if(Check_save_config()) Pop.SaveConfig(path+"_Dip_Config.xyz");
    ecrire(Gamma,SommeV,SommeH,path+"Main_Polar.txt");
    Gamma.clear();
    SommeH.clear();
    SommeV.clear();
}

void Setup::PolarPattern(Eliptic_Electric_Field& Ew, Population& Pop, double dGamma, string path)
{
    RunParallelSimulation<Population, Eliptic_Electric_Field>(
        Pop, Ew, dGamma, 45.0, path,
        // Update Function
        [](Population& p) {
            p.Randomize_Orientation();
            p.Place_Element_in_Lab_Frame();
        },
        // Calculate Function
        [this](Population& p, Eliptic_Electric_Field& ew, vector<Electric_Field>& e2w, int idx, double angle) {
            Eliptic_Electric_Field Ew_rot = ew;
            Ew_rot.Rotate_Field(0., angle);

            if(Treat_eme()) {
                Electric_Field* ptr = Full_p_developement(p, Ew_rot);
                e2w[idx] += *ptr;
                delete ptr;
            } else {
                Electric_Field* ptr = Get_Field_Amplitude(p, Ew_rot);
                e2w[idx] += *ptr;
                delete ptr;
            }
            
            if(Treat_MEE()) {
                Electric_Field* ptr = Get_M_Contribution(p, Ew_rot);
                e2w[idx] -= *ptr;
                delete ptr;
            }
        }
    );
    
    if(Check_save_config()) Pop.SaveConfig(path+"_Dip_Config.xyz");
    ecrire(Gamma,SommeV,SommeH,path+"Main_Polar.txt");
    Gamma.clear();
    SommeH.clear();
    SommeV.clear();
}

void Setup::PolarPattern(Electric_Field& Ew, vector<Population> Pop, double dGamma, string path)
{
    RunParallelSimulation<vector<Population>, Electric_Field>(
        Pop, Ew, dGamma, 0.0, path,
        // Update Function
        [](vector<Population>& pops) {
            for(unsigned int cpt=0; cpt<pops.size(); cpt++) {
                pops[cpt].Move_Pop(enoise(-250.,250), enoise(-250.,250), enoise(-250.,250));
                pops[cpt].Place_Element_in_Lab_Frame();
            }
        },
        // Calculate Function
        [this](vector<Population>& pops, Electric_Field& ew, vector<Electric_Field>& e2w, int idx, double angle) {
            Electric_Field Ew_rot = ew;
            Ew_rot.Rotate_Field(angle);
            
            for(unsigned int cpt=0; cpt<pops.size(); cpt++) {
                Electric_Field* ptr = Get_Field_Amplitude(pops[cpt], Ew_rot);
                e2w[idx] += *ptr;
                delete ptr;
            }
        }
    );
    
    if(Check_save_config()) Save_Config(Pop, path+"_Dip_Config.xyz");
    ecrire(Gamma,SommeV,SommeH,path+"_Main_Polar.txt");
    Gamma.clear();
    SommeH.clear();
    SommeV.clear();
}

    void Setup::ecrire(vector<double> t1, vector<double> t2, vector<double> t3, string name) const
    {
        ofstream fichier(name.c_str(), ios::out |ios::trunc);
        if (!fichier.is_open()) {
            cerr << "Error: Failed to open file " << name << endl;
            return;
        }
        fichier<<"PolAngle\tPolarV\tPolarH"<<endl;
        for(unsigned int i=0;i<t1.size();i++)
            fichier << t1[i]<<"\t"<< t2[i]<<"\t"<<t3[i] <<endl;
        fichier.close();
    }


    void Dipole::Save_Beta(string file)
    {
        ofstream fichier(file.c_str(),ios::out|ios::trunc);
        fichier<<"Beta EEE"<<endl<<Beta_Dip()<<endl<<endl<<endl;
        fichier<<"Beta EEM"<<endl<<Beta_EME()<<endl<<endl<<endl;
        fichier<<"Beta MEE"<<endl<<Beta_MEE()<<endl<<endl<<endl;
        fichier<<"Beta QEE"<<endl<<Beta_QEE(0)<<endl<<endl<<Beta_QEE(1)<<endl<<endl<<Beta_QEE(2);

    }

// Performance optimized version of PolarPattern
void Setup::PolarPattern_Optimized(Electric_Field& Ew, Population& Pop, double dGamma, string path)
{
    const int nGamma = int(1 + 360.0/dGamma);
    const int nFrames = static_cast<int>(Frame_);
    
    // Preallocate arrays
    SommeV.assign(nGamma, 0.0);
    SommeH.assign(nGamma, 0.0);
    
    // Precompute angles and trig functions
    vector<double> gamma_angles(nGamma);
    const double cos_col = cos(Ang_Col_u * M_PI / 180.0);
    const double sin_col = sin(Ang_Col_u * M_PI / 180.0);
    
    for(int i = 0; i < nGamma; ++i) {
        gamma_angles[i] = i * dGamma;
    }
    
    if(Gamma.empty()) {
        Gamma = gamma_angles; // Only set once
    }
    
    auto start_time = chrono::high_resolution_clock::now();
    
    // Standard sequential version with optimizations
    for(int f = 0; f < nFrames; ++f) {
        Pop.Randomize_Orientation();
        Pop.Place_Element_in_Lab_Frame();
        
        // Process all gamma angles for this frame
        for(int g = 0; g < nGamma; ++g) {
            // Set field orientation
            Electric_Field local_Ew = Ew;
            local_Ew.Rotate_Field(gamma_angles[g]);
            
            // Compute field contributions
            Electric_Field total_field;
            
            if(Treat_eme()) {
                Electric_Field* contrib = Full_p_developement(Pop, local_Ew);
                total_field += *contrib;
                delete contrib;
            } else {
                Electric_Field* contrib = Get_Field_Amplitude(Pop, local_Ew);
                total_field += *contrib;
                delete contrib;
            }
            
            if(Treat_MEE()) {
                Electric_Field* contrib = Get_M_Contribution(Pop, local_Ew);
                total_field -= *contrib;
                delete contrib;
            }
            
            if(Treat_QEE()) {
                Electric_Field* contrib = Get_RetardationContrib(Pop, local_Ew);
                total_field -= *contrib;
                delete contrib;
            }
            
            // Accumulate results
            const auto& comp = total_field.Comp();
            SommeV[g] += norm(comp[0]);
            
            const complex<double> h_field = -comp[1] * cos_col + comp[2] * sin_col;
            SommeH[g] += norm(h_field);
        }
        
        // Progress reporting
        if(f % 10 == 0) {
            cout << "Progress: " << setprecision(3) << setw(6) 
                 << 100.0 * (f + 1) / nFrames << "%\r" << flush;
        }
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    
    cout << "\nComputation completed in " << duration.count() << " ms" << endl;
    
    // Save configuration if requested
    if(Check_save_config()) {
        Pop.SaveConfig(path + "_Dip_Config.xyz");
    }
    
    // Write results with metadata
    ecrire_with_metadata(Gamma, SommeV, SommeH, path + "Main_Polar.txt", 
                        nFrames, dGamma, duration.count());
}

// Enhanced output function with metadata
void Setup::ecrire_with_metadata(const vector<double>& t1, const vector<double>& t2, 
                                 const vector<double>& t3, string name, 
                                 int frames, double dGamma, long duration_ms) const
{
    ofstream fichier(name.c_str(), ios::out | ios::trunc);
    if(!fichier) {
        cerr << "Error: Cannot open file " << name << endl;
        return;
    }
    
    // Write metadata header
    fichier << "# HRS-f Polar Pattern Data\n";
    fichier << "# Generated: " << __DATE__ << " " << __TIME__ << "\n";
    fichier << "# Frames: " << frames << "\n";
    fichier << "# dGamma: " << dGamma << "\n";
    fichier << "# Computation time: " << duration_ms << " ms\n";
    fichier << "# Collection angles: u=" << Ang_Col_u << ", v=" << Ang_Col_v << "\n";
    fichier << "# Physics flags: EME=" << Treat_eme() << ", MEE=" << Treat_MEE() 
            << ", QEE=" << Treat_QEE() << "\n";
    fichier << "# Columns: Gamma(deg), SumV, SumH\n";
    
    fichier << scientific << setprecision(6);
    for(size_t i = 0; i < t1.size(); ++i) {
        fichier << t1[i] << "\t" << t2[i] << "\t" << t3[i] << "\n";
    }
    fichier.close();
}

