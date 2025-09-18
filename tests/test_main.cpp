#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

// Simple Catch2-like testing framework (basic implementation)
// For production, use actual Catch2 library

#include "../include/class.h"
#include <cmath>

TEST_CASE("Position class tests", "[position]") {
    SECTION("Constructor and getters") {
        Position p(1.0, 2.0, 3.0);
        REQUIRE(p.x() == Approx(1.0));
        REQUIRE(p.y() == Approx(2.0));
        REQUIRE(p.z() == Approx(3.0));
    }
    
    SECTION("Norm calculation") {
        Position p(3.0, 4.0, 0.0);
        REQUIRE(p.norme() == Approx(5.0));
    }
    
    SECTION("Translation") {
        Position p(1.0, 1.0, 1.0);
        p.Deplacer_de_(1.0, 2.0, 3.0);
        REQUIRE(p.x() == Approx(2.0));
        REQUIRE(p.y() == Approx(3.0));
        REQUIRE(p.z() == Approx(4.0));
    }
}

TEST_CASE("Euler_Angles class tests", "[euler]") {
    SECTION("Constructor and getters") {
        Euler_Angles angles(30.0, 60.0, 90.0);
        REQUIRE(angles.Psi() == Approx(30.0));
        REQUIRE(angles.Theta() == Approx(60.0));
        REQUIRE(angles.Phi() == Approx(90.0));
    }
    
    SECTION("Rotation matrix orthogonality") {
        Euler_Angles angles(45.0, 30.0, 60.0);
        auto matrix = angles.Matrice_Passage();
        
        // Check that it's 3x3
        REQUIRE(matrix.size() == 3);
        REQUIRE(matrix[0].size() == 3);
        
        // Check orthogonality: M * M^T should be identity
        vector<vector<double>> identity(3, vector<double>(3, 0.0));
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                for(int k = 0; k < 3; ++k) {
                    identity[i][j] += matrix[i][k] * matrix[j][k];
                }
            }
        }
        
        // Check diagonal elements are 1
        for(int i = 0; i < 3; ++i) {
            REQUIRE(identity[i][i] == Approx(1.0).epsilon(1e-10));
        }
        
        // Check off-diagonal elements are 0
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                if(i != j) {
                    REQUIRE(identity[i][j] == Approx(0.0).epsilon(1e-10));
                }
            }
        }
    }
}

TEST_CASE("Vector operations tests", "[vector_ops]") {
    SECTION("Vector cross product") {
        vector<complex<double>> a = {complex<double>(1,0), complex<double>(0,0), complex<double>(0,0)};
        vector<complex<double>> b = {complex<double>(0,0), complex<double>(1,0), complex<double>(0,0)};
        
        auto cross = VV_Cross_Product(a, b);
        
        REQUIRE(cross.size() == 3);
        REQUIRE(real(cross[0]) == Approx(0.0));
        REQUIRE(real(cross[1]) == Approx(0.0));
        REQUIRE(real(cross[2]) == Approx(1.0));
    }
    
    SECTION("Vector scalar product") {
        vector<complex<double>> a = {complex<double>(1,0), complex<double>(2,0), complex<double>(3,0)};
        Position b(2.0, 3.0, 4.0);
        
        auto result = VV_Scal_Prod(a, b);
        
        // Expected: 1*2 + 2*3 + 3*4 = 2 + 6 + 12 = 20
        REQUIRE(real(result) == Approx(20.0));
        REQUIRE(imag(result) == Approx(0.0));
    }
    
    SECTION("Vector norm") {
        vector<double> v = {3.0, 4.0, 0.0};
        REQUIRE(norm(v) == Approx(5.0));
    }
}

TEST_CASE("Electric_Field class tests", "[electric_field]") {
    SECTION("Constructor with polarization angle") {
        Electric_Field field(45.0, complex<double>(1.0, 0.0), 800.0);
        
        auto comp = field.Comp();
        REQUIRE(comp.size() == 3);
        
        // At 45 degrees, x and y components should be equal
        REQUIRE(real(comp[0]) == Approx(real(comp[1])).epsilon(1e-10));
        REQUIRE(real(comp[2]) == Approx(0.0));
        
        // Magnitude should be preserved
        double mag = sqrt(norm(comp[0]) + norm(comp[1]) + norm(comp[2]));
        REQUIRE(mag == Approx(1.0).epsilon(1e-10));
    }
    
    SECTION("Field rotation") {
        Electric_Field field(0.0, complex<double>(1.0, 0.0), 800.0);
        
        // Initially polarized in x
        auto comp_initial = field.Comp();
        REQUIRE(real(comp_initial[0]) == Approx(1.0));
        REQUIRE(real(comp_initial[1]) == Approx(0.0));
        
        // Rotate by 90 degrees
        field.Rotate_Field(90.0);
        auto comp_rotated = field.Comp();
        
        // Should now be polarized in y
        REQUIRE(real(comp_rotated[0]) == Approx(0.0).epsilon(1e-10));
        REQUIRE(real(comp_rotated[1]) == Approx(1.0).epsilon(1e-10));
    }
}

TEST_CASE("Setup class tests", "[setup]") {
    SECTION("Constructor and initialization") {
        Setup setup(100, 45.0, 60.0);
        
        REQUIRE(setup.get_nFrame() == 100);
        
        // Check that boolean flags are properly initialized
        REQUIRE(setup.Treat_eme() == false);
        REQUIRE(setup.Treat_MEE() == false);
        REQUIRE(setup.Treat_QEE() == false);
    }
    
    SECTION("Physics flags") {
        Setup setup;
        
        setup.Select_eme(true);
        setup.Select_MEE(true);
        setup.Select_QEE(false);
        
        REQUIRE(setup.Treat_eme() == true);
        REQUIRE(setup.Treat_MEE() == true);
        REQUIRE(setup.Treat_QEE() == false);
    }
}

TEST_CASE("Population class tests", "[population]") {
    SECTION("Adding dipoles") {
        Population pop;
        
        pop.Add_Dipole(0, 0, 0, "test", 1.0, 2.0, 3.0);
        pop.Add_Dipole(45, 30, 60, "test2", 4.0, 5.0, 6.0);
        
        REQUIRE(pop.Get_Nb_Dip() == 2);
        
        auto* dipole1 = pop.Get_Dip(0);
        REQUIRE(dipole1 != nullptr);
        REQUIRE(dipole1->getPosition().x() == Approx(1.0));
        REQUIRE(dipole1->getPosition().y() == Approx(2.0));
        REQUIRE(dipole1->getPosition().z() == Approx(3.0));
    }
}

// Physics validation tests
TEST_CASE("Physics validation tests", "[physics][validation]") {
    SECTION("Single dipole Rayleigh scattering") {
        // Test that a single dipole produces the expected angular pattern
        Population pop;
        pop.Add_Dipole(0, 0, 0, "rayleigh", 0.0, 0.0, 0.0);
        
        Electric_Field field(0.0, complex<double>(1.33, 0.0), 800.0);
        Setup setup(1, 90.0, 90.0); // Single frame, perpendicular observation
        
        // This would require implementing analytical comparison
        // For now, just check that computation doesn't crash
        auto* result = setup.Get_Field_Amplitude(pop, field);
        REQUIRE(result != nullptr);
        
        auto comp = result->Comp();
        REQUIRE(comp.size() == 3);
        
        delete result;
    }
    
    SECTION("Field amplitude conservation") {
        // Test that total field amplitude behaves as expected
        Population pop;
        
        // Add dipoles in a symmetric configuration
        const int N = 4;
        const double radius = 10.0;
        for(int i = 0; i < N; ++i) {
            double angle = 2.0 * M_PI * i / N;
            pop.Add_Dipole(0, 0, 0, "sym", 
                          radius * cos(angle), 
                          radius * sin(angle), 
                          0.0);
        }
        
        Electric_Field field(0.0, complex<double>(1.33, 0.0), 800.0);
        Setup setup(1, 90.0, 90.0);
        
        auto* result = setup.Get_Field_Amplitude(pop, field);
        REQUIRE(result != nullptr);
        
        // Check that result has reasonable magnitude
        auto comp = result->Comp();
        double total_intensity = norm(comp[0]) + norm(comp[1]) + norm(comp[2]);
        REQUIRE(total_intensity > 0.0);
        
        delete result;
    }
}


