#include "../../include/class.h"
#include "../../include/config.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

class RayleighValidator {
public:
    struct ValidationResult {
        bool passed;
        double max_error;
        double rms_error;
        std::string error_message;
        std::vector<double> angles;
        std::vector<double> analytical;
        std::vector<double> numerical;
    };
    
    // Analytical Rayleigh scattering for a single dipole
    static double rayleigh_intensity(double theta_deg, double phi_deg, 
                                   double observation_theta_deg) {
        // Convert to radians
        double theta = theta_deg * M_PI / 180.0;
        double obs_theta = observation_theta_deg * M_PI / 180.0;
        
        // For perpendicular observation (obs_theta = 90°), 
        // intensity ∝ sin²(θ) where θ is angle from dipole axis
        double cos_angle = cos(theta);
        double sin_angle = sin(theta);
        
        // Rayleigh scattering intensity pattern
        // I(θ) = I₀ * sin²(θ) for observation perpendicular to dipole
        return sin_angle * sin_angle;
    }
    
    // Validate single dipole against analytical Rayleigh pattern
    static ValidationResult validate_single_dipole(double tolerance = 1e-3) {
        ValidationResult result;
        result.passed = false;
        result.max_error = 0.0;
        result.rms_error = 0.0;
        
        try {
            // Create single dipole at origin
            Population pop;
            pop.Add_Dipole(0, 0, 0, "rayleigh", 0.0, 0.0, 0.0);
            
            // Create field (x-polarized)
            Electric_Field field(0.0, complex<double>(1.33, 0.0), 800.0);
            
            // Setup with single frame, perpendicular observation
            Setup setup(1, 90.0, 90.0);
            setup.Select_eme(false);
            setup.Select_MEE(false);
            setup.Select_QEE(false);
            
            // Test various field orientations
            const double dGamma = 5.0;
            const int nAngles = int(360.0 / dGamma) + 1;
            
            result.angles.reserve(nAngles);
            result.analytical.reserve(nAngles);
            result.numerical.reserve(nAngles);
            
            double sum_sq_error = 0.0;
            
            for (int i = 0; i < nAngles; ++i) {
                double gamma = i * dGamma;
                if (gamma > 360.0) break;
                
                // Set field orientation
                Electric_Field test_field = field;
                test_field.Rotate_Field(gamma);
                
                // Get numerical result
                auto* field_result = setup.Get_Field_Amplitude(pop, test_field);
                if (!field_result) {
                    result.error_message = "Failed to compute field amplitude";
                    return result;
                }
                
                auto comp = field_result->Comp();
                double numerical_intensity = norm(comp[0]) + norm(comp[1]) + norm(comp[2]);
                
                // Get analytical result
                double analytical_intensity = rayleigh_intensity(gamma, 0.0, 90.0);
                
                // Normalize (since we're comparing patterns, not absolute intensities)
                if (i == 0) {
                    // Store normalization factors
                    if (numerical_intensity > 0 && analytical_intensity > 0) {
                        // Continue with normalized comparison
                    }
                }
                
                result.angles.push_back(gamma);
                result.analytical.push_back(analytical_intensity);
                result.numerical.push_back(numerical_intensity);
                
                // Calculate error (relative)
                double error = 0.0;
                if (analytical_intensity > 1e-10) {
                    error = std::abs(numerical_intensity - analytical_intensity) / analytical_intensity;
                } else {
                    error = std::abs(numerical_intensity - analytical_intensity);
                }
                
                result.max_error = std::max(result.max_error, error);
                sum_sq_error += error * error;
                
                delete field_result;
            }
            
            result.rms_error = sqrt(sum_sq_error / nAngles);
            result.passed = (result.max_error < tolerance);
            
            if (!result.passed) {
                result.error_message = "Maximum error " + std::to_string(result.max_error) + 
                                     " exceeds tolerance " + std::to_string(tolerance);
            }
            
        } catch (const std::exception& e) {
            result.error_message = "Exception during validation: " + std::string(e.what());
        }
        
        return result;
    }
    
    // Validate dipole ring symmetry
    static ValidationResult validate_ring_symmetry(int n_dipoles = 8, double tolerance = 1e-6) {
        ValidationResult result;
        result.passed = false;
        
        try {
            // Create symmetric ring of dipoles
            Population pop;
            const double radius = 10.0;
            
            for (int i = 0; i < n_dipoles; ++i) {
                double angle = 2.0 * M_PI * i / n_dipoles;
                pop.Add_Dipole(0, 0, 0, "ring", 
                              radius * cos(angle), 
                              radius * sin(angle), 
                              0.0);
            }
            
            Electric_Field field(0.0, complex<double>(1.33, 0.0), 800.0);
            Setup setup(1, 90.0, 90.0);
            
            // Test rotational symmetry
            const double symmetry_angle = 360.0 / n_dipoles;
            const int n_tests = 10;
            
            std::vector<double> intensities;
            
            for (int test = 0; test < n_tests; ++test) {
                double gamma = test * symmetry_angle;
                
                Electric_Field test_field = field;
                test_field.Rotate_Field(gamma);
                
                auto* field_result = setup.Get_Field_Amplitude(pop, test_field);
                if (!field_result) {
                    result.error_message = "Failed to compute field amplitude";
                    return result;
                }
                
                auto comp = field_result->Comp();
                double intensity = norm(comp[0]) + norm(comp[1]) + norm(comp[2]);
                intensities.push_back(intensity);
                
                delete field_result;
            }
            
            // Check that all intensities are approximately equal
            double mean_intensity = 0.0;
            for (double I : intensities) {
                mean_intensity += I;
            }
            mean_intensity /= intensities.size();
            
            double max_deviation = 0.0;
            for (double I : intensities) {
                double deviation = std::abs(I - mean_intensity) / mean_intensity;
                max_deviation = std::max(max_deviation, deviation);
            }
            
            result.max_error = max_deviation;
            result.passed = (max_deviation < tolerance);
            
            if (!result.passed) {
                result.error_message = "Ring symmetry broken: max deviation " + 
                                     std::to_string(max_deviation) + 
                                     " exceeds tolerance " + std::to_string(tolerance);
            }
            
        } catch (const std::exception& e) {
            result.error_message = "Exception during ring symmetry validation: " + std::string(e.what());
        }
        
        return result;
    }
    
    // Save validation results to file
    static void save_validation_report(const std::vector<ValidationResult>& results, 
                                     const std::string& filename) {
        std::ofstream file(filename);
        
        file << "# HRS-f Validation Report\n";
        file << "# Generated: " << __DATE__ << " " << __TIME__ << "\n\n";
        
        int passed = 0;
        int total = results.size();
        
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& result = results[i];
            
            file << "Test " << (i+1) << ": ";
            if (result.passed) {
                file << "PASSED";
                ++passed;
            } else {
                file << "FAILED";
            }
            file << "\n";
            
            file << "  Max Error: " << result.max_error << "\n";
            file << "  RMS Error: " << result.rms_error << "\n";
            
            if (!result.error_message.empty()) {
                file << "  Message: " << result.error_message << "\n";
            }
            
            // Save detailed data if available
            if (!result.angles.empty() && !result.analytical.empty() && !result.numerical.empty()) {
                file << "  Data points: " << result.angles.size() << "\n";
                file << "  # Angle(deg) Analytical Numerical\n";
                for (size_t j = 0; j < result.angles.size(); ++j) {
                    file << "  " << result.angles[j] << " " 
                         << result.analytical[j] << " " 
                         << result.numerical[j] << "\n";
                }
            }
            file << "\n";
        }
        
        file << "Summary: " << passed << "/" << total << " tests passed\n";
        
        if (passed == total) {
            file << "All validation tests PASSED!\n";
        } else {
            file << "Some validation tests FAILED!\n";
        }
    }
};

// Main validation runner
int main() {
    std::cout << "Running HRS-f validation tests...\n\n";
    
    std::vector<RayleighValidator::ValidationResult> results;
    
    // Test 1: Single dipole Rayleigh scattering
    std::cout << "Test 1: Single dipole Rayleigh scattering pattern... ";
    auto result1 = RayleighValidator::validate_single_dipole(0.1); // 10% tolerance
    results.push_back(result1);
    
    if (result1.passed) {
        std::cout << "PASSED (max error: " << result1.max_error << ")\n";
    } else {
        std::cout << "FAILED: " << result1.error_message << "\n";
    }
    
    // Test 2: Ring symmetry
    std::cout << "Test 2: Ring symmetry validation... ";
    auto result2 = RayleighValidator::validate_ring_symmetry(8, 0.01); // 1% tolerance
    results.push_back(result2);
    
    if (result2.passed) {
        std::cout << "PASSED (max error: " << result2.max_error << ")\n";
    } else {
        std::cout << "FAILED: " << result2.error_message << "\n";
    }
    
    // Save detailed report
    RayleighValidator::save_validation_report(results, "validation/validation_report.txt");
    
    std::cout << "\nValidation report saved to validation/validation_report.txt\n";
    
    // Return non-zero if any test failed
    for (const auto& result : results) {
        if (!result.passed) {
            return 1;
        }
    }
    
    std::cout << "All validation tests passed!\n";
    return 0;
}
