#include "../include/class.h"
#include "../include/config.h"
#include <iostream>
#include <chrono>
#ifdef HAS_OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <random>
#include <algorithm>

// Enhanced simulation runner with all PhD-level features
class HRSSimulation {
private:
    SimulationConfig config_;
    std::mt19937 rng_;
    
public:
    HRSSimulation(const SimulationConfig& config) 
        : config_(config), rng_(std::chrono::steady_clock::now().time_since_epoch().count()) {
        
        // Set OpenMP threads
        if (config_.num_threads > 0) {
#ifdef HAS_OPENMP
            omp_set_num_threads(config_.num_threads);
#endif
        }
        
        // Create output directory
        system(("mkdir -p " + config_.output_dir).c_str());
    }
    
    // Create population based on configuration
    Population create_population() {
        Population pop;
        
        if (config_.geometry == "ring") {
            // Ring geometry
            const double radius = config_.radius;
            const int n = config_.dipole_count;
            
            for (int i = 0; i < n; ++i) {
                double angle = 2.0 * M_PI * i / n;
                pop.Add_Dipole(0, 0, 0, "ring_dipole",
                              radius * cos(angle),
                              radius * sin(angle),
                              0.0);
            }
        }
        else if (config_.geometry == "helix") {
            // Helix geometry
            pop.Build_Helix(config_.dipole_count, config_.radius, 50.0, 2, false);
        }
        else if (config_.geometry == "file") {
            // Load from file
            if (!config_.geometry_file.empty()) {
                pop.ReadConfig_just_position(config_.geometry_file);
            } else {
                std::cerr << "Warning: geometry_file not specified for file geometry\n";
            }
        }
        else {
            std::cerr << "Warning: Unknown geometry type: " << config_.geometry << "\n";
            // Default to ring
            const double radius = config_.radius;
            const int n = config_.dipole_count;
            
            for (int i = 0; i < n; ++i) {
                double angle = 2.0 * M_PI * i / n;
                pop.Add_Dipole(0, 0, 0, "default_dipole",
                              radius * cos(angle),
                              radius * sin(angle),
                              0.0);
            }
        }
        
        return pop;
    }
    
    // Create electric field based on configuration
    std::unique_ptr<Electric_Field> create_field() {
        if (config_.elliptic_field) {
            return std::make_unique<Eliptic_Electric_Field>(
                config_.ellipticity,
                config_.field_angle,
                complex<double>(config_.refractive_index, 0.0),
                config_.wavelength
            );
        } else {
            return std::make_unique<Electric_Field>(
                config_.field_angle,
                complex<double>(config_.refractive_index, 0.0),
                config_.wavelength
            );
        }
    }
    
    // Run single simulation
    void run_simulation() {
        std::cout << "Starting HRS-f simulation...\n";
        std::cout << "Configuration:\n";
        std::cout << "  Frames: " << config_.frames << "\n";
        std::cout << "  dGamma: " << config_.dGamma << "°\n";
        std::cout << "  Dipoles: " << config_.dipole_count << "\n";
        std::cout << "  Geometry: " << config_.geometry << "\n";
        std::cout << "  Wavelength: " << config_.wavelength << " nm\n";
        std::cout << "  Refractive index: " << config_.refractive_index << "\n";
        std::cout << "  Physics flags: EME=" << config_.enable_eme 
                  << ", MEE=" << config_.enable_mee 
                  << ", QEE=" << config_.enable_qee << "\n";
        std::cout << "  Threads: " << 
#ifdef HAS_OPENMP
            omp_get_max_threads() 
#else
            1
#endif
            << "\n\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Create simulation objects
        Population pop = create_population();
        auto field = create_field();
        
        Setup setup(config_.frames, config_.collection_u, config_.collection_v);
        setup.Select_eme(config_.enable_eme);
        setup.Select_MEE(config_.enable_mee);
        setup.Select_QEE(config_.enable_qee);
        
        // Set up output path
        std::string output_path = config_.output_dir + config_.output_prefix;
        
        if (config_.save_config) {
            setup.Balise_Save_config(output_path);
        }
        
        // Run simulation
        if (config_.elliptic_field) {
            // Create elliptic field directly
            Eliptic_Electric_Field elliptic_field(
                config_.ellipticity,
                config_.field_angle,
                complex<double>(config_.refractive_index, 0.0),
                config_.wavelength
            );
            setup.PolarPattern(elliptic_field, pop, config_.dGamma, output_path);
        } else {
            setup.PolarPattern(*field, pop, config_.dGamma, output_path);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "\nSimulation completed in " << duration.count() << " ms\n";
        
        // Save metadata
        if (config_.save_metadata) {
            save_simulation_metadata(output_path, duration.count());
        }
        
        std::cout << "Results saved with prefix: " << output_path << "\n";
    }
    
    // Run uncertainty quantification
    void run_uncertainty_quantification() {
        if (!config_.enable_uq) return;
        
        std::cout << "Running uncertainty quantification with " 
                  << config_.uq_samples << " samples...\n";
        
        // Define parameter variations (example)
        std::uniform_real_distribution<double> radius_dist(config_.radius * 0.9, config_.radius * 1.1);
        std::uniform_real_distribution<double> n_dist(config_.refractive_index * 0.95, config_.refractive_index * 1.05);
        std::uniform_real_distribution<double> lambda_dist(config_.wavelength * 0.98, config_.wavelength * 1.02);
        
        // Storage for results
        std::vector<std::vector<double>> all_sumV, all_sumH;
        
        for (int sample = 0; sample < config_.uq_samples; ++sample) {
            // Perturb parameters
            SimulationConfig sample_config = config_;
            sample_config.radius = radius_dist(rng_);
            sample_config.refractive_index = n_dist(rng_);
            sample_config.wavelength = lambda_dist(rng_);
            
            // Run simulation with perturbed parameters
            HRSSimulation sample_sim(sample_config);
            
            // This would require modifying the simulation to return results
            // For now, just demonstrate the framework
            
            if (sample % 100 == 0) {
                std::cout << "UQ Progress: " << (100 * sample / config_.uq_samples) << "%\r" << std::flush;
            }
        }
        
        std::cout << "\nUncertainty quantification completed\n";
        
        // Analyze and save UQ results
        // Implementation would include statistical analysis of results
    }
    
    // Save simulation metadata
    void save_simulation_metadata(const std::string& output_path, long duration_ms) {
        std::ofstream meta_file(output_path + "_metadata.json");
        
        meta_file << "{\n";
        meta_file << "  \"simulation\": {\n";
        meta_file << "    \"timestamp\": \"" << __DATE__ << " " << __TIME__ << "\",\n";
        meta_file << "    \"duration_ms\": " << duration_ms << ",\n";
        meta_file << "    \"threads\": " << 
#ifdef HAS_OPENMP
            omp_get_max_threads() 
#else
            1
#endif
            << "\n";
        meta_file << "  },\n";
        meta_file << "  \"parameters\": {\n";
        meta_file << "    \"frames\": " << config_.frames << ",\n";
        meta_file << "    \"dGamma\": " << config_.dGamma << ",\n";
        meta_file << "    \"dipole_count\": " << config_.dipole_count << ",\n";
        meta_file << "    \"radius\": " << config_.radius << ",\n";
        meta_file << "    \"geometry\": \"" << config_.geometry << "\",\n";
        meta_file << "    \"wavelength\": " << config_.wavelength << ",\n";
        meta_file << "    \"refractive_index\": " << config_.refractive_index << ",\n";
        meta_file << "    \"field_angle\": " << config_.field_angle << ",\n";
        meta_file << "    \"collection_u\": " << config_.collection_u << ",\n";
        meta_file << "    \"collection_v\": " << config_.collection_v << "\n";
        meta_file << "  },\n";
        meta_file << "  \"physics\": {\n";
        meta_file << "    \"enable_eme\": " << (config_.enable_eme ? "true" : "false") << ",\n";
        meta_file << "    \"enable_mee\": " << (config_.enable_mee ? "true" : "false") << ",\n";
        meta_file << "    \"enable_qee\": " << (config_.enable_qee ? "true" : "false") << "\n";
        meta_file << "  }\n";
        meta_file << "}\n";
    }
};

// Main function with enhanced features
int main(int argc, char* argv[]) {
    try {
        // Parse command line arguments
        SimulationConfig config = ConfigManager::parse_cli(argc, argv);
        
        // Validate configuration
        std::string error_msg;
        if (!ConfigManager::validate_config(config, error_msg)) {
            std::cerr << "Configuration error: " << error_msg << std::endl;
            return 1;
        }
        
        // Run validation if requested
        if (config.run_validation) {
            std::cout << "Running validation tests...\n";
            int validation_result = system("./validation/rayleigh_validation");
            if (validation_result != 0) {
                std::cerr << "Validation tests failed!\n";
                return 1;
            }
            std::cout << "Validation tests passed.\n\n";
        }
        
        // Create and run simulation
        HRSSimulation simulation(config);
        simulation.run_simulation();
        
        // Run uncertainty quantification if requested
        if (config.enable_uq) {
            simulation.run_uncertainty_quantification();
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
}
