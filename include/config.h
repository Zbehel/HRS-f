#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <vector>
#include <map>

struct SimulationConfig {
    // Basic parameters
    int frames = 100;
    double dGamma = 5.0;
    double collection_u = 90.0;
    double collection_v = 90.0;
    
    // Field parameters
    double field_angle = 0.0;
    double refractive_index = 1.33;
    double wavelength = 800.0;
    bool elliptic_field = false;
    double ellipticity = 0.0;
    
    // Population parameters
    int dipole_count = 12;
    double radius = 10.0;
    std::string geometry = "ring"; // ring, helix, sphere, file
    std::string geometry_file = "";
    
    // Physics flags
    bool enable_eme = false;
    bool enable_mee = false;
    bool enable_qee = false;
    
    // Output parameters
    std::string output_prefix = "simulation";
    std::string output_dir = "output/data/";
    bool save_config = true;
    bool save_metadata = true;
    
    // Performance parameters
    int num_threads = 0; // 0 = auto-detect
    bool use_optimization = true;
    
    // Validation parameters
    bool run_validation = false;
    std::string validation_target = "rayleigh_single";
    double validation_tolerance = 1e-6;
    
    // Uncertainty quantification
    bool enable_uq = false;
    int uq_samples = 1000;
    std::vector<std::string> uq_parameters;
    std::vector<double> uq_ranges;
};

struct ExperimentConfig {
    std::string name;
    std::string description;
    SimulationConfig base_config;
    
    // Parameter sweep definitions
    std::map<std::string, std::vector<double>> parameter_sweeps;
    
    // Batch processing
    bool parallel_experiments = true;
    int max_concurrent = 4;
};

// Configuration management
class ConfigManager {
public:
    static SimulationConfig load_yaml(const std::string& filename);
    static void save_yaml(const SimulationConfig& config, const std::string& filename);
    static ExperimentConfig load_experiment_yaml(const std::string& filename);
    static void save_experiment_yaml(const ExperimentConfig& config, const std::string& filename);
    
    // CLI parsing
    static SimulationConfig parse_cli(int argc, char* argv[]);
    static void print_help();
    
    // Validation
    static bool validate_config(const SimulationConfig& config, std::string& error_msg);
    
    // Template generation
    static void generate_template_config(const std::string& filename);
    static void generate_experiment_template(const std::string& filename);
};

#endif // CONFIG_H_INCLUDED
