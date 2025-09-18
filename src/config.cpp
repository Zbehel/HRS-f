#include "../include/config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstring>

// Simple YAML-like parser (basic implementation)
// For production use, consider yaml-cpp library

SimulationConfig ConfigManager::load_yaml(const std::string& filename) {
    SimulationConfig config;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open config file " << filename 
                  << ". Using defaults." << std::endl;
        return config;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        // Simple key: value parsing
        size_t colon = line.find(':');
        if (colon == std::string::npos) continue;
        
        std::string key = line.substr(0, colon);
        std::string value = line.substr(colon + 1);
        
        // Trim whitespace
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        
        // Parse values
        if (key == "frames") config.frames = std::stoi(value);
        else if (key == "dGamma") config.dGamma = std::stod(value);
        else if (key == "collection_u") config.collection_u = std::stod(value);
        else if (key == "collection_v") config.collection_v = std::stod(value);
        else if (key == "field_angle") config.field_angle = std::stod(value);
        else if (key == "refractive_index") config.refractive_index = std::stod(value);
        else if (key == "wavelength") config.wavelength = std::stod(value);
        else if (key == "elliptic_field") config.elliptic_field = (value == "true");
        else if (key == "ellipticity") config.ellipticity = std::stod(value);
        else if (key == "dipole_count") config.dipole_count = std::stoi(value);
        else if (key == "radius") config.radius = std::stod(value);
        else if (key == "geometry") config.geometry = value;
        else if (key == "geometry_file") config.geometry_file = value;
        else if (key == "enable_eme") config.enable_eme = (value == "true");
        else if (key == "enable_mee") config.enable_mee = (value == "true");
        else if (key == "enable_qee") config.enable_qee = (value == "true");
        else if (key == "output_prefix") config.output_prefix = value;
        else if (key == "output_dir") config.output_dir = value;
        else if (key == "save_config") config.save_config = (value == "true");
        else if (key == "save_metadata") config.save_metadata = (value == "true");
        else if (key == "num_threads") config.num_threads = std::stoi(value);
        else if (key == "use_optimization") config.use_optimization = (value == "true");
        else if (key == "run_validation") config.run_validation = (value == "true");
        else if (key == "validation_target") config.validation_target = value;
        else if (key == "validation_tolerance") config.validation_tolerance = std::stod(value);
        else if (key == "enable_uq") config.enable_uq = (value == "true");
        else if (key == "uq_samples") config.uq_samples = std::stoi(value);
    }
    
    return config;
}

void ConfigManager::save_yaml(const SimulationConfig& config, const std::string& filename) {
    std::ofstream file(filename);
    
    file << "# HRS-f Simulation Configuration\n";
    file << "# Generated automatically\n\n";
    
    file << "# Basic parameters\n";
    file << "frames: " << config.frames << "\n";
    file << "dGamma: " << config.dGamma << "\n";
    file << "collection_u: " << config.collection_u << "\n";
    file << "collection_v: " << config.collection_v << "\n\n";
    
    file << "# Field parameters\n";
    file << "field_angle: " << config.field_angle << "\n";
    file << "refractive_index: " << config.refractive_index << "\n";
    file << "wavelength: " << config.wavelength << "\n";
    file << "elliptic_field: " << (config.elliptic_field ? "true" : "false") << "\n";
    file << "ellipticity: " << config.ellipticity << "\n\n";
    
    file << "# Population parameters\n";
    file << "dipole_count: " << config.dipole_count << "\n";
    file << "radius: " << config.radius << "\n";
    file << "geometry: " << config.geometry << "\n";
    file << "geometry_file: " << config.geometry_file << "\n\n";
    
    file << "# Physics flags\n";
    file << "enable_eme: " << (config.enable_eme ? "true" : "false") << "\n";
    file << "enable_mee: " << (config.enable_mee ? "true" : "false") << "\n";
    file << "enable_qee: " << (config.enable_qee ? "true" : "false") << "\n\n";
    
    file << "# Output parameters\n";
    file << "output_prefix: " << config.output_prefix << "\n";
    file << "output_dir: " << config.output_dir << "\n";
    file << "save_config: " << (config.save_config ? "true" : "false") << "\n";
    file << "save_metadata: " << (config.save_metadata ? "true" : "false") << "\n\n";
    
    file << "# Performance parameters\n";
    file << "num_threads: " << config.num_threads << "\n";
    file << "use_optimization: " << (config.use_optimization ? "true" : "false") << "\n\n";
    
    file << "# Validation parameters\n";
    file << "run_validation: " << (config.run_validation ? "true" : "false") << "\n";
    file << "validation_target: " << config.validation_target << "\n";
    file << "validation_tolerance: " << config.validation_tolerance << "\n\n";
    
    file << "# Uncertainty quantification\n";
    file << "enable_uq: " << (config.enable_uq ? "true" : "false") << "\n";
    file << "uq_samples: " << config.uq_samples << "\n";
}

SimulationConfig ConfigManager::parse_cli(int argc, char* argv[]) {
    SimulationConfig config;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            print_help();
            exit(0);
        }
        else if (arg == "--config" && i + 1 < argc) {
            config = load_yaml(argv[++i]);
        }
        else if (arg == "--frames" && i + 1 < argc) {
            config.frames = std::stoi(argv[++i]);
        }
        else if (arg == "--dGamma" && i + 1 < argc) {
            config.dGamma = std::stod(argv[++i]);
        }
        else if (arg == "--radius" && i + 1 < argc) {
            config.radius = std::stod(argv[++i]);
        }
        else if (arg == "--dipoles" && i + 1 < argc) {
            config.dipole_count = std::stoi(argv[++i]);
        }
        else if (arg == "--wavelength" && i + 1 < argc) {
            config.wavelength = std::stod(argv[++i]);
        }
        else if (arg == "--refractive-index" && i + 1 < argc) {
            config.refractive_index = std::stod(argv[++i]);
        }
        else if (arg == "--output-prefix" && i + 1 < argc) {
            config.output_prefix = argv[++i];
        }
        else if (arg == "--output-dir" && i + 1 < argc) {
            config.output_dir = argv[++i];
        }
        else if (arg == "--enable-eme") {
            config.enable_eme = true;
        }
        else if (arg == "--enable-mee") {
            config.enable_mee = true;
        }
        else if (arg == "--enable-qee") {
            config.enable_qee = true;
        }
        else if (arg == "--threads" && i + 1 < argc) {
            config.num_threads = std::stoi(argv[++i]);
        }
        else if (arg == "--validation") {
            config.run_validation = true;
        }
        else if (arg == "--uq" && i + 1 < argc) {
            config.enable_uq = true;
            config.uq_samples = std::stoi(argv[++i]);
        }
        else if (arg == "--generate-config" && i + 1 < argc) {
            generate_template_config(argv[++i]);
            exit(0);
        }
    }
    
    return config;
}

void ConfigManager::print_help() {
    std::cout << "HRS-f: Harmonic Rayleigh Scattering Simulation\n\n";
    std::cout << "Usage: HRS-f [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  --help, -h                 Show this help message\n";
    std::cout << "  --config FILE              Load configuration from YAML file\n";
    std::cout << "  --frames N                 Number of simulation frames (default: 100)\n";
    std::cout << "  --dGamma ANGLE             Angular resolution in degrees (default: 5.0)\n";
    std::cout << "  --radius R                 Dipole arrangement radius (default: 10.0)\n";
    std::cout << "  --dipoles N                Number of dipoles (default: 12)\n";
    std::cout << "  --wavelength WL            Wavelength in nm (default: 800.0)\n";
    std::cout << "  --refractive-index N       Refractive index (default: 1.33)\n";
    std::cout << "  --output-prefix PREFIX     Output file prefix (default: simulation)\n";
    std::cout << "  --output-dir DIR           Output directory (default: output/data/)\n";
    std::cout << "  --enable-eme               Enable EME contributions\n";
    std::cout << "  --enable-mee               Enable MEE contributions\n";
    std::cout << "  --enable-qee               Enable QEE contributions\n";
    std::cout << "  --threads N                Number of OpenMP threads (default: auto)\n";
    std::cout << "  --validation               Run validation tests\n";
    std::cout << "  --uq N                     Enable uncertainty quantification with N samples\n";
    std::cout << "  --generate-config FILE     Generate template configuration file\n\n";
    std::cout << "Examples:\n";
    std::cout << "  HRS-f --frames 1000 --dGamma 2.5 --enable-eme\n";
    std::cout << "  HRS-f --config experiments/high_precision.yaml\n";
    std::cout << "  HRS-f --validation --uq 500\n";
}

bool ConfigManager::validate_config(const SimulationConfig& config, std::string& error_msg) {
    if (config.frames <= 0) {
        error_msg = "frames must be positive";
        return false;
    }
    if (config.dGamma <= 0 || config.dGamma > 360) {
        error_msg = "dGamma must be between 0 and 360 degrees";
        return false;
    }
    if (config.dipole_count <= 0) {
        error_msg = "dipole_count must be positive";
        return false;
    }
    if (config.wavelength <= 0) {
        error_msg = "wavelength must be positive";
        return false;
    }
    if (config.refractive_index <= 0) {
        error_msg = "refractive_index must be positive";
        return false;
    }
    if (config.uq_samples <= 0 && config.enable_uq) {
        error_msg = "uq_samples must be positive when UQ is enabled";
        return false;
    }
    
    return true;
}

void ConfigManager::generate_template_config(const std::string& filename) {
    SimulationConfig default_config;
    save_yaml(default_config, filename);
    std::cout << "Generated template configuration file: " << filename << std::endl;
}


