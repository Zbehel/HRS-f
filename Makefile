# Enhanced HRS-f Project Makefile with PhD-level features
# Unified version with all latest improvements

# Compiler and flags
CXX = clang++
CXXFLAGS = -std=gnu++14 -fcolor-diagnostics -fansi-escape-codes -g -Wall -Wextra -O2
CXXFLAGS_DEBUG = -std=gnu++14 -fcolor-diagnostics -fansi-escape-codes -g -Wall -Wextra -O0 -DDEBUG
CXXFLAGS_RELEASE = -std=gnu++14 -fcolor-diagnostics -fansi-escape-codes -O3 -DNDEBUG -march=native
INCLUDES = -I include

# Try to detect OpenMP support automatically
OPENMP_FLAGS := $(shell echo | $(CXX) -fopenmp -E -dM - 2>/dev/null | grep -q _OPENMP && echo "-fopenmp -DHAS_OPENMP")
LIBS = $(OPENMP_FLAGS)

# Directories
SRCDIR = src
BUILDDIR = build
OBJDIR = $(BUILDDIR)/obj
OUTPUTDIR = output
TESTDIR = tests
VALIDATIONDIR = validation
CONFIGDIR = config
DOCSDIR = docs
PROFILINGDIR = profiling
BENCHMARKDIR = benchmarks

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# Exclude main files from library objects
LIB_SOURCES = $(filter-out $(SRCDIR)/main.cpp $(SRCDIR)/main_enhanced.cpp, $(SOURCES))
LIB_OBJECTS = $(LIB_SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# Test sources
TEST_SOURCES = $(wildcard $(TESTDIR)/*.cpp)
TEST_OBJECTS = $(TEST_SOURCES:$(TESTDIR)/%.cpp=$(OBJDIR)/test_%.o)

# Validation sources
VALIDATION_SOURCES = $(wildcard $(VALIDATIONDIR)/**/*.cpp)
VALIDATION_OBJECTS = $(VALIDATION_SOURCES:$(VALIDATIONDIR)/%.cpp=$(OBJDIR)/validation_%.o)

# Targets
TARGET = $(BUILDDIR)/HRS-f
TARGET_ENHANCED = $(BUILDDIR)/HRS-f-enhanced
TARGET_TESTS = $(BUILDDIR)/run_tests
TARGET_VALIDATION = $(BUILDDIR)/rayleigh_validation

# Default target
all: $(TARGET_ENHANCED)

# Legacy target
legacy: $(TARGET)

# Create build directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Build the enhanced executable
$(TARGET_ENHANCED): $(LIB_OBJECTS) $(OBJDIR)/main_enhanced.o | $(BUILDDIR)
	$(CXX) $(LIB_OBJECTS) $(OBJDIR)/main_enhanced.o $(LIBS) -o $(TARGET_ENHANCED)

# Build the legacy executable
$(TARGET): $(LIB_OBJECTS) $(OBJDIR)/main.o | $(BUILDDIR)
	$(CXX) $(LIB_OBJECTS) $(OBJDIR)/main.o $(LIBS) -o $(TARGET)

# Build test executable
$(TARGET_TESTS): $(LIB_OBJECTS) $(TEST_OBJECTS) | $(BUILDDIR)
	$(CXX) $(LIB_OBJECTS) $(TEST_OBJECTS) $(LIBS) -o $(TARGET_TESTS)

# Build validation executable
$(TARGET_VALIDATION): $(LIB_OBJECTS) $(VALIDATION_OBJECTS) | $(BUILDDIR)
	$(CXX) $(LIB_OBJECTS) $(VALIDATION_OBJECTS) $(LIBS) -o $(TARGET_VALIDATION)

# Compile source files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile test files
$(OBJDIR)/test_%.o: $(TESTDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile validation files (handle nested directories)
$(OBJDIR)/validation_%.o: $(VALIDATIONDIR)/*/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create output directories
$(OUTPUTDIR):
	mkdir -p $(OUTPUTDIR)/data $(OUTPUTDIR)/configs $(OUTPUTDIR)/plots

# Create profiling directory
$(PROFILINGDIR):
	mkdir -p $(PROFILINGDIR)

# Create benchmark directory
$(BENCHMARKDIR):
	mkdir -p $(BENCHMARKDIR)

# Debug build
debug: CXXFLAGS = $(CXXFLAGS_DEBUG)
debug: $(TARGET_ENHANCED)

# Release build
release: CXXFLAGS = $(CXXFLAGS_RELEASE)
release: $(TARGET_ENHANCED)

# Run the enhanced program with default config
run: $(TARGET_ENHANCED) $(OUTPUTDIR)
	./$(TARGET_ENHANCED) --config $(CONFIGDIR)/experiments/default.yaml

# Run with high precision config
run-precision: $(TARGET_ENHANCED) $(OUTPUTDIR)
	./$(TARGET_ENHANCED) --config $(CONFIGDIR)/experiments/high_precision.yaml

# Run with custom parameters
run-custom: $(TARGET_ENHANCED) $(OUTPUTDIR)
	@echo "Running with custom parameters..."
	./$(TARGET_ENHANCED) --frames 1000 --dGamma 2.5 --enable-eme --enable-mee

# Run unit tests
test: $(TARGET_TESTS)
	./$(TARGET_TESTS)

# Run validation tests
validate: $(TARGET_VALIDATION) $(OUTPUTDIR)
	./$(TARGET_VALIDATION)

# Run all tests (unit + validation)
test-all: test validate

# Performance profiling
profile: CXXFLAGS += -pg
profile: $(TARGET_ENHANCED) $(PROFILINGDIR)
	./$(TARGET_ENHANCED) --config $(CONFIGDIR)/experiments/default.yaml
	gprof $(TARGET_ENHANCED) gmon.out > $(PROFILINGDIR)/profile_report.txt
	@echo "Profile report saved to $(PROFILINGDIR)/profile_report.txt"

# Memory checking with valgrind (if available)
memcheck: $(TARGET_ENHANCED)
	@echo "Running memory check with valgrind..."
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all \
		./$(TARGET_ENHANCED) --config $(CONFIGDIR)/experiments/default.yaml

# Benchmark performance
benchmark: $(TARGET_ENHANCED) $(OUTPUTDIR) $(BENCHMARKDIR)
	@echo "Running performance benchmarks..."
	@for frames in 100 500 1000 5000; do \
		echo "Testing with $$frames frames..."; \
		time ./$(TARGET_ENHANCED) --frames $$frames --output-prefix benchmark_$$frames; \
	done
	@echo "Benchmark results saved to $(OUTPUTDIR)/data/"

# Generate documentation
docs: $(DOCSDIR)
	@echo "Generating documentation..."
	@mkdir -p $(DOCSDIR)/generated
	@echo "# HRS-f API Documentation" > $(DOCSDIR)/generated/API.md
	@echo "Generated on: $$(date)" >> $(DOCSDIR)/generated/API.md
	@echo "" >> $(DOCSDIR)/generated/API.md
	@grep -n "class\|struct" include/*.h >> $(DOCSDIR)/generated/API.md || true
	@echo "Documentation generated in $(DOCSDIR)/generated/"

# Static analysis
analyze:
	@echo "Running static analysis..."
	@if command -v clang-tidy >/dev/null 2>&1; then \
		clang-tidy $(SOURCES) -- $(INCLUDES); \
	else \
		echo "clang-tidy not found. Install clang-tools for static analysis."; \
	fi

# Format code
format:
	@echo "Formatting code..."
	@if command -v clang-format >/dev/null 2>&1; then \
		clang-format -i $(SOURCES) include/*.h; \
		echo "Code formatted successfully."; \
	else \
		echo "clang-format not found. Install clang-tools for code formatting."; \
	fi

# Generate configuration templates
generate-configs: $(TARGET_ENHANCED)
	./$(TARGET_ENHANCED) --generate-config $(CONFIGDIR)/template.yaml
	@echo "Template configuration generated: $(CONFIGDIR)/template.yaml"

# Quick development cycle
dev: clean debug test
	@echo "Development build completed successfully!"

# Full development cycle with validation
dev-full: clean debug test-all validate
	@echo "Full development cycle completed successfully!"

# Clean build files
clean:
	rm -rf $(BUILDDIR)
	rm -rf $(OUTPUTDIR)
	rm -f gmon.out

# Clean only object files
clean-obj:
	rm -rf $(OBJDIR)

# Clean everything including generated files
clean-all: clean
	rm -rf $(DOCSDIR)/generated
	rm -rf $(PROFILINGDIR)/*.txt
	rm -rf $(VALIDATIONDIR)/*.txt
	rm -rf $(BENCHMARKDIR)/*

# Install (copy executable to system path)
install: $(TARGET_ENHANCED)
	@echo "Installing HRS-f to system..."
	cp $(TARGET_ENHANCED) /usr/local/bin/hrs-f
	@echo "HRS-f installed successfully. Run with 'hrs-f --help'"

# Package for distribution
package: release
	@echo "Creating distribution package..."
	@mkdir -p dist/hrs-f
	@cp $(TARGET_ENHANCED) dist/hrs-f/
	@cp README_ENHANCED.md dist/hrs-f/README.md
	@cp -r $(CONFIGDIR) dist/hrs-f/
	@cp -r $(DOCSDIR) dist/hrs-f/
	@tar -czf dist/hrs-f.tar.gz -C dist hrs-f
	@echo "Package created: dist/hrs-f.tar.gz"

# Continuous integration target
ci: clean test-all benchmark
	@echo "Continuous integration completed successfully!"

# Show build information
info:
	@echo "HRS-f Build Information:"
	@echo "  Compiler: $(CXX)"
	@echo "  Flags: $(CXXFLAGS)"
	@echo "  OpenMP: $(if $(OPENMP_FLAGS),Enabled,Disabled)"
	@echo "  Source files: $(words $(SOURCES))"
	@echo "  Test files: $(words $(TEST_SOURCES))"
	@echo "  Validation files: $(words $(VALIDATION_SOURCES))"

# Help
help:
	@echo "HRS-f Enhanced: PhD-Level Harmonic Rayleigh Scattering Simulation"
	@echo ""
	@echo "Available targets:"
	@echo "  all          - Build enhanced version (default)"
	@echo "  legacy       - Build legacy version"
	@echo "  debug        - Build debug version"
	@echo "  release      - Build optimized release version"
	@echo "  run          - Build and run with default config"
	@echo "  run-precision- Build and run with high precision config"
	@echo "  run-custom   - Build and run with custom parameters"
	@echo "  test         - Run unit tests"
	@echo "  validate     - Run validation tests"
	@echo "  test-all     - Run all tests"
	@echo "  profile      - Build and run with profiling"
	@echo "  memcheck     - Run with memory checking"
	@echo "  benchmark    - Run performance benchmarks"
	@echo "  docs         - Generate documentation"
	@echo "  analyze      - Run static analysis"
	@echo "  format       - Format source code"
	@echo "  generate-configs - Generate template configuration"
	@echo "  dev          - Quick development cycle (clean+debug+test)"
	@echo "  dev-full     - Full development cycle with validation"
	@echo "  clean        - Remove build files"
	@echo "  clean-obj    - Remove only object files"
	@echo "  clean-all    - Remove all generated files"
	@echo "  install      - Install to system"
	@echo "  package      - Create distribution package"
	@echo "  ci           - Continuous integration target"
	@echo "  info         - Show build information"
	@echo "  help         - Show this help"
	@echo ""
	@echo "Examples:"
	@echo "  make all                    # Build enhanced version"
	@echo "  make run                    # Run with default settings"
	@echo "  make run-precision          # Run high-precision simulation"
	@echo "  make test-all               # Run all tests"
	@echo "  make benchmark              # Performance benchmarks"
	@echo "  make dev                    # Quick development cycle"

.PHONY: all legacy debug release run run-precision run-custom test validate test-all profile memcheck benchmark docs analyze format generate-configs dev dev-full clean clean-obj clean-all install package ci info help