#ifndef CATCH_HPP_INCLUDED
#define CATCH_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include <vector>
#include <functional>
#include <cmath>
#include <string>

// Simple testing framework inspired by Catch2
// For production use, replace with actual Catch2

namespace Catch {
    class Approx {
        double value_;
        double epsilon_;
        
    public:
        explicit Approx(double value) : value_(value), epsilon_(1e-5) {}
        
        Approx& epsilon(double eps) {
            epsilon_ = eps;
            return *this;
        }
        
        bool operator==(double other) const {
            return std::abs(value_ - other) < epsilon_;
        }
        
        friend bool operator==(double lhs, const Approx& rhs) {
            return rhs == lhs;
        }
    };
}

struct TestCase {
    std::string name;
    std::string tags;
    std::function<void()> test_func;
};

class TestRegistry {
public:
    static TestRegistry& instance() {
        static TestRegistry registry;
        return registry;
    }
    
    void add_test(const std::string& name, const std::string& tags, 
                  std::function<void()> test_func) {
        tests_.push_back({name, tags, test_func});
    }
    
    int run_all() {
        int passed = 0;
        int failed = 0;
        
        std::cout << "Running " << tests_.size() << " test cases...\n\n";
        
        for (const auto& test : tests_) {
            try {
                std::cout << "Running: " << test.name << " " << test.tags << " ... ";
                test.test_func();
                std::cout << "PASSED\n";
                ++passed;
            } catch (const std::exception& e) {
                std::cout << "FAILED\n";
                std::cout << "  Error: " << e.what() << "\n";
                ++failed;
            }
        }
        
        std::cout << "\nResults: " << passed << " passed, " << failed << " failed\n";
        return failed;
    }
    
private:
    std::vector<TestCase> tests_;
};

class TestCaseRegistrar {
public:
    TestCaseRegistrar(const std::string& name, const std::string& tags,
                     std::function<void()> test_func) {
        TestRegistry::instance().add_test(name, tags, test_func);
    }
};

class AssertionFailed : public std::exception {
    std::string message_;
public:
    explicit AssertionFailed(const std::string& msg) : message_(msg) {}
    const char* what() const noexcept override { return message_.c_str(); }
};

#define TEST_CASE(name, tags) \
    static void test_case_##__LINE__##_func(); \
    static TestCaseRegistrar registrar_##__LINE__##_reg(name, tags, test_case_##__LINE__##_func); \
    static void test_case_##__LINE__##_func()

#define SECTION(name) \
    std::cout << "\n    Section: " << name << " ... ";

#define REQUIRE(condition) \
    do { \
        if (!(condition)) { \
            std::ostringstream oss; \
            oss << "REQUIRE failed at " << __FILE__ << ":" << __LINE__ \
                << "\n  Expression: " << #condition; \
            throw AssertionFailed(oss.str()); \
        } \
    } while(0)

using Catch::Approx;

// Main function for test executable
#ifdef CATCH_CONFIG_MAIN
int main() {
    return TestRegistry::instance().run_all();
}
#endif

#endif // CATCH_HPP_INCLUDED
