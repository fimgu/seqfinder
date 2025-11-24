/**
 * @file main.cpp
 * @brief Sequence Pattern Recognition CLI Application
 * 
 * This application analyzes numerical sequences and identifies patterns using
 * various specialized solvers. It supports multiple input formats including
 * fractions (e.g., "1/2"), powers

(e.g., "2^3"), and decimals.
 * 
 * @author Sequence Solver Team
 * @version 2.0
 */

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "SequenceSolvers.h"

// ============================================================================
// JSON OUTPUT UTILITIES
// ============================================================================

/**
 * @brief Escapes special characters in a string for JSON output
 * @param s The input string to escape
 * @return JSON-safe escaped string
 */
std::string escapeJSON(const std::string& s) {
    std::ostringstream o;
    for (char c : s) {
        switch (c) {
            case '"':  o << "\\\""; break;
            case '\\': o << "\\\\"; break;
            case '\b': o << "\\b";  break;
            case '\f': o << "\\f";  break;
            case '\n': o << "\\n";  break;
            case '\r': o << "\\r";  break;
            case '\t': o << "\\t";  break;
            default:
                // Skip control characters (0x00-0x1F)
                if (c >= 0 && c <= 0x1f) {
                    // Ignore
                } else {
                    o << c;
                }
        }
    }
    return o.str();
}

// ============================================================================
// POWER NOTATION UTILITIES
// ============================================================================

/**
 * @brief Represents a number in power notation (base^exponent)
 */
struct PowerNotation {
    bool isPower;    ///< Whether the number was successfully parsed as a power
    double base;     ///< The base of the power
    double exp;      ///< The exponent
    double value;    ///< The calculated value (base^exp)
};

/**
 * @brief Parses a string in power notation (e.g., "2^3" -> 8)
 * @param s The input string to parse
 * @return PowerNotation struct with parsed values
 */
PowerNotation parsePowerNotation(const std::string& s) {
    PowerNotation pn;
    pn.isPower = false;
   
    size_t caretPos = s.find('^');
    if (caretPos != std::string::npos) {
        try {
            pn.base = std::stod(s.substr(0, caretPos));
            pn.exp = std::stod(s.substr(caretPos + 1));
            pn.value = std::pow(pn.base, pn.exp);
            pn.isPower = true;
            return pn;
        } catch (...) {
            // Parsing failed, return isPower = false
        }
    }
   
    return pn;
}

/**
 * @brief Converts a decimal value to power notation if possible
 * @param val The value to convert
 * @return Power notation string (e.g., "2^3") or empty string if not a clean power
 */
std::string toPowerString(double val) {
    if (val <= 1) return "";
    
    // Try bases from 2 to 100
    for (int b = 2; b <= 100; ++b) {
        double e = std::log(val) / std::log(b);
        if (std::abs(e - std::round(e)) < 1e-9) {
            int exp = std::round(e);
            if (exp >= 2) {
                return std::to_string(b) + "^" + std::to_string(exp);
            }
        }
    }
    return "";
}

// ============================================================================
// INPUT PARSING
// ============================================================================

/**
 * @brief Parses a number from string, supporting fractions and powers
 * @param s The input string (e.g., "1/2", "2^3", "3.14")
 * @return The parsed numerical value
 */
double parseNumber(const std::string& s) {
    // Check for fraction notation (e.g., "1/2")
    size_t slashPos = s.find('/');
    if (slashPos != std::string::npos) {
        double num = std::stod(s.substr(0, slashPos));
        double den = std::stod(s.substr(slashPos + 1));
        return num / den;
    }
    
    // Check for power notation (e.g., "2^3")
    size_t caretPos = s.find('^');
    if (caretPos != std::string::npos) {
        double base = std::stod(s.substr(0, caretPos));
        double exp = std::stod(s.substr(caretPos + 1));
        return std::pow(base, exp);
    }
    
    // Regular decimal number
    return std::stod(s);
}

// ============================================================================
// SOLVER INITIALIZATION
// ============================================================================

/**
 * @brief Initializes all sequence solvers in priority order
 * @return Vector of solver instances
 * 
 * Solvers are ordered from most specific to most general to ensure
 * the most accurate pattern is identified first.
 */
std::vector<std::unique_ptr<SequenceSolver>> initializeSolvers() {
    std::vector<std::unique_ptr<SequenceSolver>> solvers;
    
    // Specialized pattern solvers (highest priority)
    solvers.push_back(std::make_unique<FactorialSolver>());
    solvers.push_back(std::make_unique<PowerNSolver>());
    solvers.push_back(std::make_unique<PrimeSolver>());

    solvers.push_back(std::make_unique<ComplexFractionSolver>());
    solvers.push_back(std::make_unique<MultiplierSolver>());
    solvers.push_back(std::make_unique<PowerSolver>());
    
    // Recurrence and difference solvers
    solvers.push_back(std::make_unique<LinearRecurrenceSolver>(2));
    solvers.push_back(std::make_unique<DifferencePatternSolver>());
    solvers.push_back(std::make_unique<PolynomialSolver>());
    
    // Interleaved solver (prioritized higher to detect fraction patterns)
    solvers.push_back(std::make_unique<InterleavedSolver>());
    
    solvers.push_back(std::make_unique<OperationCycleSolver>());
    solvers.push_back(std::make_unique<LinearRecurrenceSolver>(4));
    
    // Position-based and interleaved solvers
    solvers.push_back(std::make_unique<SkipPositionSolver>());
    solvers.push_back(std::make_unique<FractionIntegerInterleavedSolver>());
    solvers.push_back(std::make_unique<PositionDependentSolver>());
    solvers.push_back(std::make_unique<ConditionalRecurrenceSolver>());
    
    // General pattern solvers
    solvers.push_back(std::make_unique<AlternatingArithmeticSolver>());
    solvers.push_back(std::make_unique<CumulativeDifferenceSolver>());
    solvers.push_back(std::make_unique<AffineRecurrenceSolver>());
    solvers.push_back(std::make_unique<RootSolver>());
    solvers.push_back(std::make_unique<DependentInterleavedSolver>());
    
    // Fallback solvers (lowest priority)
    solvers.push_back(std::make_unique<DigitSolver>());
    solvers.push_back(std::make_unique<FractionSolver>());
    
    return solvers;
}

// ============================================================================
// OUTPUT FORMATTING
// ============================================================================

/**
 * @brief Formats the output value based on input format preferences
 * @param val The numerical value to format
 * @param formattedOutput Solver-provided formatted output (may be empty)
 * @param inputHasFractions Whether input contained fractions
 * @param inputHasPowers Whether input contained powers
 * @return Formatted string representation
 */
std::string formatOutput(double val, const std::string& formattedOutput, 
                        bool inputHasFractions, bool inputHasPowers) {
    std::string displayStr;
    
    // Prefer fraction output if input had fractions
    if (inputHasFractions) {
        if (!formattedOutput.empty() && formattedOutput.find('/') != std::string::npos) {
            displayStr = formattedOutput;
        } else {
            displayStr = toFractionString(val);
        }
    } 
    // Prefer power output if input had powers
    else if (inputHasPowers) {
        if (!formattedOutput.empty() && formattedOutput.find('^') != std::string::npos) {
            displayStr = formattedOutput;
        } else {
            displayStr = toPowerString(val);
        }
    }
    
    // Fallback to decimal/integer representation
    if (displayStr.empty()) {
        if (std::abs(val - std::round(val)) < 1e-9) {
            // Integer value
            displayStr = std::to_string((long long)std::round(val));
        } else {
            // Decimal value with high precision
            std::stringstream ss;
            ss << std::setprecision(15) << val;
            displayStr = ss.str();
        }
    }
    
    return displayStr;
}

/**
 * @brief Outputs a JSON result to stdout
 * @param pattern Pattern name
 * @param nextValue Predicted next value
 * @param formula Pattern formula/description
 * @param formattedNumber Formatted representation of the next value
 */
void outputResult(const std::string& pattern, double nextValue, 
                 const std::string& formula, const std::string& formattedNumber) {
    std::cout << "{";
    std::cout << "\"pattern\": \"" << escapeJSON(pattern) << "\", ";
    std::cout << "\"next_number\": " << nextValue << ", ";
    std::cout << "\"formula\": \"" << escapeJSON(formula) << "\", ";
    std::cout << "\"formatted_number\": \"" << escapeJSON(formattedNumber) << "\"";
    std::cout << "}" << std::endl;
}

/**
 * @brief Outputs a JSON error message to stdout
 * @param message Error message
 */
void outputError(const std::string& message) {
    std::cout << "{ \"error\": \"" << escapeJSON(message) << "\" }" << std::endl;
}

// ============================================================================
// MAIN APPLICATION
// ============================================================================

/**
 * @brief Main entry point for the sequence solver CLI
 * @param argc Argument count
 * @param argv Argument values (sequence numbers)
 * @return Exit code (0 for success, 1 for error)
 */
int main(int argc, char* argv[]) {
    // Set high precision for floating-point output
    std::cout << std::setprecision(15);
   
    // Validate command-line arguments
    if (argc < 2) {
        std::cerr << "Usage: seqfinder_cli.exe <num1> <num2> ...\n";
        return 1;
    }

    // Parse and validate input sequence
    std::vector<double> series;
    std::vector<PowerNotation> powerData;
    bool allPowers = true;
    bool inputHasFractions = false;
    bool inputHasPowers = false;
   
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        // Track input format types
        if (arg.find('/') != std::string::npos) inputHasFractions = true;
        if (arg.find('^') != std::string::npos) inputHasPowers = true;

        try {
            // Try parsing as power notation first
            PowerNotation pn = parsePowerNotation(arg);
            if (pn.isPower) {
                powerData.push_back(pn);
                series.push_back(pn.value);
            } else {
                allPowers = false;
                powerData.clear();
                series.push_back(parseNumber(arg));
            }
        } catch (...) {
            // Silently ignore invalid inputs
        }
    }

    // Validate that we have at least some valid numbers
    if (series.empty()) {
        outputError("No valid numbers provided");
        return 0;
    }

    // Special case: If all inputs were in power notation, try PowerNotationSolver
    if (allPowers && !powerData.empty()) {
        std::vector<double> bases, exps;
        for (const auto& pn : powerData) {
            bases.push_back(pn.base);
            exps.push_back(pn.exp);
        }
       
        PolynomialSolver polyBase, polyExp;
        auto resBase = polyBase.solve(bases);
        auto resExp = polyExp.solve(exps);
       
        if (!resBase.formula.empty() && !resExp.formula.empty()) {
            double nextBase = std::round(resBase.nextValue);
            double nextExp = std::round(resExp.nextValue);
            double result = std::pow(nextBase, nextExp);
           
            std::string formula = std::to_string((int)nextBase) + "^" + std::to_string((int)nextExp)
                                + " (Base: " + resBase.formula + ", Exp: " + resExp.formula + ")";
            std::string formatted = std::to_string((int)nextBase) + "^" + std::to_string((int)nextExp);
            
            outputResult("Power Notation", result, formula, formatted);
            return 0;
        }
    }

    // Initialize all solvers in priority order
    auto solvers = initializeSolvers();

    // Try each solver until one finds a pattern
    for (const auto& solver : solvers) {
        auto result = solver->solve(series);
        
        if (!result.formula.empty()) {
            // Pattern found! Format and output the result
            std::string formattedNumber = formatOutput(
                result.nextValue, 
                result.formattedOutput,
                inputHasFractions, 
                inputHasPowers
            );
            
            outputResult(
                solver->getName(),
                result.nextValue,
                result.formula,
                formattedNumber
            );
            
            return 0;
        }
    }

    // No pattern found
    outputError("Pattern not recognized");
    return 0;
}
