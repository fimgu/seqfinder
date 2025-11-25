#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "SequenceSolvers.h"

std::string escapeJSON(const std::string& s) {
    std::ostringstream o;
    for (char c : s) {
        if (c == '"') o << "\\\"";
        else if (c == '\\') o << "\\\\";
        else if (c < 32) {} 
        else o << c;
    }
    return o.str();
}

BigRat parseInput(const std::string& arg) {
    std::string s = arg;
    s.erase(std::remove(s.begin(), s.end(), ','), s.end());
    size_t decPos = s.find('.');
    if (decPos != std::string::npos) {
        std::string numStr = s.substr(0, decPos) + s.substr(decPos + 1);
        std::string denStr = "1";
        for (size_t i = 0; i < s.length() - decPos - 1; ++i) denStr += "0";
        return BigRat(BigInt(numStr), BigInt(denStr));
    }
    if (s.find('/') != std::string::npos) return BigRat(s);
    try { return BigRat(s); } catch (...) { return BigRat(0); }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "{ \"error\": \"No input provided\" }" << std::endl;
        return 1;
    }

    std::vector<BigRat> series;
    bool hasFloatInput = false;

    try {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg.find('.') != std::string::npos) hasFloatInput = true;
            series.push_back(parseInput(arg));
        }
    } catch (...) {
        std::cout << "{ \"error\": \"Invalid numerical input\" }" << std::endl;
        return 1;
    }

    if (series.empty()) {
        std::cout << "{ \"error\": \"Empty sequence\" }" << std::endl;
        return 1;
    }

    GrandUnifiedSolver solver;
    std::vector<SolverResult> results = solver.solve(series);

    std::cout << "{\n";
    if (!results.empty() && results[0].confidence > 0.0) {
        SolverResult best = results[0];
        
        std::cout << "  \"success\": true,\n";
        
        // Output next_number as string "num/den" if fraction, else number
        if (best.nextValue.get_den() != 1) {
             std::cout << "  \"next_number\": \"" << best.nextValue.get_str() << "\",\n";
        } else {
             std::cout << "  \"next_number\": " << best.nextValue.get_str() << ",\n";
        }

        if (!best.unsimplifiedNum.empty() && !best.unsimplifiedDen.empty()) {
             std::cout << "  \"numerator\": \"" << best.unsimplifiedNum << "\",\n";
             std::cout << "  \"denominator\": \"" << best.unsimplifiedDen << "\",\n";
        } else {
             std::cout << "  \"numerator\": \"" << best.nextValue.get_num().get_str() << "\",\n";
             std::cout << "  \"denominator\": \"" << best.nextValue.get_den().get_str() << "\",\n";
        }

        std::cout << "  \"formula\": \"" << escapeJSON(best.formula) << "\",\n";
        std::cout << "  \"pattern\": \"" << escapeJSON(best.formula) << "\",\n";
        std::cout << "  \"explanation\": \"" << escapeJSON(best.explanation) << "\",\n";
        std::cout << "  \"confidence\": " << best.confidence << ",\n";
        
        // Print other candidates
        std::cout << "  \"candidates\": [\n";
        for (size_t i = 0; i < std::min((size_t)5, results.size()); ++i) {
            std::cout << "    { ";
            if (results[i].nextValue.get_den() != 1) {
                std::cout << "\"next\": \"" << results[i].nextValue.get_str() << "\"";
            } else {
                std::cout << "\"next\": " << results[i].nextValue.get_str();
            }
            std::cout << ", \"pattern\": \"" << escapeJSON(results[i].formula) << "\""
                      << ", \"explanation\": \"" << escapeJSON(results[i].explanation) << "\""
                      << ", \"confidence\": " << results[i].confidence << " }";
            if (i < std::min((size_t)5, results.size()) - 1) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "  ]\n";

    } else {
        std::cout << "  \"success\": false,\n";
        std::cout << "  \"error\": \"No pattern found\"\n";
    }
    std::cout << "}" << std::endl;

    return 0;
}