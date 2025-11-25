#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include "SequenceSolvers.h"

using namespace emscripten;

// Helper to escape JSON strings (copied from main.cpp)
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

// Helper to parse input (copied from main.cpp)
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

std::string solve_sequence_json(std::string input_str) {
    std::vector<BigRat> series;
    std::stringstream ss(input_str);
    std::string segment;
    
    // Split by space (input is expected to be space-separated numbers)
    while(std::getline(ss, segment, ' ')) {
        if (!segment.empty()) {
             try {
                series.push_back(parseInput(segment));
             } catch(...) {
                 return "{ \"error\": \"Invalid numerical input\" }";
             }
        }
    }

    if (series.empty()) {
        return "{ \"error\": \"Empty sequence\" }";
    }

    GrandUnifiedSolver solver;
    std::vector<SolverResult> results = solver.solve(series);

    std::stringstream out;
    out << "{\n";
    if (!results.empty() && results[0].confidence > 0.0) {
        SolverResult best = results[0];
        
        out << "  \"success\": true,\n";
        
        if (best.nextValue.get_den() != 1) {
             out << "  \"next_number\": \"" << best.nextValue.get_str() << "\",\n";
        } else {
             out << "  \"next_number\": " << best.nextValue.get_str() << ",\n";
        }

        if (!best.unsimplifiedNum.empty() && !best.unsimplifiedDen.empty()) {
             out << "  \"numerator\": \"" << best.unsimplifiedNum << "\",\n";
             out << "  \"denominator\": \"" << best.unsimplifiedDen << "\",\n";
        } else {
             out << "  \"numerator\": \"" << best.nextValue.get_num().get_str() << "\",\n";
             out << "  \"denominator\": \"" << best.nextValue.get_den().get_str() << "\",\n";
        }

        out << "  \"formula\": \"" << escapeJSON(best.formula) << "\",\n";
        out << "  \"pattern\": \"" << escapeJSON(best.formula) << "\",\n";
        out << "  \"explanation\": \"" << escapeJSON(best.explanation) << "\",\n";
        out << "  \"confidence\": " << best.confidence << ",\n";
        
        out << "  \"candidates\": [\n";
        for (size_t i = 0; i < std::min((size_t)5, results.size()); ++i) {
            out << "    { ";
            if (results[i].nextValue.get_den() != 1) {
                out << "\"next\": \"" << results[i].nextValue.get_str() << "\"";
            } else {
                out << "\"next\": " << results[i].nextValue.get_str();
            }
            out << ", \"pattern\": \"" << escapeJSON(results[i].formula) << "\""
                      << ", \"explanation\": \"" << escapeJSON(results[i].explanation) << "\""
                      << ", \"confidence\": " << results[i].confidence << " }";
            if (i < std::min((size_t)5, results.size()) - 1) out << ",";
            out << "\n";
        }
        out << "  ]\n";

    } else {
        out << "  \"success\": false,\n";
        out << "  \"error\": \"No pattern found\"\n";
    }
    out << "}";
    
    return out.str();
}

EMSCRIPTEN_BINDINGS(my_module) {
    function("solve_sequence_json", &solve_sequence_json);
}
