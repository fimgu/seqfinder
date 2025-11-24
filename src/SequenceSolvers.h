#ifndef SEQUENCE_SOLVERS_H
#define SEQUENCE_SOLVERS_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <gmpxx.h>

// ============================================================================
// TYPES & UTILS
// ============================================================================

using BigRat = mpq_class;
using BigInt = mpz_class;

inline std::string toString(const BigInt& n) { return n.get_str(); }
inline std::string toString(const BigRat& n) { return n.get_str(); }

inline BigRat absRat(const BigRat& n) {
    return (n < 0) ? -n : n;
}

// --- PRIME UTILITIES ---
inline bool isPrime(const BigInt& n) {
    if (n <= 1) return false;
    return mpz_probab_prime_p(n.get_mpz_t(), 25) > 0;
}

inline BigInt getNthPrime(long n) {
    if (n < 1) return 2;
    BigInt p = 2;
    for (long i = 1; i < n; ++i) mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    return p;
}

inline long getPrimeIndex(const BigInt& n) {
    if (!isPrime(n)) return -1;
    BigInt p = 2;
    long idx = 1;
    while (p < n) {
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        idx++;
    }
    return (p == n) ? idx : -1;
}

// --- LINEAR ALGEBRA UTILS ---
inline std::vector<BigRat> solveHomogeneousSystem(std::vector<std::vector<BigRat>>& mat) {
    if (mat.empty()) return {};
    size_t rows = mat.size();
    size_t cols = mat[0].size();
    size_t pivot_row = 0;
    std::vector<size_t> pivot_cols;
    
    for (size_t col = 0; col < cols && pivot_row < rows; ++col) {
        size_t sel = pivot_row;
        while (sel < rows && mat[sel][col] == 0) sel++;
        if (sel == rows) continue;
        std::swap(mat[sel], mat[pivot_row]);
        pivot_cols.push_back(col);
        BigRat pivot_val = mat[pivot_row][col];
        for (size_t j = col; j < cols; ++j) mat[pivot_row][j] /= pivot_val;
        for (size_t i = 0; i < rows; ++i) {
            if (i != pivot_row && mat[i][col] != 0) {
                BigRat factor = mat[i][col];
                for (size_t j = col; j < cols; ++j) mat[i][j] -= factor * mat[pivot_row][j];
            }
        }
        pivot_row++;
    }

    if (pivot_cols.size() < cols) {
        std::vector<BigRat> sol(cols, 0);
        int free_col = -1;
        for(int c = (int)cols-1; c >= 0; --c) {
             bool is_pivot = false;
             for(size_t p : pivot_cols) if((int)p == c) is_pivot = true;
             if(!is_pivot) { free_col = c; break; }
        }
        if(free_col == -1) return {};
        sol[free_col] = 1;
        for (int i = (int)pivot_row - 1; i >= 0; --i) {
            int p_col = -1;
            for(size_t j=0; j<cols; ++j) {
                if(mat[i][j] != 0) { p_col = j; break; }
            }
            if (p_col == -1 || p_col > free_col) continue;
            BigRat sum = 0;
            for (size_t j = p_col + 1; j < cols; ++j) sum += mat[i][j] * sol[j];
            sol[p_col] = -sum;
        }
        return sol;
    }
    return {};
}

struct SolverResult {
    BigRat nextValue;
    std::string formula;
    std::string formattedOutput;
    double confidence; 
    std::string unsimplifiedNum = "";
    std::string unsimplifiedDen = "";
    std::string explanation = "";
};

inline bool compareResults(const SolverResult& a, const SolverResult& b) {
    return a.confidence > b.confidence;
}

class SequenceSolver {
public:
    virtual ~SequenceSolver() = default;
    virtual std::string getName() const = 0;
    virtual std::vector<SolverResult> solve(const std::vector<BigRat>& series) = 0;
};

// ============================================================================
// 1. BASE SOLVERS
// ============================================================================

class PolynomialSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Polynomial"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        
        // Check Constant
        bool allSame = true;
        for(size_t i=1; i<series.size(); ++i) if(series[i] != series[0]) allSame = false;
        if(allSame) return {{series[0], "Constant", "", 1.0, "", "", "All terms are equal to " + toString(series[0].get_num())}};

        std::vector<std::vector<BigRat>> table;
        table.push_back(series);
        int depth = 0;
        size_t max_depth = series.size() - 1;

        while (table.back().size() > 1 && depth < (int)max_depth) {
            const auto& row = table.back();
            bool constantRow = true;
            for (size_t i = 1; i < row.size(); ++i) {
                if (row[i] != row[0]) { constantRow = false; break; }
            }
            if (constantRow) {
                BigRat nextVal = row[0]; 
                for (int i = depth - 1; i >= 0; --i) nextVal = nextVal + table[i].back();
                
                double conf = 0.97;
                // Penalize overfitting: if Degree is close to N
                if (depth > 1 && (size_t)depth >= series.size() - 1) conf = 0.5;
                
                std::string expl = "Polynomial of degree " + std::to_string(depth) + ". Constant difference found at depth " + std::to_string(depth) + ".";
                return {{nextVal, "Polynomial (Degree " + std::to_string(depth) + ")", "", conf, "", "", expl}};
            }
            std::vector<BigRat> nextRow;
            for(size_t i = 0; i < row.size() - 1; ++i) nextRow.push_back(row[i+1] - row[i]);
            table.push_back(nextRow);
            depth++;
        }
        return {};
    }
};

class PowerSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Geometric"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3 || series[0] == 0) return {};
        BigRat ratio = series[1] / series[0];
        if (ratio == 1) return {}; 
        for(size_t i=1; i<series.size()-1; ++i) {
            if (series[i] == 0 || series[i+1] / series[i] != ratio) return {};
        }
        return {{series.back() * ratio, "Geometric", "", 1.0, "", "", "Geometric progression with common ratio " + toString(ratio)}};
    }
};

class StandardMathSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Standard Math"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        int offset = (series[0] == 1 && series[1] == 4) ? 1 : 0;
        bool isPower = true;
        for(size_t i=0; i<series.size(); ++i) {
            unsigned long idx = i + offset;
            BigInt p; mpz_ui_pow_ui(p.get_mpz_t(), idx, idx);
            if (series[i].get_num() != p) { isPower = false; break; }
        }
        if (isPower) {
            unsigned long nextIdx = series.size() + offset;
            BigInt res; mpz_ui_pow_ui(res.get_mpz_t(), nextIdx, nextIdx);
            return {{BigRat(res), "n^n Sequence", "", 1.0, "", "", "Each term is n raised to the power of n (n^n)."}};
        }
        return {};
    }
};

class DigitSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Digit Patterns"; }
    std::string lookAndSayNext(std::string s) {
        if (s.empty()) return "";
        std::stringstream ss;
        for (size_t i = 0; i < s.length(); i++) {
            int count = 1;
            while (i + 1 < s.length() && s[i] == s[i+1]) { i++; count++; }
            ss << count << s[i];
        }
        return ss.str();
    }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        bool isLAS = true;
        for(size_t i=0; i<series.size()-1; ++i) {
             if (series[i].get_den() != 1) { isLAS = false; break; }
             if (lookAndSayNext(series[i].get_num().get_str()) != series[i+1].get_num().get_str()) { isLAS = false; break; }
        }
        if (isLAS && series.size() > 2) {
            return {{BigRat(BigInt(lookAndSayNext(series.back().get_num().get_str()))), "Look-and-Say Sequence", "", 1.0, "", "", "Each term describes the digits of the previous term (e.g. '1211' is read as 'one 1, one 2, two 1s')."}};
        }
        return {};
    }
};

class BerlekampMasseySolver : public SequenceSolver {
public:
    std::string getName() const override { return "Linear Recurrence"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 4) return {};
        int N = series.size();
        std::vector<BigRat> C = {1}, B = {1};
        int L = 0, m = -1;
        BigRat b(1);

        for (int n = 0; n < N; ++n) {
            BigRat d(0);
            for (size_t i = 0; i < C.size(); ++i) 
                if (n >= (int)i) d += C[i] * series[n - i];
            if (d == 0) continue;
            std::vector<BigRat> T = C;
            BigRat coef = d / b;
            int shift = n - m;
            std::vector<BigRat> scaledB(shift, 0);
            for(const auto& val : B) scaledB.push_back(val * coef);
            if (scaledB.size() > T.size()) T.resize(scaledB.size(), 0);
            for(size_t i=0; i<scaledB.size(); ++i) T[i] -= scaledB[i];
            if (2 * L <= n) { L = n + 1 - L; m = n; B = C; b = d; }
            C = T;
        }
        
        if (L == 0 || 2 * L > N + 1) return {}; 

        BigRat nextVal(0);
        BigRat c0 = C[0];
        for (size_t i = 1; i < C.size(); ++i) {
            if (series.size() >= i) nextVal += ((C[i] * -1) / c0) * series[series.size() - i];
        }
        return {{nextVal, "Linear Recurrence", "", 0.85, "", "", "Linear recurrence relation found using Berlekamp-Massey algorithm."}}; 
    }
};

class RecursiveSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Recursive Non-Linear"; }
    
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        std::vector<SolverResult> results;
        size_t N = series.size();
        if (N < 3) return {};
        
        // 1. Pure Square
        bool sq = true;
        for(size_t i=1; i<N; ++i) if(series[i] != series[i-1]*series[i-1]) sq=false;
        if (sq) results.push_back({series.back()*series.back(), "x_n = x_{n-1}^2", "", 1.0, "", "", "Each term is the square of the previous term."});

        // 2. Sylvester
        bool sylv = true;
        for(size_t i=1; i<N; ++i) if(series[i] != (series[i-1]*series[i-1]) + series[i-1]) sylv=false;
        if (sylv) results.push_back({(series.back()*series.back())+series.back(), "x_n = x_{n-1}^2 + x_{n-1}", "", 1.0, "", "", "Each term is the square of the previous term plus the previous term (Sylvester's sequence logic)."});
        
        // 3. Lag 2 Relations
        if (N >= 3) {
            bool lag1 = true, lag2 = true;
            for(size_t i=2; i<N; ++i) {
                if (series[i] != (series[i-1] * series[i-2]) - series[i-2]) lag1 = false;
                if (series[i] != (series[i-1] * series[i-2]) - series[i-1]) lag2 = false;
            }
            if (lag1) results.push_back({(series.back() * series[N-2]) - series[N-2], "x_n = x_{n-1}x_{n-2} - x_{n-2}", "", 1.0, "", "", "Term equals product of two previous terms minus the second previous term."});
            if (lag2) results.push_back({(series.back() * series[N-2]) - series.back(), "x_n = x_{n-1}x_{n-2} - x_{n-1}", "", 1.0, "", "", "Term equals product of two previous terms minus the previous term."});
        }

        // 4. Fibonacci / Subtraction Logic
        if (N >= 3) {
            // Check for AP first
            bool isAP = true;
            BigRat diff = series[1] - series[0];
            for(size_t i=1; i<N; ++i) if(series[i] - series[i-1] != diff) isAP = false;

            if (!isAP) {
                bool fibSum = true, fibSub1 = true, fibSub2 = true;
                for(size_t i=2; i<N; ++i) {
                    if (series[i] != series[i-1] + series[i-2]) fibSum = false;
                    if (series[i] != series[i-1] - series[i-2]) fibSub1 = false;
                    if (series[i] != series[i-2] - series[i-1]) fibSub2 = false;
                }
                if (fibSum) results.push_back({series.back() + series[N-2], "x_n = x_{n-1} + x_{n-2}", "", 0.99, "", "", "Fibonacci-like sequence: each term is the sum of the two preceding terms."});
                if (fibSub1) results.push_back({series.back() - series[N-2], "x_n = x_{n-1} - x_{n-2}", "", 0.99, "", "", "Each term is the difference between the two preceding terms (x_{n-1} - x_{n-2})."});
                if (fibSub2) results.push_back({series[N-2] - series.back(), "x_n = x_{n-2} - x_{n-1}", "", 0.99, "", "", "Each term is the difference between the two preceding terms (x_{n-2} - x_{n-1})."});
            }
        }

        // 5. Variable Increment/Multiplier
        if (N >= 3) {
            bool varInc = true;
            BigRat k = (series[1] - series[0]) - 1;
            for(size_t i=1; i<N; ++i) {
                if (series[i] != series[i-1] + BigRat(static_cast<unsigned long>(i)) + k) { varInc = false; break; }
            }
            if (varInc) results.push_back({series.back() + BigRat(static_cast<unsigned long>(N)) + k, "x_n = x_{n-1} + n + K", "", 1.0, "", "", "The difference between terms increases linearly."});
        }

        // 6. Coupled Step Recurrence (Solves "Variable Increment Sum")
        // x_n = x_{n-1} + x_{n-2} if n odd, x_{n-1} + x_{n-3} if n even?
        // Checking: 4, 5, 9, 13, 22, 31
        // 9 = 5+4. 13 = 9+4. 22 = 13+9. 31 = 22+9.
        // Rule: x_n = x_{n-1} + x_{n-2} (if n is even index 2, 4, 6)
        // Rule: x_n = x_{n-1} + x_{n-3} (if n is odd index 3, 5)
        // Actually: Indices 2,3,4,5. 
        // 2(9)=1(5)+0(4).  3(13)=2(9)+0(4).  4(22)=3(13)+2(9).  5(31)=4(22)+2(9).
        // Parity Pattern: x_n = x_{n-1} + x_{n - (2 + (n%2))}
        if (N >= 4) {
            bool coupled = true;
            for(size_t i=2; i<N; ++i) {
                size_t offset = 2 + (i % 2); // if i=2 offset=2. i=3 offset=3.
                if (i < offset) continue; 
                if (series[i] != series[i-1] + series[i-offset]) { coupled = false; break; }
            }
            if (coupled) {
                size_t offset = 2 + (N % 2);
                results.push_back({series.back() + series[N-offset], "Coupled Step Recurrence", "", 1.0, "", "", "Coupled recurrence relation depending on the parity of the index."});
            }
        }

        return results;
    }
};

class HolonomicSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Holonomic"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        int N = series.size();
        if (N < 5) return {}; 
        std::vector<std::pair<int, int>> configs = {
            {1, 1}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {3, 0}, {3, 1}
        };
        for (auto [R, D] : configs) {
            int num_coeffs = (R + 1) * (D + 1);
            if (N - R < num_coeffs) continue;
            std::vector<std::vector<BigRat>> matrix;
            for (int k = 0; k < N - R; ++k) {
                std::vector<BigRat> row;
                BigInt n_val = k; 
                for (int r = 0; r <= R; ++r) {
                    BigRat a_term = series[k + r];
                    BigInt n_pow = 1;
                    for (int d = 0; d <= D; ++d) { row.push_back(a_term * BigRat(n_pow)); n_pow *= n_val; }
                }
                matrix.push_back(row);
            }
            auto coeffs = solveHomogeneousSystem(matrix);
            if (!coeffs.empty()) {
                BigInt n_target = N - R;
                BigRat rhs = 0, lhs_coeff = 0; 
                int coeff_idx = 0;
                for (int r = 0; r <= R; ++r) {
                    BigRat p_r_val = 0; BigInt n_pow = 1;
                    for (int d = 0; d <= D; ++d) { p_r_val += coeffs[coeff_idx++] * BigRat(n_pow); n_pow *= n_target; }
                    if (r == R) lhs_coeff = p_r_val;
                    else rhs -= p_r_val * series[n_target.get_si() + r];
                }
                if (lhs_coeff != 0) return {{rhs / lhs_coeff, "Holonomic (Order " + std::to_string(R) + ")", "", 0.90, "", "", "Holonomic sequence (D-finite) solution found."}};
            }
        }
        return {};
    }
};

// ============================================================================
// 2. ADVANCED SOLVERS
// ============================================================================

class PrimeSolver : public SequenceSolver {
    std::unique_ptr<PolynomialSolver> poly;
public:
    PrimeSolver() { poly = std::make_unique<PolynomialSolver>(); }
    std::string getName() const override { return "Primes"; }
    
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        for(auto& s : series) if(s.get_den() != 1) return {};

        std::vector<BigRat> indices;
        bool allPrime = true;
        for(const auto& s : series) {
            BigInt val = s.get_num();
            if (val <= 0) { allPrime = false; break; }
            long idx = getPrimeIndex(val);
            if (idx == -1) { allPrime = false; break; }
            indices.push_back(BigRat(idx));
        }

        if (allPrime) {
            auto results = poly->solve(indices);
            if (!results.empty() && results[0].confidence > 0.5) {
                BigInt nextIdx = results[0].nextValue.get_num();
                if (nextIdx > 0) return {{BigRat(getNthPrime(nextIdx.get_si())), "Indexed Primes", "", 0.99, "", "", "Sequence consists of primes at indices following a pattern: " + results[0].formula}};
            }
        }
        
        std::vector<BigRat> diffs;
        bool diffsArePrime = true;
        bool decreasing = true; 
        
        for(size_t i=0; i<series.size()-1; ++i) {
            BigRat d = series[i+1] - series[i];
            if (d > 0) decreasing = false;
            BigInt num = d.get_num();
            if (num < 0) num = -num; 
            if (!isPrime(num)) { diffsArePrime = false; break; }
            diffs.push_back(d);
        }
        
        if (diffsArePrime && !diffs.empty()) {
            std::vector<BigRat> absDiffs;
            for(auto& d : diffs) absDiffs.push_back(absRat(d));
            
            std::vector<BigRat> diffIndices;
            for(auto& d : absDiffs) diffIndices.push_back(BigRat(getPrimeIndex(d.get_num())));
            
            auto results = poly->solve(diffIndices);
            
            if(!results.empty() && results[0].confidence > 0.5) {
                BigInt nextIdx = results[0].nextValue.get_num();
                BigInt nextPrimeVal = getNthPrime(nextIdx.get_si());
                
                BigRat nextDiff = BigRat(nextPrimeVal);
                if (decreasing && diffs.back() < 0) nextDiff = -nextDiff;
                else if (!decreasing && diffs.back() > 0) nextDiff = nextDiff;
                else if (diffs.back() < 0) nextDiff = -nextDiff;

                return {{series.back() + nextDiff, "Prime Differences", "", 0.94, "", "", "The differences between terms are prime numbers following a pattern."}}; 
            }
        }
        return {};
    }
};

class CyclicOpSolver : public SequenceSolver {
    struct Op { char type; BigRat val; }; 

    void generateCombinations(
        size_t L, 
        size_t currentStep, 
        std::vector<Op>& currentOps, 
        const std::vector<std::vector<Op>>& possibleOpsPerStep, 
        std::vector<std::vector<Op>>& allCombinations
    ) {
        if (currentStep == L) {
            allCombinations.push_back(currentOps);
            return;
        }

        for (const auto& op : possibleOpsPerStep[currentStep]) {
            currentOps[currentStep] = op;
            generateCombinations(L, currentStep + 1, currentOps, possibleOpsPerStep, allCombinations);
        }
    }

public:
    std::string getName() const override { return "Cyclic Operations"; }
    
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        size_t N = series.size();
        if (N < 4) return {};

        std::vector<SolverResult> results;

        for (size_t L = 1; L < N; ++L) {
            std::vector<std::vector<Op>> possibleOpsPerStep(L);
            bool cyclePossible = true;

            for (size_t k = 0; k < L; ++k) {
                if (k + 1 >= N) { cyclePossible = false; break; }

                // 1. Check Addition
                BigRat diff = series[k+1] - series[k];
                bool isAdd = true;
                for (size_t i = k; i < N - 1; i += L) {
                    if (series[i+1] != series[i] + diff) { isAdd = false; break; }
                }
                if (isAdd) possibleOpsPerStep[k].push_back({'a', diff});

                // 2. Check Multiplication
                bool isMult = true;
                BigRat div(0);
                if (series[k] == 0) {
                    isMult = false; 
                } else {
                    div = series[k+1] / series[k];
                    for (size_t i = k; i < N - 1; i += L) {
                        if (series[i] == 0 || series[i+1] != series[i] * div) { isMult = false; break; }
                    }
                }
                if (isMult) possibleOpsPerStep[k].push_back({'m', div});

                // 3. Check Power (Square)
                bool isSq = true;
                for (size_t i = k; i < N - 1; i += L) {
                    if (series[i+1] != series[i] * series[i]) { isSq = false; break; }
                }
                if (isSq) possibleOpsPerStep[k].push_back({'p', 2});

                if (possibleOpsPerStep[k].empty()) {
                    cyclePossible = false;
                    break;
                }
            }

            if (cyclePossible) {
                // Check if we have enough data to confirm cycle
                bool confirmed = false;
                if (L <= N/2) confirmed = true;
                else if (N > L + 1) confirmed = true; 

                if (confirmed) {
                    std::vector<std::vector<Op>> allCombinations;
                    std::vector<Op> currentOps(L);
                    generateCombinations(L, 0, currentOps, possibleOpsPerStep, allCombinations);

                    for (const auto& cycleOps : allCombinations) {
                        size_t transIdx = N - 1;
                        size_t opIdx = transIdx % L;
                        Op op = cycleOps[opIdx];
                        BigRat nextVal;
                        
                        if (op.type == 'a') nextVal = series.back() + op.val;
                        else if (op.type == 'm') nextVal = series.back() * op.val;
                        else if (op.type == 'p') nextVal = series.back() * series.back();
                        
                        std::stringstream opDesc;
                        opDesc << "Pattern: ";
                        for (size_t i = 0; i < L && i < N - 1; ++i) {
                            if (i > 0) opDesc << ", ";
                            opDesc << toString(series[i]) << "→" << toString(series[i+1]);
                            const Op& o = cycleOps[i];
                            if (o.type == 'a') {
                                if (o.val >= 0) opDesc << " (+" << toString(o.val) << ")";
                                else opDesc << " (" << toString(o.val) << ")";
                            } else if (o.type == 'm') {
                                opDesc << " (×" << toString(o.val) << ")";
                            } else if (o.type == 'p') {
                                opDesc << " (^2)";
                            }
                        }
                        opDesc << ", repeating every " << L << " steps";
                        
                        results.push_back({nextVal, "Cyclic Operations (Length " + std::to_string(L) + ")", "", 0.99, "", "", opDesc.str()});
                    }
                }
            }
        }
        return results;
    }
};

class RobustSolverWrapper : public SequenceSolver {
    std::vector<std::unique_ptr<SequenceSolver>> layers;
public:
    RobustSolverWrapper() {
        layers.push_back(std::make_unique<PolynomialSolver>());
        layers.push_back(std::make_unique<PowerSolver>());
        layers.push_back(std::make_unique<RecursiveSolver>()); 
        layers.push_back(std::make_unique<PrimeSolver>()); 
        layers.push_back(std::make_unique<BerlekampMasseySolver>());
    }
    std::string getName() const override { return "RobustWrapper"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        std::vector<SolverResult> all;
        for (const auto& solver : layers) {
            auto res = solver->solve(series);
            all.insert(all.end(), res.begin(), res.end());
        }
        std::sort(all.begin(), all.end(), compareResults);
        return all; // Return all candidates
    }
};

class InterleavedSolver : public SequenceSolver {
    std::unique_ptr<RobustSolverWrapper> sub;
public:
    InterleavedSolver() { sub = std::make_unique<RobustSolverWrapper>(); }
    std::string getName() const override { return "Interleaved"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 4) return {}; 
        std::vector<BigRat> evens, odds;
        for(size_t i=0; i<series.size(); ++i) {
            if (i % 2 == 0) evens.push_back(series[i]); else odds.push_back(series[i]);
        }
        auto rEven = sub->solve(evens);
        auto rOdd = sub->solve(odds);
        
        if (!rEven.empty() && !rOdd.empty()) {
            auto bestEven = rEven[0];
            auto bestOdd = rOdd[0];
            if (bestEven.confidence > 0.5 && bestOdd.confidence > 0.5) {
                double avgConf = (bestEven.confidence + bestOdd.confidence) / 2.0;
                if (avgConf > 0.9) avgConf = 0.98; 
                
                // [PARITY FIX]
                // If size is 6 (Indices 0..5), next index is 6 (Even). We want rEven.
                bool nextIsEven = (series.size() % 2 == 0);
                return {{nextIsEven ? bestEven.nextValue : bestOdd.nextValue, "Interleaved Sequence", "", avgConf, "", "", "Two interleaved sequences found. Even indices: " + bestEven.formula + ". Odd indices: " + bestOdd.formula}};
            }
        }
        return {};
    }
};

class WeightedSolver : public SequenceSolver {
    std::unique_ptr<RobustSolverWrapper> sub;
public:
    WeightedSolver() { sub = std::make_unique<RobustSolverWrapper>(); }
    std::string getName() const override { return "Weighted Solver"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 3) return {};
        
        for (int p = 1; p <= 3; ++p) {
             std::vector<BigRat> wSeries;
             for (size_t i = 0; i < series.size(); ++i) {
                 BigInt n = static_cast<unsigned long>(i + 1); 
                 BigInt w = 1;
                 for(int k=0; k<p; ++k) w *= n;
                 wSeries.push_back(series[i] * BigRat(w));
             }
             auto res = sub->solve(wSeries);
             if (!res.empty() && res[0].confidence > 0.8) {
                  BigInt nextN = static_cast<unsigned long>(series.size() + 1); 
                  BigInt w = 1;
                  for(int k=0; k<p; ++k) w *= nextN;
                  return {{res[0].nextValue / BigRat(w), "Weighted Pattern (n^" + std::to_string(p) + ")", "", res[0].confidence, "", "", "Sequence weighted by n^" + std::to_string(p) + ". Underlying pattern: " + res[0].formula}};
             }
        }
        return {};
    }
};

class DifferenceSolver : public SequenceSolver {
    std::unique_ptr<InterleavedSolver> inter;
    std::unique_ptr<PolynomialSolver> poly;
public:
    DifferenceSolver() { 
        inter = std::make_unique<InterleavedSolver>(); 
        poly = std::make_unique<PolynomialSolver>();
    }
    std::string getName() const override { return "Difference Solver"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.size() < 4) return {};
        std::vector<BigRat> diffs;
        for(size_t i=0; i<series.size()-1; ++i) diffs.push_back(series[i+1]-series[i]);

        auto diffRes = inter->solve(diffs);
        if (!diffRes.empty() && diffRes[0].confidence > 0.8) {
             return {{series.back() + diffRes[0].nextValue, "Interleaved Differences", "", 0.98, "", "", "Differences between terms follow an interleaved pattern."}};
        }
        
        auto polyRes = poly->solve(diffs);
        if (!polyRes.empty() && polyRes[0].confidence > 0.8) {
             return {{series.back() + polyRes[0].nextValue, "Polynomial Differences", "", 0.98, "", "", "Differences between terms follow a polynomial pattern."}};
        }

        return {};
    }
};

class FractionSolver : public SequenceSolver {
    std::unique_ptr<RobustSolverWrapper> sub;
    std::unique_ptr<InterleavedSolver> inter;
    std::unique_ptr<WeightedSolver> weight;
public:
    FractionSolver() { 
        sub = std::make_unique<RobustSolverWrapper>();
        inter = std::make_unique<InterleavedSolver>();
        weight = std::make_unique<WeightedSolver>();
    }
    std::string getName() const override { return "Fraction Splitter"; }
    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        bool hasFrac = false;
        std::vector<BigRat> nums, dens;
        for(const auto& s : series) { 
            if(s.get_den() != 1) hasFrac = true;
            nums.push_back(BigRat(s.get_num())); 
            dens.push_back(BigRat(s.get_den())); 
        }

        if (!hasFrac) return {};

        auto solvePart = [&](const std::vector<BigRat>& p) {
            auto r = sub->solve(p);
            if (!r.empty() && r[0].confidence > 0.5) return r;
            return inter->solve(p); 
        };
        auto rNum = solvePart(nums);
        auto rDen = solvePart(dens);
        
        std::vector<SolverResult> results;

        if (!rNum.empty() && !rDen.empty()) {
            double conf = (rNum[0].confidence + rDen[0].confidence)/2;
            if (hasFrac) conf = std::min(0.99, conf + 0.1);
            
            SolverResult finalRes = {rNum[0].nextValue / rDen[0].nextValue, "Fraction Pattern", "", conf, "", "", "Numerator pattern: " + rNum[0].formula + ". Denominator pattern: " + rDen[0].formula};
            finalRes.unsimplifiedNum = toString(rNum[0].nextValue.get_num());
            finalRes.unsimplifiedDen = toString(rDen[0].nextValue.get_num());
            results.push_back(finalRes);
        }
        
        auto wRes = weight->solve(series);
        if (!wRes.empty()) {
            if (results.empty() || wRes[0].confidence > results[0].confidence) {
                results.insert(results.begin(), wRes[0]);
            } else {
                results.push_back(wRes[0]);
            }
        }
        
        return results;
    }
};

// ============================================================================
// 4. ORCHESTRATOR
// ============================================================================

class GrandUnifiedSolver : public SequenceSolver {
    std::vector<std::unique_ptr<SequenceSolver>> layers;
public:
    GrandUnifiedSolver() {
        layers.push_back(std::make_unique<DigitSolver>());
        layers.push_back(std::make_unique<StandardMathSolver>());
        layers.push_back(std::make_unique<RecursiveSolver>()); 

        layers.push_back(std::make_unique<InterleavedSolver>()); 
        layers.push_back(std::make_unique<DifferenceSolver>()); 
        layers.push_back(std::make_unique<CyclicOpSolver>()); 
        
        layers.push_back(std::make_unique<FractionSolver>()); 

        layers.push_back(std::make_unique<PolynomialSolver>()); 
        layers.push_back(std::make_unique<PowerSolver>());
        layers.push_back(std::make_unique<PrimeSolver>()); 
        layers.push_back(std::make_unique<HolonomicSolver>()); 
        layers.push_back(std::make_unique<BerlekampMasseySolver>());
    }

    std::string getName() const override { return "GUSS Orchestrator"; }

    std::vector<SolverResult> solve(const std::vector<BigRat>& series) override {
        if (series.empty()) return {};
        
        std::vector<SolverResult> candidates;
        
        for (const auto& solver : layers) {
            try {
                std::vector<SolverResult> res = solver->solve(series);
                candidates.insert(candidates.end(), res.begin(), res.end());
            } catch (...) { continue; }
        }
        
        std::sort(candidates.begin(), candidates.end(), compareResults);
        
        return candidates;
    }
};

#endif // SEQUENCE_SOLVERS_H