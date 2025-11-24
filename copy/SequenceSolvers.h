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

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

inline bool isInteger(double d) {
    return std::abs(d - std::round(d)) < 1e-9;
}

inline bool areEqual(double a, double b) {
    return std::abs(a - b) < 1e-9;
}

inline int gcd(int a, int b) {
    while (b) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}

inline std::string toFractionString(double d) {
    for(int den=1; den<=10000; ++den) {
        int num = std::round(d * den);
        if (std::abs(d - (double)num/den) < 1e-9) {
            int g = gcd(std::abs(num), den);
            if (den/g == 1) return std::to_string(num/g);
            return std::to_string(num/g) + "/" + std::to_string(den/g);
        }
    }
    return "";
}

inline std::string addFractions(std::string f1, std::string f2) {
    if (f1.empty() || f2.empty()) return "";
    long long n1, d1, n2, d2;
    char slash;
    std::stringstream ss1(f1), ss2(f2);
    ss1 >> n1 >> slash >> d1;
    ss2 >> n2 >> slash >> d2;
    
    long long num = n1 * d2 + n2 * d1;
    long long den = d1 * d2;
    long long g = gcd(std::abs(num), den);
    return std::to_string(num/g) + "/" + std::to_string(den/g);
}

// ============================================================================
// RESULT STRUCTURE
// ============================================================================

struct SolverResult {
    double nextValue;
    std::string formula;
    std::string formattedOutput;
};

// ============================================================================
// BASE SOLVER CLASS
// ============================================================================

class SequenceSolver {
public:
    virtual ~SequenceSolver() = default;
    virtual std::string getName() const = 0;
    virtual SolverResult solve(const std::vector<double>& series) = 0;
};

// ============================================================================
// POLYNOMIAL SOLVER
// ============================================================================

class PolynomialSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Polynomial Series"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};

        std::vector<std::vector<double>> differences;
        differences.push_back(series);

        int depth = 0;
        while (true) {
            std::vector<double>& current = differences.back();
            if (current.empty()) break;

            bool allEqual = true;
            if (current.size() < 2) {
                allEqual = false;
            } else {
                for (size_t i = 1; i < current.size(); ++i) {
                    if (!areEqual(current[i], current[0])) {
                        allEqual = false;
                        break;
                    }
                }
            }

            if (allEqual) {
                double nextVal = current[0];
                for (int i = differences.size() - 2; i >= 0; --i) {
                    nextVal += differences[i].back();
                }
               
                std::string formula;
                if (depth == 1) formula = "Linear (Arithmetic)";
                else if (depth == 2) formula = "Quadratic";
                else if (depth == 3) formula = "Cubic";
                else formula = "Polynomial of degree " + std::to_string(depth);

                return {nextVal, formula, ""};
            }
            
            if (current.size() < 2) break;

            std::vector<double> nextDiff;
            for (size_t i = 0; i < current.size() - 1; ++i) {
                nextDiff.push_back(current[i+1] - current[i]);
            }
            differences.push_back(nextDiff);
            depth++;
            
            if (depth > 10) break;
        }

        return {0, "", ""};
    }
};

// ============================================================================
// LINEAR RECURRENCE SOLVER
// ============================================================================

class LinearRecurrenceSolver : public SequenceSolver {
    int maxOrder;
    int minExcess;

public:
    LinearRecurrenceSolver(int maxOrder = 4, int minExcess = 0) : maxOrder(maxOrder), minExcess(minExcess) {}

    std::string getName() const override {
        return maxOrder <= 2 ? "Linear Recurrence (Simple)" : "Linear Recurrence";
    }

    SolverResult solve(const std::vector<double>& series) override {
        int n = series.size();
        for (int k = 1; k <= maxOrder && k < n; ++k) {
            if (n < 2 * k + minExcess) continue;

            std::vector<std::vector<double>> A(k, std::vector<double>(k));
            std::vector<double> B(k);

            for (int i = 0; i < k; ++i) {
                B[i] = series[k + i];
                for (int j = 0; j < k; ++j) {
                    A[i][j] = series[k + i - 1 - j];
                }
            }

            std::vector<double> C = gaussianElimination(A, B);
            if (C.empty()) continue;

            // Require integer coefficients for integer inputs
            bool inputsAreIntegers = true;
            for(double d : series) if(!isInteger(d)) { inputsAreIntegers = false; break; }
            
            if (inputsAreIntegers) {
                 bool coeffsAreIntegers = true;
                 for (double c : C) {
                    if (std::abs(c - std::round(c)) > 1e-4) {
                        coeffsAreIntegers = false;
                        break;
                    }
                }
                if (!coeffsAreIntegers) continue;
            }

            bool matches = true;
            for (int i = k; i < n; ++i) {
                double predicted = 0;
                for (int j = 0; j < k; ++j) {
                    predicted += C[j] * series[i - 1 - j];
                }
                if (!areEqual(predicted, series[i])) {
                    matches = false;
                    break;
                }
            }

            if (matches) {
                double nextVal = 0;
                std::string formula = "a(n) = ";
                for (int j = 0; j < k; ++j) {
                    nextVal += C[j] * series[n - 1 - j];
                    if (!areEqual(C[j], 0)) {
                        if (j > 0 && C[j] > 0) formula += " + ";
                        if (C[j] < 0) formula += " - ";
                        if (!areEqual(std::abs(C[j]), 1)) formula += std::to_string(std::abs(C[j])) + "*";
                        formula += "a(n-" + std::to_string(j + 1) + ")";
                    }
                }
                return {nextVal, formula, ""};
            }
        }
        return {0, "", ""};
    }

private:
    std::vector<double> gaussianElimination(std::vector<std::vector<double>> A, std::vector<double> B) {
        int n = A.size();
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(A[j][i]) > std::abs(A[pivot][i])) pivot = j;
            }
            std::swap(A[i], A[pivot]);
            std::swap(B[i], B[pivot]);

            if (std::abs(A[i][i]) < 1e-9) return {};

            for (int j = i + 1; j < n; ++j) {
                double factor = A[j][i] / A[i][i];
                B[j] -= factor * B[i];
                for (int k = i; k < n; ++k) {
                    A[j][k] -= factor * A[i][k];
                }
            }
        }

        std::vector<double> X(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < n; ++j) {
                sum += A[i][j] * X[j];
            }
            X[i] = (B[i] - sum) / A[i][i];
        }
        return X;
    }
};

// ============================================================================
// DIGIT SOLVER
// ============================================================================

class DigitSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Digit Pattern"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        
        for(double d : series) if(!isInteger(d)) return {0, "", ""};

        std::vector<double> reversedSeries;
        for (double num : series) {
            reversedSeries.push_back(reverseInt((int)std::round(num)));
        }
        
        PolynomialSolver polySolver;
        auto revResult = polySolver.solve(reversedSeries);
        if (!revResult.formula.empty() && isInteger(revResult.nextValue)) {
            return {(double)reverseInt((int)std::round(revResult.nextValue)), "Reverse of " + revResult.formula, ""};
        }

        return {0, "", ""};
    }

private:
    int reverseInt(int n) {
        int r = 0;
        while (n > 0) {
            r = r * 10 + n % 10;
            n /= 10;
        }
        return r;
    }
};

// ============================================================================
// PRIME SOLVER
// ============================================================================

class PrimeSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Prime Numbers"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};

        std::vector<int> ints;
        std::vector<int> signs;
        bool allNeg = true;
        bool allPos = true;
        bool alternating = true;
        
        for(size_t i=0; i<series.size(); ++i) {
            double d = series[i];
            if (!isInteger(d)) return {0, "", ""};
            int n = std::round(d);
            if (!isPrime(n)) return {0, "", ""};
            
            int absN = std::abs(n);
            ints.push_back(absN);
            
            int sign = (n >= 0) ? 1 : -1;
            signs.push_back(sign);
            
            if (sign > 0) allNeg = false;
            if (sign < 0) allPos = false;
            if (i > 0 && signs[i] == signs[i-1]) alternating = false;
        }
        
        int maxVal = *std::max_element(ints.begin(), ints.end());
        std::vector<int> primes = generatePrimes(maxVal + 100);
        
        std::vector<int> indices;
        for(int n : ints) {
            auto it = std::find(primes.begin(), primes.end(), n);
            if (it == primes.end()) return {0, "", ""};
            indices.push_back(std::distance(primes.begin(), it));
        }
        
        PolynomialSolver poly;
        std::vector<double> indicesD(indices.begin(), indices.end());
        auto res = poly.solve(indicesD);
        if (!res.formula.empty()) {
            int nextIndex = std::round(res.nextValue);
            if (nextIndex < 0) return {0, "", ""};
            if (nextIndex >= primes.size()) primes = generatePrimes(maxVal * 2 + nextIndex * 10);
            if (nextIndex >= primes.size()) return {0, "", ""};
            
            double nextVal = (double)primes[nextIndex];
            
            if (allNeg) nextVal = -nextVal;
            else if (alternating) {
                int nextSign = (signs.back() == 1) ? -1 : 1;
                nextVal = nextVal * nextSign;
            } else if (!allPos) {
                return {0, "", ""};
            }
            
            return {nextVal, "Primes with index pattern: " + res.formula, ""};
        }
        
        return {0, "", ""};
    }

private:
    bool isPrime(int n) {
        n = std::abs(n);
        if (n <= 1) return false;
        for (int i = 2; i * i <= n; ++i) if (n % i == 0) return false;
        return true;
    }

    std::vector<int> generatePrimes(int limit) {
        std::vector<int> p;
        if (limit < 2) return p;
        std::vector<bool> sieve(limit + 1, true);
        sieve[0] = sieve[1] = false;
        for (int i = 2; i * i <= limit; ++i) {
            if (sieve[i]) {
                for (int j = i * i; j <= limit; j += i) sieve[j] = false;
            }
        }
        for (int i = 0; i <= limit; ++i) if (sieve[i]) p.push_back(i);
        return p;
    }
};

// ============================================================================
// FRACTION SOLVER
// ============================================================================

class FractionSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Fraction Pattern"; }

    SolverResult solve(const std::vector<double>& series) override {
        std::vector<int> nums, dens;
        for (double d : series) {
            int bestN = 0, bestD = 1;
            bool found = false;
            for(int den=1; den<=1000; ++den) {
                int num = std::round(d * den);
                if (std::abs(d - (double)num/den) < 1e-6) {
                    int g = gcd(num, den);
                    nums.push_back(num/g);
                    dens.push_back(den/g);
                    found = true;
                    break;
                }
            }
            if (!found) return {0, "", ""};
        }

        bool allIntegers = true;
        for(int d : dens) if(d != 1) { allIntegers = false; break; }
        if (allIntegers) return {0, "", ""};

        PolynomialSolver poly;
        LinearRecurrenceSolver lin;
        
        std::vector<double> numsD(nums.begin(), nums.end());
        std::vector<double> densD(dens.begin(), dens.end());

        auto resNum = poly.solve(numsD);
        if (resNum.formula.empty()) resNum = lin.solve(numsD);

        auto resDen = poly.solve(densD);
        if (resDen.formula.empty()) resDen = lin.solve(densD);

        if (!resNum.formula.empty() && !resDen.formula.empty()) {
            double nextVal = resNum.nextValue / resDen.nextValue;
            long long numInt = std::round(resNum.nextValue);
            long long denInt = std::round(resDen.nextValue);
            std::string formatted;
            if (denInt == 1 || (denInt != 0 && numInt % denInt == 0)) {
                formatted = std::to_string(numInt / denInt);
            } else {
                formatted = std::to_string(numInt) + "/" + std::to_string(denInt);
            }
            return {nextVal, "Fraction: Num(" + resNum.formula + ") / Den(" + resDen.formula + ")", formatted};
        }

        return {0, "", ""};
    }
};

// Forward declarations for complex solvers
class PositionDependentSolver;
class AlternatingArithmeticSolver;
class CumulativeDifferenceSolver;
class DifferencePatternSolver;
class OperationCycleSolver;
class ConditionalRecurrenceSolver;
class DependentInterleavedSolver;
class MultiplierSolver;
class AffineRecurrenceSolver;
class FractionIntegerInterleavedSolver;
class PowerSolver;
class ComplexFractionSolver;
class RootSolver;
class SkipPositionSolver;
class InterleavedSolver;

// ============================================================================
// POSITION DEPENDENT SOLVER
// ============================================================================

class PositionDependentSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Position Dependent"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        int n = series.size();
        
        bool match = true;
        for (size_t i = 1; i < series.size(); ++i) {
            double pos_prev = (double)i;
            double expected = series[i-1] * (pos_prev / 2.0) + (pos_prev - 1.0);
            if (!areEqual(series[i], expected)) {
                match = false;
                break;
            }
        }

        if (match) {
            double pos_last = (double)n;
            double nextVal = series.back() * (pos_last / 2.0) + (pos_last - 1.0);
            return {nextVal, "a(n) = a(n-1) * (n-1)/2 + (n-2)", ""};
        }

        return {0, "", ""};
    }
};

// ============================================================================
// ALTERNATING ARITHMETIC SOLVER
// ============================================================================

class AlternatingArithmeticSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Alternating Arithmetic"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        
        std::vector<double> absValues;
        for (double d : series) absValues.push_back(std::abs(d));
        
        PolynomialSolver poly;
        auto res = poly.solve(absValues);
        
        if (!res.formula.empty()) {
            bool alternating = true;
            for (size_t i = 0; i < series.size() - 1; ++i) {
                if (series[i] * series[i+1] >= 0) {
                    alternating = false;
                    break;
                }
            }
            
            if (alternating) {
                double nextAbs = res.nextValue;
                double nextVal = (series.back() > 0) ? -nextAbs : nextAbs;
                return {nextVal, "Alternating " + res.formula, ""};
            }
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// CUMULATIVE DIFFERENCE SOLVER
// ============================================================================

class CumulativeDifferenceSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Cumulative Diffs"; }
   
    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        
        std::vector<double> diffs;
        for (size_t i = 0; i < series.size() - 1; ++i) {
            diffs.push_back(series[i+1] - series[i]);
        }
        
        if (diffs.size() < 2) return {0, "", ""};
        
        bool alternating = true;
        for (size_t i = 0; i < diffs.size() - 1; ++i) {
            if (diffs[i] * diffs[i+1] >= 0) {
                alternating = false;
                break;
            }
        }
        
        if (alternating) {
            bool matchesNOverNPlus1 = true;
            for (size_t i = 0; i < diffs.size(); ++i) {
                double expected = (double)(i + 1) / (double)(i + 2);
                if (!areEqual(std::abs(diffs[i]), expected)) {
                    matchesNOverNPlus1 = false;
                    break;
                }
            }
            
            if (matchesNOverNPlus1) {
                size_t n = diffs.size() + 1;
                double nextAbsDiff = (double)n / (double)(n + 1);
                double nextDiff = (diffs.back() > 0) ? -nextAbsDiff : nextAbsDiff;
                double nextVal = series.back() + nextDiff;
                
                std::string lastValStr = toFractionString(series.back());
                std::string nextDiffStr = toFractionString(nextDiff);
                std::string formatted = addFractions(lastValStr, nextDiffStr);
                
                return {nextVal, "Cumulative Diffs (n/(n+1) alternating)", formatted};
            }
        }
        
        std::vector<double> absDiffs;
        for (double d : diffs) absDiffs.push_back(std::abs(d));
        
        PolynomialSolver poly;
        auto res = poly.solve(absDiffs);
        
        if (res.formula.empty() && absDiffs.size() == 2) {
            double diff_of_diffs = absDiffs[1] - absDiffs[0];
            double nextAbsDiff = absDiffs.back() + diff_of_diffs;
            res = {nextAbsDiff, "Linear (Arithmetic)", ""};
        }
        
        if (!res.formula.empty()) {
            bool alternating = true;
            for (size_t i = 0; i < diffs.size() - 1; ++i) {
                if (diffs[i] * diffs[i+1] >= 0) {
                    alternating = false;
                    break;
                }
            }
            
            if (alternating && diffs.size() > 0) {
                double nextAbsDiff = res.nextValue;
                double nextDiff = (diffs.back() > 0) ? -nextAbsDiff : nextAbsDiff;
                double nextVal = series.back() + nextDiff;
                return {nextVal, "Cumulative Diffs (Alternating " + res.formula + ")", ""};
            }
        }
        
        auto resPoly = poly.solve(diffs);
        if (!resPoly.formula.empty()) {
            return {series.back() + resPoly.nextValue, "Cumulative Diffs (" + resPoly.formula + ")", ""};
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// DIFFERENCE PATTERN SOLVER
// ============================================================================

class DifferencePatternSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Difference Pattern"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 4) return {0, "", ""};
        
        std::vector<double> diffs;
        for (size_t i = 0; i < series.size() - 1; ++i) {
            diffs.push_back(series[i+1] - series[i]);
        }

        std::vector<std::unique_ptr<SequenceSolver>> subSolvers;
        subSolvers.push_back(std::make_unique<PrimeSolver>());
        subSolvers.push_back(std::make_unique<CumulativeDifferenceSolver>());
        subSolvers.push_back(std::make_unique<AlternatingArithmeticSolver>());
        subSolvers.push_back(std::make_unique<FractionSolver>());
        subSolvers.push_back(std::make_unique<LinearRecurrenceSolver>(2, 1));
        subSolvers.push_back(std::make_unique<PolynomialSolver>());
        
        for(const auto& solver : subSolvers) {
            auto res = solver->solve(diffs);
            if (!res.formula.empty()) {
                double nextDiff = res.nextValue;
                return {series.back() + nextDiff, "Differences follow: " + res.formula, ""};
            }
        }

        return {0, "", ""};
    }
};

// ============================================================================
// OPERATION CYCLE SOLVER
// ============================================================================

class OperationCycleSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Operation Cycle"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 4) return {0, "", ""};
        
        for (int k = 2; k <= 4; ++k) {
            if (series.size() <= k) continue;
            
            struct Op { char type; double val; };
            std::vector<Op> finalCycle(k);
            bool cycleFound = true;
            
            for (int j = 0; j < k; ++j) {
                double diff = series[j+1] - series[j];
                double ratio = (std::abs(series[j]) > 1e-9) ? series[j+1] / series[j] : 0;
                
                bool addConsistent = true;
                double addVal = series[j+1] - series[j];
                for (size_t i = j; i < series.size() - 1; i += k) {
                    if (!areEqual(series[i+1] - series[i], addVal)) {
                        addConsistent = false;
                        break;
                    }
                }
                
                bool mulConsistent = true;
                double mulVal = (std::abs(series[j]) > 1e-9) ? series[j+1] / series[j] : 0;
                for (size_t i = j; i < series.size() - 1; i += k) {
                    if (std::abs(series[i]) < 1e-9 || !areEqual(series[i+1] / series[i], mulVal)) {
                        mulConsistent = false;
                        break;
                    }
                }
                
                if (addConsistent && mulConsistent) {
                    auto getScore = [](char type, double val) -> double {
                        double score = 0;
                        if (type == '+') {
                            score = std::abs(val);
                            if (!isInteger(val)) score += 1000;
                        } else {
                            if (isInteger(val)) {
                                score = std::abs(val);
                            } else if (std::abs(val) > 1e-9 && isInteger(1.0/val)) {
                                score = std::abs(1.0/val);
                            } else {
                                score = std::abs(val) + 1000;
                            }
                        }
                        return score;
                    };
                    
                    double addScore = getScore('+', addVal);
                    double mulScore = getScore('*', mulVal);
                    
                    if (mulScore < addScore) {
                         finalCycle[j] = {'*', mulVal};
                    } else {
                         finalCycle[j] = {'+', addVal};
                    }
                } else if (addConsistent) {
                    finalCycle[j] = {'+', addVal};
                } else if (mulConsistent) {
                    finalCycle[j] = {'*', mulVal};
                } else {
                    cycleFound = false;
                    break;
                }
            }
            
            if (cycleFound) {
                int nextOpIdx = (series.size() - 1) % k;
                Op op = finalCycle[nextOpIdx];
                double nextVal = (op.type == '+') ? series.back() + op.val : series.back() * op.val;
                
                std::string formula = "Cycle(" + std::to_string(k) + "): ";
                for(int j=0; j<k; ++j) {
                    if (finalCycle[j].type == '+') {
                        if (finalCycle[j].val < 0) formula += "-" + toFractionString(std::abs(finalCycle[j].val));
                        else formula += "+" + toFractionString(finalCycle[j].val);
                    } else {
                        double inv = 1.0 / finalCycle[j].val;
                        if (isInteger(inv) && std::abs(inv) > 1) {
                             formula += "/" + toFractionString(inv);
                        } else {
                             formula += "*" + toFractionString(finalCycle[j].val);
                        }
                    }
                    formula += (j < k-1 ? ", " : "");
                }
                return {nextVal, formula, ""};
            }
        }

        return {0, "", ""};
    }
};

// ============================================================================
// CONDITIONAL RECURRENCE SOLVER
// ============================================================================

class ConditionalRecurrenceSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Conditional Recurrence"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 6) return {0, "", ""};
        
        bool match1 = true;
        for (size_t i = 3; i < series.size(); ++i) {
            if (i % 2 != 0) {
                if (!areEqual(series[i], series[i-1] + series[i-3])) {
                    match1 = false;
                    break;
                }
            } else {
                if (!areEqual(series[i], series[i-1] + series[i-2])) {
                    match1 = false;
                    break;
                }
            }
        }
        
        if (match1) {
            double nextVal;
            size_t nextIdx = series.size();
            if (nextIdx % 2 != 0) {
                 nextVal = series[nextIdx-1] + series[nextIdx-3];
            } else {
                 nextVal = series[nextIdx-1] + series[nextIdx-2];
            }
            return {nextVal, "Even pos: a(n-1)+a(n-3), Odd pos: a(n-1)+a(n-2)", ""};
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// DEPENDENT INTERLEAVED SOLVER
// ============================================================================

class DependentInterleavedSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Dependent Interleaved"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 4) return {0, "", ""};
        
        bool matchSub = true;
        for (size_t i = 2; i < series.size(); ++i) {
            if (!areEqual(series[i], series[i-2] - series[i-1])) {
                matchSub = false;
                break;
            }
        }
        if (matchSub) {
            return {series[series.size()-2] - series.back(), "a(n) = a(n-2) - a(n-1)", ""};
        }

        bool matchAdd = true;
        for (size_t i = 2; i < series.size(); ++i) {
            if (!areEqual(series[i], series[i-2] + series[i-1])) {
                matchAdd = false;
                break;
            }
        }
        if (matchAdd) {
            return {series[series.size()-2] + series.back(), "a(n) = a(n-2) + a(n-1)", ""};
        }

        bool matchSkipSub = true;
        for (size_t i = 2; i < series.size(); i += 2) {
            if (!areEqual(series[i], series[i-2] - series[i-1])) {
                matchSkipSub = false;
                break;
            }
        }
        
        if (matchSkipSub) {
            if (series.size() % 2 == 0) {
                return {series[series.size()-2] - series.back(), "a(2k) = a(2k-2) - a(2k-1)", ""};
            }
        }
        return {0, "", ""};
    }
};

// ============================================================================
// MULTIPLIER SOLVER
// ============================================================================

class MultiplierSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Multiplier Pattern"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};

        std::vector<double> ratios;
        for(size_t i=0; i<series.size()-1; ++i) {
            if (std::abs(series[i]) < 1e-9) return {0, "", ""};
            ratios.push_back(series[i+1] / series[i]);
        }

        bool allEqual = true;
        for(size_t i=1; i<ratios.size(); ++i) {
            if (!areEqual(ratios[i], ratios[0])) { allEqual = false; break; }
        }
        
        if (allEqual) {
            double ratio = ratios[0];
            return {series.back() * ratio, "a(n) = a(n-1) * " + std::to_string(ratio), ""};
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// AFFINE RECURRENCE SOLVER
// ============================================================================

class AffineRecurrenceSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Affine Recurrence"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 4) return {0, "", ""};
        
        double num = series[2] - series[1];
        double den = series[1] - series[0];
        if (std::abs(den) < 1e-9) return {0, "", ""};
        
        double m = num / den;
        double c = series[1] - m * series[0];
        
        bool match = true;
        for (size_t i = 1; i < series.size(); ++i) {
            if (!areEqual(series[i], m * series[i-1] + c)) {
                match = false;
                break;
            }
        }
        
        if (match) {
            return {m * series.back() + c, "a(n) = " + std::to_string(m) + "*a(n-1) + " + std::to_string(c), ""};
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// FRACTION-INTEGER INTERLEAVED SOLVER
// ============================================================================

class FractionIntegerInterleavedSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Fraction-Integer Interleaved"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 4) return {0, "", ""};
        
        std::vector<double> integers, nonIntegers;
        std::vector<size_t> intPositions, nonIntPositions;
        
        for (size_t i = 0; i < series.size(); ++i) {
            if (isInteger(series[i])) {
                integers.push_back(series[i]);
                intPositions.push_back(i);
            } else {
                nonIntegers.push_back(series[i]);
                nonIntPositions.push_back(i);
            }
        }
        
        if (integers.size() < 2 || nonIntegers.size() < 2) return {0, "", ""};
        
        bool alternating = true;
        for (size_t i = 0; i < intPositions.size() - 1; ++i) {
            if (intPositions[i+1] - intPositions[i] != 2) {
                alternating = false;
                break;
            }
        }
        if (alternating) {
            for (size_t i = 0; i < nonIntPositions.size() - 1; ++i) {
                if (nonIntPositions[i+1] - nonIntPositions[i] != 2) {
                    alternating = false;
                    break;
                }
            }
        }
        
        if (!alternating) return {0, "", ""};
        
        PolynomialSolver polySolver;
        auto intResult = polySolver.solve(integers);
        if (intResult.formula.empty()) return {0, "", ""};
        
        std::vector<double> numerators, denominators;
        for (double val : nonIntegers) {
            bool found = false;
            for (int den = 2; den <= 1000; ++den) {
                int num = std::round(val * den);
                if (std::abs(val - (double)num / den) < 1e-6) {
                    numerators.push_back(num);
                    denominators.push_back(den);
                    found = true;
                    break;
                }
            }
            if (!found) return {0, "", ""};
        }
        
        auto numResult = polySolver.solve(numerators);
        auto denResult = polySolver.solve(denominators);
        
        if (numResult.formula.empty() || denResult.formula.empty()) return {0, "", ""};
        
        size_t nextPos = series.size();
        bool nextIsInteger = (intPositions[0] % 2 == nextPos % 2);
        
        if (nextIsInteger) {
            return {intResult.nextValue,
                    "Integers: " + intResult.formula + "; Fractions: Num(" + numResult.formula + ")/Den(" + denResult.formula + ")",
                    ""};
        } else {
            double nextNum = numResult.nextValue;
            double nextDen = denResult.nextValue;
            double nextVal = nextNum / nextDen;
            
            long long numInt = std::round(nextNum);
            long long denInt = std::round(nextDen);
            std::string formatted;
            if (denInt == 1 || (denInt != 0 && numInt % denInt == 0)) {
                formatted = std::to_string(numInt / denInt);
            } else {
                formatted = std::to_string(numInt) + "/" + std::to_string(denInt);
            }
            
            return {nextVal,
                    "Integers: " + intResult.formula + "; Fractions: Num(" + numResult.formula + ")/Den(" + denResult.formula + ")",
                    formatted};
        }
    }
};

// ============================================================================
// POWER SOLVER
// ============================================================================

class PowerSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Power Sequence"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        
        std::vector<double> bases;
        std::vector<double> exps;
        
        bool allDecomposed = true;
        for (double val : series) {
            bool found = false;
            for (int b = 2; b <= 100; ++b) {
                double e = std::log(val) / std::log(b);
                if (std::abs(e - std::round(e)) < 1e-4) {
                    bases.push_back((double)b);
                    exps.push_back(std::round(e));
                    found = true;
                    break;
                }
            }
            if (!found) {
                allDecomposed = false;
                break;
            }
        }
        
        if (allDecomposed) {
            PolynomialSolver poly;
            MultiplierSolver mult;
            
            auto resBase = poly.solve(bases);
            if (resBase.formula.empty()) resBase = mult.solve(bases);

            auto resExp = poly.solve(exps);
            if (resExp.formula.empty()) resExp = mult.solve(exps);
            
            if (!resBase.formula.empty() && !resExp.formula.empty()) {
                double nextBase = resBase.nextValue;
                double nextExp = resExp.nextValue;
                std::string baseStr = std::to_string((long long)std::round(nextBase));
                std::string expStr = std::to_string((long long)std::round(nextExp));
                return {std::pow(nextBase, nextExp), "Base(" + resBase.formula + ")^Exp(" + resExp.formula + ")", baseStr + "^" + expStr};
            }
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// COMPLEX FRACTION SOLVER (Reconstructs unsimplified fractions)
// ============================================================================

class ComplexFractionSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Complex Fraction Pattern"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};

        std::vector<int> nums, dens;
        std::vector<int> knownIndices;
        std::vector<double> knownDens;
        
        for (size_t i = 0; i < series.size(); ++i) {
            double d = series[i];
            int bestN = 0, bestD = 1;
            bool found = false;
            for(int den=1; den<=1000; ++den) {
                int num = std::round(d * den);
                if (std::abs(d - (double)num/den) < 1e-6) {
                    int g = gcd(std::abs(num), den);
                    bestN = num/g;
                    bestD = den/g;
                    found = true;
                    break;
                }
            }
            if (!found) return {0, "", ""};
            nums.push_back(bestN);
            dens.push_back(bestD);
            if (bestD > 1) {
                knownIndices.push_back((int)i);
                knownDens.push_back((double)bestD);
            }
        }

        if (knownIndices.size() < 2) return {0, "", ""};

        std::vector<double> predictedDens(series.size());
        double predictedNextDen = 0;
        bool denPatternFound = false;
        std::string denFormula = "";

        // Check for Power Pattern in denominators
        if (knownDens.size() >= 3) {
            std::vector<double> bases;
            std::vector<double> exps;
            bool allPowers = true;
            for(double val : knownDens) {
                bool isPow = false;
                for(int b=2; b<=100; ++b) {
                    double e = std::log(val)/std::log(b);
                    if(std::abs(e - std::round(e)) < 1e-4) {
                        bases.push_back(b);
                        exps.push_back(std::round(e));
                        isPow = true;
                        break;
                    }
                }
                if(!isPow) { allPowers = false; break; }
            }
            if (allPowers) {
                PolynomialSolver poly;
                auto resBase = poly.solve(bases);
                bool constExp = true;
                for(size_t k=1; k<exps.size(); ++k) if(exps[k] != exps[0]) constExp = false;
                if (!resBase.formula.empty() && constExp) {
                    bool linearShift = true;
                    int shift = bases[0] - knownIndices[0];
                    for(size_t k=1; k<bases.size(); ++k) {
                        if (bases[k] - knownIndices[k] != shift) { linearShift = false; break; }
                    }
                    
                    if (linearShift) {
                        double exp = exps[0];
                        for(size_t i=0; i<series.size(); ++i) {
                            predictedDens[i] = std::pow(i + shift, exp);
                        }
                        predictedNextDen = std::pow(series.size() + shift, exp);
                        denPatternFound = true;
                        denFormula = "(n + " + std::to_string(shift) + ")^" + std::to_string((int)exp);
                    }
                }
            }
        }

        if (!denPatternFound) {
            for (size_t startIdx = 0; startIdx <= knownDens.size() - 2; ++startIdx) {
                std::vector<double> currentAnchors;
                std::vector<int> currentIndices;
                for(size_t i=startIdx; i<knownDens.size(); ++i) {
                    currentAnchors.push_back(knownDens[i]);
                    currentIndices.push_back(knownIndices[i]);
                }
                
                if (currentAnchors.size() < 2) break;
                bool equidistant = true;
                int step = currentIndices[1] - currentIndices[0];
                for (size_t i = 1; i < currentIndices.size() - 1; ++i) {
                    if (currentIndices[i+1] - currentIndices[i] != step) {
                        equidistant = false;
                        break;
                    }
                }
                
                if (!equidistant) continue;
                PolynomialSolver poly;
                LinearRecurrenceSolver linRec(2);
                DifferencePatternSolver diffSolver;
                
                auto resPoly = poly.solve(currentAnchors);
                auto resLin = linRec.solve(currentAnchors);
                auto resDiff = diffSolver.solve(currentAnchors);
                
                bool isFib = false;
                if (currentAnchors.size() >= 3) {
                     bool fibCheck = true;
                     for(size_t k=2; k<currentAnchors.size(); ++k) {
                         if (!areEqual(currentAnchors[k], currentAnchors[k-1] + currentAnchors[k-2])) {
                             fibCheck = false; break;
                         }
                     }
                     if (fibCheck) isFib = true;
                }
                bool usePoly = !resPoly.formula.empty();
                bool useLin = !resLin.formula.empty() || isFib;
                bool useDiff = !resDiff.formula.empty();
                
                if (usePoly || useLin || useDiff) {
                    int firstAnchorIdx = currentIndices[0];
                    
                    if (useLin || useDiff) {
                        if (currentAnchors.size() >= 2) {
                            if (step == 1) {
                                std::vector<double> extended = currentAnchors;
                                int currentStart = firstAnchorIdx;
                                
                                while (currentStart > 0) {
                                    double prev;
                                    if (isFib) prev = extended[1] - extended[0];
                                    else if (usePoly && resPoly.formula.find("Linear") != std::string::npos) prev = extended[0] - (extended[1]-extended[0]);
                                    else if (useDiff) break;
                                    else break;
                                    
                                    extended.insert(extended.begin(), prev);
                                    currentStart--;
                                }
                                
                                while (extended.size() <= series.size()) {
                                    double next;
                                    if (isFib) next = extended.back() + extended[extended.size()-2];
                                    else if (usePoly) {
                                        auto r = poly.solve(extended);
                                        next = r.nextValue;
                                    }
                                    else if (useDiff) {
                                        auto r = diffSolver.solve(extended);
                                        next = r.nextValue;
                                    }
                                    else break;
                                    extended.push_back(next);
                                }
                                
                                if (extended.size() > series.size()) {
                                    predictedDens = std::vector<double>(extended.begin(), extended.begin() + series.size());
                                    predictedNextDen = extended[series.size()];
                                    denPatternFound = true;
                                    if (isFib) denFormula = "Fibonacci-like";
                                    else if (useDiff) denFormula = resDiff.formula;
                                    else denFormula = resPoly.formula;
                                }
                            }
                        }
                    }
                    
                    if (!denPatternFound && usePoly) {
                        if (step == 1) {
                            std::vector<double> extended = currentAnchors;
                            int currentStart = firstAnchorIdx;
                            
                            bool possible = true;
                            while (currentStart > 0) {
                                if (extended.size() < 2) { possible = false; break; }
                                
                                std::vector<std::vector<double>> diffs;
                                diffs.push_back(extended);
                                int depth = 0;
                                while(diffs.back().size() > 1 && depth < 10) {
                                    std::vector<double> nextD;
                                    bool allEq = true;
                                    for(size_t k=0; k<diffs.back().size()-1; ++k) {
                                        nextD.push_back(diffs.back()[k+1] - diffs.back()[k]);
                                        if (k > 0 && !areEqual(nextD.back(), nextD[k-1])) allEq = false;
                                    }
                                    diffs.push_back(nextD);
                                    if (allEq) break;
                                    depth++;
                                }
                                
                                double prev = diffs.back()[0];
                                for(int d=diffs.size()-2; d>=0; --d) {
                                    prev = diffs[d][0] - prev;
                                }
                                extended.insert(extended.begin(), prev);
                                currentStart--;
                            }
                            
                            if (possible) {
                                 while (extended.size() <= series.size()) {
                                    auto r = poly.solve(extended);
                                    extended.push_back(r.nextValue);
                                 }
                                 
                                 predictedDens = std::vector<double>(extended.begin(), extended.begin() + series.size());
                                 predictedNextDen = extended[series.size()];
                                 denPatternFound = true;
                                 denFormula = resPoly.formula;
                            }
                        }
                    }
                }
                
                if (denPatternFound) {
                    bool match = true;
                    for(size_t i=0; i<series.size(); ++i) {
                        double val = series[i] * predictedDens[i];
                        if (!isInteger(val)) {
                            match = false;
                            break;
                        }
                    }
                    
                    if (match) break;
                    else {
                        denPatternFound = false;
                    }
                }
            }
        }

        if (!denPatternFound) return {0, "", ""};

        std::vector<double> reconstructedNums;
        for(size_t i=0; i<series.size(); ++i) {
            double num = series[i] * predictedDens[i];
            reconstructedNums.push_back(std::round(num));
        }

        std::vector<std::unique_ptr<SequenceSolver>> numSolvers;
        numSolvers.push_back(std::make_unique<OperationCycleSolver>());
        numSolvers.push_back(std::make_unique<PolynomialSolver>());
        numSolvers.push_back(std::make_unique<LinearRecurrenceSolver>());
        numSolvers.push_back(std::make_unique<AlternatingArithmeticSolver>());
        numSolvers.push_back(std::make_unique<PowerSolver>());
        numSolvers.push_back(std::make_unique<DifferencePatternSolver>());
        
        for(const auto& solver : numSolvers) {
            auto res = solver->solve(reconstructedNums);
            if (!res.formula.empty()) {
                double nextNum = res.nextValue;
                double nextDen = predictedNextDen;
                
                double nextVal = nextNum / nextDen;
                
                long long numInt = std::round(nextNum);
                long long denInt = std::round(nextDen);
                std::string formatted;
                if (denInt == 1 || (denInt != 0 && numInt % denInt == 0)) {
                    formatted = std::to_string(numInt / denInt);
                } else {
                    formatted = std::to_string(numInt) + "/" + std::to_string(denInt);
                }
                
                return {nextVal, "Complex Fraction: Num(" + res.formula + ") / Den(" + denFormula + ")", formatted};
            }
        }
        return {0, "", ""};
    }
};

// ============================================================================
// ROOT SOLVER
// ============================================================================

class RootSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Root Sequence"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        
        bool match = true;
        for (size_t i = 1; i < series.size(); ++i) {
            if (series[i-1] < 0) { match = false; break; }
            if (!areEqual(series[i], std::sqrt(series[i-1]))) {
                match = false;
                break;
            }
        }
        if (match) {
            return {std::sqrt(series.back()), "a(n) = sqrt(a(n-1))", ""};
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// SKIP POSITION SOLVER
// ============================================================================

class SkipPositionSolver : public SequenceSolver {
private:
    std::vector<std::unique_ptr<SequenceSolver>> internalSolvers;

public:
    SkipPositionSolver() {
        internalSolvers.push_back(std::make_unique<MultiplierSolver>());
        internalSolvers.push_back(std::make_unique<DifferencePatternSolver>());
        internalSolvers.push_back(std::make_unique<OperationCycleSolver>());
        internalSolvers.push_back(std::make_unique<PolynomialSolver>());
        internalSolvers.push_back(std::make_unique<AffineRecurrenceSolver>());
    }

    std::string getName() const override { return "Skip Position Pattern"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 4) return {0, "", ""};
        
        std::vector<double> evenPos;
        for (size_t i = 0; i < series.size(); i += 2) {
            evenPos.push_back(series[i]);
        }
        
        if (evenPos.size() >= 3) {
            for (const auto& solver : internalSolvers) {
                auto result = solver->solve(evenPos);
                if (!result.formula.empty()) {
                    if (evenPos.size() >= 3) {
                        if (series.size() % 2 == 0) {
                            return {result.nextValue, "Even positions only: " + result.formula, result.formattedOutput};
                        }
                    }
                }
            }
        }
        
        std::vector<double> oddPos;
        for (size_t i = 1; i < series.size(); i += 2) {
            oddPos.push_back(series[i]);
        }
        
        if (oddPos.size() >= 3) {
            for (const auto& solver : internalSolvers) {
                auto result = solver->solve(oddPos);
                if (!result.formula.empty()) {
                    if (oddPos.size() >= 3) {
                        if (series.size() % 2 == 1) {
                            return {result.nextValue, "Odd positions only: " + result.formula, result.formattedOutput};
                        }
                    }
                }
            }
        }
        
        return {0, "", ""};
    }
};

// ============================================================================
// INTERLEAVED SOLVER
// ============================================================================

class InterleavedSolver : public SequenceSolver {
private:
    std::vector<std::unique_ptr<SequenceSolver>> internalSolvers;

public:
    InterleavedSolver() {
        internalSolvers.push_back(std::make_unique<MultiplierSolver>());
        internalSolvers.push_back(std::make_unique<DifferencePatternSolver>());
        internalSolvers.push_back(std::make_unique<PositionDependentSolver>());
        internalSolvers.push_back(std::make_unique<OperationCycleSolver>());
        internalSolvers.push_back(std::make_unique<AlternatingArithmeticSolver>());
        internalSolvers.push_back(std::make_unique<PowerSolver>());
        internalSolvers.push_back(std::make_unique<CumulativeDifferenceSolver>());
        internalSolvers.push_back(std::make_unique<AffineRecurrenceSolver>());
        internalSolvers.push_back(std::make_unique<RootSolver>());
        internalSolvers.push_back(std::make_unique<DependentInterleavedSolver>());
        internalSolvers.push_back(std::make_unique<FractionSolver>());
        internalSolvers.push_back(std::make_unique<ComplexFractionSolver>());
        internalSolvers.push_back(std::make_unique<PolynomialSolver>());
        internalSolvers.push_back(std::make_unique<LinearRecurrenceSolver>());
    }

    std::string getName() const override { return "Interleaved Series"; }

    SolverResult solve(const std::vector<double>& series) override {
        int n = series.size();
        for (int k = 2; k <= 4; ++k) {
            if (n < 2 * k) continue;

            std::vector<std::vector<double>> subSeries(k);
            for (int i = 0; i < n; ++i) {
                subSeries[i % k].push_back(series[i]);
            }

            bool allSolved = true;
            std::vector<double> nextValues(k);
            std::vector<std::string> nextFormatted(k);
            std::string combinedFormula = "Interleaved(" + std::to_string(k) + "): ";

            for (int i = 0; i < k; ++i) {
                bool subSolved = false;
                for (const auto& solver : internalSolvers) {
                    auto result = solver->solve(subSeries[i]);
                    if (!result.formula.empty()) {
                        if (subSeries[i].size() < 3) {
                            if (result.formula.find("Linear") != std::string::npos ||
                                result.formula.find("Arithmetic") != std::string::npos) {
                                continue;
                            }
                        }
                        
                        nextValues[i] = result.nextValue;
                        nextFormatted[i] = result.formattedOutput;
                        combinedFormula += "[" + std::to_string(i+1) + "]: " + result.formula + "; ";
                        subSolved = true;
                        break;
                    }
                }
                if (!subSolved) {
                    allSolved = false;
                    break;
                }
            }

            if (allSolved) {
                return {nextValues[n % k], combinedFormula, nextFormatted[n % k]};
            }
        }
        return {0, "", ""};
    }
};


// ============================================================================
// FACTORIAL SOLVER
// ============================================================================

class FactorialSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Factorial"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        
        // Check all values are positive integers
        for (double d : series) {
            if (!isInteger(d) || d < 1) return {0, "", ""};
        }
        
        // Check if series matches factorial pattern
        // Try different starting positions (0!, 1!, 2!, etc.)
        for (int startN = 0; startN <= 2; ++startN) {
            bool matches = true;
            for (size_t i = 0; i < series.size(); ++i) {
                int n = startN + i;
                double expectedFactorial = factorial(n);
                if (!areEqual(series[i], expectedFactorial)) {
                    matches = false;
                    break;
                }
            }
            
            if (matches) {
                int nextN = startN + series.size();
                double nextVal = factorial(nextN);
                std::string formula = "n! starting at " + std::to_string(startN) + "!";
                return {nextVal, formula, ""};
            }
        }
        
        return {0, "", ""};
    }

private:
    double factorial(int n) {
        if (n < 0) return 0;
        if (n == 0 || n == 1) return 1;
        double result = 1;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }
};

// ============================================================================
// POWER N SOLVER (n^n pattern)
// ============================================================================

class PowerNSolver : public SequenceSolver {
public:
    std::string getName() const override { return "Power n^n"; }

    SolverResult solve(const std::vector<double>& series) override {
        if (series.size() < 3) return {0, "", ""};
        
        // Try different starting positions (1^1, 2^2, etc.)
        for (int startN = 1; startN <= 3; ++startN) {
            bool matches = true;
            for (size_t i = 0; i < series.size(); ++i) {
                int n = startN + i;
                double expected = std::pow(n, n);
                if (!areEqual(series[i], expected)) {
                    matches = false;
                    break;
                }
            }
            
            if (matches) {
                int nextN = startN + series.size();
                double nextVal = std::pow(nextN, nextN);
                std::string formula = "n^n starting at " + std::to_string(startN) + "^" + std::to_string(startN);
                return {nextVal, formula, ""};
            }
        }
        
        return {0, "", ""};
    }
};


#endif // SEQUENCE_SOLVERS_H