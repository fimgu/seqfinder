# Examples

This directory contains example sequences and test cases for the sequence solver.

## Files

### `seq.txt`
Collection of challenging sequence puzzles with solutions and explanations.

Contains sequences like:
- Position-dependent sequences
- Interleaved sequences
- Operation cycles
- Complex fraction patterns
- Power notation patterns
- And more!

Use these as test cases or inspiration for new solver implementations.

---

## Using These Examples

The sequences in `seq.txt` are useful for:

1. **Testing** - Verify solver behavior on specific patterns
2. **Benchmarking** - Evaluate solver capabilities
3. **Inspiration** - Find new sequence types to support
4. **Documentation** - Examples of what the solver can handle

## Creating Your Own Test Program

Create a new `.cpp` file to test specific sequences:

```cpp
#include <iostream>
#include <vector>
#include "../src/SequenceSolvers.h"

int main() {
    std::vector<double> series = {/* your sequence */};
    
    // Try different solvers
    YourSolver solver;
    auto result = solver.solve(series);
    
    if (!result.formula.empty()) {
        std::cout << "Pattern: " << result.formula << std::endl;
        std::cout << "Next: " << result.nextValue << std::endl;
    }
    
    return 0;
}
```

Compile with: `g++ -std=c++17 -I.. yourfile.cpp -o yourfile.exe`

