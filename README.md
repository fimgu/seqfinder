# Sequence Solver

A comprehensive sequence pattern recognition and solving tool with support for multiple pattern types.

## Quick Start

### Running the CLI
```powershell
# Compile main application
g++ -std=c++17 src\main.cpp -o seqfinder_cli.exe

# Test a sequence
.\seqfinder_cli.exe 2 3 5 7 11
```

### Running Tests
```powershell
# Quick test
.\run_tests.ps1

# Or using Make
cd tests
make test
```

## Project Structure

```
seqfinder/
├── src/
│   ├── main.cpp              # CLI application
│   └── SequenceSolvers.h     # All solver implementations
├── tests/
│   ├── test_framework.h      # Test framework
│   ├── test_suite.cpp        # 31 comprehensive tests
│   └── Makefile              # Test build configuration
├── gui/
│   └── app.py                # GUI application
├── examples/
│   ├── test_prime_diffs.cpp  # Example: testing prime differences
│   ├── reproduce_issue.cpp   # Example: debugging specific sequences
│   └── seq.txt               # Sample sequences and explanations
├── run_tests.ps1             # Test runner script
├── setup_hooks.ps1           # Pre-commit hook installer
├── TESTING.md                # Complete testing guide
└── README.md                 # This file
```

## Supported Sequence Types

The solver can recognize and solve:

- **Prime Sequences** - Prime numbers with patterns
- **Polynomial** - Linear, quadratic, cubic sequences
- **Fractions** - Simple and complex fraction patterns
- **Powers** - Exponential growth patterns
- **Multipliers** - Constant ratio sequences
- **Interleaved** - Multiple sequences interleaved
- **Operation Cycles** - Repeating operations (+, -, *, /)
- **Fibonacci/Recurrence** - Linear recurrence relations
- **Difference Patterns** - Arithmetic and changing differences
- **Digit Patterns** - Palindromes, repeating digits
- **And more!** - See [src/SequenceSolvers.h](src/SequenceSolvers.h) for all solvers

## Testing

### Running Tests
```powershell
.\run_tests.ps1              # Recommended
cd tests && make test        # Using Make
```

### Adding Tests
Open `tests/test_suite.cpp` and add:

```cpp
TEST_CASE(MyTest, "Description") {
    std::vector<double> series = {1, 2, 3, 4};
    PolynomialSolver solver;
    auto result = solver.solve(series);
    
    EXPECT_FALSE(result.formula.empty());
    EXPECT_APPROX(result.nextValue, 5, 0.01);
    PASS();
}
```

See [TESTING.md](TESTING.md) for complete guide.

### Pre-Commit Testing
Automatically run tests before commits:
```powershell
.\setup_hooks.ps1            # Install
.\setup_hooks.ps1 -Remove    # Uninstall
```

## GUI Application

Start the web interface:
```powershell
python gui\app.py
```

Then open http://localhost:5000 in your browser.

## Examples

See the [examples/](examples/) directory for:
- `test_prime_diffs.cpp` - Testing sequences with prime differences
- `reproduce_issue.cpp` - Debugging specific issues
- `seq.txt` - Collection of interesting sequences with explanations

## Development

### Building
```powershell
# CLI application
g++ -std=c++17 src\main.cpp -o seqfinder_cli.exe

# Test suite
cd tests
g++ -std=c++17 -Wall -Wextra test_suite.cpp -o test_suite.exe
```

### Before Committing
1. Run tests: `.\run_tests.ps1`
2. Verify all pass before pushing

## Documentation

- **[TESTING.md](TESTING.md)** - Complete testing guide
- **[examples/seq.txt](examples/seq.txt)** - Sample sequences
- **[src/SequenceSolvers.h](src/SequenceSolvers.h)** - Solver implementations

## License

[Add your license here]
