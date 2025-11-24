# PowerShell script to set up Git pre-commit hooks for running tests

param(
    [switch]$Help = $false,
    [switch]$Remove = $false
)

# Show help
if ($Help) {
    Write-Host ""
    Write-Host "Git Pre-Commit Hook Setup" -ForegroundColor Cyan
    Write-Host "==========================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "This script installs a Git pre-commit hook that automatically runs"
    Write-Host "the test suite before each commit. If any tests fail, the commit"
    Write-Host "will be blocked."
    Write-Host ""
    Write-Host "Usage: .\setup_hooks.ps1 [options]"
    Write-Host ""
    Write-Host "Options:"
    Write-Host "  -Remove   Remove the pre-commit hook"
    Write-Host "  -Help     Show this help message"
    Write-Host ""
    Write-Host "Bypass hook: To skip the hook for a specific commit, use:"
    Write-Host "  git commit --no-verify"
    Write-Host ""
    exit 0
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$hookPath = Join-Path $scriptDir ".git\hooks\pre-commit"

# Remove hook if requested
if ($Remove) {
    if (Test-Path $hookPath) {
        Remove-Item $hookPath -Force
        Write-Host "Pre-commit hook removed successfully." -ForegroundColor Green
    } else {
        Write-Host "No pre-commit hook found." -ForegroundColor Yellow
    }
    exit 0
}

# Check if .git directory exists
if (-not (Test-Path (Join-Path $scriptDir ".git"))) {
    Write-Host "Error: Not a Git repository!" -ForegroundColor Red
    Write-Host "Please run this script from the root of your Git repository." -ForegroundColor Yellow
    exit 1
}

# Create hooks directory if it doesn't exist
$hooksDir = Join-Path $scriptDir ".git\hooks"
if (-not (Test-Path $hooksDir)) {
    New-Item -ItemType Directory -Path $hooksDir -Force | Out-Null
}

# Create pre-commit hook script
$hookContent = @'
#!/bin/sh
# Pre-commit hook to run test suite

echo ""
echo "Running pre-commit tests..."
echo "=============================="

# Run PowerShell test runner
powershell.exe -ExecutionPolicy Bypass -File "./run_tests.ps1"

if [ $? -ne 0 ]; then
    echo ""
    echo "Tests failed! Commit blocked."
    echo "To commit anyway, use: git commit --no-verify"
    exit 1
fi

echo ""
echo "All tests passed! Proceeding with commit..."
exit 0
'@

# Write the hook file
Set-Content -Path $hookPath -Value $hookContent -Encoding UTF8

Write-Host ""
Write-Host "Pre-commit hook installed successfully!" -ForegroundColor Green
Write-Host ""
Write-Host "The test suite will now run automatically before each commit." -ForegroundColor Cyan
Write-Host ""
Write-Host "To bypass the hook for a specific commit, use:" -ForegroundColor Yellow
Write-Host "  git commit --no-verify" -ForegroundColor Yellow
Write-Host ""
Write-Host "To remove the hook later, run:" -ForegroundColor Yellow
Write-Host "  .\setup_hooks.ps1 -Remove" -ForegroundColor Yellow
Write-Host ""
