#!/bin/bash
set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Set environment variables
export CARGO_TERM_COLOR=always

# Function to run standard CI checks
run_standard_checks() {
    echo -e "${YELLOW}Starting CI checks locally...${NC}\n"

    # Check
    echo -e "${YELLOW}==> Running cargo check${NC}"
    cargo check --all-targets
    echo -e "${GREEN}✓ Check passed${NC}\n"

    # Format
    echo -e "${YELLOW}==> Running cargo fmt${NC}"
    cargo fmt --check
    echo -e "${GREEN}✓ Format check passed${NC}\n"

    # Clippy
    echo -e "${YELLOW}==> Running cargo clippy${NC}"
    cargo clippy --all-targets -- -D warnings
    echo -e "${GREEN}✓ Clippy passed${NC}\n"

    # Tests
    echo -e "${YELLOW}==> Running tests${NC}"
    cargo test --all-targets
    echo -e "${GREEN}✓ Tests passed${NC}\n"

    echo -e "${YELLOW}==> Running doc tests${NC}"
    cargo test --doc
    echo -e "${GREEN}✓ Doc tests passed${NC}\n"

    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}All CI checks passed! ✓${NC}"
    echo -e "${GREEN}========================================${NC}"
}

# Function to run coverage checks
run_coverage() {
    echo -e "${YELLOW}Running coverage checks (matching GitHub CI)...${NC}\n"

    # Clean previous coverage data
    echo -e "${YELLOW}==> Cleaning previous coverage data${NC}"
    cargo +nightly llvm-cov clean --workspace
    rm -rf target/llvm-cov-target/*.profraw
    echo -e "${GREEN}✓ Coverage data cleaned${NC}\n"

    # Run tests and collect coverage (no report yet)
    echo -e "${YELLOW}==> Running tests with coverage instrumentation${NC}"
    cargo +nightly llvm-cov --workspace --all-features --doctests --no-report
    echo -e "${GREEN}✓ Coverage data collected${NC}\n"

    # Generate Cobertura XML report
    echo -e "${YELLOW}==> Generating Cobertura XML report${NC}"
    cargo +nightly llvm-cov report --cobertura --output-path coverage.xml
    echo -e "${GREEN}✓ coverage.xml generated${NC}\n"

    # Generate HTML report
    echo -e "${YELLOW}==> Generating HTML report${NC}"
    cargo +nightly llvm-cov report --html --output-dir coverage_html
    echo -e "${GREEN}✓ HTML report generated in coverage_html/${NC}\n"

    # Check coverage threshold
    echo -e "${YELLOW}==> Checking coverage threshold (≥90% lines)${NC}"
    cargo +nightly llvm-cov --workspace --all-features --doctests --no-run --fail-under-lines 90
    echo -e "${GREEN}✓ Coverage threshold met${NC}\n"

    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Coverage checks passed! ✓${NC}"
    echo -e "${GREEN}View HTML report: coverage_html/index.html${NC}"
    echo -e "${GREEN}========================================${NC}"
}

# Parse command line arguments
case "${1:-}" in
    coverage)
        run_coverage
        ;;
    all)
        run_standard_checks
        run_coverage
        ;;
    "")
        run_standard_checks
        ;;
    *)
        echo -e "${RED}Usage: $0 [coverage|all]${NC}"
        echo -e "  (no args)  - Run standard CI checks (check, fmt, clippy, test)"
        echo -e "  coverage   - Run coverage checks only"
        echo -e "  all        - Run both standard checks and coverage"
        exit 1
        ;;
esac
