#!/bin/bash
set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Set environment variables
export CARGO_TERM_COLOR=always
export SIDERUST_STUBS=1

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
