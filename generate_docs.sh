#!/bin/bash
# Script to generate Doxygen documentation for Cancer Cell Agent Simulation

echo "Generating Doxygen documentation for Cancer Cell Agent Simulation..."

# Check if Doxygen is installed
if ! command -v doxygen &> /dev/null; then
    echo "Error: Doxygen is not installed."
    echo "Please install Doxygen to generate documentation."
    echo "On Ubuntu/Debian: sudo apt-get install doxygen"
    echo "On Fedora/RHEL: sudo dnf install doxygen"
    echo "On macOS: brew install doxygen"
    exit 1
fi

# Check if Graphviz is installed (for generating graphs)
if ! command -v dot &> /dev/null; then
    echo "Warning: Graphviz is not installed."
    echo "Some diagrams may not be generated correctly."
    echo "On Ubuntu/Debian: sudo apt-get install graphviz"
    echo "On Fedora/RHEL: sudo dnf install graphviz"
    echo "On macOS: brew install graphviz"
fi

# Create docs directory if it doesn't exist
mkdir -p docs/doxygen

# Clean the output directory to ensure old files are removed
echo "Cleaning previous documentation output..."
rm -rf docs/doxygen/html

# Run Doxygen
doxygen Doxyfile

if [ $? -eq 0 ]; then
    echo "Documentation successfully generated!"
    echo "HTML documentation is available at: docs/doxygen/html/index.html"
else
    echo "Error: Documentation generation failed."
    exit 1
fi

# Open documentation in browser if on a system with a GUI
if [ -n "$DISPLAY" ]; then
    if command -v xdg-open &> /dev/null; then
        echo "Opening documentation in browser..."
        xdg-open docs/doxygen/html/index.html
    elif command -v open &> /dev/null; then
        echo "Opening documentation in browser..."
        open docs/doxygen/html/index.html
    fi
fi

exit 0