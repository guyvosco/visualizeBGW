#!/usr/bin/env bash
set -e

# Install the package (editable) with dev dependencies if you defined them
if [ -f "pyproject.toml" ]; then
  pip install -e ".[dev]"
fi

# Run the test suite
pytest -q
