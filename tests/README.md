# Tests

This directory contains tests for the gnomad-toolbox package.

## Running Tests

To run all tests:

```bash
python -m pytest
```

To run specific test files:

```bash
python -m pytest tests/test_clinvar_loading.py
python -m pytest tests/test_clinvar_filtering.py
```

To run tests with verbose output:

```bash
python -m pytest -v
```

To run tests and show coverage:

```bash
python -m pytest --cov=gnomad_toolbox
```

## Test Structure

- `conftest.py`: Pytest configuration and fixtures
- `test_clinvar_loading.py`: Tests for ClinVar resource loading functionality
- `test_clinvar_filtering.py`: Tests for ClinVar filtering functionality
- `test_placeholder.py`: Placeholder test to ensure test suite works

## Test Conventions

- Tests are organized in classes that group related functionality
- Each test method has a descriptive name and docstring
- Tests use pytest fixtures and parametrization where appropriate
- Tests include both positive and negative test cases
- Tests check for expected exceptions and error conditions
