# Tests Directory

Unit and integration tests for `src/` modules.

## Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_preprocessing.py

# Run with coverage
pytest --cov=src/

# Verbose output
pytest -v
```

## Structure

```
tests/
├── test_preprocessing.py   # Tests for src/preprocessing.py
├── test_features.py        # Tests for src/features.py
├── conftest.py             # Fixtures & setup
└── README.md
```

## Example Test

```python
import pytest
from src.preprocessing import clean_data

def test_clean_data_removes_nulls():
    """Test that clean_data removes null values."""
    input_data = [1, None, 3]
    result = clean_data(input_data)
    assert None not in result

def test_clean_data_empty_input():
    """Test error handling on empty input."""
    with pytest.raises(ValueError):
        clean_data([])
```

## Best Practices

✅ Test one thing per test  
✅ Use descriptive names  
✅ Use fixtures for setup  
✅ Test edge cases  

❌ Don't test external APIs  
❌ Don't test data files (use small fixtures)  

---
Before pushing: Run `pytest` to ensure quality.
