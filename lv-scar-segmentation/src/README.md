# Source Code Directory

This directory contains production-ready Python modules and core logic.

## Usage

All pipeline logic goes here, not in notebooks.

```python
# Example structure (add as needed)
src/
├── __init__.py
├── preprocessing.py    # Data cleaning
├── features.py         # Feature extraction
├── models.py          # Model definitions
├── utils.py           # Helper functions
└── config.py          # Configuration loading
```

## Imports

To use your modules from notebooks:

```python
from src.preprocessing import clean_data
from src.features import extract_features
```

Or from tests:

```bash
python -m pytest tests/
```

## Best Practices

✅ Keep modular (one function = one job)  
✅ Add docstrings  
✅ Use type hints  
✅ Test with pytest  

---
Update as your project grows.
