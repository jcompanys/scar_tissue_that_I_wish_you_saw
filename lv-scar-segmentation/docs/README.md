# Documentation Directory

Project documentation, methodology, and guides.

## Contents (as needed)

```
docs/
├── README.md                    # Overview
├── METHODOLOGY.md              # Methods & algorithms
├── DATA_DICTIONARY.md          # Variable definitions
├── SETUP.md                    # Installation & setup
└── API_REFERENCE.md            # Function documentation
```

## Building with Sphinx (Optional)

When ready to generate professional docs:

```bash
pip install sphinx sphinx-rtd-theme
sphinx-quickstart docs/
make -C docs/ html
open docs/_build/html/index.html
```

## Best Practices

✅ Keep README at project root  
✅ Link to specific modules/functions  
✅ Include examples  
✅ Update when code changes  

---
Documentation should be as valuable as code.
