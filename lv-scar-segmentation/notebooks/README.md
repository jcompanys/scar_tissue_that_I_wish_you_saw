# Notebooks Directory

Use Jupyter notebooks for **exploration and visualization only**.

⚠️ **Critical**: Notebooks are temporary. Core logic must be in `src/`.

## Naming Convention

Use numerical prefixes to show workflow order:

```
01_data_exploration.ipynb
02_preprocessing_validation.ipynb
03_analysis_results.ipynb
```

## Best Practices

✅ Import from `src/` modules  
✅ Use relative imports  
✅ Clear outputs before committing  
✅ Document findings  

❌ Don't hardcode data paths  
❌ Don't bury analysis logic in notebooks  
❌ Don't commit huge outputs  

## Example

```python
# Good: Use modular import
from src.preprocessing import clean_data

df = clean_data(input_file)
print(df.describe())
```

```python
# Bad: Copy-paste code
def clean_data(x):
    # 50 lines of complex logic
    pass
```

## Clean Before Committing

```bash
# Remove output cells
jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True *.ipynb

# Then commit
git add 01_*.ipynb
```

---
Notebooks are exploration. Code is in `src/`.
