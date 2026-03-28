# Contributing

## Before You Commit - CRITICAL DATA CHECK

Before pushing, run:

```bash
# Check for accidentally staged data files
git diff --cached --name-only | grep -E '\.(nii|nii\.gz|dcm|h5|db|npy|npz)$'

# If any results, STOP and remove them:
git reset HEAD data/raw/*.nii.gz
```

## Workflow

1. **Create a branch** for your work
   ```bash
   git checkout -b feature/my-feature
   ```

2. **Work locally**
   - Add code to `src/`
   - Explore in `notebooks/`
   - Test with `pytest tests/`

3. **Before pushing - verify NO DATA**
   ```bash
   git status | grep data/
   # Should be empty!
   ```

4. **Commit & push**
   ```bash
   git add src/ tests/ docs/
   git commit -m "Add feature X"
   git push origin feature/my-feature
   ```

5. **Open a pull request**

## Code Guidelines

- Keep `src/` modular
- Write docstrings
- Use type hints
- Test your code: `pytest tests/`
- Keep notebooks clean (no outputs on commit)

## Data Protection Rules

✅ **DO**
- Keep sensitive data local only
- Use `.gitignore` to protect data
- Store paths in config files
- Back up data separately

❌ **DON'T**
- Add data files to git
- Commit `.nii.gz`, `.dcm`, `.h5` files
- Hardcode personal identifiers
- Store passwords/API keys in code

## Need Help?

See [README.md](README.md) and [data/README.md](data/README.md)
