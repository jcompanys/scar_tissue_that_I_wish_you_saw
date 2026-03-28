# Research Project Template

**A minimal, scalable repository structure for research with strict data protection.**

> ⚠️ **CRITICAL**: This repo contains sensitive data protection measures. See [data/README.md](data/README.md) for details.

## Quick Start

```bash
# Clone and setup
git clone <your-repo-url>
cd research-template

# View structure
tree -L 2

# Check data protection
cat .gitignore
```

## Directory Structure

```
research-template/
├── data/                 # ⚠️  SENSITIVE - GITIGNORED
│   ├── raw/             # Original data (never committed)
│   ├── processed/       # Cleaned data (never committed)
│   └── interim/         # Working files (never committed)
│
├── src/                 # Source code modules
├── notebooks/           # Jupyter exploratory analysis
├── tests/               # Test files
├── configs/             # Configuration files
├── results/             # Output & results (gitignored)
├── docs/                # Documentation
│
├── .gitignore           # Data protection rules
├── README.md            # This file
├── LICENSE              # License
└── .gitkeep files       # Keep empty folders in git
```

### Folder Purposes

| Folder | Purpose | Committed? |
|--------|---------|-----------|
| `data/` | Raw & processed sensitive data | ❌ NO (gitignored) |
| `src/` | Production code & modules | ✅ YES |
| `notebooks/` | Jupyter exploration | ✅ YES (outputs cleared) |
| `tests/` | Unit & integration tests | ✅ YES |
| `configs/` | YAML/JSON configs, hyperparams | ✅ YES |
| `results/` | Pipeline outputs, metrics | ❌ NO (gitignored) |
| `docs/` | Project documentation | ✅ YES |

## Data Protection

### ⚠️ Critical Rules

1. **NEVER** commit to `data/` folders
2. **NEVER** push sensitive formats (`.nii.gz`, `.dcm`, `.h5`, etc.)
3. **ALWAYS** verify with `.gitignore` before `git add`
4. **ALWAYS** use `.gitkeep` placeholder files

### Check Before Committing

```bash
# See what's staged
git diff --cached

# See all untracked files in data/
git clean -n data/

# Remove accidentally added data file
git rm --cached data/raw/patient_01.nii.gz
echo "data/raw/*.nii.gz" >> .gitignore
git commit "Remove sensitive data file"
```

## Setup Instructions

### 1. Initialize Git

```bash
cd research-template
git init
git add .
git commit -m "initial: scalable research template"
git branch -M main
git remote add origin <github-url>
git push -u origin main
```

### 2. Add Your Data

```bash
# Copy raw data locally (not tracked)
cp /path/to/data/* data/raw/

# Verify it's ignored
git status  # Should NOT show data/ files
```

### 3. Create Project Config

```bash
# Create configs for your project
cp configs/README.md configs/default.yaml
vim configs/default.yaml
```

### 4. Start Analysis

```bash
# Add your code to src/
# Use notebooks/ for exploration
# Run tests with: pytest tests/
```

## Workflow

```
raw data (local)
    ↓
notebooks/ (exploratory)
    ↓
src/ (modularize good code)
    ↓
tests/ (validate)
    ↓
results/ (output, gitignored)
```

## Key Files

- **[.gitignore](.gitignore)** - Strict data protection rules
- **[data/README.md](data/README.md)** - Data governance guidelines
- **[CONTRIBUTING.md](CONTRIBUTING.md)** - How to contribute

## Best Practices

✅ **DO**
- Keep data local: `data/` folders are gitignored
- Version code: All `src/` changes go to git
- Document: Use `docs/` and comments
- Test: Run `pytest` before pushing
- Config externalize: Store settings in `configs/`

❌ **DON'T**
- Commit data files (`.nii.gz`, `.csv`, `.h5`, `.db`)
- Hard-code paths (use config files)
- Skip tests before push
- Commit API keys / passwords (.env in `.gitignore`)
- Push large outputs to git (use results/ which is gitignored)

## Scaling Up

As your project grows, add:

- **Continuous Integration**: `.github/workflows/` for testing
- **Pre-commit hooks**: `.pre-commit-config.yaml` to prevent data commits
- **Docker**: `Dockerfile` for reproducible environments
- **Package management**: `requirements.txt` or `pyproject.toml`
- **Docs build**: `docs/` → Sphinx documentation

## Support Files

- Read full data policy: [data/README.md](data/README.md)
- See contributing guide: [CONTRIBUTING.md](CONTRIBUTING.md)
- Check license: [LICENSE](LICENSE)

---

**Template Version**: 1.0  
**Last Updated**: March 2025
