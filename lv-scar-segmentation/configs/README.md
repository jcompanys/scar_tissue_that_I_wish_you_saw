# Configuration Directory

Store all configuration as YAML or JSON files.

## Best Practices

✅ **Externalize settings**: No hardcoded paths/parameters in code  
✅ **Version control**: Store configs in git  
✅ **Override easily**: Different configs per environment  

## Example Structure

```yaml
# configs/default.yaml
data:
  raw_path: "data/raw/"
  processed_path: "data/processed/"

preprocessing:
  normalize: true
  method: "minmax"

model:
  type: "random_forest"
  params:
    n_estimators: 100
    max_depth: 10
```

## Usage in Code

```python
import yaml
from pathlib import Path

def load_config(config_file="configs/default.yaml"):
    """Load config from YAML."""
    with open(config_file) as f:
        return yaml.safe_load(f)

config = load_config()
raw_data = config["data"]["raw_path"]
```

## Naming

```
configs/
├── default.yaml           # Default settings
├── experiment_01.yaml     # Experiment override
├── experiment_02.yaml     # Different experiment
└── local.yaml             # Local overrides (gitignored if sensitive)
```

⚠️ **Never put secrets here** (API keys, passwords)  
Use environment variables or `.env` instead.

---
All non-sensitive config should be versioned in git.
