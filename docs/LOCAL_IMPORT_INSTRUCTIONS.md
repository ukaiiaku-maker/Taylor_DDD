# Local full-code import instructions

The GitHub repository has been initialized with documentation and project metadata. The full code bundle should be imported from the local export directory or ZIP created during the ChatGPT session.

Recommended local workflow:

```bash
cd /Users/sdillon

git clone https://github.com/ukaiiaku-maker/Taylor_DDD.git
cd Taylor_DDD

# Copy the exported files into this clone. For example, if the local export folder is available:
rsync -avh --exclude '.git/' /path/to/taylor_ddd_repo/ ./

# Or, if using the ZIP bundle:
# unzip taylor_ddd_repo_initial_import.zip
# rsync -avh taylor_ddd_repo/ ./

python3 -m venv .venv-opendis
source .venv-opendis/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 -m py_compile clean_arrhenius_taylor_ddd_v17.py analyze_v6_results.py analyze_depin_burst_statistics.py plot_depin_burst_ccdfs.py

git status --short
git add .
git commit -m 'import Arrhenius Taylor DDD drivers and analysis scripts'
git push origin main
```

After this, Codex and VS Code should be able to work directly against the browseable Python and shell files in the repository.
