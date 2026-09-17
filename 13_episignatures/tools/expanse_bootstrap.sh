#!/usr/bin/env bash
# Sync this repo to Expanse and (first time) build the conda envs there.
# Run from Git Bash on the Windows dev box, from the repo root. SSH goes through WSL, where the
# authorized Expanse key lives:   ssh from WSL  ->  login.expanse.sdsc.edu   (alias "expanse").
#
#   tools/expanse_bootstrap.sh            # push current main, pull on Expanse, ensure envs
#   tools/expanse_bootstrap.sh --no-envs  # skip the micromamba step
#
# Layout on Expanse:  $DEST.git (bare remote)   $DEST (working checkout, results live here)
set -euo pipefail
DEST="${EXPANSE_PROJECT_DIR:-/expanse/lustre/projects/ddp195/jsebat/dnmt3a-episignature}"
WSSH=(wsl.exe -e ssh -o BatchMode=yes expanse)
NOISE='CreateProcessCommon|Lmod has detected|cannot be loaded|module spider|^$'
# the remote command's exit status must survive the noise filter: with `|| true` every failed git clone / micromamba create
# was hidden and the script printed "done" under set -e (X25)
rssh() { local out rc; out=$("${WSSH[@]}" "$@" 2>&1); rc=$?; printf '%s\n' "$out" | tr -d '\0' | grep -vE "$NOISE" || true; return $rc; }

echo ">> remote: $DEST"
rssh "mkdir -p '$DEST' && { [ -d '$DEST.git' ] || git init -q --bare '$DEST.git'; }"
git remote get-url expanse >/dev/null 2>&1 || git remote add expanse "expanse:$DEST.git"
GIT_SSH_COMMAND="wsl.exe -e ssh" git push -q expanse main
rssh "cd '$DEST' && { [ -d .git ] || git clone -q '$DEST.git' . ; } && git pull -q && git log --oneline | head -1"

if [ "${1:-}" != "--no-envs" ]; then
  echo ">> conda envs via micromamba (skipped if present)"
  rssh "cd '$DEST' && export MAMBA_ROOT_PREFIX=\${MAMBA_ROOT_PREFIX:-\$HOME/micromamba} && \
    { micromamba env list | grep -q dnmt3a-py || micromamba create -y -q -f envs/py.yaml; } && \
    { micromamba env list | grep -q dnmt3a-r  || micromamba create -y -q -f envs/r.yaml; } && \
    micromamba run -n dnmt3a-py snakemake -n --snakefile workflow/Snakefile results/01_variant/annotation.md | tail -3"
fi
echo ">> done. On Expanse: cd $DEST && micromamba activate dnmt3a-py"
