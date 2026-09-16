#!/usr/bin/env python3
"""Build transfer inventory from clean annotated campaign worktrees."""
from pathlib import Path
import hashlib, subprocess, tomllib
out = Path(__file__).resolve().parent
prior = out.parent / "v16-deployment"
remote = "/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v17"
packages = [
    ("FLOWPanel", "flowpanel", "FLOWPanel.jl", Path("/private/tmp/flowpanel-p021-r4-counters-v17"), "campaign/p021-r4-counters-source-20260915-v17"),
    ("FastMultipole", "fastmultipole", "FastMultipole", Path("/private/tmp/fastmultipole-p021-r4-activity-v11"), "campaign/p021-r4-activity-source-20260915-v11"),
]
def git(w, *args):
    return subprocess.check_output(["git", "-C", str(w), *args]).decode().strip()
pins = []
for name, short, directory, worktree, tag in packages:
    assert not git(worktree, "status", "--porcelain"), f"Dirty {worktree}"
    assert git(worktree, "cat-file", "-t", tag) == "tag"
    sha = git(worktree, "rev-parse", "HEAD")
    assert sha == git(worktree, "rev-parse", tag + "^{commit}")
    files = subprocess.check_output(["git", "-C", str(worktree), "ls-files", "-z"]).decode().split("\0")
    files = sorted(f for f in files if f and not f.startswith("data/") and (worktree/f).is_file() and not (worktree/f).is_symlink())
    (out/f"{short}.files0").write_bytes("\0".join(files).encode()+b"\0")
    manifest = "".join(hashlib.sha256((worktree/f).read_bytes()).hexdigest()+"  "+f+"\n" for f in files)
    (out/f"{short}.sha256").write_text(manifest)
    digest = hashlib.sha256(manifest.encode()).hexdigest()
    pins.append(f'[packages.{name}]\npath = "{remote}/{directory}"\ntag = "{tag}"\nsha = "{sha}"\ndeployment = "rsync"\ncontent_manifest = "{remote}/{short}.sha256"\ncontent_manifest_sha256 = "{digest}"\n')
    print(name, sha, len(files), sum((worktree/f).stat().st_size for f in files), digest)
vpm = tomllib.loads((prior/"pins.toml").read_text())["packages"]["FLOWVPM"]
pins.append("[packages.FLOWVPM]\n"+"".join(f'{k} = "{v}"\n' for k,v in vpm.items()))
(out/"pins.toml").write_text("\n".join(pins))
for name in ("Project.toml", "Manifest.toml"):
    (out/name).write_text((prior/name).read_text().replace("/counters-v16/", "/counters-v17/"))
