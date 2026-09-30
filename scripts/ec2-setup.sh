#!/usr/bin/env bash
# Run on an Ubuntu 22.04/24.04 x86_64 EC2 instance as your SSH user.
set -euo pipefail

if [[ $(uname -s) != Linux || $(uname -m) != x86_64 ]]; then
    echo 'Use an x86_64 Linux EC2 instance for the DIA-NN containers.' >&2
    exit 1
fi
source /etc/os-release
if [[ $ID != ubuntu || ! $VERSION_ID =~ ^(22\.04|24\.04)$ ]]; then
    echo 'This installer supports Ubuntu 22.04 and 24.04.' >&2
    exit 1
fi
if (( EUID == 0 )); then
    echo 'Run this script as your SSH user, without sudo; it calls sudo where needed.' >&2
    exit 1
fi

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
sudo apt-get update
sudo apt-get install -y ca-certificates curl openjdk-17-jre-headless python3 tmux unzip

if ! command -v docker >/dev/null 2>&1; then
    sudo install -m 0755 -d /etc/apt/keyrings
    sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
    sudo chmod a+r /etc/apt/keyrings/docker.asc
    sudo tee /etc/apt/sources.list.d/docker.sources >/dev/null <<EOF
Types: deb
URIs: https://download.docker.com/linux/ubuntu
Suites: ${UBUNTU_CODENAME:-$VERSION_CODENAME}
Components: stable
Architectures: amd64
Signed-By: /etc/apt/keyrings/docker.asc
EOF
    sudo apt-get update
    sudo apt-get install -y docker-ce docker-ce-cli containerd.io
fi
sudo systemctl enable --now docker
# Docker group membership grants root-level access; appropriate for the trusted SSH user.
sudo usermod -aG docker "$(id -un)"

# Pin the version seen in this checkout's previous run; its minimum is 26.04.1.
export NXF_VER=26.04.6
installer=$(mktemp)
trap 'rm -f -- "$installer"' EXIT
curl -fsSL https://get.nextflow.io -o "$installer"
mkdir -p "$HOME/.local/bin"
(
    cd "$HOME/.local/bin"
    bash "$installer"
    chmod +x nextflow
)

# Windows archives may lose executable permissions or introduce CRLF endings.
python3 - "$repo_dir" <<'PY'
from pathlib import Path
import sys
for path in (Path(sys.argv[1]) / 'bin').iterdir():
    if path.is_file():
        content = path.read_bytes()
        if b'\r\n' in content:
            path.write_bytes(content.replace(b'\r\n', b'\n'))
        path.chmod(path.stat().st_mode | 0o111)
PY

"$HOME/.local/bin/nextflow" -version
echo 'Setup complete. Exit SSH and reconnect to activate Docker group membership.'
echo 'After reconnecting, verify: docker run --rm hello-world'
