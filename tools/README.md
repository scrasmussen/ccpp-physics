# Preview the documentation locally

Run the documentation builds from `.github/workflows/doxygen.yml` using
[`act`](https://nektosact.com/usage/index.html), then view the generated site
in your browser. All workflow checkouts and generated pages stay inside this
repository under `tools/.build/doxygen/`, which Git ignores.

## Requirements

- Git, GNU Make, and Python 3.8 or newer.
- A recent installation of [`act`](https://nektosact.com/installation/index.html).
- A running Docker engine with Linux container support. Docker Desktop works
  on macOS and Windows; on Windows, run these commands from WSL with Docker
  Desktop's WSL integration enabled.
- Internet access to download the runner image, GitHub actions, repositories,
  and documentation packages.
- A GitHub token that can read the repositories checked out by the workflow.
  `act` reads `GITHUB_TOKEN` from your environment or prompts for it. Publishing
  permissions are not needed for this preview.

Doxygen, Graphviz, and LaTeX are installed inside the runner container by the
workflow. They do not need to be installed on your computer. The first build
downloads the container image and documentation dependencies and can take a while.

## Build and view the pages

From the repository root:

```sh
make -C tools build
make -C tools pages
```

Or run `make build` and `make pages` after changing into `tools/`.

Open <http://127.0.0.1:8000> to view the landing page and its HSD, GWv1, LandDA,
V7, and SRWv3 documentation. Press Ctrl-C to stop the server.

`make pages` serves the existing build locally; it does not publish to GitHub
Pages or rebuild the site. Run `make build` again when you want a fresh build.
The entry point is `tools/.build/doxygen/site/index.html`.
Each build clears the previous generated site before rebuilding it.

## How the preview works

`make build` first runs `prepare_doxygen.py`, a Python standard-library helper
that derives `tools/.build/doxygen/workflow.yml` from the deployment workflow.
The preview retains the build steps through **Build landing page**, removes
the upstream-repository restriction and deployment environment, and omits
the Pages setup, artifact upload, and deployment steps. It then runs `act`
with `--bind` so the generated files remain available locally.
The current repository is already mounted in the container, so the preview
reads its stored documentation configuration directly from that working tree.

The workflow file in your working tree is used, including uncommitted workflow
changes. Physics sources still follow the workflow's configured branches and
tags: `main`, v7, and SRWv3. Uncommitted changes to local physics source files or
documentation are not overlaid onto those checkouts.

Build directories:

```text
tools/.build/doxygen/
  workflow.yml   Generated preview workflow
  scm_repo/      SCM, framework, and physics checkouts used for the builds
  site/         Complete generated website
```

`act` and Docker also maintain their own action and image caches outside this
directory. The source deployment workflow is not changed by the helper.

## Configuration

Override Make variables on the command line:

```sh
make -C tools pages PAGES_PORT=8080
make -C tools build ACT_ARGS='--verbose'
make -C tools build ACT_ARGS='--pull=false'
make -C tools workflow
```

`make workflow` generates the preview YAML without starting a container.

| Variable | Default | Purpose |
| --- | --- | --- |
| `ACT` | `act` | Path to the `act` executable |
| `PYTHON` | `python3` | Path to Python |
| `ACT_IMAGE` | `catthehacker/ubuntu:act-latest` | Ubuntu runner image |
| `ACT_ARCH` | `linux/amd64` | Container architecture |
| `ACT_ARGS` | empty | Additional arguments passed to `act` |
| `PAGES_BIND` | `127.0.0.1` | Local server bind address |
| `PAGES_PORT` | `8000` | Local server port |

The default architecture also works on Apple Silicon using Docker's amd64
emulation. You can override `ACT_ARCH` and `ACT_IMAGE` together when using a
runner image with native support for another architecture. See the
[`act` runner documentation](https://nektosact.com/usage/runners.html).

If you build on a remote computer, forward the chosen server port over SSH
to view it in your local browser.
