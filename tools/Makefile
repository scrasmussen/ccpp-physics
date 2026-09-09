.DEFAULT_GOAL := help

ACT ?= act
PYTHON ?= python3
ACT_IMAGE ?= catthehacker/ubuntu:act-latest
ACT_ARCH ?= linux/amd64
ACT_ARGS ?=
PAGES_BIND ?= 127.0.0.1
PAGES_PORT ?= 8000

.PHONY: help workflow build pages

help:
	@printf '%s\n' \
	  'make build     Build the documentation locally with act.' \
	  'make pages     Serve the built site at http://$(PAGES_BIND):$(PAGES_PORT).' \
	  'make workflow  Generate the preview workflow without running act.' \
	  'See README.md for prerequisites and configuration.'

workflow:
	"$(PYTHON)" prepare_doxygen.py

build: workflow
	@command -v "$(ACT)" >/dev/null 2>&1 || { \
	  printf '%s\n' 'act was not found. Install act and start Docker; see tools/README.md.' >&2; \
	  exit 1; \
	}
	cd .. && "$(ACT)" workflow_dispatch \
	  --workflows tools/.build/doxygen/workflow.yml \
	  --job preview --bind \
	  --platform "ubuntu-latest=$(ACT_IMAGE)" \
	  --container-architecture "$(ACT_ARCH)" \
	  --secret GITHUB_TOKEN $(ACT_ARGS)
	@test -f .build/doxygen/site/index.html || { \
	  printf '%s\n' 'The build did not produce tools/.build/doxygen/site/index.html.' >&2; \
	  exit 1; \
	}
	@printf '%s\n' 'Site built in tools/.build/doxygen/site. Run make pages to preview it.'

pages:
	@test -f .build/doxygen/site/index.html || { \
	  printf '%s\n' 'No built site found. Run make build first.' >&2; \
	  exit 1; \
	}
	@printf '%s\n' 'Serving http://$(PAGES_BIND):$(PAGES_PORT) (Ctrl-C to stop).'
	"$(PYTHON)" -m http.server "$(PAGES_PORT)" --bind "$(PAGES_BIND)" --directory .build/doxygen/site
