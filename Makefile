.PHONY: help test fixture run clean

SERVICE_NAME ?= processor-cytometry-umap
PYTHON ?= python3

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test     - pytest: synthetic fixture through the processor, both input layouts, Lambda handler (no Docker)"
	@echo "make fixture  - write the synthetic fixture to data/fixture/{linked,flat}/"
	@echo "make run      - build the image and run it in Docker over data/fixture/linked/, writing to data/output/"
	@echo "make clean    - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	$(PYTHON) -m tests.make_fixture data/fixture

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

clean:
	rm -rf data/output
