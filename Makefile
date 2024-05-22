# Define variables
PYTHON = python3
VENV_DIR = python/bindings/venv
ACTIVATE = . $(VENV_DIR)/bin/activate
REQUIREMENTS = python/bindings/requirements.txt
TEST_DIR = python/tests
INPUT_DIR = python/bindings/build/inputs
OUTPUT_DIR = python/bindings/build/outputs
DATASET ?= $(INPUT_DIR)

# Default target
all: clean setup build install test

# Setup target: Create virtual environment and install dependencies
setup: $(VENV_DIR)/bin/activate

$(VENV_DIR)/bin/activate: $(REQUIREMENTS)
	if [ ! -d "$(VENV_DIR)" ]; then \
		$(PYTHON) -m venv $(VENV_DIR); \
	fi
	$(ACTIVATE) && pip install -U pip && pip install -r $(REQUIREMENTS)
	#touch $(VENV_DIR)/bin/activate

# Build target
build: setup
	$(ACTIVATE) && cd python/bindings && $(PYTHON) setup.py build_ext --inplace

# Install target
install: build
	$(ACTIVATE) && cd python/bindings && pip install .

# Generate dataset if necessary
dataset:
	@if [ ! -d "$(DATASET)" ]; then \
		echo "Creating a random dataset in $(DATASET)"; \
		mkdir -p $(DATASET); \
		$(ACTIVATE) && $(PYTHON) python/create_test_dataset.py $(DATASET); \
	fi

# Test target
test: setup dataset
	$(ACTIVATE) && for test_file in $(TEST_DIR)/*.py; do \
		echo "Running $$test_file with dataset $(DATASET)"; \
		mkdir -p $(OUTPUT_DIR); \
		$(PYTHON) $$test_file $(DATASET) $(OUTPUT_DIR) || exit 1; \
	done

# Clean target
clean:
	rm -rf $(VENV_DIR) python/bindings/build python/bindings/*.so python/bindings/*.egg-info

# Phony targets
.PHONY: all setup build install dataset test clean
