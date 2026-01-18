PROJECT_NAME := atomization_scorer
IMAGE_NAME := atomization_scorer

.PHONY: venv install test clean docker-build docker-run

# -------------------------
# Local (uv)
# -------------------------

# Create venv if missing
venv:
	uv venv

# Install deps + your package (includes dev extras)
install: venv
	uv sync --all-extras

# Run tests inside the uv environment
test:
	uv run pytest tests/

# Optional: format/lint/typecheck shortcuts if you use them
lint:
	uv run ruff check .

format:
	uv run ruff format .

typecheck:
	uv run mypy .

clean:
	rm -rf .pytest_cache
	rm -rf dist build *.egg-info
	rm -rf output

# -------------------------
# Docker
# -------------------------
docker-build:
	docker build -t $(IMAGE_NAME) .

# Runs the container image (assumes image ENTRYPOINT is your CLI)
docker-run:
	docker run --rm -it \
		-v $(PWD):/data \
		$(IMAGE_NAME) \
		/data/mini.fa /data/mini.geese /data/output
