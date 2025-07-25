.PHONY: clean clean-build clean-pyc clean-test coverage dist docs help install lint lint/flake8 cmodules

.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts, rm -f {TARGET}? but removing cutils.c
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +
	rm -f csolver.so

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

lint/flake8: ## check style with flake8
	python -m flake8 src/cube_solver --count --select=E9,F63,F7,F82 --show-source --statistics
	python -m flake8 src/cube_solver --count --exit-zero --max-complexity=23 --max-line-length=127 --statistics

lint: lint/flake8 ## check style

test: ## run tests quickly with the default Python
	python -m pytest
	python -m pytest --doctest-modules src/cube_solver/

coverage: ## check code coverage quickly with the default Python
	python -m coverage run -m pytest
	python -m coverage report -m
	python -m coverage html
	$(BROWSER) htmlcov/index.html

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ src/cube_solver -e -M -T
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

dist: clean ## builds source and wheel package
	python -m build
	ls -l dist

release: dist ## package and upload a release
	python -m twine upload dist/*

install: clean ## install the package to the active Python's site-packages
	python -m pip install --upgrade pip
	python -m pip install -r requirements_dev.txt
