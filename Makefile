test:
	python -m pytest -qx tests
debug:
	python -m pytest -qx tests --pdb
coverage:
	python -m pytest -v tests --cov zmat.py --cov molconvert.py --cov-report html --cov-report term

