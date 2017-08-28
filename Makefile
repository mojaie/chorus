
benchmark:
	@python3 ./scripts/benchmark.py

build:
	cd chorus/cython; python3 cython_setup.py build_ext --inplace

builddocs:
	cd docs; make html

test:
	@python3 -m unittest discover -s chorus/test
