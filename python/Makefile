../figures/%.pdf: %.py build_ext
	python $<

build_ext:
	python setup.py build_ext --inplace

clean:
	rm -rf build
	rm -f DSolveGC.cpp DSolveParallel.cpp findmax.c *.so *.pyc
