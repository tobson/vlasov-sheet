../figures/%.pdf: %.py build_ext
	python $<

build_ext:
	python setup.py build_ext --inplace

clean:
	rm -rf build *.so *.cpp *.c *.pyc
