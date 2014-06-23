figures/%.pdf: python/*.py
	cd python; $(MAKE) ../$@

paper.pdf: paper.tex preamble.tex stix.tex references.bib figures
	latexmk paper

figures = drel guiding-center vf-mri vf-mri-2D
figures: $(addprefix figures/, $(addsuffix .pdf, $(figures)))

clean:
	latexmk -C paper
	cd python; make clean
