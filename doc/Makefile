FILE = Manual
all: pdf
dvi:
	latex $(FILE).tex
	latex $(FILE).tex
	latex $(FILE).tex
xdvi: dvi
	xdvi $(FILE).dvi
ps: dvi
	dvips -o $(FILE).ps $(FILE).dvi 
pdf:
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex
gv: ps
	gv $(FILE).ps&
clean:
	(rm -f *~ *.aux *.log *.dvi *.toc core)

