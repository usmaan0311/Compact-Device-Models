CXX=g++ -Wall -O3 -o
TAR=mos.exe
SRC=mos.cc Poisson.cc Charge.cc


$(TAR): $(SRC)
	$(CXX) $(TAR) $(SRC)


plot:
	python PoissonPlot

cplot:
	python ChargePlot.py

clean:
	$(RM) *.exe *.dat ~* *~ *.swp
