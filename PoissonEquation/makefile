CXX=g++ -Wall -O3 -o
TAR=mos.exe
SRC=mos.cc Poisson.cc


$(TAR): $(SRC)
	$(CXX) $(TAR) $(SRC)


plot:
	python PoissonPlot



clean:
	$(RM) *.exe *.dat ~* *~ *.swp
