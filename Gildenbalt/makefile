SRC = mos.cc GildenBalt.cc 
TAR = Gildenbalt.exe
CXX = g++ -Wall -O3 -o

$(TAR) : $(SRC)
	$(CXX) $(TAR) $(SRC)
mos:mos1.cc GildenBalt.cc
	$(CXX) mos1.exe mos1.cc GildenBalt.cc

plot: 
	python PlotPs

clean:
	$(RM) *.exe ~* *~
