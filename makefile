ga_epx_nkland : global.o file_man.o ga.o selection.o statistics.o transformation.o aux_functions.o
	g++ -Wall global.o file_man.o ga.o selection.o statistics.o transformation.o aux_functions.o -o ga_epx_nkland

global.o : global.cpp	
	g++ -Wall -o global.o -c global.cpp

file_man.o : file_man.cpp	
	g++ -Wall -o file_man.o -c file_man.cpp

ga.o : ga.cpp	
	g++ -Wall -o ga.o -c ga.cpp

selection.o : selection.cpp	
	g++ -Wall -o selection.o -c selection.cpp

statistics.o : statistics.cpp	
	g++ -Wall -o statistics.o -c statistics.cpp

transformation.o : transformation.cpp	
	g++ -Wall -o transformation.o -c transformation.cpp

aux_functions.o : aux_functions.cpp	
	g++ -Wall -o aux_functions.o -c aux_functions.cpp

