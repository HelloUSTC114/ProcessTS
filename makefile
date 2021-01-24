CXXFLAGS = -g -Wall -I${DIR_INC} -I. `root-config --cflags`

Compare:  Compare.cpp 
	@echo Generating $@
	@`root-config --cxx `	-o	$@	$<	`root-config --libs` $(CXXFLAGS) $(CXXFLAGS) -lSpectrum

Match:  Match.cpp 
	@echo Generating $@
	@`root-config --cxx `	-o	$@	$<	`root-config --libs` $(CXXFLAGS) $(CXXFLAGS) -lSpectrum


