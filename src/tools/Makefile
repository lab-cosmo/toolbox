LDFLAGS=-L../libs -ltoolbox
-include ../../make.in

SOURCES=$(shell ls *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
DEPENDS=$(SOURCES:.cpp=.d)
TOOLS=autocorr crosscorr gaussmix histogram ndhistogram trajworks fourier
.PHONY: all clean $(TOOLS)

CXXFLAGS:=$(CXXFLAGS) -I../libs

all: $(TOOLS)

-include $(DEPENDS)
%.d: %.cpp
	$(MAKEDEPEND) $< > $@.$$$$;                      \
           sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
           rm -f $@.$$$$

trajworks.o: trajworks.cpp trajworks.d
	$(CXX) -c $(DREGEX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(OBJECTS) $(DEPENDS)

# generates -- and targets -- executables in another directory
$(TOOLS): % : ../../bin/%
$(TOOLS:%=../../bin/%): ../../bin/%: %.o 
	$(CXX) $(CXXFLAGS) -o $@ $< $(TBLDFLAGS)

../../bin/trajworks: trajworks.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(TBLDFLAGS) 
