-include ../make.in
.PHONY: all clean

CXXINC:=$(CXXINC) -I../libs
TOOLS=autocorr crosscorr gaussmix histogram ndhistogram trajworks

all: $(TOOLS)

$(TOOLS): % : %.cpp
	$(CXX) $(CXXFLAGS) $(CXXINC) -o $@ $< $(LDFLAGS)

trajworks: trajworks.cpp
	$(CXX) -D__USEREGEX $(CXXFLAGS) $(CXXINC) -o $@ $< $(LDFLAGS) -lboost_regex
 
clean:
	rm -f $(TOOLS)
