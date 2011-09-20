-include $(TOPDIR)/make.in
-include ../make.in
.PHONY: all clean 
all: libtoolbox.a autocorr crosscorr trajworks histogram ndhistogram dimred dimproj dimdist dimlandmark
libinterpol.o: libinterpol.cpp interpol.hpp tbdefs.hpp tensor-full.hpp
	$(CXX) -c -o libinterpol.o $(CCINC) libinterpol.cpp $(CCOPTS)
libdimred.o: libdimred.cpp dimreduce.hpp tbdefs.hpp
	$(CXX) -c -o libdimred.o $(CCINC) libdimred.cpp $(CCOPTS) -DNOARPACK
libioatoms.o: libioatoms.cpp ioatoms.hpp tbdefs.hpp
	$(CXX) -c -o libioatoms.o $(CCINC) libioatoms.cpp $(CCOPTS) 
libtb.o: libtb.cpp tbdefs.hpp
	$(CXX) -c -o libtb.o $(CCINC) libtb.cpp $(CCOPTS) 
libminsearch.o: libminsearch.cpp tbdefs.hpp minsearch.hpp
	$(CXX) -c -o libminsearch.o $(CCINC) libminsearch.cpp $(CCOPTS) 
libioparser.o: libioparser.cpp ioparser.hpp tbdefs.hpp
	$(CXX) -c -o libioparser.o $(CCINC) libioparser.cpp $(CCOPTS)
libclparser.o: libclparser.cpp clparser.hpp tbdefs.hpp
	$(CXX) -c -o libclparser.o $(CCINC) libclparser.cpp $(CCOPTS)
librndgen.o: librndgen.cpp rndgen.hpp
	$(CXX) -c -o librndgen.o $(CCINC) librndgen.cpp $(CCOPTS) 
libfmblas.o: matrix-full-blas.hpp libfmblas.cpp
	$(CXX) -c -o libfmblas.o $(CCINC) libfmblas.cpp $(CCOPTS) 
liblinalg.o: liblinalg.cpp linalg.hpp
	$(CXX) -c -o liblinalg.o $(CCINC) liblinalg.cpp $(CCOPTS) 
libunits.o: libunits.cpp conv-units.hpp
	$(CXX) -c -o libunits.o $(CCINC) libunits.cpp $(CCOPTS) 
libtoolbox.a: libtb.o libminsearch.o libioparser.o libclparser.o librndgen.o libfmblas.o liblinalg.o libioatoms.o libunits.o libdimred.o libinterpol.o
	ar rcs libtoolbox.a libtb.o libminsearch.o libioparser.o libclparser.o librndgen.o libfmblas.o liblinalg.o libioatoms.o libunits.o libdimred.o libinterpol.o
#gendist: tbdefs.hpp gendist.cpp librndgen.o
#	$(CXX) -o gendist $(CCINC) gendist.cpp $(CCOPTS) $(LDOPTS) librndgen.o
#matrix-test: tbdefs.hpp matrices.hpp matrix-full.hpp matrix-crs.hpp matrix-fun.hpp matrix-test.cpp matrix-coord.hpp matrix-conv.hpp ioparser.hpp libioparser.o libclparser.o librndgen.o libtb.o matrix-pcrs.hpp
#	$(CXX) -o matrix-test $(CCINC) matrix-test.cpp $(CCOPTS) $(LDOPTS) libioparser.o librndgen.o libtb.o
#matconv: tbdefs.hpp matrices.hpp matrix-full.hpp matrix-crs.hpp matrix-coord.hpp matrix-conv.hpp matconv.cpp libioparser.o libclparser.o
#	$(CXX) -o matconv $(CCINC) matconv.cpp $(CCOPTS) $(LDOPTS) libioparser.o libclparser.o
#io-test: libioparser.o io-test.cpp 
#	$(CXX) -o io-test $(CCINC) $(CCOPTS) $(LDOPTS) io-test.cpp  libioparser.o
autocorr: autocorr.cpp tools-autocorr.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o autocorr $(CCINC) $(CCOPTS) $(LDOPTS) autocorr.cpp libtoolbox.a -lfftw3 -lm 
	cp autocorr ../bin
crosscorr: crosscorr.cpp tools-autocorr.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o crosscorr $(CCINC) $(CCOPTS) $(LDOPTS) crosscorr.cpp libtoolbox.a -lfftw3 -lm 
	cp crosscorr ../bin
histogram: histogram.cpp tools-histogram.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o histogram $(CCINC) $(CCOPTS) $(LDOPTS) histogram.cpp libtoolbox.a -lm
ndhistogram: ndhistogram.cpp tools-histogram.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o ndhistogram $(CCINC) $(CCOPTS) $(LDOPTS) ndhistogram.cpp libtoolbox.a -lm
trajworks.o: trajworks.cpp ioatoms.hpp libtoolbox.a
	$(CXX) -c $(CCINC) $(CCOPTS) -I/usr/include trajworks.cpp -D__USEREGEX
trajworks: trajworks.o libtoolbox.a
	$(CXX) -o trajworks $(CCINC) $(CCOPTS) $(LDOPTS) -I/usr/include trajworks.o libtoolbox.a -lm -lfftw3 -llapack -lboost_regex
dimred: dimred.cpp libtoolbox.a
	$(CXX) -o dimred  $(CCINC) $(CCOPTS) $(LDOPTS) dimred.cpp libtoolbox.a -llapack 
dimproj: dimproj.cpp libtoolbox.a
	$(CXX) -o dimproj  $(CCINC) $(CCOPTS) $(LDOPTS) dimproj.cpp libtoolbox.a -llapack 
dimlandmark: dimlandmark.cpp libtoolbox.a
	$(CXX) -o dimlandmark $(CCINC) $(CCOPTS) $(LDOPTS) dimlandmark.cpp libtoolbox.a -llapack 
dimdist: dimdist.cpp libtoolbox.a
	$(CXX) -o dimdist $(CCINC) $(CCOPTS) $(LDOPTS) dimdist.cpp libtoolbox.a -llapack 
#pmult: libtoolbox.a pmult.cpp matrix-pcrs.hpp
#	$(CXX) -DTB_MPI -o pmult $(CCINC) $(CCOPTS) $(LDOPTS) pmult.cpp libtoolbox.a -lm 
#gaussmix: libtoolbox.a gaussmix.cpp
#	$(CXX) -o gaussmix $(CCINC) $(CCOPTS) $(LDOPTS) -I/usr/include gaussmix.cpp libtoolbox.a -lm -llapack
clean:
	rm -f *~ *.o
