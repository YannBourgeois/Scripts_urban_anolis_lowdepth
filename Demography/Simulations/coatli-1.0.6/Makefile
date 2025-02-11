CXXFLAGS+=-Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -Wno-unused-result -O3
PROG:=colmean matmean m2diff m3diff m2elsum m3elsum fold nt2mom nt3mom nt3mom_slow fs3mom mscond mscondfs ms2loci msfs msM2fs M2fs Smom testbias mscut mssub ms2fastphase ms2SF ms2vcf ms2tikzhaplo msTsub ms2tikztree
PARSER=msScan.cpp
#PARSER=libsequenceMsScan.cpp
#LIBS=-lz -lsequence
prefix:=/usr/local/bin

.PHONY: all install lib clean

all:	$(PROG)

colmean: colmean.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

m2diff: m2diff.cpp m2IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

m3diff: m3diff.cpp m3IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

m2elsum: m2elsum.cpp m2IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

m3elsum: m3elsum.cpp m3IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

fold: fold.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

matmean: matmean.cpp m2IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

fs3mom: fs3mom.cpp m3IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

nt2mom: nt2mom.cpp xit1.cpp xit2.cpp m2IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

nt3mom: nt3mom.cpp xit1.cpp xit2.cpp xit3.cpp xi3mom.cpp m3IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

nt3mom_slow: nt3mom.cpp xit1.cpp xit2.cpp xit3_slow.cpp xi3mom.cpp m3IO.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

mscondfs: mscondfs.cpp mscondfsCalc.cpp condfsOutput.cpp msScan.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

Smom: Smom.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

testbias: testbias.cpp xi3mom.cpp xit1.cpp xit2.cpp xit3.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

M2fs: M2fs.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

msfs: msfs.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

msM2fs: msM2fs.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

mscond: mscond.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

ms2loci: ms2loci.cpp mscondfsCalc.cpp msScan.cpp
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

mscut: mscut.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

mssub: mssub.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

ms2fastphase: ms2fastphase.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

ms2SF: ms2SF.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

ms2vcf: ms2vcf.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

ms2tikzhaplo: ms2tikzhaplo.cpp $(PARSER)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

msTsub: newicksubtree.cpp msTsub.cpp msScan.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

ms2tikztree: parseNewick.cpp ms2tikztree.cpp tree2tikz.cpp msScan.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

	
install: all
	for prog in $(PROG); do \
		install -m 0755 $$prog $(prefix); \
	done

lib: xi3mom.cpp xit1.cpp xit2.cpp xit3.cpp
	$(CXX) $(CXXFLAGS) -c xi3mom.cpp xit1.cpp xit2.cpp xit3.cpp
	ar rcs libxi3mom.a xi3mom.o xit1.o xit2.o xit3.o
	rm xi3mom.o xit1.o xit2.o xit3.o

clean:
	@- $(RM) $(PROG)
