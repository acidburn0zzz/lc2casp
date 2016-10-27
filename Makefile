-include FLAGS

CLINGO_ROOT?=${HOME}/gringo-claspoutput
CLINGCON_ROOT?=${HOME}/clingcon-3.0
CLINGO_BUILD?=debug
CXXFLAGS?=-W -Wall

_CXXFLAGS=$(CXXFLAGS) -std=c++11 -I$(CLINGO_ROOT)/liblp -I$(CLINGO_ROOT)/libprogram_opts -I$(CLINGO_ROOT)/libgringo -I.
_LDFLAGS=-L$(CLINGO_ROOT)/build/$(CLINGO_BUILD) -llp -lprogram_opts -lgringo $(LDFLAGS)

TARGET=lc2casp
OBJECTS=main.o translator.o printer.o

all: $(TARGET)

test: $(TARGET)
	./test.sh $(CLINGO_ROOT)/build/$(CLINGO_BUILD)/gringo ./$(TARGET) $(CLINGCON_ROOT)/build/bin/clingcon

%.o: %.cc FLAGS
	$(CXX) $(_CXXFLAGS) -c -o $@ $<

$(TARGET): $(OBJECTS)
	$(CXX) $(_CXXFLAGS) -o $@ $^ $(_LDFLAGS)
ifdef STRIP
	$(STRIP) $@
endif

clean:
	rm -f $(OBJECTS) $(TARGET)

translator.o: translator.hh
printer.o: printer.hh
main.o: translator.hh printer.hh

FLAGS:
	echo 'CLINGO_ROOT=$(CLINGO_ROOT)' > FLAGS
	echo 'CLINGCON_ROOT=$(CLINGCON_ROOT)' >> FLAGS
	echo 'CLINGO_BUILD=$(CLINGO_BUILD)' >> FLAGS
	echo 'CXX=$(CXX)' >> FLAGS
	echo 'CXXFLAGS=$(CXXFLAGS)' >> FLAGS
	echo 'LDFLAGS=$(LDFLAGS)' >> FLAGS
	echo 'STRIP=' >> FLAGS

.PHONY: all test clean
