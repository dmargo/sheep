include ../Makefile.config

BIN = efennel graph2adj read_partition tree2adj tree2dot vfennel

all: $(BIN)

clean: force
	rm -f $(BIN)
force:



DEPPATH = -I../lib -I../../llama/llama/include
DEPCPP = ../lib/jnode.cpp ../lib/jtree.cpp ../lib/partition.cpp
DEPH 	 = ../lib/defs.h ../lib/graph_wrapper.h ../lib/jdata.h ../lib/jnode.h ../lib/jtree.h \
			   ../lib/merge.h ../lib/partition.h ../lib/readerwriter.h ../lib/sequence.h ../lib/unionfind.h

efennel: efennel.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o efennel efennel.cpp $(LDFLAGS) $(LIBS)

graph2adj: graph2adj.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o graph2adj graph2adj.cpp $(LDFLAGS) $(LIBS)
	
read_partition: read_partition.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o read_partition read_partition.cpp $(LDFLAGS) $(LIBS)
	
tree2adj: tree2adj.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o tree2adj tree2adj.cpp $(LDFLAGS) $(LIBS)

tree2dot: tree2dot.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o tree2dot tree2dot.cpp $(LDFLAGS) $(LIBS)
	
vfennel: vfennel.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o vfennel vfennel.cpp $(LDFLAGS) $(LIBS)

