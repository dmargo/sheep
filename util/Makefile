include ../Makefile.config

BIN = efennel graph2adj read_partition tree2adj tree2dot vfennel

all: $(BIN)

clean: force
	rm -f $(BIN)
force:



DEPPATH = -I../lib -I../../llama/llama/include
DEPCPP =
DEPH 	 = ../lib/defs.h ../lib/graph_wrapper.h ../lib/jdata.h ../lib/jnode.h ../lib/jtree.h \
			   ../lib/merge.h ../lib/partition.h ../lib/sequence.h ../lib/unionfind.h

efennel: efennel.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o efennel efennel.cpp $(DEPCPP) $(LDFLAGS) $(LIBS)

graph2adj: graph2adj.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o graph2adj graph2adj.cpp $(DEPCPP) $(LDFLAGS) $(LIBS)
	
read_partition: read_partition.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o read_partition read_partition.cpp $(DEPCPP) $(LDFLAGS) $(LIBS)
	
tree2adj: tree2adj.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o tree2adj tree2adj.cpp $(DEPCPP) $(LDFLAGS) $(LIBS)

tree2dot: tree2dot.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o tree2dot tree2dot.cpp $(DEPCPP) $(LDFLAGS) $(LIBS)
	
vfennel: vfennel.cpp $(DEPCPP) $(DEPH) 
	$(CC) $(CXXFLAGS) $(DEPPATH) -o vfennel vfennel.cpp $(DEPCPP) $(LDFLAGS) $(LIBS)

