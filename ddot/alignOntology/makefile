OBJS = alignOntology.o alignmentUtils.o
CC = g++
LIBS = -I ./boost_1_47_0
CFLAGS = -Wall -c -O4 -ffast-math -fomit-frame-pointer -fforce-addr -static
LFLAGS = -Wall -O4 -ffast-math -fomit-frame-pointer -fforce-addr -static

all: alignOntology calculateFDRsForAlignment collapseRedundantNodes ontologyTermStats

alignOntology: $(OBJS)
	$(CC) $(LFLAGS) $(LIBS) $(OBJS) -o alignOntology

alignOntology.o: alignOntology.cpp graph.h util.h alignmentUtils.h removalList.h
	$(CC) $(CFLAGS) $(LIBS) alignOntology.cpp

alignmentUtils.o: alignmentUtils.cpp alignmentUtils.h graph.h removalList.h util.h
	$(CC) $(CFLAGS) $(LIBS) alignmentUtils.cpp

calculateFDRsForAlignment: calculateFDRsForAlignment.o
	$(CC) $(LFLAGS) $(LIBS) calculateFDRsForAlignment.o -o calculateFDRsForAlignment

calculateFDRsForAlignment.o: calculateFDRsForAlignment.cpp util.h graph.h calculateFDRsForAlignment.h
	$(CC) $(CFLAGS) $(LIBS) calculateFDRsForAlignment.cpp

collapseRedundantNodes: collapseRedundantNodes.o collapseRedundantNodes.cpp collapseRedundantNodes.h
	$(CC) $(LFLAGS) $(LIBS) collapseRedundantNodes.o -o collapseRedundantNodes

collapseRedundantNodes.o: collapseRedundantNodes.cpp collapseRedundantNodes.h util.h graph.h
	$(CC) $(CFLAGS) $(LIBS) collapseRedundantNodes.cpp

ontologyTermStats: ontologyTermStats.o
	$(CC) $(LFLAGS) $(LIBS) ontologyTermStats.o -o ontologyTermStats

ontologyTermStats.o: ontologyTermStats.cpp util.h graph.h
	$(CC) $(CFLAGS) $(LIBS) ontologyTermStats.cpp

clean:
	rm *.o
