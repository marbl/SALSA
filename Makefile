COMPILER = g++
CFLAGS = -O3 -Wall -Wextra -std=c++11

ALL = break_contigs break_contigs_start correct_links

all: $(ALL)

break_contigs: break_contigs.cpp
	$(COMPILER) $(CFLAGS) -o break_contigs break_contigs.cpp

break_contigs_start: break_contigs_start.cpp
	$(COMPILER) $(CFLAGS) -o break_contigs_start break_contigs_start.cpp

correct_links: correct_links.cpp
	$(COMPILER) $(CFLAGS) -o correct_links correct_links.cpp


clean:
	rm -f $(ALL)
