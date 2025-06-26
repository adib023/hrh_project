
CC = c++
CC_INCLUDE = -I/home/adib023/eigen_3_4_0/


TARGETS = main

all: $(TARGETS)

main: main.cpp $(wildcard *.h)
	$(CC) -g $< -O3 -o $@ -O3 $(CC_INCLUDE) -fopenmp -Wall -Wextra -pedantic


clean:
	rm -f $(TARGETS)
  
again: clean $(TARGETS)main
