SOURCE=exp.cpp test_word_embedding.cpp 
MYPROGRAM=exp
MYPROGRAM2=test_word_embedding
CC=g++
CFLAGS = -lm -pthread -Wall -funroll-loops -Wno-unused-result -fopenmp -lpthread

all: $(MYPROGRAM) $(MYPROGRAM2)

$(MYPROGRAM): $(SOURCE)
	$(CC) ../src/$(MYPROGRAM).cpp -O3 -std=c++11 -o ../bin/$(MYPROGRAM) $(CFLAGS)


$(MYPROGRAM2): test_word_embedding.cpp word_embedding.h
	$(CC) ../src/$(MYPROGRAM2).cpp -O3 -std=c++11 -o ../bin/$(MYPROGRAM2) $(CFLAGS)

clean:

	rm -f ../bin/$(MYPROGRAM) ../bin/$(MYPROGRAM2)

test_word_embedding.o:word_embedding.h
