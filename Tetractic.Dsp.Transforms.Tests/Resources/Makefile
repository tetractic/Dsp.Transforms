EXTRA_CFLAGS = $(shell pkg-config --cflags fftw3l)
EXTRA_LDLIBS = -lm $(shell pkg-config --libs fftw3l)

all: genref1d

genref1d.o: genref1d.c
	$(CC) -c -o $@ $^ $(EXTRA_CFLAGS)

genref1d: genref1d.o
	$(CC) -o $@ $^ $(EXTRA_LDLIBS)

clean:
	$(RM) genref1d genref1d.o
