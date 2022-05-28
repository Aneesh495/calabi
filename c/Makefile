CC ?= clang
CFLAGS ?= -O3 -Wall -Wextra -std=c11 -Iinclude -fPIC
LDFLAGS ?= 

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    SO_EXT := dylib
    SHARED_FLAG := -dynamiclib
else
    SO_EXT := so
    SHARED_FLAG := -shared
endif

all: bin/bench_c lib/libcalabi_simd.$(SO_EXT)

bin/bench_c: src/bench_c.c src/calabi_simd.c
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ src/bench_c.c src/calabi_simd.c -lm

lib/libcalabi_simd.$(SO_EXT): src/calabi_simd.c
	@mkdir -p lib
	$(CC) $(CFLAGS) $(SHARED_FLAG) -o $@ src/calabi_simd.c -lm

clean:
	rm -rf bin lib *.o

.PHONY: all clean
