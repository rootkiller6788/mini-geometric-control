CC = gcc
CFLAGS = -std=c11 -O2 -Iinclude -w
LDFLAGS = -lm

SRCS = src/polynomial.c src/flatness_core.c src/lie_algebra.c src/trajectory_generation.c src/dynamic_feedback.c src/flatness_applications.c src/flatness_algorithms.c src/flatness_analysis.c src/flatness_advanced.c src/flatness_models.c src/flatness_numerical.c src/flatness_observer.c src/flatness_control.c src/flatness_advanced_control.c
OBJS = $(SRCS:.c=.o)

.PHONY: all clean test examples

all: libflatness.a

libflatness.a: $(OBJS)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: libflatness.a tests/test_runner.c
	$(CC) $(CFLAGS) tests/test_runner.c -L. -lflatness $(LDFLAGS) -o tests/test_runner
	./tests/test_runner

examples: libflatness.a examples/example_main.c
	$(CC) $(CFLAGS) examples/example_main.c -L. -lflatness $(LDFLAGS) -o examples/example_main
	./examples/example_main

clean:
	rm -f src/*.o
	rm -f libflatness.a
	rm -f tests/test_runner examples/example_main
