CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libnonholonomic_planning.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_nonholonomic.c
TEST_BIN = build/test_nonholonomic
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
DEMO_SRC = demos/demo_overview.c
DEMO_BIN = build/demo_overview
BENCH_SRC = benches/bench_core.c
BENCH_BIN = build/bench_core

.PHONY: all builddir clean test examples demo bench

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lnonholonomic_planning $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/example_%: examples/example_%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lnonholonomic_planning $(LDFLAGS)

demo: examples
	@echo "=== Running Unicycle Demo ==="
	./build/example_unicycle
	@echo ""
	@echo "=== Running Car Parking Demo ==="
	./build/example_car_parking
	@echo ""
	@echo "=== Running Trailer Demo ==="
	./build/example_trailer

bench: all $(BENCH_BIN)
	./$(BENCH_BIN)

build/bench_%: benches/bench_%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lnonholonomic_planning $(LDFLAGS)

clean:
	rm -rf build
