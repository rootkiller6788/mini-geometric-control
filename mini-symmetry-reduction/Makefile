CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libsymmetry_reduction.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_symmetry_reduction.c
TEST_BIN = build/test_symmetry_reduction
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
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lsymmetry_reduction $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lsymmetry_reduction $(LDFLAGS)

demo: examples
	@echo "=== Running example1_rigid_body ==="
	./build/example1_rigid_body
	@echo ""
	@echo "=== Running example2_heavy_top ==="
	./build/example2_heavy_top
	@echo ""
	@echo "=== Running example3_spacecraft_control ==="
	./build/example3_spacecraft_control

bench: all
	$(CC) $(CFLAGS) -o $(BENCH_BIN) $(BENCH_SRC) -Lbuild -lsymmetry_reduction $(LDFLAGS)
	./$(BENCH_BIN)

clean:
	rm -rf build
