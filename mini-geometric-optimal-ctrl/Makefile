CC = gcc
CFLAGS = -std=c11 -Iinclude -Wall -Wextra -O2
LDFLAGS = -lm
LIB = build/libgeo_optctrl.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_geo_optctrl.c
TEST_BIN = build/test_geo_optctrl
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir clean test examples demo bench

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lgeo_optctrl $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo All examples built.

build/%: examples/%.c
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lgeo_optctrl $(LDFLAGS)

demo: examples
	@echo =========================================
	@echo   Geometric Optimal Control Demo
	@echo =========================================
	./build/example1_so3_attitude
	@echo ---
	./build/example2_subriemannian_car
	@echo ---
	./build/example3_brockett_integrator

bench: all
	$(CC) $(CFLAGS) -o build/bench_geo benches/bench_geo.c -Lbuild -lgeo_optctrl $(LDFLAGS)
	./build/bench_geo

clean:
	rm -rf build
