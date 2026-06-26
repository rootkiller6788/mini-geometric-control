CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libfeedback_linearization.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_feedback_linearization.c
TEST_BIN = build/test_feedback_linearization
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir clean test examples demo

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lfeedback_linearization $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lfeedback_linearization $(LDFLAGS)

demo: examples
	@echo "=== Running example_pendulum ==="
	./build/example_pendulum
	@echo ""
	@echo "=== Running example_ball_beam ==="
	./build/example_ball_beam
	@echo ""
	@echo "=== Running example_robot_arm ==="
	./build/example_robot_arm

clean:
	rm -rf build
