CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libconnection_principal_bundle.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_bundle.c
TEST_BIN = build/test_bundle
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
DEMO_SRC = demos/demo_connection.c
DEMO_BIN = build/demo_connection

.PHONY: all builddir clean test examples demo

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lconnection_principal_bundle $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lconnection_principal_bundle $(LDFLAGS)

demo: examples
	@for exe in $(EXAMPLE_BINS); do \
		echo ""; \
		echo "=== Running $$(basename $$exe) ==="; \
		./$$exe; \
	done

clean:
	rm -rf build
