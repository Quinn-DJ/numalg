CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
INCLUDES := -Iinclude

SRC := src/gauss.cpp src/solve.cpp
OBJ_DIR := build/obj
OBJ := $(SRC:src/%.cpp=$(OBJ_DIR)/%.o)

TEST_SRC := test/test.cpp
TIME_SRC := test/timeConsumingTest.cpp

TEST_BIN := build/test
TIME_BIN := build/timeConsumingTest

.PHONY: all test timeConsumingTest run-test run-timeConsumingTest clean

all: test timeConsumingTest

test: $(TEST_BIN)

timeConsumingTest: $(TIME_BIN)

run-test: $(TEST_BIN)
	./$(TEST_BIN)

run-timeConsumingTest: $(TIME_BIN)
	./$(TIME_BIN) $(START_N) $(END_N) $(STEP)

# Defaults for timeConsumingTest arguments (override via make START_N=... END_N=... STEP=...)
START_N ?= 10
END_N ?= 500
STEP ?= 10

$(TEST_BIN): $(OBJ) $(TEST_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(TEST_SRC) -o $@

$(TIME_BIN): $(OBJ) $(TIME_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(TIME_SRC) -o $@

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build
