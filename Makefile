CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
INCLUDES := -Iinclude

SRC := src/gauss.cpp src/solve.cpp src/inverse_norm_estimate.cpp
OBJ_DIR := build/obj
OBJ := $(SRC:src/%.cpp=$(OBJ_DIR)/%.o)

TEST_SRC := test/test.cpp
TIME_SRC := test/timeConsumingTest.cpp
HOMEWORK_SRC := test/homework03.cpp
NORM_SRC := test/test_inverse_norm.cpp

TEST_BIN := build/test
TIME_BIN := build/timeConsumingTest
HOMEWORK_BIN := build/homework03
NORM_BIN := build/test_inverse_norm

.PHONY: all test timeConsumingTest homework inverseNorm run-test run-timeConsumingTest run-homework run-inverseNorm clean

all: test timeConsumingTest homework inverseNorm

test: $(TEST_BIN)

timeConsumingTest: $(TIME_BIN)

homework: $(HOMEWORK_BIN)

inverseNorm: $(NORM_BIN)

run-test: $(TEST_BIN)
	./$(TEST_BIN)

run-timeConsumingTest: $(TIME_BIN)
	./$(TIME_BIN) $(START_N) $(END_N) $(STEP)

run-homework: $(HOMEWORK_BIN)
	./$(HOMEWORK_BIN)

run-inverseNorm: $(NORM_BIN)
	./$(NORM_BIN)

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

$(HOMEWORK_BIN): $(OBJ) $(HOMEWORK_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HOMEWORK_SRC) -o $@

$(NORM_BIN): $(OBJ) $(NORM_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(NORM_SRC) -o $@

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build
