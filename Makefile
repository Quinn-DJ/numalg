CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
INCLUDES := -Iinclude

SRC := src/gauss.cpp src/solve.cpp src/inverse_norm_estimate.cpp
OBJ_DIR := build/obj
OBJ := $(SRC:src/%.cpp=$(OBJ_DIR)/%.o)

TEST_SRC := test/test.cpp
TIME_SRC := test/timeConsumingTest.cpp
HW2_SRC := test/homework02.cpp
HW3_SRC := test/homework03.cpp
HW5_SRC := test/homework05.cpp
NORM_SRC := test/test_inverse_norm.cpp

TEST_BIN := build/test
TIME_BIN := build/timeConsumingTest
HW2_BIN := build/homework02
HW3_BIN := build/homework03
HW5_BIN := build/homework05
NORM_BIN := build/test_inverse_norm

.PHONY: all test timeConsumingTest hw2 hw3 hw5 inverseNorm run-test run-timeConsumingTest run-hw2 run-hw3 run-hw5 run-inverseNorm clean

all: test timeConsumingTest hw2 hw3 hw5 inverseNorm

test: $(TEST_BIN)

timeConsumingTest: $(TIME_BIN)

hw2: $(HW2_BIN)

hw3: $(HW3_BIN)

hw5: $(HW5_BIN)

inverseNorm: $(NORM_BIN)

run-test: $(TEST_BIN)
	./$(TEST_BIN)

run-timeConsumingTest: $(TIME_BIN)
	./$(TIME_BIN) $(START_N) $(END_N) $(STEP)

run-hw2: $(HW2_BIN)
	./$(HW2_BIN)

run-hw3: $(HW3_BIN)
	./$(HW3_BIN)

run-hw5: $(HW5_BIN)
	./$(HW5_BIN)

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

$(HW2_BIN): $(OBJ) $(HW2_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW2_SRC) -o $@

$(HW3_BIN): $(OBJ) $(HW3_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW3_SRC) -o $@

$(HW5_BIN): $(OBJ) $(HW5_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW5_SRC) -o $@

$(NORM_BIN): $(OBJ) $(NORM_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(NORM_SRC) -o $@

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build
