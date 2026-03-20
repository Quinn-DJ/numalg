CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
INCLUDES := -Iinclude

SRC := src/gauss.cpp src/solve.cpp
OBJ_DIR := build/obj
OBJ := $(SRC:src/%.cpp=$(OBJ_DIR)/%.o)

TEST_SRC := test/test.cpp
HW_SRC := test/homework02.cpp
TIME_SRC := test/timeConsumingTest.cpp

TEST_BIN := build/test
HW_BIN := build/homework02
TIME_BIN := build/timeConsumingTest

.PHONY: all test homework02 timeConsumingTest run-test run-homework02 run-timeConsumingTest clean

all: test homework02 timeConsumingTest

test: $(TEST_BIN)

homework02: $(HW_BIN)

timeConsumingTest: $(TIME_BIN)

run-test: $(TEST_BIN)
	./$(TEST_BIN)

run-homework02: $(HW_BIN)
	./$(HW_BIN)

run-timeConsumingTest: $(TIME_BIN)
	./$(TIME_BIN)

$(TEST_BIN): $(OBJ) $(TEST_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(TEST_SRC) -o $@

$(HW_BIN): $(OBJ) $(HW_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(HW_SRC) -o $@

$(TIME_BIN): $(OBJ) $(TIME_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJ) $(TIME_SRC) -o $@

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build
