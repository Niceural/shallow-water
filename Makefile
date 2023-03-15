CXX      := g++ # -c++
CXXFLAGS := -std=c++11 -pedantic-errors -Wall -Wextra -Werror -Wshadow
BUILD    := ./build
INCLUDE  := -Iinclude/

# src variables
SRC_TARGET   := shallowWater
SRC_LDFLAGS  := -L/usr/lib -lstdc++ -lm -lblas -lboost_program_options
SRC_OBJ_DIR  := $(BUILD)/srcObjects
SRC_BIN_DIR  := $(BUILD)/srcBin
SRCS      :=                      \
   $(wildcard src/matrices/*.cpp) \
   $(wildcard src/*.cpp)
SRC_OBJECTS  := $(SRCS:%.cpp=$(SRC_OBJ_DIR)/%.o)

# tests variables
TEST_LDFLAGS  := -L/usr/lib -lstdc++ -lm -lblas -lboost_unit_test_framework
TEST_OBJ_DIR  := $(BUILD)/testObjects
TEST_BIN_DIR  := $(BUILD)/testBin
TESTS := $(wildcard tests/*.cpp)
TEST_OBJECTS  := $(TESTS:%.cpp=$(TEST_OBJ_DIR)/%.o)
TEST_TARGETS :=  $(TESTS:%.cpp=$(TEST_BIN_DIR)/%)

# src targets

# 
build: init $(SRC_BIN_DIR)/$(SRC_TARGET)

$(SRC_OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(SRC_BIN_DIR)/$(SRC_TARGET): $(SRC_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(SRC_BIN_DIR)/$(TARGET) $^ $(SRC_LDFLAGS)


.PHONY: build clean debug release init test

debug: CXXFLAGS += -DDEBUG -g
debug: build

release: CXXFLAGS += -O2
release: build

# test targets
test: $(TEST_TARGETS)

$(TEST_OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(TEST_TARGETS): $(TEST_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGETS) $^ $(TEST_LDFLAGS)

init:
	@mkdir -p $(SRC_BIN_DIR)
	@mkdir -p $(SRC_OBJ_DIR)
	@mkdir -p $(TEST_BIN_DIR)
	@mkdir -p $(TEST_OBJ_DIR)

clean:
	-@rm -rvf $(SRC_BIN_DIR)/*
	-@rm -rvf $(SRC_OBJ_DIR)/*
	-@rm -rvf $(TEST_BIN_DIR)/*
	-@rm -rvf $(TEST_OBJ_DIR)/*

# CC = g++
# CFLAGS = -Wall -Wshadow
# LIBS = -lblas -lboost_unit_test_framework

# test1:
# 	${CC} ${CFLAGS} -o target/test1 src/main.cpp -lboost_program_options

# GeneralMatrix.o:
# 	${CC} ${CFLAGS} -o target/GeneralMatrix.o -c src/matrices/GeneralMatrix.cpp

# SquareBandedMatrix.o:
# 	${CC} ${CFLAGS} -o target/SquareBandedMatrix.o -c src/matrices/SquareBandedMatrix.cpp

# matrices: GeneralMatrix.o SquareBandedMatrix.o

# ShallowWater.o:
# 	${CC} ${CFLAGS} -o target/ShallowWater.o -c src/ShallowWater.cpp

# build: matrices ShallowWater.o

# initialCondition: matrices ShallowWater.o
# 	${CC} ${CFLAGS} -o target/initialCondition.o -c tests/initialCondition.cpp
# 	${CC} ${CFLAGS} -o target/initialCondition target/GeneralMatrix.o target/SquareBandedMatrix.o target/ShallowWater.o target/initialCondition.o ${LIBS}
# 	./target/initialCondition

# centralDifference: matrices ShallowWater.o
# 	${CC} ${CFLAGS} -o target/centralDifference.o -c tests/centralDifference.cpp
# 	${CC} ${CFLAGS} -o target/centralDifference target/GeneralMatrix.o target/SquareBandedMatrix.o target/ShallowWater.o target/centralDifference.o ${LIBS}
# 	./target/centralDifference

# test: initialCondition centralDifference
