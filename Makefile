# inspired from https://www.partow.net/programming/makefile/index.html
CXX      := g++ # -c++
CXXFLAGS := -std=c++11 -pedantic-errors -Wall -Wextra -Werror -Wshadow
LDFLAGS  := -L/usr/lib -lstdc++ -lm -lblas -lboost_program_options
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
TARGET   := shallowWater
INCLUDE  := -Iinclude/
SRC      :=                      \
   $(wildcard src/matrices/*.cpp) \
   $(wildcard src/*.cpp)         \

TESTS := $(wildcard tests/*.cpp)

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"

test1: release
	./build/apps/shallowWater --dt 0.1 --T 1.0 --Nx 100 --Ny 100 --ic 1

test3: release
	./build/apps/shallowWater --dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 3

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
