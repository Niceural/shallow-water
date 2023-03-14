CC = g++
CFLAGS = -O2 -Wall -Wshadow
LIBS = -lblas -lboost_unit_test_framework

test1:
	${CC} ${CFLAGS} -o target/test1 src/main.cpp -lboost_program_options

GeneralMatrix.o:
	${CC} ${CFLAGS} -o target/GeneralMatrix.o -c src/matrices/GeneralMatrix.cpp

SquareBandedMatrix.o:
	${CC} ${CFLAGS} -o target/SquareBandedMatrix.o -c src/matrices/SquareBandedMatrix.cpp

matrices: GeneralMatrix.o SquareBandedMatrix.o

ShallowWater.o:
	${CC} ${CFLAGS} -o target/ShallowWater.o -c src/ShallowWater.cpp

build: matrices ShallowWater.o

initialCondition: matrices ShallowWater.o
	${CC} ${CFLAGS} -o target/initialCondition.o -c tests/initialCondition.cpp
	${CC} ${CFLAGS} -o target/initialCondition target/GeneralMatrix.o target/SquareBandedMatrix.o target/ShallowWater.o target/initialCondition.o ${LIBS}
	./target/initialCondition

centralDifference: matrices ShallowWater.o
	${CC} ${CFLAGS} -o target/centralDifference.o -c tests/centralDifference.cpp
	${CC} ${CFLAGS} -o target/centralDifference target/GeneralMatrix.o target/SquareBandedMatrix.o target/ShallowWater.o target/centralDifference.o ${LIBS}
	./target/centralDifference

test: initialCondition centralDifference
