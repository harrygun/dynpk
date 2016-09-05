include ./include/ElVars

all: FisherMatrix WF1D WF1Dmean, WF1Dlensingmatrix WF1Dmulti WF2D WF1Dvectors ExampleHPD

FisherMatrix: FisherMatrix.o FisherHelperFunctions.o MatrixSetUpFunctions.o VectorFunctions.o
	$(CXX) $(EL_LINK_FLAGS) -o FisherMatrix FisherMatrix.o FisherHelperFunctions.o MatrixSetUpFunctions.o VectorFunctions.o $(EL_LIBS)

FisherMatrix.o: FisherMatrix.cpp ./include/FisherHelperFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c FisherMatrix.cpp

FisherHelperFunctions.o: ./src/FisherHelperFunctions.cpp ./include/MatrixSetUpFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/FisherHelperFunctions.cpp

WF1Dlensingmatrix: WF1Dlensingmatrix.o WienerFilter1Dlensematrix.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o
	$(CXX) $(EL_LINK_FLAGS) -o WF1Dlensingmatrix WF1Dlensingmatrix.o WienerFilter1Dlensematrix.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o $(EL_LIBS)

WF1Dlensingmatrix.o: WF1Dlensingmatrix.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp ./include/WienerFilter1Dlensematrix.hpp ./include/Statistics.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c WF1Dlensingmatrix.cpp

WF1Dmean: WF1Dmean.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o
	$(CXX) $(EL_LINK_FLAGS) -o WF1Dmean WF1Dmean.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o $(EL_LIBS)

WF1Dmean.o: WF1Dmean.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp ./include/WienerFilter1D.hpp ./include/Statistics.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c WF1Dmean.cpp

WF2D: WF2D.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o
	$(CXX) $(EL_LINK_FLAGS) -o WF2D WF2D.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o $(EL_LIBS)

WF2D.o: WF2D.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp ./include/WienerFilter1D.hpp ./include/Statistics.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c WF2D.cpp

WF1Dmulti: WF1Dmulti.o WienerFilter1Dmulti.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o
	$(CXX) $(EL_LINK_FLAGS) -o WF1Dmulti WF1Dmulti.o WienerFilter1Dmulti.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o $(EL_LIBS)

WF1Dmulti.o: WF1Dmulti.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp ./include/WienerFilter1Dmulti.hpp ./include/Statistics.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c WF1Dmulti.cpp

ExampleHPD: ExampleHPD.o
	$(CXX) $(EL_LINK_FLAGS) -o ExampleHPD ExampleHPD.o $(EL_LIBS)

ExampleHPD.o: ExampleHPD.cpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ExampleHPD.cpp

WF1Dvectors: WF1Dvectors.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o
	$(CXX) $(EL_LINK_FLAGS) -o WF1Dvectors WF1Dvectors.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o $(EL_LIBS)

WF1Dvectors.o: WF1Dvectors.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp ./include/WienerFilter1D.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c WF1Dvectors.cpp

ScalingTest: ScalingTest.o
	$(CXX) $(EL_LINK_FLAGS) -o ScalingTest ScalingTest.o $(EL_LIBS)

ScalingTest.o: ScalingTest.cpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ScalingTest.cpp

WF1D: WF1D.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o
	$(CXX) $(EL_LINK_FLAGS) -o WF1D WF1D.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o Statistics.o $(EL_LIBS)

WF1D.o: WF1D.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp ./include/WienerFilter1D.hpp ./include/Statistics.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c WF1D.cpp

Statistics.o: ./src/Statistics.cpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/Statistics.cpp

WienerFilter1Dmulti.o: ./src/WienerFilter1Dmulti.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/WienerFilter1Dmulti.cpp

WienerFilter1D.o: ./src/WienerFilter1D.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/WienerFilter1D.cpp

WienerFilter1Dlensematrix.o: ./src/WienerFilter1Dlensematrix.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/WienerFilter1Dlensematrix.cpp

VectorFunctions.o: ./src/VectorFunctions.cpp ./include/VectorFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/VectorFunctions.cpp

MatrixSetUpFunctions.o: ./src/MatrixSetUpFunctions.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/MatrixSetUpFunctions.cpp

clean:
	rm *o WF1D ScalingTest