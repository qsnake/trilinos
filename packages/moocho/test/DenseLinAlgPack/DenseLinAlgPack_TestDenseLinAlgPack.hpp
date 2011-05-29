// ///////////////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_TestDenseLinAlgPack.hpp

#ifndef TEST_LIN_ALG_PACK_H
#define TEST_LIN_ALG_PACK_H

#include <iosfwd>

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack{
  namespace TestingPack {

    // Tests entire package
    bool TestDenseLinAlgPack(std::ostream* out);

    // Full automated testing with optional output
    bool TestVectorClass(std::ostream* out);
    bool TestVectorOp(std::ostream* out);
    bool TestGenMatrixClass(std::ostream* out);
    bool TestGenMatrixOp(std::ostream* out);

//		// Just output
//		void TestVectorBasicOp(std::ostream& out);
//		void TestGenMatrixBasicOp(std::ostream& out);
    void TestDenseLinAlgPackIO(std::istream& in, std::ostream& out);
//		void TestPivotVecMat(std::ostream& out);
  }
}

#endif	// TEST_LIN_ALG_PACK_H
