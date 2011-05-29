//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceIQI_VecVecBLF_UiVj_Cell.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

IQI_VecVecBLF_UiVj_Cell::IQI_VecVecBLF_UiVj_Cell( int spatialDim ,
						      const CellType & maxCellType ,
						      const BasisFamily &testBasis ,
						      int test_component ,
						      const BasisFamily &unkBasis ,
						      int unk_component ,
						      const QuadratureFamily &quad ,
						      const ParameterList& verbParams ):
  ElementIntegralBilinearFormCell( spatialDim ,
				   maxCellType ,
				   testBasis ,
				   unkBasis ,
				   quad ,
				   verbParams ),
  j_test_( test_component ) , i_trial_( unk_component )
{
  TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_VecVecBLF_UiVj_Cell is not implemented" );
}

void IQI_VecVecBLF_UiVj_Cell::evaluate( CellJacobianBatch& JTrans,
					  const double* const coeff,
					  RefCountPtr<Array<double> >& A) const
{
  TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_VecVecBLF_UiVj_Cell is not implemented" );
}

