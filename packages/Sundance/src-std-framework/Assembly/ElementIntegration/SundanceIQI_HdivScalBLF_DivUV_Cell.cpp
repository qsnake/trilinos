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

#include "SundanceIQI_HdivScalBLF_DivUV_Cell.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

IQI_HdivScalBLF_DivUV_Cell::IQI_HdivScalBLF_DivUV_Cell( int spatialDim ,
						      const CellType & maxCellType ,
						      const BasisFamily &testBasis ,
						      const BasisFamily &unkBasis ,
						      const QuadratureFamily &quad ,
						      const ParameterList& verbParams ):
  ElementIntegralBilinearFormCell( spatialDim ,
				   maxCellType ,
				   testBasis ,
				   unkBasis ,
				   quad ,
				   verbParams )
{
  TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_HdivScalBLF_DivUV_Cell is not implemented" );
}

void IQI_HdivScalBLF_DivUV_Cell::evaluate( CellJacobianBatch& JTrans,
					  const double* const coeff,
					  RefCountPtr<Array<double> >& A) const
{
  TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_HdivScalBLF_DivUV_Cell is not implemented" );
}

