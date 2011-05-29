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

#include "SundanceIQI_ScalVecBLF_UVj_Facet.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

IQI_ScalVecBLF_UVj_Facet::IQI_ScalVecBLF_UVj_Facet( int spatialDim ,
						    const CellType & maxCellType ,
						    int dim ,
						    const CellType &cellType ,
						    const BasisFamily &testBasis ,
						    int test_component ,
						    const BasisFamily &unkBasis ,
						    const QuadratureFamily &quad ,
						    const ParameterList& verbParams ):
  ElementIntegralBilinearFormFacet( spatialDim ,
				    maxCellType ,
				    dim ,
				    cellType ,
				    testBasis ,
				    unkBasis ,
				    quad ,
				    verbParams ),
  j_test_( test_component ) 
{
  TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_ScalVecBLF_UVj_Cell is not implemented" );
}

void IQI_ScalVecBLF_UVj_Facet::evaluate( CellJacobianBatch& JTrans,
					CellJacobianBatch& JVol ,
					const Array<int> &facetIndex ,
					const double* const coeff,
					RefCountPtr<Array<double> >& A) const
{
  TEST_FOR_EXCEPTION(true,
		     InternalError,
		     "IQI_ScalVecBLF_UVj_Facet is not implemented" );
}

