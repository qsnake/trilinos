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

#include "SundanceIQI_HdivLF_DivV_Cell.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HDIV_TET_I1_FEM.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

 
IQI_HdivLF_DivV_Cell::IQI_HdivLF_DivV_Cell( int spatialDim,
					    const CellType& maxCellType,
					    const BasisFamily& testBasis,
					    const QuadratureFamily& quad,
					    const ParameterList& verbParams):
  ElementIntegralLinearFormCell( spatialDim ,
				 maxCellType ,
				 testBasis ,
				 quad ,
				 verbParams ),
  DivV_( testBasis.nReferenceDOFs(maxCellType,maxCellType) , quad.getNumPoints(maxCellType) ),
  QP_( quad.getNumPoints( maxCellType ) , spatialDim ),
  QW_( quad.getNumPoints( maxCellType ) )
{
  // bypass testBasis.refEval and use Intrepid to fill DivV
  // warning: only works on tets right now.
  Intrepid::Basis_HDIV_TET_I1_FEM<double,Intrepid::FieldContainer<double> > myBasis;

  // move quadrature points into a field container
  Array<Point> qpSundance;
  Array<double> qwSundance;
  quad.getPoints( maxCellType , qpSundance , qwSundance );

  for (int i=0;i<(int)qpSundance.size();i++) {
    for (int j=0;j<3;j++) {
      QP_(i,j) = qpSundance[i][j];
    }
    QW_(i)=qwSundance[i];
  }

  // now tabulate the divergences 
  myBasis.getValues( DivV_ , QP_ , Intrepid::OPERATOR_DIV );

}

void IQI_HdivLF_DivV_Cell::evaluate( CellJacobianBatch& JTrans,
				     const double* const coeff,
				     RefCountPtr<Array<double> >& A) const
{
  const int nqp = quad().getNumPoints( maxCellType() );
  const int ncell = JTrans.numCells();
  const int nbf = testBasis().nReferenceDOFs(maxCellType(),maxCellType());

  // wrap A into a rank 2 field container
  Teuchos::Array<int> Aindices(2);
  Aindices[0] = nqp;
  Aindices[1] = nbf;
  Intrepid::FieldContainer<double> AFC(Aindices,*A);

  // wrap coeff into another field container.
  // by way of a Teuchos array
  // by way of discarding the const

  /* this surprisingly doesn't depend on the Jacobian ! */

  for (int c=0;c<ncell;c++) {
    for (int bf=0;bf<nbf;bf++) {
      AFC(c,bf) = 0.0;
      for (int q=0;q<nqp;q++) {
	AFC(c,bf) += QW_(q) * coeff[c*nqp+q] * DivV_(bf,q );
      }
    }
  }

  return;
}
