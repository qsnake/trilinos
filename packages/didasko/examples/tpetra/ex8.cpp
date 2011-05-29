// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
// 
// ***********************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

/* A simple matlab script to convert pdetool-generated meshes to this format:
 
[p,e,t] = initmesh('lshapeg');
fid = fopen('mesh.dat', 'w');
fprintf(fid, 'NumVertices %d\n', size(p, 2));
fprintf(fid, 'NumElements %d\n', size(t, 2));
fprintf(fid, 'NumUnknownsPerElement 3\n');
fprintf(fid, 'NumDimensions 2\n');
fprintf(fid, 'ElementData\n');
for i=1:size(t,2)
  fprintf(fid, '%d %d %d\n', t(1, i) - 1, t(2, i) - 1, t(3, i) - 1);
end
fprintf(fid, 'EndElementData\n');
fprintf(fid, 'VertexData\n');
for i=1:size(p, 2)
  fprintf(fid, '%e %e\n', p(1, i), p(2, i));
end
fprintf(fid, 'EndVertexData\n');
fclose(fid);
*/
 
// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 28-Nov-05

const int ElementInfoSize = 3;
const int ElementInfoPackSize = ElementInfoSize * sizeof(int) / sizeof(char);

class ElementInfo
{
public:
  inline ElementInfo(const int value = 0)
  {
    for (int i = 0 ; i < ElementInfoSize ; ++i)
      data[i] = value;
  }

  inline ElementInfo(const vector<int> rhs)
  {
    assert (ElementInfoSize == (signed)rhs.size());
    for (int i = 0 ; i < ElementInfoSize ; ++i)
      data[i] = rhs[i];
  }

  inline ElementInfo& operator=(const int rhs)
  {
    for (int i = 0 ; i < ElementInfoSize ; ++i)
      data[i] = rhs;
    return *this;
  }

  inline ElementInfo& operator=(const ElementInfo rhs)
  {
    for (int i = 0 ; i < ElementInfoSize ; ++i)
      data[i] = rhs.data[i];
    return *this;
  }

  inline ElementInfo& operator+=(const ElementInfo rhs)
  {
    for (int i = 0 ; i < ElementInfoSize ; ++i)
      data[i] += rhs.data[i];
    return *this;
  }

  inline ElementInfo operator*(const int i) const
  {
    return(ElementInfo(i));
  }

  inline int operator[](const int i) const
  {
    return(data[i]);
  }

  inline int Magnitude() const
  {
    throw(1);
  }

  int data[ElementInfoSize];
};

inline ostream& operator<<(ostream& os, const ElementInfo& obj)
{
  for (int i = 0 ; i < ElementInfoSize ; ++i)
    os << obj.data[i] << " ";
  return(os);
}

namespace Teuchos {
  template<>
  struct ScalarTraits<ElementInfo>
  {
    typedef int magnitudeType;
    //! Returns representation of zero for this scalar type.
    static inline ElementInfo zero()                     
    { 
      return(ElementInfo(0));
    }

    //! Returns representation of one for this scalar type.
    static inline ElementInfo one()                
    { 
      return(ElementInfo(1));
    }
  };

} // namespace Teuchos

#ifdef HAVE_MPI
namespace Tpetra 
{
  template<>
  struct MpiTraits<ElementInfo> 
  {
    static inline MPI_Datatype datatype() {return(MPI_CHAR);};
    static inline int count(int userCount) 
    {
      return(userCount * (ElementInfoPackSize));
    }

    static MPI_Op sumOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)ElementInfoSumOp, true, &myOp);
      return(myOp);
    };
    static MPI_Op maxOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)ElementInfoMaxOp, true, &myOp);
      return(myOp);
    };
    static MPI_Op minOp() {
      MPI_Op myOp;
      MPI_Op_create((MPI_User_function*)ElementInfoMinOp, true, &myOp);
      return(myOp);
    };

    static void ElementInfoSumOp(ElementInfo* in, ElementInfo* inout, 
                                  int* len, MPI_Datatype* dptr) 
    {
      for(int i = 0; i < *len / ElementInfoPackSize ; i++) 
      {
        inout[i] += in[i];
      }
    };

    static void ElementInfoMaxOp(ElementInfo* in, ElementInfo* inout, 
                                  int* len, MPI_Datatype* dptr) 
    {
      for(int i = 0; i < *len / ElementInfoPackSize ; i++) 
      {
        if (in[i].Magnitude() > inout[i].Magnitude())
          inout[i] = in[i];
      }
    };

    static void ElementInfoMinOp(ElementInfo* in, ElementInfo* inout, 
                                  int* len, MPI_Datatype* dptr) 
    {
      for(int i = 0; i < *len / ElementInfoPackSize ; i++) 
      {
        if (in[i].Magnitude() < inout[i].Magnitude())
          inout[i] = in[i];
      }
    };
  };
} // namespace Tpetra
#endif

template<class OrdinalType, class ScalarType>
class BaseGrid
{
  public:
    virtual ~BaseGrid() {}

    virtual Tpetra::ElementSpace<OrdinalType>& getElementSpace() = 0;

    virtual Tpetra::VectorSpace<OrdinalType, ScalarType>& getVectorElementSpace() = 0;

    virtual Tpetra::ElementSpace<OrdinalType>& getVertexSpace() = 0;

    virtual Tpetra::VectorSpace<OrdinalType, ScalarType>& getVectorVertexSpace() = 0;

    virtual Tpetra::ElementSpace<OrdinalType>& getOverlappingVertexSpace() = 0;

    virtual Tpetra::VectorSpace<OrdinalType, ScalarType>& getVectorOverlappingVertexSpace() = 0;

    virtual Tpetra::Vector<OrdinalType, ElementInfo>& getGridInfo() = 0;

    virtual double getCoord(const int dim, const int vertex) = 0;

    virtual int getNumUnknownsPerFiniteElement() = 0;

    virtual void write(Tpetra::Vector<OrdinalType, ScalarType>& v) = 0;
};

template<class OrdinalType, class ScalarType>
class MATLABGrid
{
  public:
    MATLABGrid(const Tpetra::Comm<OrdinalType, ScalarType>& CommST, 
               const Tpetra::Comm<OrdinalType, ElementInfo>& CommEI,
               char* FileName) :
      FileName_(FileName),
      CommST_(CommST),
      CommEI_(CommEI)
    {
#ifdef HAVE_MPI
      platformOT_ = new Tpetra::MpiPlatform <OrdinalType, OrdinalType>(MPI_COMM_WORLD);
      platformST_ = new Tpetra::MpiPlatform <OrdinalType, ScalarType>(MPI_COMM_WORLD);
      platformEI_ = new Tpetra::MpiPlatform <OrdinalType, ElementInfo>(MPI_COMM_WORLD);
#else
      platformOT_ = new Tpetra::SerialPlatform <OrdinalType, OrdinalType>;
      platformST_ = new Tpetra::SerialPlatform <OrdinalType, ScalarType>;
      platformEI_ = new Tpetra::SerialPlatform <OrdinalType, ElementInfo>;
#endif

      FILE* fp;

      OrdinalType NumGlobalVertices;
      OrdinalType NumGlobalFiniteElements;
      OrdinalType NumDimensions;
      char FiniteElementType[10];
      char what[80];

      if ((fp = fopen(FileName, "r")) == 0)
      {
        cerr << "Error opening file." << endl;
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        exit(EXIT_FAILURE);
      }
      fscanf(fp, "%s %d", &what, &NumGlobalVertices);
      fscanf(fp, "%s %d", &what, &NumGlobalFiniteElements);
      fscanf(fp, "%s %d", &what, &NumUnknownsPerFiniteElement_);
      fscanf(fp, "%s %d", &what, &NumDimensions);

      if (CommST_.getMyImageID() == 0)
      {
        cout << "NumGlobalVertices           = " << NumGlobalVertices << endl;
        cout << "NumGlobalFiniteElements     = " << NumGlobalFiniteElements << endl;
        cout << "NumUnknownsPerFiniteElement = " << NumUnknownsPerFiniteElement_ << endl;
        cout << "NumDimensions               = " << NumDimensions << endl;
      }

      // build a linear map for the vertices

      VertexSpace_ = new Tpetra::ElementSpace<OrdinalType>(NumGlobalVertices, 0, *platformOT_);
      VectorVertexSpace_ = new Tpetra::VectorSpace<OrdinalType, ScalarType>(*VertexSpace_, *platformST_);

      // read through the file and count how many local elements we have

      fscanf(fp, "%s", what);

      int NumMyFiniteElements;
      vector<OrdinalType> MyGlobalFiniteElements;

      int offset;
      offset = ftell(fp);

      NumMyFiniteElements = 0;
      for (int i = 0 ; i < NumGlobalFiniteElements ; ++i)
      {
        int unk;
        bool ok = false;
        for (int j = 0 ; j < NumUnknownsPerFiniteElement_ ; ++j)
        {
          fscanf(fp, "%d", &unk);
          if (VertexSpace_->isMyGID(unk))
            ok = true;
        }
        if (ok)
        {
          MyGlobalFiniteElements.push_back(i);
          ++NumMyFiniteElements;
        }
      }

      cout << "Local finite elements = " << NumMyFiniteElements << endl;

      // build the map corresponding to elements and prepare the
      // map for all the vertices information required to assemble
      // the local matrices. This means all the vertices 
      // included in the local elements.

      FiniteElementSpace_ = new Tpetra::ElementSpace<OrdinalType>(-1, NumMyFiniteElements, MyGlobalFiniteElements, 0, *platformOT_);
      VectorFiniteElementSpace_ = new Tpetra::VectorSpace<OrdinalType, ElementInfo>(*FiniteElementSpace_, *platformEI_);

      GridInfo_ = new Tpetra::Vector<OrdinalType, ElementInfo>(*VectorFiniteElementSpace_);

      map<int, bool> list;

      fseek(fp, offset, SEEK_SET);

      vector<int> nodes(NumUnknownsPerFiniteElement_);

      NumMyFiniteElements = 0;
      for (int i = 0 ; i < NumGlobalFiniteElements ; ++i)
      {
        int node;
        bool ok = false;
        for (int j = 0 ; j < NumUnknownsPerFiniteElement_ ; ++j)
        {
          fscanf(fp, "%d", &(nodes[j]));
          if (VertexSpace_->isMyGID(nodes[j]))
            ok = true;
        }
        if (ok)
        {
          (*GridInfo_)[NumMyFiniteElements] = ElementInfo(nodes);
          ++NumMyFiniteElements;
          for (int j = 0 ; j < NumUnknownsPerFiniteElement_ ; ++j)
            list[nodes[j]] = true;
        }
      }

      int NumMyOverlappingVertices = list.size();
      vector<OrdinalType> MyGlobalOverlappingVertices(NumMyOverlappingVertices);

      map<int, bool>::iterator iter;

      int count = 0;
      for (iter = list.begin() ; iter != list.end() ; ++iter)
      {
        MyGlobalOverlappingVertices[count++] = iter->first;
      }

      cout << "NumMyOverlappingVertices = " << NumMyOverlappingVertices << endl;

      // now reads the coordinates

      // FIXME: the following is `double'
      OverlappingVertexSpace_ = new Tpetra::ElementSpace<OrdinalType>(-1, NumMyOverlappingVertices, MyGlobalOverlappingVertices, 0, *platformOT_);
      VectorOverlappingVertexSpace_ = new Tpetra::VectorSpace<OrdinalType, ScalarType>(*OverlappingVertexSpace_, *platformST_);

      XCoord_ = new Tpetra::Vector<OrdinalType, ScalarType>(*VectorOverlappingVertexSpace_);
      YCoord_ = new Tpetra::Vector<OrdinalType, ScalarType>(*VectorOverlappingVertexSpace_);

      fscanf(fp, "%s", what); // EndElementData
      fscanf(fp, "%s", what); // VertexData
      vector<double> coord(NumDimensions);

      NumMyOverlappingVertices = 0;
      for (int i = 0 ; i < NumGlobalVertices ; ++i)
      {
        // should skip data if not owned
        for (int j = 0 ; j < NumDimensions ; ++j)
        {
          fscanf(fp, "%lf", &(coord[j]));
        }

        if (OverlappingVertexSpace_->isMyGID(i))
        {
          (*XCoord_)[NumMyOverlappingVertices] = coord[0];
          (*YCoord_)[NumMyOverlappingVertices] = coord[1];
          ++NumMyOverlappingVertices;
        }
      }

      fscanf(fp, "%s", what); // VertexData
      fclose(fp);

      //cout << *GridInfo_;
      //cout << *XCoord_;
      //cout << *YCoord_;
    }

    virtual ~MATLABGrid() {}

    virtual Tpetra::ElementSpace<OrdinalType>& getFiniteElementSpace()
    {
      return(*FiniteElementSpace_);
    }

    virtual Tpetra::VectorSpace<OrdinalType, ElementInfo>& getVectorFiniteElementSpace()
    {
      return(*VectorFiniteElementSpace_);
    }

    virtual Tpetra::ElementSpace<OrdinalType>& getVertexSpace()
    {
      return(*VertexSpace_);
    }

    virtual Tpetra::VectorSpace<OrdinalType, ScalarType>& getVectorVertexSpace()
    {
      return(*VectorVertexSpace_);
    }

    virtual Tpetra::ElementSpace<OrdinalType>& getOverlappingVertexSpace()
    {
      return(*OverlappingVertexSpace_);
    }

    virtual Tpetra::VectorSpace<OrdinalType, ScalarType>& getVectorOverlappingVertexSpace() 
    {
      return(*VectorOverlappingVertexSpace_);
    }

    virtual Tpetra::Vector<OrdinalType, ElementInfo>& getGridInfo()
    {
      return(*GridInfo_);
    }

    virtual double getCoord(const int dim, const int vertex)
    {
      if (dim == 0)
        return((*XCoord_)[vertex]);
      else
        return((*YCoord_)[vertex]);
    }

    virtual int getNumUnknownsPerFiniteElement()
    {
      return(NumUnknownsPerFiniteElement_);
    }

    virtual void write(Tpetra::Vector<OrdinalType, ScalarType>& v)
    {
      if (v.vectorSpace() != *VectorOverlappingVertexSpace_)
      {
        cerr << "Input vector is not based on OverlappingVertexSpace" << endl;
        return;
      }

      // define a space that contains all the elements on 
      // processor 0, then gather the vector v on processor 0,
      // and write it on file in the specified format.

      int NumOutputElements = 0;
      if (CommST_.getMyImageID() == 0)
        NumOutputElements = VertexSpace_->getNumGlobalElements();

      Tpetra::ElementSpace<OrdinalType> OutputSpace(-1, NumOutputElements, 0, *platformOT_);
      Tpetra::VectorSpace<OrdinalType, ScalarType> VectorOutputSpace(OutputSpace, *platformST_);
      Tpetra::Vector<OrdinalType, ScalarType> Output(VectorOutputSpace);

      Tpetra::Import<OrdinalType> Importer(*OverlappingVertexSpace_, OutputSpace);

      Output.doImport(v, Importer, Tpetra::Insert);

      // at this point I have everything I need on processor 0,
      // because I can read the mesh information from the input file.

      if (CommST_.getMyImageID() == 0)
      {
        FILE* fp;
        FILE* fp_output;

        OrdinalType NumGlobalVertices;
        OrdinalType NumGlobalFiniteElements;
        OrdinalType NumDimensions;
        char FiniteElementType[10];
        char what[80];

        if ((fp = fopen(FileName_.c_str(), "r")) == 0)
        {
          cerr << "Error opening file." << endl;
#ifdef HAVE_MPI
          MPI_Finalize();
#endif
          exit(EXIT_FAILURE);
        }

        if ((fp_output = fopen("medit.mesh", "w")) == 0)
        {
          cerr << "Error opening file." << endl;
#ifdef HAVE_MPI
          MPI_Finalize();
#endif
          exit(EXIT_FAILURE);
        }

        fscanf(fp, "%s %d", &what, &NumGlobalVertices);
        fscanf(fp, "%s %d", &what, &NumGlobalFiniteElements);
        fscanf(fp, "%s %d", &what, &NumUnknownsPerFiniteElement_);
        fscanf(fp, "%s %d", &what, &NumDimensions);

        fscanf(fp, "%s", what);

        fprintf(fp_output, "MeshVersionFormatted 1\n");
        fprintf(fp_output, "Dimension 3\n");
        fprintf(fp_output, "# mesh from Tpetra/Trilinos\n");
        fprintf(fp_output, "Triangles %d\n", NumGlobalFiniteElements);

        int unk;
        for (int i = 0 ; i < NumGlobalFiniteElements ; ++i)
        {
          for (int j = 0 ; j < NumUnknownsPerFiniteElement_ ; ++j)
          {
            fscanf(fp, "%d", &unk);
            fprintf(fp_output, "%d ", unk + 1);
          }
          fprintf(fp_output, "0\n");
        }

        fscanf(fp, "%s", what); // EndElementData
        fscanf(fp, "%s", what); // VertexData

        fprintf(fp_output, "Vertices %d\n", NumGlobalVertices);

        assert (NumDimensions == 2); // code to be changed is below
        double coord;
        for (int i = 0 ; i < NumGlobalVertices ; ++i)
        {
          for (int j = 0 ; j < NumDimensions ; ++j)
          {
            fscanf(fp, "%lf", &coord);
            fprintf(fp_output, "%lf ", coord);
          }
          fprintf(fp_output, "0.0 1\n");
        }

        fscanf(fp, "%s", what); // VertexData
        fclose(fp);

        fprintf(fp_output, "End\n");
        fclose(fp_output);

        // ======== //
        // .bb file //
        // ======== //

        fp_output = fopen("medit.bb", "w");

        fprintf(fp_output, "3 1 %d 2\n", NumGlobalVertices);
                
        for (int i = 0 ; i < NumGlobalVertices ; ++i)
        {
          fprintf(fp_output, "%e\n", Output[i]);
        }

        fclose(fp_output);
      }
      
      // synchronize the processors and return
      CommST_.barrier();
    }

  private:
    string FileName_;
    
    const Tpetra::Comm<OrdinalType, ScalarType>& CommST_;
    const Tpetra::Comm<OrdinalType, ElementInfo>& CommEI_;
    Tpetra::Platform<OrdinalType, OrdinalType>* platformOT_;
    Tpetra::Platform<OrdinalType, ScalarType>*  platformST_;
    Tpetra::Platform<OrdinalType, ElementInfo>* platformEI_;

    Tpetra::ElementSpace<OrdinalType>* FiniteElementSpace_;
    Tpetra::VectorSpace<OrdinalType, ElementInfo>* VectorFiniteElementSpace_;

    Tpetra::ElementSpace<OrdinalType>* VertexSpace_;
    Tpetra::VectorSpace<OrdinalType, ScalarType>* VectorVertexSpace_;

    Tpetra::ElementSpace<OrdinalType>* OverlappingVertexSpace_;
    Tpetra::VectorSpace<OrdinalType, ScalarType>* VectorOverlappingVertexSpace_;

    Tpetra::Vector<OrdinalType, ElementInfo>* GridInfo_;
    Tpetra::Vector<OrdinalType, ScalarType>* XCoord_;
    Tpetra::Vector<OrdinalType, ScalarType>* YCoord_;

    int NumUnknownsPerFiniteElement_;
};



typedef int OrdinalType;
typedef double ScalarType;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
  Tpetra::MpiComm<OrdinalType, ElementInfo> CommEI(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
  Tpetra::SerialComm<OrdinalType, ElementInfo> CommEI;
#endif

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformOT(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformST(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ElementInfo> platformEI(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformOT;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformST;
  const Tpetra::SerialPlatform <OrdinalType, ElementInfo> platformEI;
#endif

  if (argc == 1)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  MATLABGrid<OrdinalType, ScalarType> LinearDistributedGrid(Comm, CommEI, argv[1]);

  // need to call ParMETIS and redistribute the grid objects here...

  // now build the RowSpace, in this case is simply VertexSpace
  Tpetra::ElementSpace<OrdinalType> RowSpace(LinearDistributedGrid.getVertexSpace());
  Tpetra::VectorSpace<OrdinalType, ScalarType> VectorRowSpace(LinearDistributedGrid.getVectorVertexSpace());

  // ================ //
  // Setup the matrix //
  // ================ //
  
  OrdinalType NumUnknownsPerFiniteElement = LinearDistributedGrid.getNumUnknownsPerFiniteElement();
  Teuchos::SerialDenseMatrix<OrdinalType, ScalarType> localMatrix(NumUnknownsPerFiniteElement, NumUnknownsPerFiniteElement);
  localMatrix.putScalar(ScalarZero);
  // diagonal now 
#if 0
  if (NumUnknownsPerFiniteElement == 4)
  {
    localMatrix(0, 0) =  2.0; localMatrix(0, 1) = -1.0; localMatrix(0, 2) =  0.0; localMatrix(0,3) = -1.0;
    localMatrix(1, 0) = -1.0; localMatrix(1, 1) =  2.0; localMatrix(1, 2) = -1.0; localMatrix(1,3) =  0.0;
    localMatrix(2, 0) =  0.0; localMatrix(2, 1) = -1.0; localMatrix(2, 2) =  2.0; localMatrix(2,3) = -1.0;
    localMatrix(3, 0) = -1.0; localMatrix(3, 1) =  0.0; localMatrix(3, 2) = -1.0; localMatrix(3,3) =  2.0;
  }
#endif
  assert (NumUnknownsPerFiniteElement == 3);

  // Allocate the matrix, based on VectorRowSpace
  
  Tpetra::CisMatrix<OrdinalType,ScalarType> matrix(VectorRowSpace);

  // add zero diagonal

  for (OrdinalType LID = OrdinalZero ; LID < RowSpace.getNumMyElements() ; ++LID)
  {
    OrdinalType GID = RowSpace.getGID(LID);
    matrix.submitEntry(Tpetra::Insert, GID, ScalarZero, GID);
  }

  Tpetra::Vector<OrdinalType, ElementInfo>& GridInfo = LinearDistributedGrid.getGridInfo();

  // Loop over all the (locally owned) elements
  
  for (OrdinalType FEID = OrdinalZero ; FEID < GridInfo.vectorSpace().getNumMyEntries() ; ++FEID)
  {
    vector<OrdinalType> GIDs(NumUnknownsPerFiniteElement);
    vector<OrdinalType> LIDs(NumUnknownsPerFiniteElement);
    // FIXME: This is `double'
    vector<ScalarType> x(NumUnknownsPerFiniteElement);
    vector<ScalarType> y(NumUnknownsPerFiniteElement);

    // load the global IDs for vertices

    for (OrdinalType i = OrdinalZero ; i < NumUnknownsPerFiniteElement ; ++i)
    {
      GIDs[i] = GridInfo[FEID][i];
      LIDs[i] = LinearDistributedGrid.getOverlappingVertexSpace().getLID(GIDs[i]);
      x[i] = LinearDistributedGrid.getCoord(0, LIDs[i]);
      y[i] = LinearDistributedGrid.getCoord(1, LIDs[i]);
    }

    // build the local matrix, in this case localMatrix
    localMatrix(0,0) = (y[2]-y[1])*(y[2]-y[1]) + (x[2]-x[1])*(x[2]-x[1]);
    localMatrix(0,1) = (y[2]-y[1])*(y[0]-y[2]) + (x[2]-x[1])*(x[0]-x[2]);
    localMatrix(0,2) = (y[1]-y[0])*(y[2]-y[1]) + (x[1]-x[0])*(x[2]-x[1]);
    
    localMatrix(1,0) = (y[2]-y[1])*(y[0]-y[2]) + (x[2]-x[1])*(x[0]-x[2]);
    localMatrix(1,1) = (y[2]-y[0])*(y[2]-y[0]) + (x[2]-x[0])*(x[2]-x[0]);
    localMatrix(1,2) = (y[0]-y[2])*(y[1]-y[0]) + (x[0]-x[2])*(x[1]-x[0]);

    localMatrix(2,0) = (y[1]-y[0])*(y[2]-y[1]) + (x[1]-x[0])*(x[2]-x[1]);
    localMatrix(2,1) = (y[0]-y[2])*(y[1]-y[0]) + (x[0]-x[2])*(x[1]-x[0]);
    localMatrix(2,2) = (y[1]-y[0])*(y[1]-y[0]) + (x[1]-x[0])*(x[1]-x[0]);

    ScalarType det_J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    det_J = 2*det_J;

    for (int ii = 0; ii < 3; ii++) 
      for (int jj = 0; jj < 3; jj++) 
        localMatrix(ii,jj) = localMatrix(ii,jj) / det_J;

    // submit entries of localMatrix into the matrix

    for (OrdinalType row = OrdinalZero ; row < NumUnknownsPerFiniteElement ; ++row)
    {
      // We skip non-locally owned vertices
      if (RowSpace.isMyGID(GIDs[row]))
      {
        for (OrdinalType col = OrdinalZero ; col < NumUnknownsPerFiniteElement ; ++col)
        {
          matrix.submitEntry(Tpetra::Add, GIDs[row], localMatrix(row, col), GIDs[col]);
        }
      }
    }
  }

  matrix.fillComplete();

  Tpetra::Vector<OrdinalType, ScalarType> X(VectorRowSpace);
  Tpetra::Vector<OrdinalType, ScalarType> Y(VectorRowSpace);

  X.setAllToScalar(ScalarOne);

  matrix.apply(X, Y, false);

  ScalarType norm2 = Y.norm2();
  if (Comm.getMyImageID() == 0)
    cout << "||A * 1||_2 = " << norm2 << endl;

  // example of visualization
  
  for (OrdinalType i = OrdinalZero ; i < RowSpace.getNumMyElements() ; ++i)
  {
    X[i] = (ScalarType) Comm.getMyImageID();
  }

  Tpetra::Vector<OrdinalType, ScalarType> OverlappingX(LinearDistributedGrid.getVectorOverlappingVertexSpace());

  Tpetra::Import<OrdinalType> Importer(RowSpace, LinearDistributedGrid.getOverlappingVertexSpace());
  OverlappingX.doImport(X, Importer, Tpetra::Insert);

  LinearDistributedGrid.write(OverlappingX);
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
}
