/*
 * SundanceCurveIntagralCalc.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: benk
 */

#include "SundanceCurveIntegralCalc.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace Sundance;

CurveIntegralCalc::CurveIntegralCalc() {
	//nothing to do
}

void CurveIntegralCalc::getCurveQuadPoints(CellType  maxCellType ,
							   int maxCellLID ,
							   const Mesh& mesh ,
							   const ParametrizedCurve& paramCurve,
							   const QuadratureFamily& quad ,
							   Array<Point>& curvePoints ,
							   Array<Point>& curveDerivs ,
							   Array<Point>& curveNormals ){

	int verb = 0;
	Tabs tabs;

    // get all the points from the maxDimCell
	int nr_point = mesh.numFacets( mesh.spatialDim() , 0 , 0);
	Array<Point> maxCellsPoints(nr_point);
	int tmp_i , point_LID;

	SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints nr points per cell: " << nr_point)
	for (int jj = 0 ; jj < nr_point ; jj++){
		point_LID = mesh.facetLID( mesh.spatialDim() , maxCellLID , 0 , jj , tmp_i );
		maxCellsPoints[jj] = mesh.nodePosition(point_LID);
		SUNDANCE_MSG3(verb, tabs << " max cell point p[" << jj << "]:"<< maxCellsPoints[jj]);
	}

	// call the simple line method
	CurveIntegralCalc::getCurveQuadPoints_line( maxCellType , maxCellsPoints , paramCurve,
							    quad , curvePoints , curveDerivs , curveNormals);

	// later here could be more sophisticated method called instead of the simple line/plane method

}

void CurveIntegralCalc::getCurveQuadPoints_line(
                                   CellType  maxCellType ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals ){

	int verb = 0;
	Tabs tabs;

    // select the dimension of the curve
	switch (paramCurve.getCurveDim()){
	  case 1: {

		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, curve has ONE dimnesion")
		  // 1D curve integral in 2D or 3D

		  // first determine the two intersection points with the edges
		  Point startPoint(0.0,0.0);
		  Point endPoint(0.0,0.0);
	      int nrPoints , total_points = 0 ;

		  switch (maxCellType){
		    case QuadCell:{
		    	// loop over each edge and detect intersection point
		    	// there can be only one
		    	TEST_FOR_EXCEPTION( cellPoints.size() != 4 ,
		    			RuntimeError ," CurveIntegralCalc::getCurveQuadPoints , QuadCell must have 4 points" );
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, on QuadCell")
		    	Array<Point> result(0);
		    	int edegIndex[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };

		    	// loop over the edges
		    	for (int jj = 0 ; jj < 4 ; jj++ ){
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints > 0){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
			    	}
					total_points += nrPoints;
					SUNDANCE_MSG3(verb, tabs << "ind:" << jj << ", nr Int points :" << nrPoints
							<< " , start:" << startPoint << ", end:"<< endPoint);
					TEST_FOR_EXCEPTION( nrPoints > 1 ,
							RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , QuadCell one edge " << jj << " , can have only one intersection point" );
		    	}
		    	// test if we have too much intersection points
		    	TEST_FOR_EXCEPTION( total_points > 2 ,RuntimeError , " CurveIntegralCalc::getCurveQuadPoints total_points > 2 : " << total_points );
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, end finding intersection points")
		    } break;
		    case TriangleCell:{
		    	TEST_FOR_EXCEPTION( cellPoints.size() != 3 , RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , TriangleCell must have 3 points" );
		    	Array<Point> result;
		    	int edegIndex[3][2] = { {0,1} , {0,2} , {1,2} };

		    	// loop over the edges
		    	for (int jj = 0 ; jj < 3 ; jj++ ){
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints > 0){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
			    	}
					total_points += nrPoints;
					TEST_FOR_EXCEPTION( nrPoints > 1 ,
							RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , TriangleCell one edge " << jj << " , can have only one intersection point" );
		    	}
		    	// test if we have too much intersection points
		    	TEST_FOR_EXCEPTION( total_points > 2 ,RuntimeError , " CurveIntegralCalc::getCurveQuadPoints total_points > 2 : " << total_points );
		    } break;
		    default : {
		    	TEST_FOR_EXCEPTION( true , RuntimeError , "CurveIntegralCalc::getCurveQuadPoints , Unknown Cell in 2D" );
		    }
		  }
		  // test for to many intersection points
		  TEST_FOR_EXCEPTION( total_points > 2 ,RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , no much intersection points: " << total_points);

		  // -> having the intersection points now we have to determine the:
		  // quad points and the gradients at the points(which norm will be in the integration)
		  // and the normalized normal vector (normalized since we need the sin and cos for Nitsche )
		  // -> as a very simple approach we consider the curve as a line between intersection points "startPoint -> endPoint"

		  // in the X direction the line should always have increasing values
		  if (startPoint[0] > endPoint[0]){
			  Point tmp = startPoint;
			  startPoint = endPoint;
			  endPoint = tmp;
		  }
		  SUNDANCE_MSG3(verb, tabs << "start end and points , start:" << startPoint << ", end:"<< endPoint)

		  // get the quadrature points for the line
		  Array<Point> quadPoints;
		  Array<double> quadWeights;
		  quad.getPoints( LineCell , quadPoints , quadWeights );
		  int nr_quad_points = quadPoints.size();


		  // The intersection points , we distribute the points along the line
		  curvePoints.resize(nr_quad_points);
		  SUNDANCE_MSG3(verb, tabs << " setting reference quadrature points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  SUNDANCE_MSG3(verb, tabs << " nr:" << ii << " Line quad point:" << quadPoints[ii] << ", line:" << (endPoint - startPoint));
			  curvePoints[ii] = startPoint + quadPoints[ii][0]*(endPoint - startPoint);
			  // we transform the intersection points to the reference cell
			  // todo: THIS MIGHT NOT WORK FOR UNSTRUCTURED CASE !!!
			  switch (maxCellType){
				  case QuadCell:{
					  curvePoints[ii][0] = -(cellPoints[0][0] - curvePoints[ii][0])/(cellPoints[3][0] - cellPoints[0][0]);
					  curvePoints[ii][1] = -(cellPoints[0][1] - curvePoints[ii][1])/(cellPoints[3][1] - cellPoints[0][1]);
				  } break;
				  case TriangleCell:{
					  curvePoints[ii][0] = -(cellPoints[0][0] - curvePoints[ii][0])/(cellPoints[1][0] - cellPoints[0][0]);
					  curvePoints[ii][1] = -(cellPoints[0][1] - curvePoints[ii][1])/(cellPoints[2][1] - cellPoints[0][1]);
				  } break;
				  default:{
					  // This should not happen
				  }
			  }
			  SUNDANCE_MSG3(verb, tabs << " quad point nr:" << ii << " = " << curvePoints[ii]);
		  }


		  // The derivatives at points, for this simple method are simple and constant over the whole quadrature
		  curveDerivs.resize(nr_quad_points);
		  Point dist_point(endPoint[0] - startPoint[0] , endPoint[1] - startPoint[1]);
		  SUNDANCE_MSG3(verb, tabs << " setting derivative values points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  curveDerivs[ii] = endPoint - startPoint;
			  SUNDANCE_MSG3(verb, tabs << " curve Derivatives point nr:" << ii << " = " << curveDerivs[ii]);
		  }


		  // calculate the norms to the curve
		  curveNormals.resize(nr_quad_points);
		  double sin_line = 0.0, cos_line = 0.0;
		  // in case we have a horizontal line
		  if ( fabs( endPoint[0] - startPoint[0])  < 1e-8 ){
			    Point tmp(startPoint[0],startPoint[1] - sqrt(dist_point*dist_point)*1e-2 );
		 	 	double inside_i =  paramCurve.curveEquation( tmp );
		 	 	if ( inside_i < 0.0)
		 	 	{
	  			    sin_line = -1.0;
				    cos_line = 0.0;
		 	 	}
			 	else
			 	{
		  		    sin_line = 1.0;
				    cos_line = 0.0;
			 	}
				SUNDANCE_MSG3(verb, tabs << "case 1.0, cos_line:" << cos_line << ", sin_line:" << sin_line);
		  } else {
		  // general case to determine the normal vector
			    Point tmp(startPoint[0] - sqrt(dist_point*dist_point)*1e-2 ,startPoint[1]);
		 	 	double inside_i =  paramCurve.curveEquation( tmp );
		 	 	double r_dev = sqrt ( dist_point * dist_point);
		 	 	// todo: check if this is OK such
		 	 	if ( (startPoint[1]-endPoint[1]) > 0.0 ){
		 	 	    // case 2.1 , positive slope
		 	 		if (inside_i < 0.0){
		 	 		  	cos_line = -fabs(endPoint[1]-startPoint[1]) / r_dev;
						sin_line = -fabs(endPoint[0]-startPoint[0]) / r_dev;
		 	 		} else {
		 	 		  	cos_line = fabs(endPoint[1]-startPoint[1]) / r_dev;
						sin_line = fabs(endPoint[0]-startPoint[0]) / r_dev;
		 	 		}
					SUNDANCE_MSG3(verb, tabs << "case 2.1, cos_line:" << cos_line << ", sin_line:" << sin_line << ", inside_i:" << inside_i);
		 	 	}else{
		 	 		// case 2.2 , negative slope
		 	 		if (inside_i < 0.0){
		 	 		  	cos_line = -fabs(endPoint[1]-startPoint[1]) / r_dev;
						sin_line = fabs(endPoint[0]-startPoint[0]) / r_dev;
		 	 		}else{
	 	 			  	cos_line = fabs(endPoint[1]-startPoint[1]) / r_dev;
						sin_line = -fabs(endPoint[0]-startPoint[0]) / r_dev;
		 	 		}
					SUNDANCE_MSG3(verb, tabs << "case 2.2, cos_line:" << cos_line << ", sin_line:" << sin_line << ", inside_i:" << inside_i);
		 		}
		  } // - end normal direction calculation

		  // set the normal direction in the corresponding way
		  SUNDANCE_MSG3(verb, tabs << " setting normal values at points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  Point tmp( 0.0 , 0.0);
			  curveNormals[ii] = tmp;
			  curveNormals[ii][0] = cos_line;
			  curveNormals[ii][1] = sin_line;
			  SUNDANCE_MSG3(verb, tabs << " curve Normals point nr:" << ii << " = " << curveNormals[ii]
                                  << ", cos_line:" << cos_line << ", sin_line:" << sin_line);
		  }

	  } break;

//========================== END 1D curve in 2D context ==========================

	  case 2: {
		  // 2D curve integral in 3D

		  // here we consider a simple surface, as it is in bilinear interpolation

		  TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints , surface integral not implemented yet ");

	  } break;

//========================== END 2D curve in 3D context ==========================

	  default: {
        // throw exception
		TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints , curve dimension must be 1 or two ");
	  }
	}

}
