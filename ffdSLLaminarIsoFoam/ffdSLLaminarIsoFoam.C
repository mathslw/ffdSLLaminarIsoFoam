/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "meshSearch.H"
#include "interpolation.H"
#include "singlePhaseTransportModel.H"
#include "volPointInterpolation.H"
#include "pisoControl.H"
#include "OFstream.H"
#include "fvOptions.H"

#define expo 2.0
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool pointInTriangle
(
    Foam::vector& pIT0,
    const Foam::Vector<double>& pIT1,
    const Foam::Vector<double>& pIT2,
    const Foam::Vector<double>& pIT3
)
{
    bool isInside = false;
    scalar dp1 = ((pIT1 - pIT2)^(pIT0 - pIT2))&((pIT1 - pIT2)^(pIT3 - pIT2));
    scalar dp2 = ((pIT2 - pIT3)^(pIT0 - pIT3))&((pIT2 - pIT3)^(pIT1 - pIT3));
    scalar dp3 = ((pIT3 - pIT1)^(pIT0 - pIT1))&((pIT3 - pIT1)^(pIT2 - pIT1));
    isInside = (dp1>=0)&&(dp2>=0)&&(dp3>=0);
    return isInside;
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    label patchInlet = mesh.boundaryMesh().findPatchID("inlet"); //Careful!!
    const surfaceVectorField& Cf  = mesh.Cf();
    const surfaceVectorField& Sf  = mesh.Sf();
    vector zero(0,0,0);
    volPointInterpolation interpolateVolPoint (mesh);
    
    while (runTime.loop())
    {
	    Info<< "Time = " << runTime.timeName() << nl << endl;
	    #include "CourantNo.H"
	    surfaceVectorField UUF = fvc::interpolate(U.oldTime());
	    pointVectorField sP = interpolateVolPoint.interpolate(U.oldTime());
	    dimensionedScalar dt("dt", dimensionSet(0,0,1,0,0,0,0), runTime.deltaT().value());
	    dictionary interpolationDict = mesh.solutionDict().subDict("interpolationSchemes");
	    autoPtr<interpolation<vector> > Uinterp = interpolation<vector>::New(interpolationDict, U.oldTime());
	    volVectorField posD(U.mesh().C());
	    volVectorField dis(dt*U.oldTime());
	    volScalarField dcg(mag(dis));
	    volVectorField posT(posD-dis);
	    volVectorField UU(U.oldTime());
	    forAll(U.internalField(),patchi)
    	    {
		vector projection = -dis.internalField()[patchi];
		bool findCell = false;
		label depCell = patchi;
		if (mag(projection) <= SMALL) //if there is no trace back!!
		{
			findCell = true;
			UU.ref()[patchi] = U.oldTime()[depCell];
		}

/*********************************************Arriving point*****/
		const labelList& facesA = mesh.cells()[depCell];
		bool inside = false;
		label priorInsideFace = 0;
		label i = 0;
		while ((!findCell)&&(i<facesA.size())&&(!inside))
		//while ((!findCell)&&(!inside))
		{
	    		label nFace = facesA[i];
			vector faceCentreTmp = Cf[nFace];
			vector normal = Sf[nFace];
  			normal /= mag(normal);
			//vector fp = faceCentreTmp - posD[depCell];//cell center to face center
			vector iPoint;
			scalar multiplierNumerator = (faceCentreTmp - posD[patchi]) & normal;
			scalar multiplierDenominator = projection & normal;
			if (mag(multiplierDenominator) > SMALL)
			{
				scalar multiplier = multiplierNumerator/multiplierDenominator;
				if(multiplier > 0)
				{
					iPoint = posD[patchi] + multiplier*projection; //locate the intersection point
					const labelList& points = mesh.faces()[nFace];
					scalar npoints = points.size();
					if (npoints == 3) //check if the ipoint is inside the face, tetra face
					{
						inside = pointInTriangle
						(
						iPoint, 
						mesh.points()[points[0]],
						mesh.points()[points[1]],
						mesh.points()[points[2]]
						);
					}
					else //check if the ipoint is inside the face, hex face
					{	
						label n = 0;
						while ((!inside) && (n<4))
						{
							label p1 = n % 4;
							label p2 = (n + 1) % 4;
							label p3 = (n + 2) % 4;
							inside = pointInTriangle
							(
							iPoint, 
							mesh.points()[points[p1]],
							mesh.points()[points[p2]],
							mesh.points()[points[p3]]
							);
							n++;
						}
					}
				}
			}
			if (inside) 
			{
				if (((iPoint-posT[patchi])&projection)>=0.0)
				{
					findCell = true;
					//#include "bfsiAlgorithm.H"
					UU.ref()[patchi] = Uinterp->interpolate(posT[patchi], depCell);
				}
				else
				{	
					if(!mesh.isInternalFace(nFace))
					{
						findCell = true;
						label patchII = mesh.boundaryMesh().whichPatch(nFace);
						if (patchII == patchInlet)
						{
							UU.ref()[patchi] = U.oldTime().boundaryField()[patchII]
							[
							mesh.boundaryMesh()[patchII].whichFace
							(
								nFace
							)
							];
						}
						else 
						{
							UU.ref()[patchi] = U.oldTime()[depCell];
						}
					}
					else
					{
						priorInsideFace = nFace;
						label ownCell = mesh.owner()[nFace];
						if (ownCell == depCell) depCell = mesh.neighbour()[nFace];
						else depCell = ownCell;
					}
				}
			}
			i++;
		}
/***************************************Arriving point finished*/

		label iii = 0;
		while ((!findCell)&&(iii<99))  //identify the depCell
		{
			const labelList& faces = mesh.cells()[depCell];
			inside = false;
			label ii = 0;
			while ((!inside) && (ii<faces.size()))
			{
		    		label nFace = faces[ii];
				vector faceCentreTmp = Cf[nFace];
				vector normal = Sf[nFace];
          			normal /= mag(normal);
				//vector fp = faceCentreTmp - posD[depCell];//cell center to face center
				vector iPoint;
				if (nFace != priorInsideFace) //choose the possible faces
				{
					scalar multiplierNumerator = (faceCentreTmp - posD[patchi]) & normal;
					scalar multiplierDenominator = projection & normal;
					if (mag(multiplierDenominator) > SMALL)
					{
						scalar multiplier = multiplierNumerator/multiplierDenominator;
						iPoint = posD[patchi] + multiplier*projection; //locate the intersection point
						const labelList& points = mesh.faces()[nFace];
						scalar npoints = points.size();
						if (npoints == 3) //check if the ipoint is inside the face, tetra face
						{
							inside = pointInTriangle
							(
							iPoint, 
							mesh.points()[points[0]],
							mesh.points()[points[1]],
							mesh.points()[points[2]]
							);
						}
						else //check if the ipoint is inside the face, hex face
						{	
							label n = 0;
							while ((!inside) && (n<4))
							{
							    label p1 = n % 4;
							    label p2 = (n + 1) % 4;
							    label p3 = (n + 2) % 4;
							    inside = pointInTriangle
							    (
								iPoint, 
								mesh.points()[points[p1]],
								mesh.points()[points[p2]],
								mesh.points()[points[p3]]
							    );
							    n++;
							}
						}
	 				}
				}
				if (inside) 
				{
					if (((iPoint-posT[patchi])&projection)>=0.0)
					{
						findCell = true;
						//#include "bfsiAlgorithm.H"
						UU.ref()[patchi] = Uinterp->interpolate(posT[patchi], depCell);
					}
					else
					{	
						if(!mesh.isInternalFace(nFace))
						{
							findCell = true;
							label patchII = mesh.boundaryMesh().whichPatch(nFace);
							if (patchII == patchInlet)
							{
								UU.ref()[patchi] = U.oldTime().boundaryField()[patchII]
								[
								mesh.boundaryMesh()[patchII].whichFace
								(
									nFace
								)
								];
							}
							else 
							{
								UU.ref()[patchi] = U.oldTime()[depCell];
							}
						}
						else
						{
							priorInsideFace = nFace;
							label ownCell = mesh.owner()[nFace];
							if (ownCell == depCell) 
							{depCell = mesh.neighbour()[nFace];}
							else 
							{depCell = ownCell;}
						}
					}
				}
				ii++;
			}
			iii++;
			if(!inside) 
			{
				iii = 100;
				findCell = true;
			}
		}

		if (iii == 100)
		{
			label deIndex = patchi;
			label deIndexNew = patchi;
			scalar dc(dcg[patchi]);
			scalar dcNew(dcg[patchi]);
			for (int ij=1; ij<100; ij++)
			{
			 	labelList adjacent = mesh.cellCells()[deIndex];
				forAll(adjacent, j)
				{
				    scalar dd(mag(posT[patchi]-posD[adjacent[j]]));
				    if (dd<dc)
				    {
					dcNew = dd;
					deIndexNew = adjacent[j];
				    }
				}
				if (deIndexNew == deIndex) ij=100;
				else
				{
				    deIndex = deIndexNew;
				    dc = dcNew;
				}
			}
			UU.ref()[patchi] = Uinterp->interpolate(posT[patchi], deIndexNew);
		//	TT.internalField()[patchi] = Tinterp->interpolate(posT[patchi], deIndexNew);
		}
	    }

//cout<<" tracing back finished "<<endl;
	U.oldTime() = UU;
//	T.oldTime() = TT;

	#include "UEqn.H"

        while (piso.correct())
        {
	    #include "pEqn.H"
        }

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
