volScalarField rUA("rUA", dt*UDEqn.A()/UDEqn.A());
surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

surfaceScalarField phiU(phi);

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
	fvm::laplacian(rUAf, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    if (piso.finalNonOrthogonalIter())
    {
	phi -= pEqn.flux();
	//p += dp;
    }
}

#include "continuityErrs.H"

U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
U.correctBoundaryConditions();
