fvVectorMatrix UDEqn
(
    fvm::ddt(U)
  - fvm::laplacian(nu, U)
 ==
    fvOptions(U)
);

UDEqn.relax();

fvOptions.constrain(UDEqn);

solve(UDEqn);

phi = fvc::interpolate(U) & mesh.Sf();
