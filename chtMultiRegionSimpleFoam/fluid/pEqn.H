{
    // Lembrete: rho é um volScalarField criado em createFields
    rho = thermo.rho();
    rho = max(rho, rhoMin[i]);
    rho = min(rho, rhoMax[i]);
    rho.relax();

    
    volScalarField rAU(1.0/UEqn().A());

    // Interpolate rAu to faces giving it "Dp" name
    surfaceScalarField rhorAUf("Dp", fvc::interpolate(rho*rAU));

    // H by A
    // Remembering that U is a volVectorField defined in setRegionFluidFields.H
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    // Libera o tmp velocidade criado em UEqn.H
    UEqn.clear();

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    // phiHbyA é um a densidade interpolada nas faces vezes
    // o produto interno do vetor HbyA interpolado e
    // os vetores de area das faces na malha
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*(fvc::interpolate(HbyA) & mesh.Sf())
    );

    // OpenFOAM 2.2.2
    // fvOptions.relativeFlux(fvc::interpolate(rho), phiHbyA);

    // OpenFOAM 2.4.0
    fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

    bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);

    phiHbyA += phig;

    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);

    // Solve pressure
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference
        (
            pRefCell,
            compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
        );

        p_rghEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    p = p_rgh + rho*gh;

    #include "continuityErrs.H"

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (closedVolume && compressible)
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /compressibility;
        p_rgh = p - rho*gh;
    }

    rho = thermo.rho();
    rho = max(rho, rhoMin[i]);
    rho = min(rho, rhoMax[i]);
    rho.relax();

    Info<< "Min/max rho:" << min(rho).value() << ' '
        << max(rho).value() << endl;
}
