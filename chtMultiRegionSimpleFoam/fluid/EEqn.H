{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      - fvm::laplacian(turb.alphaEff(), he)
     ==
        rad.Sh(thermo)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    solverPerformance sp = EEqn.solve();
//    Info << " --- ADDED: sp.finalResidual = " << sp.finalResidual() << endl;
    
    fvOptions.correct(he);

    thermo.correct();
    rad.correct();

    Info<< "(Energy) Min/max T:" << min(thermo.T()).value() << ' '
	<< max(thermo.T()).value() << endl;
}
