    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turb.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve
    (
        UEqn()
     ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
// IMPORTANTE:
//
// o termo ghf*fvc::snGrad(rho) surge na equação devido ao algoritmo
// SIMPLE, no ajuste da pressão. p_rgh já leva uma compensação do rho*g.
// Quando o gradiente é calculado para p_rgh surge um termo espúrio
// (não físico) que é compensado pela subtração do termo em questão.	    
        )
    );

    fvOptions.correct(U);
