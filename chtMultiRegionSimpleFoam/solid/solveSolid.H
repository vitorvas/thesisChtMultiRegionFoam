{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix hEqn
        (
           - fvm::laplacian(betav*alpha, h, "laplacian(alpha,h)")
           - Q
//           + fvOptions(rho, h)
//          Retorna um tmp<fvMatrix<type> > depois de fazer operações
//          na matriz do campo e na matriz criada internamente.
//          implementado em fvOptionListTemplates.C
//          Ainda não sei o que o addSup faz nem como o this->operator[](i) funciona
        );

        hEqn.relax();

//        fvOptions.constrain(hEqn);
//        faz um source.setValue(na equacao e no campo). O campo vem de eqn.psi().name()
        hEqn.solve();

//        fvOptions.correct(h);
//        chama o correct do campo, que no caso de h aparentemente está implementado
//        em temperatureLimitsConstraint

// ------------ Resta saber se estão todas reimplementadas neste caso, para este solver
    }
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()) << ' ' << max(thermo.T()) << endl;

