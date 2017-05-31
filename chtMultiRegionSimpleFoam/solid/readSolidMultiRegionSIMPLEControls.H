/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
  chtMultiRegionSimpleFoam

  Description
  Steady-state version of chtMultiRegionFoam

  ---------------------------------------------------------------------------
  ---
  --- This file is modified version of the original OpenFOAM code by
  --- Vitor Vasconcelos Araújo Silva Ph.D thesis work as student of 
  --- Universidade Federal de Minas Gerais - UFMG, Brazil and worker
  --- of Centro de Desenvolvimento da Tecnologia Nuclear, CDTN/CNEN, Brazil
  ---
  --- Copyright 2016 Vitor Silva
  ---
  --- contact: vitors@cdtn.br
  ---
  --- UFMG, CDTN/CNEN and hereby disclaims all copyright in the present
  --- program written by Vitor Vasconcelos Araújo Silva.
  ---
  --- Vitor Vasconcelos Araújo Silva, 19 december 2016

  Application
  thesisCoupledFoam

  Description
  Steady-state version of chtMultiRegionFoam adapted to get the source-term for 
  the energy equation from an external source. The source term is defined as a 
  volScalarField and read from a C array structure. This version is implemented 
  and tested to run in parallel and sequentialy.

  IMPORTANT:
  ----------
  OpenFOAM and Milonga communication is based in POSIX semaphores.
  Milonga should be started before OpenFOAM. Otherwise OpenFOAM will 
  hold soon after initialization waiting for the first semaphore.
  (A call to sem_wait(calcOf) in createCouplingFields.H)

  ---------------------------------------------------------------------------

  \*---------------------------------------------------------------------------*/
    const dictionary& simple = solidRegions[i].solutionDict().subDict("SIMPLE");

    const int nNonOrthCorr =
        simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

