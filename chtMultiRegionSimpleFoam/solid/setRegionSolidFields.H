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
  --- This file is a modified version of the original OpenFOAM code by
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

    solidThermo& thermo = thermos[i];

    tmp<volScalarField> trho = thermo.rho();
    tmp<volScalarField> tcp = thermo.Cp();
    tmp<volScalarField> talpha = thermo.alpha();

    const volScalarField& alpha = talpha();

    tmp<volScalarField> tkappa = thermo.kappa();

    volScalarField& h = thermo.he();

    const volScalarField& betav = betavSolid[i];

    // OpenFoam's original IOoptionsList introduced in version 2.2.2 is removed
    // fv::IOoptionList& fvOptions = solidHeatSources[i];

    // Added to use my own source-term
    // When reading from another file, it should be done at
    // this point, before atributing the qVol to const Q field

    // 13/04/2015
    // ------------------------------------------------------------------
    // CHANGE from former implementation.
    // Q (former qvol) is only read as a default file in OpenFOAM format.
    // if it is not present, Q is assumed 0 in the first step.
    // ------------------------------------------------------------------

    // First of all, only do that if the solid is the fuel
if (solidRegions[i].name() == "fuel")
  {
    // Create a scalar field from the values read from Q file in
    // solids directories
    scalarField sF = scalarField(qVol[i].internalField());
    scalarField completeSF(qVol[i].internalField().size()*Pstream::nProcs(), 0.0);
    
    // Values in qVol[i] are dimensioned values.
    // Pay attention on this before assigning it to the vector.
    // the code below only assing 0 to the vector
    // v.assign(qVol[i].size(), 0);
      
    // This is only necessary if running in parallel
    if(Pstream::parRun())
      {
	fileName cellCorr(runTime.rootPath()+"/"+runTime.caseName()+
			  "/constant/fuel/polyMesh");
      }

    // First step: get data from processors and write it to a vector in the main processor
    // label n = Pstream::nProcs();
    scalarList localDataQ(qVol[i].internalField().size(), 0.0);
    
  }
// If the file is not present, use current qVol[i] (the volScalarField)
const volScalarField& Q = qVol[i];



 
