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
    const wordList solidsNames(rp["solid"]);
    PtrList<fvMesh> solidRegions(solidsNames.size());
    List<labelList> solidRegionsLists(solidsNames.size());

    forAll(solidsNames, i)
    {
        Info<< "Create solid mesh for region " << solidsNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

	    solidRegions.set
		(
		    i,
		    new fvMesh
		    (
			IOobject
			(
			    solidsNames[i],
			    runTime.timeName(),
			    runTime,
			    IOobject::MUST_READ
			    )
			)
		    );

	Info << "Create correspondence to the complete mesh for region "
	     << solidsNames[i] << " for time = " << runTime.timeName()
	     << nl << endl;

	if(Pstream::master())
	{
	  std::string caseConst(runTime.caseConstant());

	  // If running in parallel, the runTime.caseConstant()
	  // adds a relative path to the sequential value.
	  // Since the directory where sets file is located is the same,
	  // the relative path is removed from the string
	  if(Pstream::parRun())
	     caseConst.erase(0,3);

	  std::string
	    regSetPath(caseConst+"/polyMesh/sets/"+solidsNames[i]);
	  std::string sedCommand("sed -i 's/cellSet/labelList/g' ");
	  std::string command(sedCommand+regSetPath);
	  
	  fileName pM("polyMesh");
	  fileName setS("sets");
	  fileName setsPath(runTime.caseConstant());

	  // Information about sets are stored in a cellSet class, which
	  // IOobject read as a hashTable. This is a problem to keep the
	  // information on data position in the coupling vector.
	  //
	  // Solution:
	  // call to sed using system, replacing cellSet in the set
	  // files of each region by labelList. Then OpenFOAM IOobject
	  // will be able to direct read data faster.
	  system(command.c_str());

	  solidRegionsLists[i] = labelIOList(IOobject
	  		 (
	  		  solidsNames[i],
	  		  setsPath/pM/setS,
	  		  runTime,
	  		  IOobject::MUST_READ,
	  		  IOobject::NO_WRITE
	  		  )
	     );
	}

    }


// Force calculation of geometric properties to prevent it being done
// later in e.g. some boundary evaluation
//(void)solidRegions[i].weights();
//(void)solidRegions[i].deltaCoeffs();

