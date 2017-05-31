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
    // Initialise solid field pointer lists
    PtrList<solidThermo> thermos(solidRegions.size());
    PtrList<radiation::radiationModel> radiations(solidRegions.size());
    PtrList<volScalarField> betavSolid(solidRegions.size());
    PtrList<volScalarField> qVol(solidRegions.size());

    List<List<labelList> > solidList(solidRegions.size(), List<labelList>(Pstream::nProcs()));
	
    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermos\n" << endl;
        thermos.set(i, solidThermo::New(solidRegions[i]));

        Info<< "    Adding to radiations\n" << endl;
        radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

	// Read Q if it exists. If not, create a null (zero)
	// volScalarField
	IOobject Qfile
	  (
	   "Q",
	   runTime.timeName(),
	   solidRegions[i],
	   IOobject::READ_IF_PRESENT,
	   IOobject::AUTO_WRITE
	   );

        // At this point, Q is zero if the file does not exist
        // or is filled with values from file "Q"

	// Must check it before creating the field
	if(Qfile.headerOk())
	  {
	    // Create a qvol field from dictionary
	    qVol.set
	      (
	       i,
	       new volScalarField (Qfile, solidRegions[i])
	       );

	  }
	else
	  {
	    // If file is not there, create a 'dummy' IOobject
	    // setting the dimensions of the field
	    qVol.set
	      (
	       i,
	       new volScalarField
	       (
		IOobject
		(
		 "Q",
		 runTime.path(),
		 solidRegions[i],
		 IOobject::NO_READ,
		 IOobject::AUTO_WRITE
		 ),
		solidRegions[i],
		dimensionedScalar("2", dimensionSet(1, -1, -3, 0, 0),
				  scalar(0.0))
		)
	       );
	  }

        IOobject betavSolidIO
        (
            "betavSolid",
            runTime.timeName(),
            solidRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (betavSolidIO.headerOk())
        {
            betavSolid.set
            (
                i,
                new volScalarField(betavSolidIO, solidRegions[i])
            );
        }
        else
        {
            betavSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betavSolid",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("1", dimless, scalar(1.0))
                )
            );
        }

	// 23/04/2015
	// This is only necessary if running in parallel
	//
	// If running in parallel, the decomposed meshes must have
	// a function to map parallel meshes cells to the region original
	// mesh. The cellProcAddressing file
	// (processor<n>/constant/<region>/polyMesh/cellProcAddressing)
	// This labelIOList will be used to map rho and T values in one processor
	// to the complete region mesh.
	//
	// Note: for solid regions, values are read and written.

	if(Pstream::parRun())
	{
	  fileName cellCorr(runTime.rootPath()+"/"+runTime.caseName()+
			    "/constant/"+solidRegions[i].name()+"/polyMesh");

	  // Read cellProcAddressing for each processor
	  solidList[i][Pstream::myProcNo()] = labelIOList(IOobject
	     (
	      "cellProcAddressing", // Filename
	      cellCorr,
	      solidRegions[i],	// Registry
	      IOobject::MUST_READ, 	// Read option
	      IOobject::NO_WRITE 	// Write Option
	      ));

	  // Create a mapping from each processor to the solidListsRegions
	}
    }

