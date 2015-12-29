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

        // Test if Q must be filled for the first iteration
//        forAll(mesh.cells(), i)
//        {
//           	Q[i] = 1.0e9;
//        }


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

