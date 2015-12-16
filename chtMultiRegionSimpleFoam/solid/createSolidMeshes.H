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

