    const wordList fluidNames(rp["fluid"]);
    PtrList<fvMesh> fluidRegions(fluidNames.size());
    List<labelList> fluidRegionsLists(fluidNames.size());

    forAll(fluidNames, i)
    {
        Info<< "Create fluid mesh for region " << fluidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

	    fluidRegions.set
		(
		    i,
		    new fvMesh
		    (
			IOobject
			(
			    fluidNames[i],
			    runTime.timeName(),
			    runTime,
			    IOobject::MUST_READ
			    )
			)
		    );

	Info << "Create correspondence to the complete mesh for region "
	     << fluidNames[i] << " for time = " << runTime.timeName()
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
	      regSetPath(caseConst+"/polyMesh/sets/"+fluidNames[i]);
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

	  fluidRegionsLists[i] =  labelIOList(IOobject
	  		 (
	  		  fluidNames[i],
	  		  setsPath/pM/setS,
	  		  runTime,
	  		  IOobject::MUST_READ,
	  		  IOobject::NO_WRITE
	  		  )
	     );
	}
    }
