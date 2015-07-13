
    fvMesh& mesh = solidRegions[i];
    solidThermo& thermo = thermos[i];
    const radiation::radiationModel& radiation = radiations[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> talpha = thermo.alpha();
    const volScalarField& alpha = talpha();
    tmp<volScalarField> tkappa = thermo.kappa();
    const volScalarField& kappa = tkappa();

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

    // Info << " --- ADDED: " << nl;
    // Info << " sF.size()= " << sF.size() << nl << " completeSF.size()= " << completeSF.size() << nl;
    // Info << " sF= " << sF << nl;
    // Info << " completeSF= " << completeSF << nl;
    
    // Values in qVol[i] are dimensioned values.
    // Pay attention on this before assigning it to the vector.
    // the code below only assing 0 to the vector
    // v.assign(qVol[i].size(), 0);
      
    //    procScalarFieldList[Pstream::myProcNo()] = sF;
    //    Info << " --- ADDED: procScalarFieldList.size()= " << procScalarFieldList.size() << nl;
    //Info << procScalarFieldList << endl;
    
    // This is only necessary if running in parallel
    if(Pstream::parRun())
      {
	fileName cellCorr(runTime.rootPath()+"/"+runTime.caseName()+
			  "/constant/fuel/polyMesh");

	// // Read cellProcAddressing for each processor
	// labelIOList procList
	//   (
	//    IOobject
	//    (
	//     "cellProcAddressing", 			   		// Name
	//     cellCorr,
	//     mesh,  			// Registry
	//     IOobject::MUST_READ, 	// Read option
	//     IOobject::NO_WRITE 		// Write Option
	//     ),
	//    procList
	//    );

	// procLabelList[Pstream::myProcNo()] = procList;

	//	Pout << " --- ADDED: procLabelList[" << Pstream::myProcNo()
	//     << "]: " << procList;
	
	//	Pstream::gatherList<scalarField>(procScalarFieldList);
	//	Pstream::scatterList<scalarField>(procScalarFieldList);

	//			Pout << " --- ADDED: procList: " << procList << endl;

      }

    // ------------------------------------------------------------------------
    //
    // As operacoes feitas aqui deverao ser feitas no armazenamento dos dados de
    // Q, sem arquivos. Por isso não foram totalmente apagadas (14/05/2015)
    //
    // ------------------------------------------------------------------------

    // First step: get data from processors and write it to a vector in the main processor
//    label n = Pstream::nProcs();

    // *************************** Corrigir
    scalar number = 29999+Pstream::myProcNo();
    scalarList localDataQ(qVol[i].internalField().size(), number);

    // If is the master process it will gather data from
    // other processes. If in sequential mode, it will not
    // enter the and jump over the gather iterations
    // if (Pstream::myProcNo() == 0)
    //   {
    // 	// Data structures
    // 	List<scalarList> allDataQ(n);
    // 	allDataQ[0] = localDataQ;

    // 	for(label j=1; j<n; j++)
    // 	  {
    // 	    // create the input stream from processor i
    // 	    //	    IPstream lStream(Pstream::blocking, j);
    // 	    //	    lStream >> allDataQ[j] >> procLabelList[j];

    // 	  }
    // 	//	Info << " --- ADDED: allData Q: " << allDataQ << endl;
    //   }
    // else
    //   {
    // 	// create the stream to send to the main proc
    // 	OPstream lStream
    // 	  (
    // 	   Pstream::blocking, 0
    // 	   );
    // 	//	lStream << localDataQ << procLabelList[Pstream::myProcNo()];
    //   }

    // In the solid 'fuel', Q must be read as input and the data from one
    // file is splitted among the processors.
    if(!Pstream::parRun())
      {
	// Since there is only one processor, the data is sequentially available
	// from procScalarFieldList[0]

	// ********************************** Erro aqui
    // 	forAll(completeSF, j)
    // 	  {
    // 	    //	    completeSF[j]=dataQ.at(j);
    // 	  }
    // 	//	Info << "* SEQ * --- ADDED: qVol[i].internalField() = completeSF: "
    // 	//	     << completeSF << nl;
    // 	qVol[i].internalField() = completeSF;
    //   }
    // else
    //   {
    // 	// If running in parallel, the Q data must be split among the
    // 	// processors respecting mesh topology. To enforce this restriction,
    // 	// cell positions are checked from procLabelList vector
    // 	forAll(sF, j)
    // 	  {
    // 	    //	    sF[j]=dataQ.at(procLabelList[Pstream::myProcNo()][j]);
    // 	  }
    // 	//	Info << "* PAR * --- ADDED: qVol[i].internalField() = sF: "
    // 	//	     << sF << nl;
    // 	qVol[i].internalField()= sF;
    //   }
      }
  }
// If the file is not present, use current qVol[i] (the volScalarField)
const volScalarField& Q = qVol[i];

//Info << " --- ADDED: Q after setRegionSolidFidls.H = " << Q << nl;

 
