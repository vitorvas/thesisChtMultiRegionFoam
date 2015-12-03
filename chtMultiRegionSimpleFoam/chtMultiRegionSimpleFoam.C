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
  --- This file was modified as part of 
  --- Vitor Vasconcelos Araújo Silva Ph.D thesis work.
  ---
  contact: vitors@cdtn.br

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
  (A call to sem_wait(semsent) in createCouplingFields.H)

  \*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulenceModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvIOoptionList.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// C++ stdlib
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstring>         // Used by memcpy when coupling

// C Posix (Shared memory communication)
#include <semaphore.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

using namespace std;

// To be able to call the Fortran function heatflux_
extern "C" void heatflux_(float*);

int main(int argc, char *argv[])
{
    
#include "setRootCase.H"

    Info << nl << nl
	 << " ---> thesisCoupledFoam. Versao no Brasil." << nl
	 << " ---> Vitor Vasconcelos - vitors@cdtn.br" << nl
	 << " ---> Built at " << __TIME__ << ", " << __DATE__ << nl << endl;

#include "createTime.H"

    // Number of CFD iterations before neutronics call
    // hardcoded:
    unsigned int nIterations = 0;
    unsigned int cFact = 100;
   
    regionProperties rp(runTime);
  
#include "createFluidMeshes.H"
#include "createSolidMeshes.H"

#include "createFluidFields.H"
#include "createSolidFields.H"

#include "initContinuityErrs.H"

#include "createCouplingFields.H"

    // Get the number of processes
    //  label nProcs = Pstream::nProcs();

    while (runTime.loop())
    {
	nIterations++;

	Info<< "Time = " << runTime.timeName() << nl << endl;

	forAll(fluidRegions, i)
	{
	    Info << "\nSolving for fluid region "
		 << fluidRegions[i].name() << endl;
	  
#include "setRegionFluidFields.H"
#include "readFluidMultiRegionSIMPLEControls.H"
#include "solveFluid.H"

	}
      
	forAll(solidRegions, i)
	{
	    Info<< "\nSolving for solid region "
		<< solidRegions[i].name() << endl;

#include "setRegionSolidFields.H"
#include "readSolidMultiRegionSIMPLEControls.H"
#include "solveSolid.H"
    
	    runTime.write();

	    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< nl << endl;
	}

	// Time to couple ----------------------------------------------
	//
	// How neutronics is called?
	//
	// - check if the defined number of iterations (defined in nIterCFD)
	//   ran using the std::fmod() function. If yes:
	// - copy all data from all processors to completeLists[]
	// - call neutronics
	// - normalize all qVol data. This means that the data returned
	//   from neutronics is not copied over qVol, but rather the data
	//   multiplies qVol ir order to have a normalized value.
	//
	// After the regions loops, neutronics is called.
	// The new Q is non-normalized. The Q data must then be multiplied
	// by the new Q obtained from neutronics.
	// - bool values updated to make
      
	// Info << " ------ %: " << nIterations%cFact
	// 	   << " time().value()= " << runTime.time().value()
	// 	   << " cFact= " << cFact
	// 	   << endl;

	if(!(nIterations%cFact))
	{
	    forAll(fluidRegions, i)
	    {
		// Structures to hold data and cell addressing
		// for all processors
		List<scalarList> dataT(Pstream::nProcs());
		List<scalarList> dataRho(Pstream::nProcs());

		scalarList localDataT(thermoFluid[i].T().size());
		scalarList localDataRho(thermoFluid[i].T().size());

		// Since T, rho and Q have the same size
		// the processing of all data is performed
		// inside the same loop
		for(int k=0; k<Pstream::nProcs(); k++)
		{
		    dataT[k].setSize(thermoFluid[i].T().size(), 0.0);
		    dataRho[k].setSize(thermoFluid[i].rho().size(), 0.0);
		}
	
		// The if structure below gathers data from all processors
		if(Pstream::master())
		{
		    // For master processor data is copied directly
		    dataT[0] = thermoFluid[i].T();
		    dataRho[0] = thermoFluid[i].rho();
		      
		    // Gather data from all processors to master
		    for(label j=1; j<Pstream::nProcs(); j++)
		    {
			IPstream inputMasterStream(Pstream::blocking, j);
			inputMasterStream >> dataT[j] >> dataRho[j];
		    }

		    // Copy from dataT and dataRho to the completeList for FLUID Regions
		    for(int k=0; k<Pstream::nProcs(); k++)
		    {
			for(int m=0; m<dataT[k].size(); m++)
			{
			    // ATTENTION:
			    // Trick mapping among the complete vector, using the fluid Regions data
			    // which maps to the data coming from processors
			    temperatureCompleteList[fluidRegionsLists[i][(k*dataT[k].size())+m]] = dataT[k][m];
			    densityCompleteList[fluidRegionsLists[i][(k*dataRho[k].size())+m]] = dataRho[k][m];
			}
		    }
		}

		else
		{
		    localDataT = thermoFluid[i].T().internalField();
		    localDataRho = thermoFluid[i].rho().internalField();
		      
		    // 0 references de master process
		    OPstream outputSlavesStream(Pstream::blocking, 0);
		    outputSlavesStream << thermoFluid[i].T().internalField()
				       << thermoFluid[i].rho().internalField()
				       << fluidList[Pstream::myProcNo()];
		}
	    } // End forAll fluid loop
	  
	    forAll(solidRegions, i)
	    {
		// Structures to hold data and cell addressing
		// for all processors
		List<scalarList> dataT(Pstream::nProcs());
		List<scalarList> dataRho(Pstream::nProcs());
		List<scalarList> dataQ(Pstream::nProcs());
	      
		scalarList localDataT(thermos[i].T().size());
		scalarList localDataRho(thermos[i].rho().size());
		scalarList localDataQ(qVol[i].internalField().size());
	      
		// Since T, rho and Q have the same size
		// the processing of all data is performed
		// inside the same loop
		for(int k=0; k<Pstream::nProcs(); k++)
		{
		    dataT[k].setSize(thermos[i].T().size(), 0.0);
		    dataRho[k].setSize(thermos[i].rho().size(), 0.0);
		    dataQ[k].setSize(qVol[i].size(), 0.0);
		}
	      
		// The if structure below gathers data from all processors
		if(Pstream::master())
		{
		    // For master data is copied directly
		    dataT[0] = thermos[i].T();
		    dataRho[0] = thermos[i].rho();
		    dataQ[0] = qVol[i].internalField();
		  
		    // Gather data from all processors to master
		    for(label j=1; j<Pstream::nProcs(); j++)
		    {
			IPstream inputMasterStream(Pstream::blocking, j);
			inputMasterStream >> dataT[j] >> dataRho[j] >> dataQ[j];
		    }
		  
		    // Copy from dataT and dataRho to the completeList for SOLID Regions
		    for(int k=0; k<Pstream::nProcs(); k++)
		    {

			// This loop invariant is the size of dataT, which is the same of dataRho and dataQ.
			// Use of any of these lists yields the same result
			for(int m=0; m<dataT[k].size(); m++)
			{
			    // ATTENTION:
			    // Trick mapping among the complete vector, using the solid Regions data
			    // which maps to the data coming from processors
			    temperatureCompleteList[solidRegionsLists[i][(k*dataT[k].size())+m]] = dataT[k][m];
			    densityCompleteList[solidRegionsLists[i][(k*dataRho[k].size())+m]] = dataRho[k][m];
			    powerCompleteList[solidRegionsLists[i][(k*dataQ[k].size())+m]] = dataQ[k][m];
			}
		    }
		}
		else
		{
		    localDataT = thermos[i].T().internalField();
		    localDataRho = thermos[i].rho().internalField();
		    localDataQ = qVol[i].internalField();
		  
		    // 0 references de master process
		    OPstream outputSlavesStream(Pstream::blocking, 0);
		    outputSlavesStream << localDataT << localDataRho << localDataQ << solidList[Pstream::myProcNo()];
		}
	      
		// -----------------------------------------------------------------------------------------
		// -----------------------------------------------------------------------------------------
		// -----------------------------------------------------------------------------------------

		// Calling neutronics
		if(solidRegions[i].name() == "fuel")
		{
		    Info << nl << "Calling NEUTRONICS..." << nl << endl;
		    for(int o=0; o<solidRegionsLists[i].size(); o++)
			powerCompleteList[solidRegionsLists[i][o]] =  0.9 + static_cast <float> (rand())
			    /( static_cast <float> (RAND_MAX/(1.1-0.9)));
		}
	      
	    } // End forAll solid loop

	  
	    // -----------------------------------------------------------------------------------------
	    // -----------------------------------------------------------------------------------------
	    // -----------------------------------------------------------------------------------------

	    forAll(solidRegions, i)
	    {
		// Initialize dataQ field with the right size
		List<scalarList> dataQ(Pstream::nProcs());

		// Used to scatter data to each processor
		scalarList localDataQ(qVol[i].internalField().size(), 0.0);
	      
		if(solidRegions[i].name() == "fuel")
		{
		    for(int k=0; k<Pstream::nProcs(); k++)
		    {
			dataQ[k].setSize(qVol[i].internalField().size(), 1.0);
		    }

		    // The if structure below scatters data to all processors
		    if(Pstream::master())
		    {
			for(int k=0; k<Pstream::nProcs(); k++)
			{
			    // This loop control is the size of dataT, which is the same of dataRho and dataQ.
			    // Use of any of these lists yields the same result
			    for(int m=0; m<dataQ[k].size(); m++)
			    {
				// ATTENTION:
				// Tricky mapping among the complete vector, using the solid Regions data
				// which maps to the data coming from processors
				dataQ[k][m] *= powerCompleteList[solidRegionsLists[i][(k*dataQ[k].size())+m]];
			    }
			}

			// Scatter data from all processors to master
			for(label j=1; j<Pstream::nProcs(); j++)
			{
			    OPstream outputMasterStream(Pstream::blocking, j);
			    outputMasterStream << dataQ[j];
			}

			// For the main process or in a sequential run
			// update qVol data from neutronics
			qVol[i].internalField() *= dataQ[0];
		    }
		    else
		    {
			localDataQ = dataQ[Pstream::myProcNo()];
		      
			// 0 references de master process
			IPstream inputSlavesStream(Pstream::blocking, 0);
			inputSlavesStream >> localDataQ;

			// Copy data to local scalarField qVol[fuel]
			qVol[i].internalField() *= localDataQ;
		      
		    }
		    Pout << " --- " << qVol[i].internalField() << nl << endl;
		    solidRegions[i].write();		  
		}
	    }
	}
      
    } // runTime.loop()

    return 0;
}




// ************************************************************************* //
