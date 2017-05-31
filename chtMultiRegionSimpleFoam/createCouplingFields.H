/*
  ---------------------------------------------------------------------------
  ---
  --- This file is original work by
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
  ---
  ---------------------------------------------------------------------------

  Application
  thesisCoupledFoam

  Description
  Steady-state version of chtMultiRegionFoam adapted to get the source-term for 
  the energy equation from an external source. The source term is defined as a 
  volScalarField and read from a C array structure. This version is implemented 
  and tested to run in parallel and sequentialy.

  ---------------------------------------------------------------------------

  This file is part of thesisCoupledFoam.

  thesisCoupledFoam is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  thesisCoupledFoam is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with thesisCoupledFoam.  If not, see <http://www.gnu.org/licenses/>.

*/

// - Includes std C function to get directory name
// - Create Lists for data interchange with neutronics
// - Total number of cells is the number of elements in the
// createSolidFields/solidRegionsLists structure and
// createFluidFields/fluidRegionsLists structure

unsigned int totalNumberOfCells = 0;

// ---------------------------------------------------------------
// Notes on shared memory coupling:
//
// - Shared memory files are created by milonga.
// - OpenFOAM seeks files. If they're not present or have problems
// reading them, runs UNCOUPLED.
// - Shared file names and semaphores are of the running directory
// case appended by T, Q, Of and Mil:
// -> T = temperature
// -> Q = power
// -> Of = OpenFOAM semaphore
// -> Mil = milonga semaphore
// ---------------------------------------------------------------

bool coupling = true;

// Three C standard arrays are created to shared data with Milonga
double *shmTarray = NULL;
double *shmQarray = NULL;

void *shmT;
void *shmQ;

int shmTfile = 0;
int shmQfile = 0;

// Posix C semaphores
sem_t *calcOf;
sem_t *calcMil;

// Directory base names allows running multiple coupled calculations
// in the same system (as long as the cases directory names are different).
char buffer[1000] = "";
getcwd(buffer, 1000);

std::string my(buffer);

// Strip string
std::string filename = my.substr(my.find_last_of("/\\")+1);

// create file strings
std::string temperatures = filename+"T";
std::string potencias = filename+"Q";
std::string semOf = filename+"Of";
std::string semMil = filename+"Mil";

// Posix structures initialization
//
// O_CREAT: create semaphore if it does not exists
// Maybe this is not the best flag. Check.
// 0666 is an optional argument mode: opens for reading
// and writing by every proccess.
//
// Check sem_open(3) man page.

calcOf = sem_open(semOf.c_str(), O_CREAT, 0666);
calcMil = sem_open(semMil.c_str(), O_CREAT, 0666);

// This is only done by the master proccess
if(Pstream::master())
{
    shmTfile = shm_open(temperatures.c_str(), O_RDWR, 0666);
    shmQfile = shm_open(potencias.c_str(), O_RDWR, 0666);
    
    // Semaphores and shared memory files are tested at this point.
    // If milonga is not running, OpenFOAM runs uncoupled.
	
    if(shmTfile == -1 || shmQfile == -1)
    {
      Info << nl << " ---: Error getting shared memory: " << strerror(errno) << ". NON-COUPLED calculations." << nl << endl;
	coupling = false;
    }
    else
      Info << nl << " ---: milonga shared memory data detected! COUPLED calculations. " << nl << endl;
}

// Get the total number of cells by summing sizes of every region
// of the mesh
if(Pstream::master())
{
  for (int i=0; i<solidsNames.size(); ++i)
    {
      totalNumberOfCells += solidRegionsLists[i].size();
    }

  for (int i=0; i<fluidNames.size(); ++i)
    {
      totalNumberOfCells += fluidRegionsLists[i].size();
    }
}

Info << " ---: totalNumberOfCells: " << totalNumberOfCells << endl;

// The last arguments are the initialization values of
// the lists. They are only *references* of physical consistent
// values. They're not used as is.
List<double> temperatureCompleteList(totalNumberOfCells, 300.0);
List<double> densityCompleteList(totalNumberOfCells, 1000.0);
List<double> powerCompleteList(totalNumberOfCells, 0.0);

if(Pstream::master())
{
    // After reading the shared memory files, they must be mapped to data.
    // Arguments:
    //            NULL:            The kernel choses the address to create the mapping
    //            totalNumber... : length bytes to be mapped
    //            PROT_WRITE:      Mapped pages can be written
    //            MAP_SHARED:      Share the mapping, making them visible to other proccesses
    //            shmXfile:        File descriptor of the mapping (already created)
    //            0:               Offset for mapping (not used).
    //
    // Further information: manpage of mmap(2).
    shmT = mmap(NULL, totalNumberOfCells*sizeof(double), PROT_WRITE, MAP_SHARED, shmTfile, 0);
    shmQ = mmap(NULL, totalNumberOfCells*sizeof(double), PROT_WRITE, MAP_SHARED, shmQfile, 0);

    // Check if all files were properly mapped
    // - coupling flag garantees not trying to run in coupled mode if mapping works but
    // files are not correctly created on the expected conditions
    if((shmT == MAP_FAILED || shmQ == MAP_FAILED) && (coupling))
    {
	Info << " --- Error mapping shared memory: " << strerror(errno) << ". Exiting..." << endl;
	exit(errno);
    }

    // Make a C++ cast of POSIX variables
    shmTarray = reinterpret_cast<double*>(shmT);
    shmQarray = reinterpret_cast<double*>(shmQ);
}
