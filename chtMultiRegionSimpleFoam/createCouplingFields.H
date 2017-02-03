// Include std C to get directory name

// Create Lists for data interchange with neutronics

// Total number of cells is the number of elements in the
// createSolidFields/solidRegionsLists structure and
// createFluidFields/fluidRegionsLists structure
unsigned int totalNumberOfCells = 0;
//unsigned int totalNumberOfRegions = solidsNames.size() + fluidNames.size();

// ---------------------------------------------------------------
// Notes on shared memory coupling:
//
// - Shared memory files are created by OpenFoam.
// - Files must be deleted by OpenFOAM.
// (changed from the last implementation)

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

// Posix structures initialization
calcOf = sem_open("calcOf", O_CREAT, 0666);
calcMil = sem_open("calcMil", O_CREAT, 0666);

if(Pstream::master())
{
// Semaphores and shared memory files are tested at this point.
// If milonga is not running, OpenFOAM runs uncoupled.

    // Names for files are derived from the current directory plus t and p
    // for temperatures and power. This allows running multiple coupled calculations
    // in the same system (as long as the cases directory names are different.
    char buffer[1000] = "";
    getcwd(buffer, 1000);

    std::string my(buffer);
    
    // Strip string
    std::string filename = my.substr(my.find_last_of("/\\")+1);

    // create file strings
    std::string temperatures = filename+"T";
    std::string potencias = filename+"Q";
	
    shmTfile = shm_open(temperatures.c_str(), O_RDWR, 0666);
    shmQfile = shm_open(potencias.c_str(), O_RDWR, 0666);
    
    if(shmTfile == -1 || shmQfile == -1)
    {
      Info << nl << " ---: Error getting shared memory: " << strerror(errno) << ". NON-COUPLED calculations." << nl << endl;
	coupling = false;
    }
    else
      Info << nl << " ---: milonga shared memory data detected! COUPLED calculations. " << nl << endl;
}

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

List<double> temperatureCompleteList(totalNumberOfCells, 300.0);
List<double> densityCompleteList(totalNumberOfCells, 1000.0);
List<double> powerCompleteList(totalNumberOfCells, 0.0);

if(Pstream::master())
{
// After reading the shared memory files, they must be mapped to data
    shmT = mmap(NULL, totalNumberOfCells*sizeof(double), PROT_WRITE, MAP_SHARED, shmTfile, 0);
    shmQ = mmap(NULL, totalNumberOfCells*sizeof(double), PROT_WRITE, MAP_SHARED, shmQfile, 0);

// Check if all files were properly mapped
    if((shmT == MAP_FAILED || shmQ == MAP_FAILED) && (coupling))
    {
	Info << " --- Error mapping shared memory: " << strerror(errno) << ". Exiting..." << endl;
	exit(errno);
    }

// Make a C++ cast
    shmTarray = reinterpret_cast<double*>(shmT);
    shmQarray = reinterpret_cast<double*>(shmQ);
}
