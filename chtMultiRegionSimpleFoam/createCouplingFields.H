// Create Lists for data interchange with neutronics

// Total number of cells is the number of elements in the
// createSolidFields/solidRegionsLists structure and
// createFluidFields/fluidRegionsLists structure
unsigned int totalNumberOfCells = 0;
//unsigned int totalNumberOfRegions = solidsNames.size() + fluidNames.size();

// ---------------------------------------------------------------
// Notes on shared memory coupling:
//
// - Shared memory files are created by Milonga.
// - Files are also deleted by Milonga.

// Three C standard arrays are created to shared data with Milonga
double *shmTarray;
double *shmDarray;
double *shmQarray;

void *shmT;
void *shmD;
void *shmQ;

int shmTfile;
int shmDfile;
int shmQfile;

// Posix C semaphores
sem_t *semsent;
sem_t *semreceived;

// Posix structures initialization
semsent = sem_open("semsent", O_CREAT, 0666);
semreceived = sem_open("semreceived", O_CREAT, 0666);

// Semaphores and shared memory files are tested at this point.
// If milonga is not running, OpenFOAM exits.
shmTfile = shm_open("temperaturas", O_RDWR, 0666);
shmDfile = shm_open("densidades", O_RDWR, 0666);
shmQfile = shm_open("potencias", O_RDWR, 0666);

if(shmTfile == -1 || shmDfile == -1 || shmQfile == -1)
{
    Info << " --- Error getting shared memory: " << strerror(errno) << ". Exiting..." << endl;
    exit(errno);
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

// After reading the shared memory files, they must be mapped to data
shmT = mmap(NULL, totalNumberOfCells*sizeof(double), PROT_WRITE, MAP_SHARED, shmTfile, 0);
shmD = mmap(NULL, totalNumberOfCells*sizeof(double), PROT_WRITE, MAP_SHARED, shmDfile, 0);
shmQ = mmap(NULL, totalNumberOfCells*sizeof(double), PROT_WRITE, MAP_SHARED, shmQfile, 0);

// Check if all files were properly mapped
if(shmT == MAP_FAILED || shmD == MAP_FAILED || shmQ == MAP_FAILED)
{
    Info << " --- Error mapping shared memory: " << strerror(errno) << ". Exiting..." << endl;
    exit(errno);
}

// Make a C++ cast
shmTarray = reinterpret_cast<double*>(shmT);
shmDarray = reinterpret_cast<double*>(shmD);
shmQarray = reinterpret_cast<double*>(shmQ);


// for(unsigned int i=0; i<totalNumberOfCells; i++)
//     powerCompleteList[i] = 0.9 + static_cast <float> (rand())
// 			      /( static_cast <float> (RAND_MAX/(1.1-0.9)));
