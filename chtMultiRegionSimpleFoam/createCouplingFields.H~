// Create Lists for data interchange with neutronics

// Total number of cells is the number of elements in the
// createSolidFields/solidRegionsLists structure and
// createFluidFields/fluidRegionsLists structure
unsigned int totalNumberOfCells = 0;
//unsigned int totalNumberOfRegions = solidsNames.size() + fluidNames.size();

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

List<float> temperatureCompleteList(totalNumberOfCells, 300.0);
List<float> densityCompleteList(totalNumberOfCells, 1000.0);
List<float> powerCompleteList(totalNumberOfCells, 0.0);

for(unsigned int i=0; i<totalNumberOfCells; i++)
    powerCompleteList[i] = 0.9 + static_cast <float> (rand())
			      /( static_cast <float> (RAND_MAX/(1.1-0.9)));
