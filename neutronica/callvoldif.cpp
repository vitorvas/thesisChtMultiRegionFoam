#include <iostream>
//#include "/home/EQN_SOFT/petsc-3.4.4/include/mpiuni/mpi.h"

extern "C" void voldif_(float*);

int main(int argc, char* argv[])
{
  float qu[10];

  for(int ik=0; ik<10; ik++)
    std::cout << " qu[" << ik << "]= " << qu[ik] << std::endl;
  std::cout << "End of C++ program" << std::endl;

  voldif_(qu);

  return 0;
}
