/*
 * Init stage binary for the Binoth Les Houches Accord 
 *
 * Author: Michael Rauch <michael.rauch@kit.edu>
 * Initial version: Jun 2015
 * Last modified: Jun 2015
 *
 */

#include <iostream>
#include <cstring>
#include "VBFNLO/utilities/BLHAinterface.h"

int main(int argc, char* argv[]) {

  if ( argc < 2 ) {
    std::cerr << "Usage: " << argv[0] << " <order-file> [<contract-file>]" << std::endl;
    return 1;
  }

  int ierr=0;
  char* contractfile;
  if ( argc > 2 ) contractfile = argv[2];
  else {
    contractfile = new char[16];
    strcpy(contractfile,"OLE_contract.lh"); // SHERPA default
  }


  OLP_Order(argv[1],contractfile,&ierr);

  if ( ierr == -1 ) return 1;
  else if ( ierr == -2 ) {
    std::cout << "Unknown option: see contract file " << contractfile << " for details" << std::endl;
    return 1;
  }

  std::cout << "Successfully generated contract file " << contractfile << std::endl;
  return 0;

}

