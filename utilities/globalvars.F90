! storage of global variables, should be used instead of globals.inc
! where possible

module globalvars

  logical lglobalprint ! switches stdout on/off, used for multiple processes in mpi
  logical ldoscales, ldoblha
  logical isggflo ! corresponds to ggflo in VBFNLO < 3.0

  integer seed ! input seed

end module

