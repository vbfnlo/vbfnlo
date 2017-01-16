module hist_stor
   type :: hist_store
      ! for more elegant structure, should add here:
      !nbin, hlo, hup, htitles, hstored
      real*8, dimension(:), allocatable :: hsum
      real*8, dimension(:), allocatable :: hsquare
      real*8, dimension(:), allocatable :: hcount
   end type hist_store

   type :: hist_store2d
      real*8, dimension(:,:), allocatable :: hsum
      real*8, dimension(:,:), allocatable :: hsquare
      real*8, dimension(:,:), allocatable :: hcount
   end type hist_store2d
   
   type :: realsave
      ! temporary storage of different kinematics that get summed up first
      ! before they are used to calculate the MC error
      ! indices: kinematic (1: RE, 2..N: subtraction), hist_call (see hist_calls)
      real*8, dimension(:,:), allocatable :: save_y
      ! store bin number corresponding to number in save_y
      integer, dimension(:,:), allocatable :: save_bin

      ! counts the number of calls to FillHist (e.g. 2 calls for pt_jet in VBF LO)
      integer :: hist_calls 
      integer :: hist_calls_max
   end type realsave

   type :: realsave2d
      real*8, dimension(:,:,:), allocatable :: save_y
      integer, dimension(:,:), allocatable :: save_bin_x
      integer, dimension(:,:), allocatable :: save_bin_y
      integer :: hist_calls
      integer :: hist_calls_max
   end type realsave2d

   type(hist_store), dimension(:), allocatable, save :: hist
   type(hist_store2d), dimension(:), allocatable, save :: hist2d

   type (realsave), dimension(:), allocatable, save :: hsave
   type (realsave2d), dimension(:), allocatable, save :: hsave2d

   integer, save :: real_kinematics


contains

   subroutine inc (a, b)
      double precision, intent(inout) :: a
      double precision, intent(in) :: b
      a = a + b
   end subroutine inc

end module hist_stor



