program read_dcd

IMPLICIT NONE

INTEGER::uno
CHARACTER*1::word_uno(20)
equivalence (word_uno, uno)
CHARACTER*6::endianness

INTEGER::i,j,ii,k,kk,needle

INTEGER(KIND=4)::NTITLE, DCD_VARS(20),NATOM,HD(2),TEST
character*6::ch
CHARACTER*4::DCD_TYPE
CHARACTER*80,ALLOCATABLE,DIMENSION(:)::TITLE
REAL(KIND=4)::delta_t
equivalence(delta_t,DCD_VARS(10))
REAL*8::cell(3,3)
REAL(KIND=4),ALLOCATABLE,DIMENSION(:)::buffer
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::coors
LOGICAL:: with_cell
integer(4), parameter :: dcd_magic_little = X'44524f43', dcd_magic_big = X'434f5244'

OPEN (unit=66,name='run.dcd',status='old',action='read',form='unformatted',access='stream')

!!!!!!!!! Checking Endianness
uno=1
IF (ichar(word_uno(1)) == 0) THEN
   !"Big Endian"
   endianness='BIG   '
ELSE
   !"Little Endian"
   endianness='LITTLE'
END IF
!!!!!!!!!

READ(66) HD(:)
IF (HD(1)==84) THEN
   IF (HD(2) == dcd_magic_little) THEN
      ! Detected standard 32-bit DCD file with little-endian data.
      IF (endianness=='BIG   ') THEN
         CLOSE(66)
         OPEN (unit=66,name='run.dcd',status='old',action='read',form='unformatted',access='stream',convert='LITTLE_ENDIAN')
         READ(66) HD(:)
      END IF
   END IF
   IF (HD(2) == dcd_magic_big) THEN
      ! Detected standard 32-bit DCD file with big-endian data.
      IF (endianness=='LITTLE') THEN
         CLOSE(66)
         OPEN (unit=66,name='run.dcd',status='old',action='read',form='unformatted',access='stream',convert='BIG_ENDIAN')
         READ(66) HD(:)
      END IF
   END IF
ELSE
   IF ((HD(1)+HD(2))==84) THEN
      ! Detected CHARMM -i8 64-bit DCD file not supported.
      print*,'# Detected CHARMM -i8 64-bit DCD file. File not supported.'
      stop
   ELSE
      print*,'# Detected Unknown DCD format. File not supported.'
      stop
   END IF
END IF

READ(66,pos=5) DCD_TYPE,DCD_VARS
READ(66,pos=97) NTITLE
ALLOCATE(TITLE(NTITLE))
READ(66) TITLE(:)
INQUIRE(66,pos=needle)
READ(66,pos=needle) HD,NATOM

IF (DCD_TYPE/='CORD') THEN
   print*,'# DCD type ',DCD_TYPE,' not supported.'
   stop
END IF

! DCD_TYPE: CORD or VELD
! DCD_VARS(1):  Number of frames in this file
! DCD_VARS(2):  Number of previous integration steps
! DCD_VARS(3):  Frequency (integration steps) to save this file
! DCD_VARS(4):  Number of integration steps in the run to create this file
! DCD_VARS(5):  Frequency of coordinate saving
! DCD_VARS(6):  
! DCD_VARS(7):
! DCD_VARS(8):  Number of degrees of freedom during the run
! DCD_VARS(9):  Number of fixed atoms
! DCD_VARS(10): Timestep in AKMA-units. Bit-copy from the 32-bit real number
! DCD_VARS(11): 1 if crystal lattice information is present in the frames
! DCD_VARS(12): 1 if this is a 4D trajectory
! DCD_VARS(13): 1 if fluctuating charges are present
! DCD_VARS(14): 1 if trajectory is the result of merge without consistency checks
! DCD_VARS(15):
! DCD_VARS(16):
! DCD_VARS(17):
! DCD_VARS(18):
! DCD_VARS(19):
! DCD_VARS(20): CHARMM version number


with_cell=.false.
IF (DCD_VARS(11)) with_cell=.true.


ALLOCATE(buffer(NATOM),coors(3,NATOM))

DO ii=1,100

   IF (with_cell) THEN
      cell=0.0d0
      READ(66) HD,cell(1,1), cell(1,2), cell(2,2), cell(1,3), cell(2,3), cell(3,3)
   END IF


   coors=0.0d0
   READ(66) HD,buffer(:)
   coors(1,:)=dble(buffer(:))
   READ(66) HD,buffer(:)
   coors(2,:)=buffer(:)
   READ(66) HD,buffer(:)
   coors(3,:)=buffer(:)

END DO

DEALLOCATE(buffer,coors,title)




end program read_dcd
