SUBROUTINE open(len_ch,file_name,option,funit)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::len_ch
  CHARACTER(80),INTENT(IN)::file_name
  CHARACTER(1),INTENT(IN)::option
  INTEGER,INTENT(OUT)::funit

  LOGICAL:: UNITOP

  INTEGER::uno
  CHARACTER*1::word_uno(20)
  equivalence (word_uno, uno)
  CHARACTER*6::endianness

  INTEGER(KIND=4)::NTITLE, DCD_VARS(20),NATOM,HD(2)
  CHARACTER*4::DCD_TYPE
  CHARACTER*80,ALLOCATABLE,DIMENSION(:)::TITLE
  REAL(KIND=4)::delta_t
  equivalence(delta_t,DCD_VARS(10))
  LOGICAL:: with_cell
  integer(4), parameter :: dcd_magic_little = X'44524f43', dcd_magic_big = X'434f5244'

  ! Checking the unit

  UNITOP=.False.

  funit=500

  do while (UNITOP)
     inquire (unit=funit,opened=UNITOP)
     if (UNITOP) funit=funit+1
  end do

  if (len_ch>80) then
     PRINT*, '# Error: Name of file too long.'
  end if

  OPEN (unit=funit,name=TRIM(file_name),status='old',action='read',form='unformatted',access='stream')

  !!!!!!!!! Checking the native Endianness

  uno=1
  IF (ichar(word_uno(1)) == 0) THEN
     !"Big Endian"
     endianness='BIG   '
  ELSE
     !"Little Endian"
     endianness='LITTLE'
  END IF
  !!!!!!!!!

  !!!!!!!!! Checking the file Endianness

  READ(funit) HD(:)
  IF (HD(1)==84) THEN
     IF (HD(2) == dcd_magic_little) THEN
        ! Detected standard 32-bit DCD file with little-endian data.
        IF (endianness=='BIG   ') THEN
           CLOSE(funit)
           OPEN (unit=funit,name='run.dcd',status='old',action='read',form='unformatted',access='stream',convert='LITTLE_ENDIAN')
           READ(funit) HD(:)
        END IF
     END IF
     IF (HD(2) == dcd_magic_big) THEN
        ! Detected standard 32-bit DCD file with big-endian data.
        IF (endianness=='LITTLE') THEN
           CLOSE(funit)
           OPEN (unit=funit,name='run.dcd',status='old',action='read',form='unformatted',access='stream',convert='BIG_ENDIAN')
           READ(funit) HD(:)
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

  !!!!!!!!!! Reading header

  READ(funit,pos=5) DCD_TYPE,DCD_VARS
  READ(funit,pos=97) NTITLE
  ALLOCATE(TITLE(NTITLE))
  READ(funit) TITLE(:)
  READ(funit) HD,NATOM
  
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

  DEALLOCATE(TITLE)

END SUBROUTINE open


SUBROUTINE read (funit,natom,with_cell,pos_i,pos_o,cell,coors)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit,natom,pos_i
  LOGICAL,INTENT(IN):: with_cell
  INTEGER,INTENT(OUT)::pos_o

  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::cell
  DOUBLE PRECISION,DIMENSION(natom,3),INTENT(OUT)::coors

  REAL*8::buffer_cell(3,3)
  REAL(KIND=4),ALLOCATABLE,DIMENSION(:)::buffer
  INTEGER(KIND=4)::NTITLE, DCD_VARS(20),NATOM,HD(2)

  cell=0.0d0
  pos_o=pos_i

  IF (with_cell) THEN
     buffer_cell=0.0d0
     READ(funit,pos=pos_o) HD,buffer_cell(1,1), buffer_cell(1,2), buffer_cell(2,2), buffer_cell(1,3), buffer_cell(2,3), buffer_cell(3,3)
     cell=dble(buffer_cell)
     INQUIRE(funit,pos=pos_o)
  END IF

  ALLOCATE(buffer(natom))

  coors=0.0d0
  READ(funit,pos=pos_o) HD,buffer(:)
  coors(:,1)=dble(buffer(:))
  READ(funit) HD,buffer(:)
  coors(:,2)=dble(buffer(:))
  READ(funit) HD,buffer(:)
  coors(:,3)=dble(buffer(:))

  DEALLOCATE(buffer)

  INQUIRE(funit,pos=pos_o)
  
END SUBROUTINE read


SUBROUTINE close(funit)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit

  CLOSE(funit)

END SUBROUTINE close
