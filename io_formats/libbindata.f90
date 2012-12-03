SUBROUTINE fopen_read(funit,file_name)

  IMPLICIT NONE
  CHARACTER(80),INTENT(IN)::file_name
  INTEGER,INTENT(IN)::funit
  
  OPEN(unit=funit,FILE=TRIM(file_name),STATUS='old',action='READ',form='unformatted',access='stream')

END SUBROUTINE fopen_read

SUBROUTINE fclose(funit)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit

  CLOSE(funit)

END SUBROUTINE fclose

SUBROUTINE read_float_frame(funit,frame,num_parts,dimensions,coors)

  INTEGER,INTENT(IN)::funit,frame,num_parts,dimensions
  DOUBLE PRECISION,DIMENSION(num_parts,dimensions),INTENT(OUT)::coors

  INTEGER(KIND=8)::needle
  
  needle=frame*num_parts*dimensions*8+1
  READ(funit,pos=needle) coors(:,:)

END SUBROUTINE read_float_frame

SUBROUTINE read_float_coor(funit,frame,particle,dim,num_parts,dimensions,coor)

  INTEGER,INTENT(IN)::funit,frame,num_parts,dimensions
  DOUBLE PRECISION,INTENT(OUT)::coor

  INTEGER(KIND=8)::needle
  
  needle=frame*num_parts*dimensions*8+particle*dimensions*8+dim*8
  READ(funit,pos=needle) coor

END SUBROUTINE read_float_coor

SUBROUTINE read_int_frame(funit,frame,num_parts,dimensions,coors)

  INTEGER,INTENT(IN)::funit,frame,num_parts,dimensions
  INTEGER,DIMENSION(num_parts,dimensions),INTENT(OUT)::coors

  INTEGER(KIND=8)::needle
  
  needle=frame*bits_frame
  READ(funit,pos=needle) coors(:,:)

END SUBROUTINE read_int_frame

SUBROUTINE read_int_coor(funit,frame,particle,dim,num_parts,dimensions,coor)

  INTEGER,INTENT(IN)::funit,frame,num_parts,dimensions
  INTEGER,INTENT(OUT)::coor

  INTEGER(KIND=8)::needle
  
  needle=frame*num_parts*dimensions*8+particle*dimensions*8+dim*8
  READ(funit,pos=needle) coor

END SUBROUTINE read_int_coor


