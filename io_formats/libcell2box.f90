!! This function is not used by gro, xtc and trr
SUBROUTINE TRICLINIC (cell,box,volume,ortho)

  IMPLICIT NONE
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::cell
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::box
  DOUBLE PRECISION,INTENT(OUT)::volume
  INTEGER,INTENT(OUT)::ortho
  DOUBLE PRECISION::alpha,beta,gamma,x,y,z

  ortho=0
  box=0.0d0
  alpha=cell(1,2)
  beta=cell(1,3)
  gamma=cell(2,3)
  x=cell(1,1)
  y=cell(2,2)
  z=cell(3,3)

  box(1,1)=x
  IF ((alpha==90.0d0).and.(beta==90.0d0).and.(gamma==90.0d0)) THEN
     box(2,2)=y
     box(3,3)=z
     volume=x*y*z
     ortho=1
  ELSE
     box(2,1)=y*cosd(gamma)
     box(2,2)=y*sind(gamma)  ! sqrt(y**2-box(2,2)**2)
     box(3,1)=z*cosd(beta)
     box(3,2)=z*(cosd(alpha)-cosd(beta)*cosd(gamma))/sind(gamma) 
     box(3,3)=sqrt(z*z-box(3,1)**2-box(3,2)**2)
     volume=box(1,1)*box(2,2)*box(3,3)
  END IF

END SUBROUTINE TRICLINIC
