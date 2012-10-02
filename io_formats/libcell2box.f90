!! This function is not used by gro, xtc and trr
SUBROUTINE CELL2BOX (cell,box,volume,ortho)

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

END SUBROUTINE CELL2BOX

!! This function is used by gro, xtc and trr
SUBROUTINE BOX2CELL (box,cell,volume,ortho)

  IMPLICIT NONE
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::cell
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
  DOUBLE PRECISION,INTENT(OUT)::volume
  INTEGER,INTENT(OUT)::ortho
  DOUBLE PRECISION::alpha,beta,gamma,x,y,z

  ortho=0
  cell=0.0d0
  volume=0.0d0

  IF ((box(2,1)==0.0d0).and.(box(3,1)==0.0d0).and.(box(3,2)==0.0d0)) THEN
     volume=box(1,1)*box(2,2)*box(3,3)
     cell(1,1)=box(1,1)
     cell(2,2)=box(2,2)
     cell(3,3)=box(3,3)
     cell(1,2)=90.0d0
     cell(1,3)=90.0d0
     cell(2,3)=90.0d0
     ortho=1
  ELSE
     volume=box(1,1)*box(2,2)*box(3,3)
     x=box(1,1)
     y=sqrt(box(2,1)**2+box(2,2)**2)
     z=sqrt(box(3,1)**2+box(3,2)**2+box(3,3)**2)
     cell(1,1)=x
     cell(2,2)=y
     cell(3,3)=z
     cell(2,3)=acosd(dot_product(box(2,:),box(1,:))/(x*y)) ! gamma
     cell(1,3)=acosd(dot_product(box(3,:),box(1,:))/(x*z)) ! beta
     cell(1,2)=acosd(dot_product(box(3,:),box(2,:))/(x*y)) ! alpha
  END IF

END SUBROUTINE BOX2CELL

