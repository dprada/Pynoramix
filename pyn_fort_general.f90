MODULE GLOB

  INTEGER,DIMENSION(:),ALLOCATABLE::cl_ind,cl_start  !! indices fortran
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::cl_val


  INTEGER,DIMENSION(:,:),ALLOCATABLE::cell_ns
  INTEGER,DIMENSION(:),ALLOCATABLE::gns_head,gns_list,gns_at_cell
  LOGICAL,DIMENSION(:),ALLOCATABLE::cell_upd
  INTEGER::mtot_c,mx_c,my_c,mz_c,my_mz_c

  !### Verlet neighbours lists should be this way, but f2py does not support derived types yet. 
  !TYPE I_NO_RECT_MAT
  !   INTEGER,DIMENSION(:),ALLOCATABLE::column
  !   INTEGER::dim
  !END TYPE I_NO_RECT_MAT
  !TYPE(I_NO_RECT_MAT),DIMENSION(:),ALLOCATABLE::ver_ic_ind,ver_oc_ind  !! indices fortran

  INTEGER,DIMENSION(:,:),ALLOCATABLE::ver_ic_ind,ver_oc_ind  !! indices fortran
  INTEGER,DIMENSION(:),ALLOCATABLE::ver_ic_dim,ver_oc_dim

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::pos_ant

CONTAINS


SUBROUTINE PBC(vector,box,ortho)
 
  IMPLICIT NONE
 
  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::vector
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
  INTEGER,INTENT(IN)::ortho
  INTEGER::i
  DOUBLE PRECISION::x,L,Lhalf
 
  IF (ortho==1) THEN
     DO i=1,3
        L=box(i,i)
        Lhalf=0.50d0*L
        x=vector(i)
        IF (abs(x)>Lhalf) THEN
           IF (x>Lhalf) THEN
              x=x-L
           ELSE
              x=x+L
           END IF
           vector(i)=x
        END IF
     END DO
  ELSE
 
     print*, 'Not implemented'
 
  END IF
  
END SUBROUTINE PBC


INTEGER FUNCTION CELL_INDEX(ix,iy,iz)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::ix,iy,iz
  cell_index=1+MOD(ix-1+mx_c,mx_c)+MOD(iy-1+my_c,my_c)*mx_c+MOD(iz-1+mz_c,mz_c)*mx_my_c

END FUNCTION CELL_INDEX


SUBROUTINE MAKE_CELL_NS (r,box,natom)

  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,INTENT(IN)::r
  double precision,DIMENSION(3,3),INTENT(IN)::box

  INTEGER::ii,jj,kk,gg,ix,iy,iz

  INTEGER::icell

  mx_c=FLOOR(box(1,1)/r)
  my_c=FLOOR(box(2,2)/r)
  mz_c=FLOOR(box(3,3)/r)
  mx_my_c=mx_c*my_c
  mtot_c=mx_my_c*mz_c

  print*,mx_c,my_c,mz_c
  print*,box(1,1)/mx_c,box(2,2)/my_c,box(3,3)/mz_c

  print*,((natom*1.0d0)/(mtot_c*1.0d0))


  IF (ALLOCATED(cell_ns)) DEALLOCATE(cell_ns)
  IF (ALLOCATED(cell_upd)) DEALLOCATE(cell_upd)

  ALLOCATE(cell_ns(mtot_c,27))
  ALLOCATE(cell_upd)

  DO ix=1,mx_c
     DO iy=1,my_c
        DO iz=1,mz_c
           icell=cell_index(ix,iy,iz)
           gg=0
           DO ii=-1,1
              DO jj=-1,1
                 DO kk=-1,1
                    gg=gg+1
                    cell_ns(icell,gg)=cell_index(ix+ii,iy+jj,iz+kk)
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
  
  IF (ALLOCATED(gns_head))    DEALLOCATE(gns_head)
  IF (ALLOCATED(gns_list))    DEALLOCATE(gns_list)
  IF (ALLOCATED(gns_at_cell)) DEALLOCATE(gns_at_cell)
  ALLOCATE(gns_head(mtot_c))
  ALLOCATE(gns_list(natom))
  ALLOCATE(gns_at_cell(natom))

END SUBROUTINE MAKE_CELL_NS


SUBROUTINE GRID_NS_LIST (coors,box,natom)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,DIMENSION(natom,3),INTENT(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box

  DOUBLE PRECISION::mx_c_L,my_c_L,mz_c_L
  INTEGER::ii,icell


  mx_c_L=(mx_c*1.0d0)/box(1,1)
  my_c_L=(my_c*1.0d0)/box(2,2)
  mz_c_L=(mz_c*1.0d0)/box(3,3)

  gns_head=0
  DO ii=1,natom
     icell=1+int(coors(ii,1)*mx_c_L)+int(coors(ii,2)*my_c_L)*mx_c+int(coors(ii,3)*mz_c_L)*mx_my_c
     gns_at_cell(ii)=icell
     gns_list(ii)=gns_head(icell)
     gns_head(icell)=ii
  END DO


END SUBROUTINE GRID_NS_LIST

SUBROUTINE MAKE_VERLET_LIST_GRID_NS (r_ic,r_oc,pbc_opt,coors,box,vol,ortho,natom)

  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN)::r_ic,r_oc
  INTEGER,INTENT(IN)::pbc_opt,ortho
  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,INTENT(IN)::vol
  DOUBLE PRECISION,DIMENSION(natom,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box

  DOUBLE PRECISION::r_ic2,r_oc2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux,pi,fact

  !!! Hasta aqui 7.2 seg 1000 frames (3 min 25000)
  INTEGER,DIMENSION(:),ALLOCATABLE::auxlist_ic,auxlist_oc
  INTEGER::ii,jj,gg,ll,kk,icell,ncell
  INTEGER::dim_ic,gg_ic
  INTEGER::dim_oc,gg_oc

  r_ic2=r_ic*r_ic
  r_oc2=r_oc*r_oc

  pi=acos(-1.0d0)
  fact=(natom/vol)*(4.0d0*pi/3.0d0)*1.50d0   !! factor 1.50 to be safe, since f2py does not support derived types
  dim_ic=INT(fact*(r_ic2*r_ic))+1
  dim_oc=INT(fact*(r_oc2*r_oc))+1

  print*,dim_ic,dim_oc

  IF (ALLOCATED(ver_ic_ind))   DEALLOCATE(ver_ic_ind)
  IF (ALLOCATED(ver_ic_dim))   DEALLOCATE(ver_ic_dim)
  IF (ALLOCATED(ver_oc_ind))   DEALLOCATE(ver_oc_ind)
  IF (ALLOCATED(ver_oc_dim))   DEALLOCATE(ver_oc_dim)
  IF (ALLOCATED(pos_ant))      DEALLOCATE(pos_ant)

  ALLOCATE(ver_ic_ind(natom,dim_ic),ver_oc_ind(natom,dim_oc))
  ALLOCATE(ver_ic_dim(natom),ver_oc_dim(natom))
  ALLOCATE(pos_ant(natom,3))
  pos_ant=coors

  ALLOCATE(auxlist_ic(dim_ic),auxlist_oc(dim_oc))

  DO ii=1,natom
     vect_aux=coors(ii,:)
     gg_ic=0
     gg_oc=0
     icell=gns_at_cell(ii)
     DO jj=1,27
        ncell=cell_ns(icell,jj)
        kk=gns_head(ncell)
        DO
           IF (kk==0) EXIT
           IF (kk/=ii) THEN
              vect=(coors(kk,:)-vect_aux)
              CALL PBC (vect,box,ortho)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=r_oc2) THEN
                 gg_oc=gg_oc+1
                 IF (gg_oc>dim_oc) THEN
                    PRINT*,"# ERROR: variable dim_oc should be greater than", dim_oc
                    STOP
                 END IF
                 auxlist_oc(gg_oc)=kk
                 IF (val_aux<=r_ic2) THEN
                    gg_ic=gg_ic+1
                    IF (gg_ic>dim_ic) THEN
                       PRINT*,"# ERROR: variable dim_ic should be greater than", dim_ic
                       STOP
                    END IF
                    auxlist_ic(gg_ic)=kk
                 END IF
              END IF
           END IF
           kk=gns_list(kk)
        END DO
     END DO
     
     ver_oc_ind(ii,1:gg_oc)=auxlist_oc(1:gg_oc)
     ver_oc_dim(ii)=gg_oc
     ver_ic_ind(ii,1:gg_ic)=auxlist_ic(1:gg_ic)
     ver_ic_dim(ii)=gg_ic

  END DO

  DEALLOCATE(auxlist_ic,auxlist_oc)

END SUBROUTINE MAKE_VERLET_LIST_GRID_NS


SUBROUTINE UPDATE_VERLET_LIST_GRID_NS (r_ic,r_oc,pbc_opt,coors,box,vol,ortho,natom)

  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN)::r_ic,r_oc
  INTEGER,INTENT(IN)::pbc_opt,ortho
  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,INTENT(IN)::vol
  DOUBLE PRECISION,DIMENSION(natom,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box

  INTEGER::ii,jj,kk,icell
  DOUBLE PRECISION::delta,val_aux
  DOUBLE PRECISION::r_diff_2,r_diff_22,r_ic2,r_oc2
  DOUBLE PRECISION,DIMENSION(3)::vect

  drneimax=0.0 
  drneimax2=0.0
  r_oc2=r_oc*r_oc
  r_ic2=r_ic*r_ic
  r_diff_2=(r_oc-r_ic)/2.0d0
  r_diff_22=r_diff_2*r_diff_2

  CALL GRID_NS_LIST (coors,box,natom)

  cell_upd=.FALSE.
  DO icell=1,mtot_c
     delta=0.0 
     kk=gns_head(icell)
     DO
        IF (kk==0) EXIT
        vect=(coors(kk,:)-pos_ant(kk,:))
        CALL PBC (vect,box,ortho)  !! Check if this can be avoid (maybe not)
        val_aux=dot_product(vect,vect)
        IF (val_aux > r_diff_22) THEN  !! It is not the standard way of doing it: first_max+second_max>r_oc-r_ic
           DO jj=1,27
              cell_upd(cell_ns(i,jj))=.TRUE.
           END DO
           EXIT
        END IF
        kk=gns_list(kk)
     END DO
  END DO

  DO icell=1,mtot_c
     IF (cell_upd(icell)==.FALSE.) THEN
        kk=gns_head(icell)
        DO
           IF (kk==0) EXIT
           

     ELSE
        kk=gns_head(icell)

        pos_ant()




SUBROUTINE MAKE_CONTACT_LIST (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2

  INTEGER,DIMENSION(:),ALLOCATABLE::ilist1,ilist2
  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,ai,aj
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val
  
  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num

  IF (ALLOCATED(cl_ind)) DEALLOCATE(cl_ind)
  IF (ALLOCATED(cl_val)) DEALLOCATE(cl_val)
  IF (ALLOCATED(cl_start)) DEALLOCATE(cl_start)

  ALLOCATE(ilist1(n1),ilist2(n2))
  ilist1(:)=list1(:)+1
  ilist2(:)=list2(:)+1

  cut_off2=cut_off*cut_off

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0


  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       IF (gg==0) THEN
                          ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                       ELSE
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(:)=pcl_ind(ii)%i1(:)
                          box_val(:)=pcl_val(ii)%d1(:)
                          DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                          ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                          pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                          pcl_val(ii)%d1(1:gg)=box_val(:)
                          DEALLOCATE(box_ind,box_val)
                       END IF
                       gg=gg+1
                       pcl_ind(ii)%i1(gg)=aj
                       pcl_val(ii)%d1(gg)=val_aux
                    END IF
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       IF (gg==0) THEN
                          ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                       ELSE
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(:)=pcl_ind(ii)%i1(:)
                          box_val(:)=pcl_val(ii)%d1(:)
                          DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                          ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                          pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                          pcl_val(ii)%d1(1:gg)=box_val(:)
                          DEALLOCATE(box_ind,box_val)
                       END IF
                       gg=gg+1
                       pcl_ind(ii)%i1(gg)=aj
                       pcl_val(ii)%d1(gg)=val_aux
                    END IF
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        END IF
     ELSE
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=pcl_num(ii)
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    ll=pcl_num(jj)
                    IF (gg==0) THEN
                       ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(:)=pcl_ind(ii)%i1(:)
                       box_val(:)=pcl_val(ii)%d1(:)
                       DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                       ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                       pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                       pcl_val(ii)%d1(1:gg)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    IF (ll==0) THEN
                       ALLOCATE(pcl_ind(jj)%i1(1),pcl_val(jj)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(ll),box_val(ll))
                       box_ind(:)=pcl_ind(jj)%i1(:)
                       box_val(:)=pcl_val(jj)%d1(:)
                       DEALLOCATE(pcl_ind(jj)%i1,pcl_val(jj)%d1)
                       ALLOCATE(pcl_ind(jj)%i1(ll+1),pcl_val(jj)%d1(ll+1))
                       pcl_ind(jj)%i1(1:ll)=box_ind(:) 
                       pcl_val(jj)%d1(1:ll)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    gg=gg+1
                    ll=ll+1
                    pcl_ind(ii)%i1(gg)=aj
                    pcl_val(ii)%d1(gg)=val_aux
                    pcl_ind(jj)%i1(ll)=ai
                    pcl_val(jj)%d1(ll)=val_aux
                    pcl_num(jj)=ll
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=pcl_num(ii)
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    ll=pcl_num(jj)
                    IF (gg==0) THEN
                       ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(:)=pcl_ind(ii)%i1(:)
                       box_val(:)=pcl_val(ii)%d1(:)
                       DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                       ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                       pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                       pcl_val(ii)%d1(1:gg)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    IF (ll==0) THEN
                       ALLOCATE(pcl_ind(jj)%i1(1),pcl_val(jj)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(ll),box_val(ll))
                       box_ind(:)=pcl_ind(jj)%i1(:)
                       box_val(:)=pcl_val(jj)%d1(:)
                       DEALLOCATE(pcl_ind(jj)%i1,pcl_val(jj)%d1)
                       ALLOCATE(pcl_ind(jj)%i1(ll+1),pcl_val(jj)%d1(ll+1))
                       pcl_ind(jj)%i1(1:ll)=box_ind(:) 
                       pcl_val(jj)%d1(1:ll)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    gg=gg+1
                    ll=ll+1
                    pcl_ind(ii)%i1(gg)=aj
                    pcl_val(ii)%d1(gg)=val_aux
                    pcl_ind(jj)%i1(ll)=ai
                    pcl_val(jj)%d1(ll)=val_aux
                    pcl_num(jj)=ll
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        END IF
     END IF
  END IF

  gg=SUM(pcl_num(:),DIM=1)
  ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))

  gg=0
  DO ii=1,natom1
     cl_start(ii)=gg
     DO jj=1,pcl_num(ii)
        gg=gg+1
        cl_ind(gg)=pcl_ind(ii)%i1(jj)
        cl_val(gg)=pcl_val(ii)%d1(jj)
     END DO
     DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
  END DO
  cl_start(natom1+1)=gg

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)
  DEALLOCATE(ilist1,ilist2)

  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF

END SUBROUTINE MAKE_CONTACT_LIST


SUBROUTINE MAKE_CONTACT_LIST2 (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2

  INTEGER,DIMENSION(:),ALLOCATABLE::ilist1,ilist2
  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,kk,ai,aj,lim
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind,box_ind2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val,box_val2
  
  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num,pcl_numi

  IF (ALLOCATED(cl_ind)) DEALLOCATE(cl_ind)
  IF (ALLOCATED(cl_val)) DEALLOCATE(cl_val)
  IF (ALLOCATED(cl_start)) DEALLOCATE(cl_start)

  ALLOCATE(ilist1(n1),ilist2(n2))
  ilist1(:)=list1(:)+1
  ilist2(:)=list2(:)+1

  cut_off2=cut_off*cut_off

  lim=30
  ALLOCATE(box_ind(lim),box_val(lim))

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0

  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       gg=gg+1
                       IF (gg>lim) THEN
                          ALLOCATE(box_ind2(lim),box_val2(lim))
                          box_ind2=box_ind
                          box_val2=box_val
                          DEALLOCATE(box_ind,box_val)
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(1:lim)=box_ind2
                          box_val(1:lim)=box_val2
                          DEALLOCATE(box_ind2,box_val2)
                          lim=gg
                       END IF
                       box_ind(gg)=aj
                       box_val(gg)=val_aux
                    END IF
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       gg=gg+1
                       IF (gg>lim) THEN
                          ALLOCATE(box_ind2(lim),box_val2(lim))
                          box_ind2=box_ind
                          box_val2=box_val
                          DEALLOCATE(box_ind,box_val)
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(1:lim)=box_ind2
                          box_val(1:lim)=box_val2
                          DEALLOCATE(box_ind2,box_val2)
                          lim=gg
                       END IF
                       box_ind(gg)=aj
                       box_val(gg)=val_aux
                    END IF
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        END IF
     ELSE
        ALLOCATE(pcl_numi(natom1))
        pcl_numi=0
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    gg=gg+1
                    IF (gg>lim) THEN
                       ALLOCATE(box_ind2(lim),box_val2(lim))
                       box_ind2=box_ind
                       box_val2=box_val
                       DEALLOCATE(box_ind,box_val)
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(1:lim)=box_ind2
                       box_val(1:lim)=box_val2
                       DEALLOCATE(box_ind2,box_val2)
                       lim=gg
                    END IF
                    box_ind(gg)=aj
                    box_val(gg)=val_aux
                    pcl_numi(aj)=pcl_numi(aj)+1
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    gg=gg+1
                    IF (gg>lim) THEN
                       ALLOCATE(box_ind2(lim),box_val2(lim))
                       box_ind2=box_ind
                       box_val2=box_val
                       DEALLOCATE(box_ind,box_val)
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(1:lim)=box_ind2
                       box_val(1:lim)=box_val2
                       DEALLOCATE(box_ind2,box_val2)
                       lim=gg
                    END IF
                    box_ind(gg)=aj
                    box_val(gg)=val_aux
                    pcl_numi(aj)=pcl_numi(aj)+1
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        END IF
     END IF
  END IF

  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        gg=SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
        END DO
        cl_start(natom1+1)=gg
        DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
     ELSE
        gg=2*SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           gg=gg+pcl_numi(ii)
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
        END DO
        cl_start(natom1+1)=gg
        pcl_numi=0
        DO ii=1,natom1
           DO jj=1,pcl_num(ii)
              kk=pcl_ind(ii)%i1(jj)
              ll=pcl_numi(kk)+1
              pcl_numi(kk)=ll
              gg=cl_start(kk)+ll
              cl_ind(gg)=ii
              cl_val(gg)=pcl_val(ii)%d1(jj)
           END DO
           DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
        END DO
        DEALLOCATE(pcl_numi)
     END IF
  END IF

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)
  DEALLOCATE(ilist1,ilist2)

  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF


END SUBROUTINE MAKE_CONTACT_LIST2


SUBROUTINE UPDATE_CONTACT_LIST2 (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2


  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,mm,kk,ai,aj,dim_caja,amigo,amigo2,lim
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind,box_ind2,caja,aux_caja
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val,box_val2
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind,pcl_ind2   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num,pcl_numi

  cut_off2=cut_off*cut_off
  dim_caja=30
  lim=30


  ALLOCATE(box_ind(lim),box_val(lim))
  ALLOCATE(filtro(natom2),caja(dim_caja))
  filtro=.FALSE.

  cl_ind=cl_ind+1
  DEALLOCATE(cl_val)

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0


  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=0
           mm=0

           filtro(ii)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 mm=mm+1
                 IF (mm>=dim_caja) THEN
                    ALLOCATE(aux_caja(dim_caja))
                    aux_caja=caja
                    DEALLOCATE(caja)
                    ALLOCATE(caja(mm))
                    caja(1:dim_caja)=aux_caja(:)
                    DEALLOCATE(aux_caja)
                    dim_caja=mm
                 END IF
                 caja(mm)=amigo
                 filtro(amigo)=.TRUE.
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo2
                    filtro(amigo2)=.TRUE.
                 END IF
              END DO
           END DO

           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 gg=gg+1
                 IF (gg>lim) THEN
                    ALLOCATE(box_ind2(lim),box_val2(lim))
                    box_ind2=box_ind
                    box_val2=box_val
                    DEALLOCATE(box_ind,box_val)
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(1:lim)=box_ind2
                    box_val(1:lim)=box_val2
                    DEALLOCATE(box_ind2,box_val2)
                    lim=gg
                 END IF
                 box_ind(gg)=aj
                 box_val(gg)=val_aux
              END IF
           END DO
           ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
           pcl_ind(ii)%i1(:)=box_ind(1:gg)
           pcl_val(ii)%d1(:)=box_val(1:gg)
           pcl_num(ii)=gg
        END DO

     ELSE

        ALLOCATE(pcl_numi(natom1))
        pcl_numi=0

        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=0
           mm=0

           filtro(ai)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 IF (amigo>ai) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo
                    filtro(amigo)=.TRUE.
                 END IF
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    IF (amigo2>ai) THEN
                       mm=mm+1
                       IF (mm>=dim_caja) THEN
                          ALLOCATE(aux_caja(dim_caja))
                          aux_caja=caja
                          DEALLOCATE(caja)
                          ALLOCATE(caja(mm))
                          caja(1:dim_caja)=aux_caja(:)
                          DEALLOCATE(aux_caja)
                          dim_caja=mm
                       END IF
                       caja(mm)=amigo2
                       filtro(amigo2)=.TRUE.
                    END IF
                 END IF
              END DO
           END DO
           
           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 gg=gg+1
                 IF (gg>lim) THEN
                    ALLOCATE(box_ind2(lim),box_val2(lim))
                    box_ind2=box_ind
                    box_val2=box_val
                    DEALLOCATE(box_ind,box_val)
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(1:lim)=box_ind2
                    box_val(1:lim)=box_val2
                    DEALLOCATE(box_ind2,box_val2)
                    lim=gg
                 END IF
                 box_ind(gg)=aj
                 box_val(gg)=val_aux
                 pcl_numi(aj)=pcl_numi(aj)+1
              END IF
           END DO
           ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
           pcl_ind(ii)%i1(:)=box_ind(1:gg)
           pcl_val(ii)%d1(:)=box_val(1:gg)
           pcl_num(ii)=gg
        END DO

     END IF
  END IF

  DEALLOCATE(filtro,caja,box_ind,box_val)
  DEALLOCATE(cl_ind,cl_start)
  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        gg=SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
        END DO
        cl_start(natom1+1)=gg
        DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
     ELSE
        gg=2*SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           gg=gg+pcl_numi(ii)
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
        END DO
        cl_start(natom1+1)=gg
        pcl_numi=0
        DO ii=1,natom1
           DO jj=1,pcl_num(ii)
              kk=pcl_ind(ii)%i1(jj)
              ll=pcl_numi(kk)+1
              pcl_numi(kk)=ll
              gg=cl_start(kk)+ll
              cl_ind(gg)=ii
              cl_val(gg)=pcl_val(ii)%d1(jj)
           END DO
           DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
        END DO
        DEALLOCATE(pcl_numi)
     END IF
  END IF

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)


  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF


END SUBROUTINE UPDATE_CONTACT_LIST2


SUBROUTINE UPDATE_CONTACT_LIST (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2


  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,mm,kk,ai,aj,dim_caja,amigo,amigo2
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind,caja,aux_caja
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind,pcl_ind2   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num,pcl_num2

  cut_off2=cut_off*cut_off
  dim_caja=10

  ALLOCATE(filtro(natom2),caja(dim_caja))
  filtro=.FALSE.

  cl_ind=cl_ind+1
  DEALLOCATE(cl_val)

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0


  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=0
           
           mm=0
           filtro(ii)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 mm=mm+1
                 IF (mm>=dim_caja) THEN
                    ALLOCATE(aux_caja(dim_caja))
                    aux_caja=caja
                    DEALLOCATE(caja)
                    ALLOCATE(caja(mm))
                    caja(1:dim_caja)=aux_caja(:)
                    DEALLOCATE(aux_caja)
                    dim_caja=mm
                 END IF
                 caja(mm)=amigo
                 filtro(amigo)=.TRUE.
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo2
                    filtro(amigo2)=.TRUE.
                 END IF
              END DO
           END DO

           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 IF (gg==0) THEN
                    ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                 ELSE
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(:)=pcl_ind(ii)%i1(:)
                    box_val(:)=pcl_val(ii)%d1(:)
                    DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                    ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                    pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                    pcl_val(ii)%d1(1:gg)=box_val(:)
                    DEALLOCATE(box_ind,box_val)
                 END IF
                 gg=gg+1
                 pcl_ind(ii)%i1(gg)=aj
                 pcl_val(ii)%d1(gg)=val_aux
              END IF
           END DO
           pcl_num(ii)=gg
        END DO

     ELSE

        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=pcl_num(ii)

           mm=0
           filtro(ai)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 IF (amigo>ai) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo
                    filtro(amigo)=.TRUE.
                 END IF
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    IF (amigo2>ai) THEN
                       mm=mm+1
                       IF (mm>=dim_caja) THEN
                          ALLOCATE(aux_caja(dim_caja))
                          aux_caja=caja
                          DEALLOCATE(caja)
                          ALLOCATE(caja(mm))
                          caja(1:dim_caja)=aux_caja(:)
                          DEALLOCATE(aux_caja)
                          dim_caja=mm
                       END IF
                       caja(mm)=amigo2
                       filtro(amigo2)=.TRUE.
                    END IF
                 END IF
              END DO
           END DO

           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 ll=pcl_num(aj)
                 IF (gg==0) THEN
                    ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                 ELSE
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(:)=pcl_ind(ii)%i1(:)
                    box_val(:)=pcl_val(ii)%d1(:)
                    DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                    ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                    pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                    pcl_val(ii)%d1(1:gg)=box_val(:)
                    DEALLOCATE(box_ind,box_val)
                 END IF
                 IF (ll==0) THEN
                    ALLOCATE(pcl_ind(aj)%i1(1),pcl_val(aj)%d1(1))
                 ELSE
                    ALLOCATE(box_ind(ll),box_val(ll))
                    box_ind(:)=pcl_ind(aj)%i1(:)
                    box_val(:)=pcl_val(aj)%d1(:)
                    DEALLOCATE(pcl_ind(aj)%i1,pcl_val(aj)%d1)
                    ALLOCATE(pcl_ind(aj)%i1(ll+1),pcl_val(aj)%d1(ll+1))
                    pcl_ind(aj)%i1(1:ll)=box_ind(:) 
                    pcl_val(aj)%d1(1:ll)=box_val(:)
                    DEALLOCATE(box_ind,box_val)
                 END IF
                 gg=gg+1
                 ll=ll+1
                 pcl_ind(ii)%i1(gg)=aj
                 pcl_val(ii)%d1(gg)=val_aux
                 pcl_ind(aj)%i1(ll)=ai
                 pcl_val(aj)%d1(ll)=val_aux
                 pcl_num(aj)=ll
              END IF
           END DO
           pcl_num(ii)=gg
        END DO

     END IF
  END IF

  gg=SUM(pcl_num(:),DIM=1)

  DEALLOCATE(filtro,caja)
  DEALLOCATE(cl_ind,cl_start)
  ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))

  gg=0
  DO ii=1,natom1
     cl_start(ii)=gg
     DO jj=1,pcl_num(ii)
        gg=gg+1
        cl_ind(gg)=pcl_ind(ii)%i1(jj)
        cl_val(gg)=pcl_val(ii)%d1(jj)
     END DO
     DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
  END DO
  cl_start(natom1+1)=gg

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)

  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF

END SUBROUTINE UPDATE_CONTACT_LIST




SUBROUTINE DISTANCE (diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,matrix)

IMPLICIT NONE

INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1
integer,intent(in)::n1,n2,natom1,natom2
INTEGER,DIMENSION(n1),INTENT(IN)::list1
INTEGER,DIMENSION(n2),INTENT(IN)::list2
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,dimension(natom2,3),intent(in)::coors2
double precision,dimension(n1,n2),intent(out)::matrix
integer::i,j,ai,aj
double precision,dimension(:),allocatable::vect,vect_aux
integer,dimension(:),allocatable::llist1,llist2
double precision::val_aux

ALLOCATE(vect(3),vect_aux(3))
ALLOCATE(llist1(n1),llist2(n2))
llist1=list1+1
llist2=list2+1

matrix=0.0d0

IF ((diff_syst==1) .or. (diff_set==1)) THEN
   IF (pbc_opt==1) THEN
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=1,n2
            aj=llist2(j)
            vect=(coors2(aj,:)-vect_aux)
            CALL PBC (vect,box1,ortho1)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
         end do
      end do
   ELSE
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=1,n2
            aj=llist2(j)
            vect=(coors2(aj,:)-vect_aux)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
         end do
      end do
   END IF
ELSE
   IF (pbc_opt==1) THEN
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=i+1,n2
            aj=llist2(j)
            vect=(coors1(aj,:)-vect_aux)
            CALL PBC (vect,box1,ortho1)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
            matrix(j,i)=val_aux
         end do
      end do
   ELSE
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=i+1,n2
            aj=llist2(j)
            vect=(coors1(aj,:)-vect_aux)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
            matrix(j,i)=val_aux
         end do
      end do
   END IF
END IF

DEALLOCATE(vect,vect_aux)
DEALLOCATE(llist1,llist2)

END SUBROUTINE DISTANCE


SUBROUTINE DISTANCE_IMAGES (diff_syst,diff_set,list1,coors1,box1,ortho1,list2,coors2,&
                            n1,n2,natom1,natom2,min_dists,ind_atoms_min,min_image)

IMPLICIT NONE
  
INTEGER,INTENT(IN)::diff_syst,diff_set,ortho1
integer,intent(in)::n1,n2,natom1,natom2
INTEGER,DIMENSION(n1),INTENT(IN)::list1
INTEGER,DIMENSION(n2),INTENT(IN)::list2
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,dimension(natom2,3),intent(in)::coors2
double precision,dimension(n1),intent(out)::min_dists
integer,dimension(n1),intent(out):: ind_atoms_min
integer,dimension(n1,3),intent(out):: min_image

integer::ii,jj,gg,kk,ai,aj,ind_aux_min,val_imin,prov_imin
double precision,dimension(:),allocatable::vect,vect_aux,vect_aux2
integer,dimension(:),allocatable::llist1,llist2
double precision::val_aux,val_aux_min,val_ref_min,prov_val_min
INTEGER,DIMENSION(:,:),ALLOCATABLE::imag
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::trans_imag

val_ref_min=10000.0d0
DO ii=1,3
   val_aux=dot_product(box1(ii,:),box1(ii,:))
   IF (val_ref_min>val_aux) THEN
      val_ref_min=val_aux
   END IF
END DO

ALLOCATE(vect(3),vect_aux(3),vect_aux2(3),imag(26,3),trans_imag(26,3))
gg=0
DO ii=-1,1
   DO jj=-1,1
      DO kk=-1,1
         IF ((ii/=0).or.(jj/=0).or.(kk/=0)) THEN
            gg=gg+1
            imag(gg,:)=(/ii,jj,kk/)
            trans_imag(gg,:)=ii*box1(1,:)+jj*box1(2,:)+kk*box1(3,:)
         END IF
      END DO
   END DO
END DO

ALLOCATE(llist1(n1),llist2(n2))
llist1=list1+1
llist2=list2+1

IF (diff_set==1) THEN
   DO ii=1,n1
      ai=llist1(ii)
      vect_aux=coors1(ai,:)
      val_aux_min=val_ref_min
      DO jj=1,n2
         aj=llist2(jj)
         vect_aux2=coors2(aj,:)-vect_aux
         DO gg=1,26
            vect=vect_aux2+trans_imag(gg,:)
            val_aux=dot_product(vect,vect)
            IF (val_aux<val_aux_min) THEN
               val_aux_min=val_aux
               val_imin=gg
               ind_aux_min=jj
            END IF
         END DO
      END DO
      min_dists(ii)=val_aux_min
      ind_atoms_min(ii)=ind_aux_min
      min_image(ii,1)=val_imin
   END DO
ELSE
   min_dists(:)=val_ref_min
   DO ii=1,n1
      ai=llist1(ii)
      vect_aux=coors1(ai,:)
      val_aux_min=min_dists(ii)
      DO jj=ii+1,n2
         aj=llist2(jj)
         vect_aux2=coors2(aj,:)-vect_aux
         prov_val_min=val_ref_min
         DO gg=1,26
            vect=vect_aux2+trans_imag(gg,:)
            val_aux=dot_product(vect,vect)
            IF (val_aux<prov_val_min) THEN
               prov_val_min=val_aux
               prov_imin=gg
            END IF
         END DO
         IF (prov_val_min<val_aux_min) THEN
            val_aux_min=prov_val_min
            val_imin=prov_imin
            ind_aux_min=jj
         END IF
         IF (prov_val_min<min_dists(jj)) THEN
            min_dists(jj)=prov_val_min
            ind_atoms_min(jj)=ii
            min_image(jj,1)=prov_imin
         END IF
      END DO
      min_dists(ii)=val_aux_min
      ind_atoms_min(ii)=ind_aux_min
      min_image(ii,1)=val_imin
   END DO
END IF

DO ii=1,n1
   min_dists(ii)=sqrt(min_dists(ii))
   ind_atoms_min(ii)=list2(ind_atoms_min(ii))
   min_image(ii,:)=imag(min_image(ii,1),:)
END DO

DEALLOCATE(vect,vect_aux,vect_aux2,imag,trans_imag)
DEALLOCATE(llist1,llist2)

END SUBROUTINE DISTANCE_IMAGES


SUBROUTINE RADIUS_GYRATION (list1,coors1,box1,ortho1,n1,natom1,val_Rg)

IMPLICIT NONE

INTEGER,INTENT(IN)::ortho1
integer,intent(in)::n1,natom1
INTEGER,DIMENSION(n1),INTENT(IN)::list1
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,INTENT(OUT)::val_Rg

integer::ii,ai
double precision,dimension(:),allocatable:: cdm,vect_aux
integer,dimension(:),allocatable::llist1

ALLOCATE(llist1(n1),cdm(3),vect_aux(3))
llist1=list1+1

val_Rg=0.0d0
cdm=0.0d0

DO ii=1,n1
   ai=llist1(ii)
   cdm=cdm+coors1(ai,:)
END DO
cdm=cdm/(n1*1.0d0)

DO ii=1,n1
   ai=llist1(ii)
   vect_aux=(coors1(ai,:)-cdm)
   val_Rg=val_Rg+dot_product(vect_aux,vect_aux)
END DO

DEALLOCATE(llist1,cdm,vect_aux)
val_Rg=sqrt(val_Rg/(1.0d0*n1))

END SUBROUTINE RADIUS_GYRATION


SUBROUTINE PRINCIPAL_INERTIA_AXIS (list1,coors1,box1,ortho1,n1,natom1,axis)

IMPLICIT NONE

INTEGER,INTENT(IN)::ortho1
integer,intent(in)::n1,natom1
INTEGER,DIMENSION(n1),INTENT(IN)::list1
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,DIMENSION(3,3),INTENT(OUT)::axis


integer::ii,ai
double precision:: dd
double precision,dimension(:),allocatable:: cdm,vect_aux,values
double precision,dimension(:,:),allocatable:: matrix
integer,dimension(:),allocatable::llist1
integer::num_val,info
INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work

ALLOCATE(work(8*3),iwork(5*3),ifail(3))
ALLOCATE(llist1(n1),cdm(3),vect_aux(3),matrix(3,3),values(3))

llist1=list1+1

cdm=0.0d0
matrix=0.0d0

DO ii=1,n1
   ai=llist1(ii)
   cdm=cdm+coors1(ai,:)
END DO
cdm=cdm/(n1*1.0d0)

DO ii=1,n1
   ai=llist1(ii)
   vect_aux=(coors1(ai,:)-cdm)
   dd=dot_product(vect_aux,vect_aux)
   matrix(1,1)=matrix(1,1)+dd-vect_aux(1)**2
   matrix(2,2)=matrix(2,2)+dd-vect_aux(2)**2
   matrix(3,3)=matrix(3,3)+dd-vect_aux(3)**2
   matrix(1,2)=matrix(1,2)-vect_aux(1)*vect_aux(2)
   matrix(1,3)=matrix(1,3)-vect_aux(1)*vect_aux(3)
   matrix(2,3)=matrix(2,3)-vect_aux(2)*vect_aux(3)
END DO
matrix(2,1)=matrix(1,2)
matrix(3,1)=matrix(1,3)
matrix(3,2)=matrix(2,3)


matrix=matrix/(n1*1.0d0)

DEALLOCATE(llist1)

values=0.0d0

CALL dsyevx ('V','I','U',3,matrix,3,0,0,1,3,0.0d0,num_val&
       &,values,axis,3,work,8*3,iwork,ifail,info)

IF (info/=0) THEN
   print*,"Error with the diagonalization."
   print*,"The array 'work' should has the dimension:",work(1)
END IF

DEALLOCATE(work,iwork,ifail)
DEALLOCATE(cdm,vect_aux,matrix,values)


END SUBROUTINE PRINCIPAL_INERTIA_AXIS


SUBROUTINE PRINCIPAL_GEOMETRIC_AXIS (list1,coors1,box1,ortho1,n1,natom1,axis)

IMPLICIT NONE

INTEGER,INTENT(IN)::ortho1
integer,intent(in)::n1,natom1
INTEGER,DIMENSION(n1),INTENT(IN)::list1
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,DIMENSION(3,3),INTENT(OUT)::axis


integer::ii,ai
double precision:: dd
double precision,dimension(:),allocatable:: cdm,vect_aux,values
double precision,dimension(:,:),allocatable:: matrix
integer,dimension(:),allocatable::llist1
integer::num_val,info
INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work

ALLOCATE(work(8*3),iwork(5*3),ifail(3))
ALLOCATE(llist1(n1),cdm(3),vect_aux(3),matrix(3,3),values(3))

llist1=list1+1

cdm=0.0d0
matrix=0.0d0

DO ii=1,n1
   ai=llist1(ii)
   cdm=cdm+coors1(ai,:)
END DO
cdm=cdm/(n1*1.0d0)

DO ii=1,n1
   ai=llist1(ii)
   vect_aux=(coors1(ai,:)-cdm)
   dd=dot_product(vect_aux,vect_aux)
   matrix(1,1)=matrix(1,1)+vect_aux(1)**2
   matrix(2,2)=matrix(2,2)+vect_aux(2)**2
   matrix(3,3)=matrix(3,3)+vect_aux(3)**2
   matrix(1,2)=matrix(1,2)+vect_aux(1)*vect_aux(2)
   matrix(1,3)=matrix(1,3)+vect_aux(1)*vect_aux(3)
   matrix(2,3)=matrix(2,3)+vect_aux(2)*vect_aux(3)
END DO
matrix(2,1)=matrix(1,2)
matrix(3,1)=matrix(1,3)
matrix(3,2)=matrix(2,3)

matrix=matrix/(n1*1.0d0)

DEALLOCATE(llist1)

values=0.0d0

CALL dsyevx ('V','I','U',3,matrix,3,0,0,1,3,0.0d0,num_val&
       &,values,axis,3,work,8*3,iwork,ifail,info)

IF (info/=0) THEN
   print*,"Error with the diagonalization."
   print*,"The array 'work' should has the dimension:",work(1)
END IF

DEALLOCATE(work,iwork,ifail)
DEALLOCATE(cdm,vect_aux,matrix,values)


END SUBROUTINE PRINCIPAL_GEOMETRIC_AXIS


SUBROUTINE neighbs_ranking (diff_syst,diff_set,pbc_opt,limit,list1,coors1,box1,ortho1,list2,coors2,&
                            n1,n2,natom1,natom2,neighb_list) !before: neighb_dist,neighb_uvect

  IMPLICIT NONE

  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1
  INTEGER,INTENT(IN)::limit
  integer,intent(in)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2
  INTEGER,dimension(n1,limit),intent(out)::neighb_list

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dist_matrix
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dist_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::neight_aux
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  INTEGER::ii,jj,gg

  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_dist
  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit,3),INTENT(OUT)::neighb_uvect

  ! The indexes in list1 and list2 not corrected yet (because of function dist)

  ALLOCATE(dist_matrix(n1,n2),dist_aux(n2),filter(n2),neight_aux(limit))
  filter=.TRUE.
  CALL distance (diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,dist_matrix)

  IF ((diff_syst==1) .or. (diff_set==1)) THEN
     DO ii=1,n1
        dist_aux=dist_matrix(ii,:)
        DO jj=1,limit
           gg=MINLOC(dist_aux(:),DIM=1,MASK=filter(:))
           neight_aux(jj)=gg
           filter(gg)=.FALSE.
        END DO
        DO jj=1,limit
           gg=neight_aux(jj)
           filter(gg)=.TRUE.
           neighb_list(ii,jj)=list2(gg) ! No correction was done on indexes
        END DO
     END DO
  ELSE
     DO ii=1,n1
        dist_aux=dist_matrix(ii,:)
        filter(ii)=.FALSE.
        DO jj=1,limit
           gg=MINLOC(dist_aux(:),DIM=1,MASK=filter(:))
           neight_aux(jj)=gg
           filter(gg)=.FALSE.
        END DO
        DO jj=1,limit
           gg=neight_aux(jj)
           filter(gg)=.TRUE.
           neighb_list(ii,jj)=list2(gg) ! No correction was done on indexes
        END DO
        filter(ii)=.TRUE.
     END DO
  END IF

  DEALLOCATE(dist_matrix,dist_aux,filter,neight_aux)

END SUBROUTINE neighbs_ranking


SUBROUTINE neighbs_dist (diff_syst,diff_set,pbc_opt,limit,list1,coors1,box1,ortho1,list2,&
                         coors2,n1,n2,natom1,natom2,contact_map,num_neighbs,dist_matrix) !before: neighb_dist,neighb_uvect

  IMPLICIT NONE

  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1
  DOUBLE PRECISION,INTENT(IN)::limit
  integer,intent(in)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2
  INTEGER,dimension(n1,n2),intent(out)::contact_map
  INTEGER,DIMENSION(n1),INTENT(OUT)::num_neighbs
  DOUBLE PRECISION,DIMENSION(n1,n2),intent(out)::dist_matrix

  INTEGER::ii,jj,gg

  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_dist
  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit,3),INTENT(OUT)::neighb_uvect

  ! The indexes in list1 and list2 not corrected yet (because of function dist)
  CALL DISTANCE (diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,dist_matrix)

  contact_map=0
  num_neighbs=0

  IF ((diff_syst==1) .or. (diff_set==1)) THEN
     DO ii=1,n1
        gg=0
        DO jj=1,n2
           IF (dist_matrix(ii,jj)<=limit) THEN
              contact_map(ii,jj)=1
              gg=gg+1
           END IF
        END DO
        num_neighbs(ii)=gg
     END DO
  ELSE
     DO ii=1,n1
        gg=num_neighbs(ii)
        DO jj=ii+1,n2
           IF (dist_matrix(ii,jj)<=limit) THEN
              contact_map(ii,jj)=1
              contact_map(jj,ii)=1
              gg=gg+1
              num_neighbs(jj)=num_neighbs(jj)+1
           END IF
        END DO
        num_neighbs(ii)=gg
     END DO
  END IF
  
END SUBROUTINE neighbs_dist

SUBROUTINE translate_list (sort,list,filter,distances,dim_out,n_list,trans_inds)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::n_list,dim_out,sort
  INTEGER,DIMENSION(n_list),INTENT(IN)::list,filter
  DOUBLE PRECISION,DIMENSION(n_list),INTENT(IN)::distances
  INTEGER,DIMENSION(dim_out),INTENT(OUT)::trans_inds

  LOGICAL,DIMENSION(:),ALLOCATABLE::ifilter
  INTEGER::ii,gg

  IF (sort==1) THEN

     ALLOCATE(ifilter(n_list))
     ifilter=filter

     DO ii=1,dim_out
        gg=MINLOC(distances(:),DIM=1,MASK=ifilter(:))
        ifilter(gg)=.FALSE.
        trans_inds(ii)=list(gg) ! The indexes were not corrected
     END DO

     DEALLOCATE(ifilter)

  ELSE

     gg=0
     DO ii=1,n_list
        IF (filter(ii)==1) THEN
           gg=gg+1
           trans_inds(gg)=list(gg)
        END IF
     END DO

  END IF

END SUBROUTINE translate_list


SUBROUTINE within (list_dists,cutoff,dim_list,ISIN)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::dim_list
  DOUBLE PRECISION,INTENT(IN)::cutoff
  DOUBLE PRECISION,DIMENSION(dim_list),INTENT(IN)::list_dists
  INTEGER,INTENT(OUT)::ISIN
  INTEGER::ii

  ISIN=0

  DO ii=1,dim_list
     IF (list_dists(ii)<=cutoff) THEN
        ISIN=1
        EXIT
     END IF
  END DO

END SUBROUTINE within

!!$SUBROUTINE min_dist_atoms (pbc_opt,eq_opt,coors,box,ortho,list_a,list_b,N_tot,N_a,N_b,ind_a,ind_b,min_dist)
!!$
!!$  IMPLICIT NONE
!!$  integer,intent(in)::N_tot,N_a,N_b,pbc_opt,eq_opt,ortho
!!$  real,dimension(N_tot,3),intent(in)::coors
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box
!!$  INTEGER,DIMENSION(N_a),INTENT(IN)::list_a
!!$  INTEGER,DIMENSION(N_b),INTENT(IN)::list_b
!!$  INTEGER,DIMENSION(N_a)::auxlist_a
!!$  INTEGER,DIMENSION(N_b)::auxlist_b
!!$  INTEGER,INTENT(OUT)::ind_a,ind_b
!!$  REAL,INTENT(OUT)::min_dist
!!$
!!$  REAL,DIMENSION(3)::vect,vect_a
!!$  REAL::aux_dist
!!$  INTEGER::i,j,ia,jb
!!$
!!$  auxlist_a=list_a+1
!!$  auxlist_b=list_b+1
!!$
!!$  vect=0.0d0
!!$  min_dist=1.0d0/0.0d0
!!$  DO i=1,N_a
!!$     ia=auxlist_a(i)
!!$     vect_a=coors(ia,:)
!!$     DO j=1,N_b
!!$        jb=auxlist_b(j)
!!$        IF ((eq_opt==0).or.(jb>ia)) THEN
!!$           vect=(coors(jb,:)-vect_a(:))
!!$           IF (pbc_opt==1) CALL PBC (vect,box,ortho)
!!$           aux_dist=sqrt(dot_product(vect,vect))
!!$           IF (aux_dist<min_dist) THEN
!!$              min_dist=aux_dist
!!$              ind_a=ia
!!$              ind_b=jb
!!$           END IF
!!$        END IF
!!$     END DO
!!$  END DO
!!$  
!!$  ind_a=ind_a-1
!!$  ind_b=ind_b-1
!!$
!!$END SUBROUTINE min_dist_atoms
!!$
!!$SUBROUTINE min_dist_atoms_ref (pbc_opt,coors,box,ortho,list_a,list_coors_b,N_tot,N_a,N_b,ind_a,ind_b,min_dist)
!!$
!!$  IMPLICIT NONE
!!$  integer,intent(in)::N_tot,N_a,N_b,pbc_opt
!!$  real,dimension(N_tot,3),intent(in)::coors
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box
!!$  INTEGER,DIMENSION(N_a),INTENT(IN)::list_a
!!$  REAL,DIMENSION(N_b,3),INTENT(IN)::list_coors_b
!!$  INTEGER,DIMENSION(N_a)::auxlist_a
!!$  INTEGER,INTENT(OUT)::ind_a,ind_b
!!$  REAL,INTENT(OUT)::min_dist
!!$
!!$  REAL,DIMENSION(3)::vect,vect_a
!!$  REAL::aux_dist
!!$  INTEGER::i,j,ia,jb
!!$
!!$  auxlist_a=list_a+1
!!$
!!$  vect=0.0d0
!!$  min_dist=1.0d0/0.0d0
!!$  do i=1,N_a
!!$     ia=auxlist_a(i)
!!$     vect_a=coors(ia,:)
!!$     do j=1,N_b
!!$        vect=(list_coors_b(j,:)-vect_a(:))
!!$        IF (pbc_opt==1) CALL PBC (vect,box,ortho)
!!$        aux_dist=sqrt(dot_product(vect,vect))
!!$        IF (aux_dist<min_dist) THEN
!!$           min_dist=aux_dist
!!$           ind_a=ia
!!$           ind_b=j
!!$        END IF
!!$     end do
!!$  end do
!!$
!!$  ind_a=ind_a-1
!!$  ind_b=ind_b-1
!!$
!!$END SUBROUTINE min_dist_atoms_ref
!!$


!!$SUBROUTINE neighbs_dist2(pbc_opt,ident,ii,dist,coors1,box1,ortho1,coors2,n_atoms2,neighb_list,neighb_dist,neighb_uvect)
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,INTENT(IN)::pbc_opt,ident,ortho1
!!$  REAL,INTENT(IN)::dist
!!$  INTEGER,INTENT(IN)::n_atoms2,ii
!!$  REAL,DIMENSION(3),INTENT(IN)::coors1
!!$  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box1
!!$
!!$  INTEGER,DIMENSION(n_atoms2),INTENT(OUT)::neighb_list
!!$  REAL,DIMENSION(n_atoms2),INTENT(OUT)::neighb_dist
!!$  REAL,DIMENSION(n_atoms2,3),INTENT(OUT)::neighb_uvect
!!$
!!$  LOGICAL::lpbc,lident
!!$  INTEGER::j,g,limit
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::list
!!$  REAL,DIMENSION(:),ALLOCATABLE::list_dists
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
!!$  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
!!$  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
!!$  REAL::norm
!!$
!!$  lident=.FALSE.
!!$  lpbc=.FALSE.
!!$  IF (ident>0) lident=.TRUE.
!!$  IF (pbc_opt>0) lpbc=.TRUE.
!!$
!!$  
!!$
!!$  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))
!!$
!!$  list=0
!!$  list_dists=0.0d0
!!$  list_vects=0.0d0
!!$  aux=0.0d0
!!$  aux2=0.0d0
!!$  filter=.false.
!!$
!!$
!!$  aux=coors1(:)
!!$  limit=0
!!$
!!$  DO j=1,n_atoms2
!!$     aux2=coors2(j,:)-aux
!!$     IF (lpbc.eqv..true.) CALL PBC (aux2,box1,ortho1)
!!$     norm=sqrt(dot_product(aux2,aux2))
!!$     IF (norm<=dist) THEN
!!$        limit=limit+1
!!$        filter(j)=.true.
!!$        list_dists(j)=norm
!!$        list_vects(j,:)=aux2
!!$     END IF
!!$  END DO
!!$  
!!$  IF (lident.eqv..true.) THEN 
!!$     filter(ii)=.false.
!!$     limit=limit-1
!!$  END IF
!!$  
!!$  
!!$  DO j=1,limit
!!$     g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
!!$     list(j)=g
!!$     norm=list_dists(g)
!!$     neighb_dist(j)=norm
!!$     neighb_uvect(j,:)=list_vects(g,:)/norm
!!$     neighb_list(j)=g
!!$     filter(g)=.false.
!!$  END DO
!!$  
!!$  DO j=1,limit
!!$     g=list(j)
!!$     filter(g)=.false.
!!$  END DO
!!$  filter(ii)=.false.
!!$
!!$
!!$  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)
!!$
!!$  neighb_list=neighb_list-1
!!$
!!$
!!$END SUBROUTINE NEIGHBS_DIST2
!!$
!!$
!!$SUBROUTINE neighbs_dist1(pbc_opt,ident,dist,coors1,box1,ortho1,coors2,n_atoms1,n_atoms2,neighb_list,neighb_dist,neighb_uvect)
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,INTENT(IN)::pbc_opt,ident,ortho1
!!$  REAL,INTENT(IN)::dist
!!$  INTEGER,INTENT(IN)::n_atoms1,n_atoms2
!!$  REAL,DIMENSION(n_atoms1,3),INTENT(IN)::coors1
!!$  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box1
!!$  INTEGER,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_list
!!$  REAL,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_dist
!!$  REAL,DIMENSION(n_atoms1,n_atoms2,3),INTENT(OUT)::neighb_uvect
!!$
!!$  LOGICAL::lpbc,lident
!!$  INTEGER::i,j,g,limit
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::list
!!$  REAL,DIMENSION(:),ALLOCATABLE::list_dists
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
!!$  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
!!$  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
!!$  REAL::norm
!!$
!!$  lident=.FALSE.
!!$  lpbc=.FALSE.
!!$  IF (ident>0) lident=.TRUE.
!!$  IF (pbc_opt>0) lpbc=.TRUE.
!!$
!!$
!!$  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))
!!$
!!$  list=0
!!$  list_dists=0.0d0
!!$  list_vects=0.0d0
!!$  aux=0.0d0
!!$  aux2=0.0d0
!!$  filter=.false.
!!$
!!$  DO i=1,n_atoms1
!!$     aux=coors1(i,:)
!!$     limit=0
!!$     DO j=1,n_atoms2
!!$        aux2=coors2(j,:)-aux
!!$        IF (lpbc.eqv..true.) CALL PBC (aux2,box1,ortho1)
!!$        norm=sqrt(dot_product(aux2,aux2))
!!$        IF (norm<=dist) THEN
!!$           limit=limit+1
!!$           filter(j)=.true.
!!$           list_dists(j)=norm
!!$           list_vects(j,:)=aux2
!!$        END IF
!!$     END DO
!!$
!!$     IF (lident.eqv..true.) THEN 
!!$        filter(i)=.false.
!!$        limit=limit-1
!!$     END IF
!!$
!!$     DO j=1,limit
!!$        g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
!!$        list(j)=g
!!$        norm=list_dists(g)
!!$        neighb_dist(i,j)=norm
!!$        neighb_uvect(i,j,:)=list_vects(g,:)/norm
!!$        neighb_list(i,j)=g
!!$        filter(g)=.false.
!!$     END DO
!!$
!!$     DO j=1,limit
!!$        g=list(j)
!!$        filter(g)=.false.
!!$     END DO
!!$     filter(i)=.false.
!!$  END DO
!!$
!!$  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)
!!$
!!$  neighb_list=neighb_list-1
!!$
!!$END SUBROUTINE NEIGHBS_DIST1


!!$subroutine min_rmsd(struct_ref,struct_2,N,U,center_ref,center_2,rmsd,g)
!!$
!!$INTEGER,INTENT(IN)::N
!!$REAL,DIMENSION(N,3),INTENT(IN)::struct_ref,struct_2
!!$
!!$DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::U
!!$DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::center_ref,center_2
!!$DOUBLE PRECISION,INTENT(OUT)::rmsd
!!$DOUBLE PRECISION,DIMENSION(N,3),INTENT(OUT)::g
!!$
!!$
!!$INTEGER::i,j
!!$DOUBLE PRECISION,DIMENSION(N,3)::x,y
!!$DOUBLE PRECISION,DIMENSION(N)::w
!!$DOUBLE PRECISION::sw,msd,x_norm,y_norm
!!$DOUBLE PRECISION,DIMENSION(3,3)::R
!!$DOUBLE PRECISION,DIMENSION(4,4)::F
!!$DOUBLE PRECISION,DIMENSION(3)::tmp
!!$
!!$!To diagonalise:
!!$DOUBLE PRECISION,DIMENSION(4,4)::CC
!!$INTEGER::num_val,info
!!$INTEGER, DIMENSION(:),ALLOCATABLE::iwork,ifail
!!$DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE::values,work
!!$DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::vectors
!!$
!!$
!!$ALLOCATE(values(4),vectors(4,4),work(8*4),iwork(5*4),ifail(4))
!!$
!!$w=0.0d0
!!$x=0.0d0
!!$y=0.0d0
!!$CC=0.0d0
!!$rmsd=0.0d0
!!$msd=0.0d0
!!$sw=0.0d0
!!$values=0.0d0
!!$vectors=0.0d0
!!$center_ref=0.0d0
!!$center_2=0.0d0
!!$U=0.0d0
!!$g=0.0d0
!!$x_norm=0.0d0
!!$y_norm=0.0d0
!!$R=0.0d0
!!$F=0.0d0
!!$tmp=0.0d0
!!$
!!$
!!$!!! copio y peso las coordenadas:
!!$w=1.0d0
!!$DO i=1,N
!!$   sw=w(i)
!!$   x(i,:)=sw*dble(struct_ref(i,:))
!!$   y(i,:)=sw*dble(struct_2(i,:))
!!$END DO
!!$
!!$!!! calculo baricentros, centroides y normas:
!!$
!!$DO i=1,3
!!$   center_ref(i)=sum(x(:,i))/dble(N)
!!$   center_2(i)=sum(y(:,i))/dble(N)
!!$   x(:,i)=x(:,i)-center_ref(i)
!!$   y(:,i)=y(:,i)-center_2(i)
!!$   x_norm=x_norm+dot_product(x(:,i),x(:,i))
!!$   y_norm=y_norm+dot_product(y(:,i),y(:,i))
!!$END DO
!!$
!!$!!! calculo la matriz R
!!$DO i=1,3
!!$   DO j=1,3
!!$      R(i,j)=dot_product(x(:,i),y(:,j))
!!$   END DO
!!$END DO
!!$
!!$!!! construimos la matriz F:
!!$
!!$F(1,1)=R(1,1)+R(2,2)+R(3,3)
!!$F(2,1)=R(2,3)-R(3,2)
!!$F(3,1)=R(3,1)-R(1,3)
!!$F(4,1)=R(1,2)-R(2,1)
!!$F(1,2)=F(2,1)
!!$F(2,2)=R(1,1)-R(2,2)-R(3,3)
!!$F(3,2)=R(1,2)+R(2,1)
!!$F(4,2)=R(1,3)+R(3,1)
!!$F(1,3)=F(3,1)
!!$F(2,3)=F(3,2)
!!$F(3,3)=-R(1,1)+R(2,2)-R(3,3)
!!$F(4,3)=R(2,3)+R(3,2)
!!$F(1,4)=F(4,1)
!!$F(2,4)=F(4,2)
!!$F(3,4)=F(4,3)
!!$F(4,4)=-R(1,1)-R(2,2)+R(3,3) 
!!$
!!$!!! calculos los autovalores y autovectores:
!!$CC=F
!!$call dsyevx ('V','I','U',4,CC,4,0,0,1,4,0.0d0,num_val&
!!$     &,values,vectors,4,work,8*4,iwork,ifail,info)
!!$
!!$!!! computo el rmsd, la matriz de rotacion y g
!!$
!!$msd=max(0.0d0,((x_norm+y_norm)-2.0d0*values(4)))/dble(N)
!!$rmsd=sqrt(msd)
!!$
!!$
!!$call rotation_matrix(vectors(:,4),U)
!!$
!!$DO i=1,N
!!$   DO j=1,3
!!$      tmp(:)=matmul(transpose(U(:,:)),y(i,:))
!!$      g(i,j)=(x(i,j)-tmp(j))/(rmsd*dble(N))
!!$   END DO
!!$END DO
!!$
!!$!!! calculo las nuevas posiciones con la traslacion y rotacion
!!$
!!$!DO i=1,N
!!$!   pos_new(i,:)=matmul(transpose(U(:,:)),struct_2(i,:)-center_2)+center_ref
!!$!END DO
!!$
!!$
!!$END subroutine min_rmsd
!!$
!!$subroutine rotation_matrix(q, U)
!!$
!!$DOUBLE PRECISION,DIMENSION(4),INTENT(in)::q
!!$DOUBLE PRECISION,DIMENSION(3,3),INTENT(out)::U
!!$DOUBLE PRECISION::q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33
!!$
!!$q0=q(1)
!!$q1=q(2)
!!$q2=q(3)
!!$q3=q(4)
!!$
!!$b0=2.0d0*q0
!!$b1=2.0d0*q1
!!$b2=2.0d0*q2
!!$b3=2.0d0*q3
!!$
!!$q00=b0*q0-1.0d0
!!$q01=b0*q1
!!$q02=b0*q2
!!$q03=b0*q3
!!$
!!$q11=b1*q1
!!$q12=b1*q2
!!$q13=b1*q3  
!!$
!!$q22=b2*q2
!!$q23=b2*q3
!!$
!!$q33=b3*q3 
!!$
!!$U(1,1)=q00+q11
!!$U(1,2)=q12-q03
!!$U(1,3)=q13+q02
!!$
!!$U(2,1)=q12+q03
!!$U(2,2)=q00+q22
!!$U(2,3)=q23-q01
!!$
!!$U(3,1)=q13-q02
!!$U(3,2)=q23+q01
!!$U(3,3)=q00+q33
!!$
!!$end subroutine rotation_matrix
!!$
!!$
!!$subroutine rot_trans(struct,rot,center_2,center_ref,N,new_struct)
!!$
!!$IMPLICIT NONE
!!$INTEGER,INTENT(IN)::N
!!$REAL,DIMENSION(N,3),INTENT(IN)::struct
!!$DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::rot
!!$DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::center_2,center_ref
!!$REAL,DIMENSION(N,3),INTENT(OUT)::new_struct
!!$INTEGER::i
!!$
!!$DO i=1,N
!!$   new_struct(i,:)=matmul(transpose(rot(:,:)),(struct(i,:)-center_2(:)))+center_ref(:)
!!$END DO
!!$
!!$END subroutine rot_trans
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$subroutine proj3d(vect1,vect2,N,val)
!!$
!!$IMPLICIT NONE
!!$INTEGER,INTENT(IN)::N
!!$REAL,DIMENSION(N,3),INTENT(IN)::vect1,vect2
!!$REAL,DIMENSION(N,3)::vect_norm
!!$REAL::norm
!!$REAL,INTENT(OUT)::val
!!$INTEGER::i
!!$
!!$vect_norm=0.0d0
!!$norm=0.0d0
!!$DO i=1,N
!!$   norm=norm+dot_product(vect2(i,:),vect2(i,:))
!!$END DO
!!$vect_norm=vect2/(sqrt(norm))
!!$
!!$val=0.0d0
!!$DO i=1,N
!!$   val=val+dot_product(vect1(i,:),vect_norm(i,:))
!!$END DO
!!$
!!$end subroutine proj3d
!!$
!!$
!!$
END MODULE GLOB



MODULE RDF

CONTAINS

SUBROUTINE rdf_frame(distances,box,segment_min,segment_max,bins,n_A,n_B,rdf)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::n_A,n_B,bins
  DOUBLE PRECISION,INTENT(IN)::segment_min,segment_max
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
  DOUBLE PRECISION,DIMENSION(n_A,n_B),INTENT(IN)::distances
  DOUBLE PRECISION,DIMENSION(bins),INTENT(OUT)::rdf

  INTEGER::ii,jj,kk,tt
  DOUBLE PRECISION::pi,factor,density,aux
  DOUBLE PRECISION::ro,delta_x

  pi=acos(-1.0d0)
  delta_x=(segment_max-segment_min)/(1.0d0*bins)
  factor=0.0d0
  factor=(4.0d0*pi*n_B*delta_x)
  density=(box(1,1)*box(2,2)*box(3,3))/(1.0d0*n_A)
  factor=density/factor
  rdf=0.0d0

  DO jj=1,n_A
     DO ii=1,n_B
        
        aux=(distances(jj,ii)-segment_min)/delta_x
        tt=FLOOR(aux)+1
        rdf(tt)=rdf(tt)+factor

     END DO
  END DO

  ro=segment_min+delta_x/2.0d0
  DO ii=1,bins
     rdf(ii)=rdf(ii)/(ro**2)
     ro=ro+delta_x
  END DO
  
END SUBROUTINE rdf_frame


END MODULE RDF




MODULE HBONDS

USE GLOB

!! Parameters
INTEGER::definition
DOUBLE PRECISION::sk_param,roh2_param,roo2_param
DOUBLE PRECISION::cos_angooh_param  ! the cosine


!!Output
INTEGER,DIMENSION(:),ALLOCATABLE::hbs_s_A,hbs_inds_A,hbs_s_B,hbs_inds_B
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::hbs_vals_A,hbs_vals_B

CONTAINS

SUBROUTINE FREE_MEMORY ()

  IF (ALLOCATED(hbs_s_A)) DEALLOCATE(hbs_s_A)
  IF (ALLOCATED(hbs_s_B)) DEALLOCATE(hbs_s_B)
  IF (ALLOCATED(hbs_inds_A)) DEALLOCATE(hbs_inds_A)
  IF (ALLOCATED(hbs_inds_B)) DEALLOCATE(hbs_inds_B)
  IF (ALLOCATED(hbs_vals_A)) DEALLOCATE(hbs_vals_A)
  IF (ALLOCATED(hbs_vals_B)) DEALLOCATE(hbs_vals_B)

END SUBROUTINE FREE_MEMORY


SUBROUTINE GET_HBONDS (effic,diff_syst,diff_set,pbc_opt,acc_A,acc_sH_A,acc_H_A,don_A,don_sH_A,don_H_A,coors1,box1,ortho1, &
     acc_B,acc_sH_B,acc_H_B,don_B,don_sH_B,don_H_B,coors2,nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, &
     nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H,natomA,natomB)

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  INTEGER,INTENT(IN)::effic,diff_syst,diff_set,pbc_opt,ortho1
  INTEGER,INTENT(IN)::nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H
  INTEGER,INTENT(IN)::nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H
  INTEGER,DIMENSION(nA_acc),INTENT(IN)    ::acc_A
  INTEGER,DIMENSION(nB_acc),INTENT(IN)    ::acc_B
  INTEGER,DIMENSION(nA_don),INTENT(IN)    ::don_A
  INTEGER,DIMENSION(nB_don),INTENT(IN)    ::don_B
  INTEGER,DIMENSION(nA_acc_sH),INTENT(IN) ::acc_sH_A
  INTEGER,DIMENSION(nB_acc_sH),INTENT(IN) ::acc_sH_B
  INTEGER,DIMENSION(nA_acc_H),INTENT(IN)  ::acc_H_A
  INTEGER,DIMENSION(nB_acc_H),INTENT(IN)  ::acc_H_B
  INTEGER,DIMENSION(nA_don_sH),INTENT(IN) ::don_sH_A
  INTEGER,DIMENSION(nB_don_sH),INTENT(IN) ::don_sH_B
  INTEGER,DIMENSION(nA_don_H),INTENT(IN)  ::don_H_A
  INTEGER,DIMENSION(nB_don_H),INTENT(IN)  ::don_H_B
  DOUBLE PRECISION,DIMENSION(natomA,3),intent(in)::coors1
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box1
  DOUBLE PRECISION,DIMENSION(natomB,3),intent(in)::coors2

  TYPE(iarray_pointer),DIMENSION(:),POINTER::hbs_a_ind,hbs_b_ind 
  TYPE(darray_pointer),DIMENSION(:),POINTER::hbs_a_val,hbs_b_val 
  INTEGER,DIMENSION(:),ALLOCATABLE::num_hbs_a,num_hbs_b

  INTEGER::ii,jj,gg
  INTEGER::don,don_h,acc
  INTEGER::lim_hbs

  DOUBLE PRECISION,DIMENSION(3)::pos_acc,pos_don,pos_h
  DOUBLE PRECISION,DIMENSION(3)::vect_don_acc,vect_h_acc,vect_don_h
  DOUBLE PRECISION,DIMENSION(3)::aux_vect_1,aux_vect_2,aux_vect_3
  DOUBLE PRECISION::dist_h_acc,dist_don_acc,dist_don_h,aux_cos,sk_val
  DOUBLE PRECISION::dist2_h_acc,dist2_don_acc,dist2_don_h

  INTEGER::acc_H1,acc_H2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_perp

  INTEGER,DIMENSION(:),ALLOCATABLE::aux_box_ind,aux2_box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::aux_box_val,aux2_box_val

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  lim_hbs=3

  ALLOCATE(aux_box_ind(lim_hbs),aux_box_val(lim_hbs),filtro(lim_hbs))
  filtro=.FALSE.

  ALLOCATE(hbs_a_ind(nA_don_H),hbs_a_val(nA_don_H))
  ALLOCATE(hbs_b_ind(nB_don_H),hbs_b_val(nB_don_H))
  ALLOCATE(num_hbs_a(nA_don_H),num_hbs_b(nB_don_H))

  num_hbs_a=0
  num_hbs_b=0


  IF (diff_syst==0) THEN

     SELECT CASE (definition)

     CASE (1) ! Skinner                                                                                  !SK 
        !!!! Source: R. Kumar, J. R. Schmidt and J. L. Skinner. J. Chem. Phys. 126, 204107 (2007)        !SK 
        ALLOCATE(vect_perp(nB_acc,3))                                                                    !SK 
        DO ii=1,nB_acc                                                                                   !SK
           acc=acc_B(jj)+1                                                                               !SK
           jj=acc_sH_B(ii)+1                                                                             !SK
           acc_H1=acc_H_B(jj)                                                                            !SK 
           acc_H2=acc_H_B(jj+1)                                                                          !SK 
           pos_acc=coors2(acc,:)                                                                         !SK 
           aux_vect_1=coors2(acc_H1,:)-pos_acc(:)                                                        !SK 
           aux_vect_2=coors2(acc_H2,:)-pos_acc(:)                                                        !SK 
           CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)                             !SK 
           vect_perp(ii,:)=aux_vect_3                                                                    !SK 
        END DO                                                                                           !SK 
                                                                                                         !SK 
        DO ii=1,nA_don                                                                                   !SK 
           DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                                           !SK 
              don_h=don_H_A(hh)+1                                                                        !SK 
              pos_h=coors1(don_h,:)                                                                      !SK 
              gg=0                                                                                       !SK 
              DO jj=1,nB_acc                                                                             !SK
                 acc=acc_B(jj)+1                                                                         !SK 
                 vect_h_acc=pos_h(:)-coors2(acc,:)                                                       !SK 
                 IF (pbc_opt) CALL PBC(vect_h_acc,box1,ortho1)                                           !SK 
                 dist_h_acc=sqrt(dot_product(vect_h_acc,vect_h_acc))                                     !SK 
                 aux_cos=dot_product(vect_perp(jj,:),vect_h_acc(:))/dist_h_acc                           !SK 
                 if (aux_cos>=1.0d0) aux_cos=1.0d0                                                       !SK 
                 if (aux_cos<=-1.0d0) aux_cos=-1.0d0                                                     !SK 
                 aux_cos=acos(aux_cos)                                                                   !SK 
                 aux_cos=aux_cos*(90/pi)                                                                 !SK 
                 IF (aux_cos>90) THEN                                                                    !SK 
                    print*,'aquiii error 3.14',aux_cos                                                   !SK 
                    STOP                                                                                 !SK 
                 END IF                                                                                  !SK 
                 IF (aux_cos<0) THEN                                                                     !SK 
                    print*,'aquiii error 3.14',aux_cos                                                   !SK 
                    STOP                                                                                 !SK 
                 END IF                                                                                  !SK 
                 sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)         !SK 
                 IF (sk_val>sk_param) THEN                                                               !SK 
                    gg=gg+1                                                                              !SK 
                    IF (gg>lim_hbs) THEN                                                                 !SK 
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                             !SK 
                       aux2_box_ind=aux_box_ind                                                          !SK 
                       aux2_box_val=aux_box_val                                                          !SK 
                       DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                        !SK 
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                              !SK 
                       aux_box_ind(1:lim_hbs)=aux2_box_ind                                               !SK 
                       aux_box_val(1:lim_hbs)=aux2_box_val                                               !SK 
                       filtro=.FALSE.
                       DEALLOCATE(aux2_box_ind,aux2_box_val)                                             !SK 
                       lim_hbs=gg                                                                        !SK 
                    END IF                                                                               !SK 
                    aux_box_ind(gg)=acc_B(jj)                                                            !SK 
                    aux_box_val(gg)=sk_val                                                               !SK 
                 END IF                                                                                  !SK 
              END DO                                                                                     !SK 
                                                                                                         !SK 
              IF (gg>0) THEN                                                                             !SK 
                 ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                                     !SK 
                 filtro(1:gg)=.TRUE.                                                                     !SK 
                 DO jj=1,gg                                                                              !SK 
                    ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                          !SK 
                    filtro(ll)=.FALSE.                                                                   !SK 
                    hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                                 !SK 
                    hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                                 !SK 
                 END DO                                                                                  !SK 
              END IF                                                                                     !SK 
                                                                                                         !SK 
              num_hbs_a(hh)=gg                                                                           !SK 
                                                                                                         !SK 
           END DO                                                                                        !SK 
        END DO                                                                                           !SK 
        DEALLOCATE(vect_perp)                                                                            !SK 
                                                                                                         !SK 
        IF (diff_set) THEN                                                                               !SK 
                                                                                                         !SK 
           DO ii=1,nA_acc                                                                                !SK 
              acc=acc_A(jj)+1                                                                            !SK
              jj=acc_sH_A(ii)+1                                                                          !SK
              acc_H1=acc_H_A(jj)                                                                         !SK 
              acc_H2=acc_H_A(jj+1)                                                                       !SK 
              pos_acc=coors1(acc,:)                                                                      !SK 
              aux_vect_1=coors1(acc_H1,:)-pos_acc(:)                                                     !SK 
              aux_vect_2=coors1(acc_H2,:)-pos_acc(:)                                                     !SK 
              CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)                          !SK 
              vect_perp(ii,:)=aux_vect_3                                                                 !SK 
           END DO                                                                                        !SK 
                                                                                                         !SK 
           DO ii=1,nB_don                                                                                !SK 
              DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)                                                        !SK 
                 don_h=don_H_B(hh)+1                                                                     !SK 
                 pos_h=coors2(don_h,:)                                                                   !SK 
                 gg=0                                                                                    !SK 
                 DO jj=1,nA_acc                                                                          !SK 
                    acc=acc_A(jj)+1                                                                      !SK 
                    vect_h_acc=pos_h(:)-coors1(acc,:)                                                    !SK 
                    IF (pbc_opt) CALL PBC(vect_h_acc,box1,ortho1)                                        !SK 
                    dist_h_acc=sqrt(dot_product(vect_h_acc,vect_h_acc))                                  !SK 
                    aux_cos=dot_product(vect_perp(jj,:),vect_h_acc(:))/dist_h_acc                        !SK 
                    if (aux_cos>=1.0d0) aux_cos=1.0d0                                                    !SK 
                    if (aux_cos<=-1.0d0) aux_cos=-1.0d0                                                  !SK 
                    aux_cos=acos(aux_cos)                                                                !SK 
                    aux_cos=aux_cos*(90/pi)                                                              !SK 
                    IF (aux_cos>90) THEN                                                                 !SK 
                       print*,'aquiii error 3.14',aux_cos                                                !SK 
                       STOP                                                                              !SK 
                    END IF                                                                               !SK 
                    IF (aux_cos<0) THEN                                                                  !SK 
                       print*,'aquiii error 3.14',aux_cos                                                !SK 
                       STOP                                                                              !SK 
                    END IF                                                                               !SK 
                    sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)      !SK 
                    IF (sk_val>sk_param) THEN                                                            !SK 
                       gg=gg+1                                                                           !SK 
                       IF (gg>lim_hbs) THEN                                                              !SK 
                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                          !SK 
                          aux2_box_ind=aux_box_ind                                                       !SK 
                          aux2_box_val=aux_box_val                                                       !SK 
                          DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                     !SK 
                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                           !SK 
                          aux_box_ind(1:lim_hbs)=aux2_box_ind                                            !SK 
                          aux_box_val(1:lim_hbs)=aux2_box_val                                            !SK 
                          filtro=.FALSE.
                          DEALLOCATE(aux2_box_ind,aux2_box_val)                                          !SK 
                          lim_hbs=gg                                                                     !SK 
                       END IF                                                                            !SK 
                       aux_box_ind(gg)=acc_A(jj)                                                         !SK 
                       aux_box_val(gg)=sk_val                                                            !SK 
                    END IF                                                                               !SK 
                 END DO                                                                                  !SK 
                                                                                                         !SK 
                 IF (gg>0) THEN                                                                          !SK 
                    ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                                  !SK 
                    filtro(1:gg)=.TRUE.                                                                  !SK 
                    DO jj=1,gg                                                                           !SK 
                       ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                       !SK 
                       filtro(ll)=.FALSE.                                                                !SK 
                       hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                                              !SK 
                       hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                                              !SK 
                    END DO                                                                               !SK 
                 END IF                                                                                  !SK 
                                                                                                         !SK 
                 num_hbs_b(hh)=gg                                                                        !SK 
                                                                                                         !SK 
              END DO                                                                                     !SK 
           END DO                                                                                        !SK 
           DEALLOCATE(vect_perp)                                                                         !SK 
                                                                                                         !SK 
        END IF                                                                                           !SK 

     CASE (2) ! R(o,h)                                                                   !ROH 
        !!!! Source: V. J. Buch. J. Chem. Phys. 96, 3814-3823 (1992)                     !ROH 
        DO ii=1,nA_don                                                                   !ROH
           DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                           !ROH 
              don_h=don_H_A(hh)+1                                                        !ROH 
              pos_h=coors1(don_h,:)                                                      !ROH 
              gg=0                                                                       !ROH 
              DO jj=1,nB_acc                                                             !ROH 
                 acc=acc_B(jj)+1                                                         !ROH 
                 vect_h_acc=coors2(acc,:)-pos_h(:)                                       !ROH 
                 IF (pbc_opt) CALL PBC(vect_h_acc,box1,ortho1)                           !ROH 
                 dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)                          !ROH 
                 IF (dist2_h_acc<roh2_param) THEN                                          !ROH 
                    gg=gg+1                                                              !ROH 
                    IF (gg>lim_hbs) THEN                                                 !ROH 
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))             !ROH 
                       aux2_box_ind=aux_box_ind                                          !ROH 
                       aux2_box_val=aux_box_val                                          !ROH 
                       DEALLOCATE(aux_box_ind,aux_box_val,filtro)                        !ROH 
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))              !ROH 
                       aux_box_ind(1:lim_hbs)=aux2_box_ind                               !ROH 
                       aux_box_val(1:lim_hbs)=aux2_box_val                               !ROH 
                       filtro=.FALSE.
                       DEALLOCATE(aux2_box_ind,aux2_box_val)                             !ROH 
                       lim_hbs=gg                                                        !ROH 
                    END IF                                                               !ROH 
                    aux_box_ind(gg)=acc_B(jj)                                            !ROH 
                    aux_box_val(gg)=dist2_h_acc                                           !ROH 
                 END IF                                                                  !ROH 
              END DO                                                                     !ROH 
                                                                                         !ROH 
              IF (gg>0) THEN                                                             !ROH 
                 ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                     !ROH 
                 filtro(1:gg)=.TRUE.                                                     !ROH 
                 DO jj=1,gg                                                              !ROH 
                    ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                          !ROH 
                    filtro(ll)=.FALSE.                                                   !ROH 
                    hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                 !ROH 
                    hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                 !ROH 
                 END DO                                                                  !ROH 
              END IF                                                                     !ROH 
                                                                                         !ROH 
              num_hbs_a(hh)=gg                                                           !ROH 
                                                                                         !ROH 
           END DO                                                                        !ROH 
        END DO                                                                           !ROH 
                                                                                         !ROH 
                                                                                         !ROH 
        IF (diff_set) THEN                                                               !ROH 
                                                                                         !ROH 
           DO ii=1,num_don_B                                                             !ROH 
              DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)                                        !ROH 
                 don_h=don_H_B(hh)+1                                                     !ROH 
                 pos_h=coors2(don_h,:)                                                   !ROH 
                 gg=0                                                                    !ROH 
                 DO jj=1,num_acc_A                                                       !ROH 
                    acc=acc_A(jj)+1                                                      !ROH 
                    vect_h_acc=coors1(acc,:)-pos_h(:)                                    !ROH 
                    IF (pbc_opt) CALL PBC(vect_h_acc,box1,ortho1)                        !ROH 
                    dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)                  !ROH 
                    IF (dist2_h_acc<roh2_param) THEN                                       !ROH 
                       gg=gg+1                                                           !ROH 
                       IF (gg>lim_hbs) THEN                                              !ROH 
                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))          !ROH 
                          aux2_box_ind=aux_box_ind                                       !ROH 
                          aux2_box_val=aux_box_val                                       !ROH 
                          DEALLOCATE(aux_box_ind,aux_box_val,filtro)                     !ROH 
                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))           !ROH 
                          aux_box_ind(1:lim_hbs)=aux2_box_ind                            !ROH 
                          aux_box_val(1:lim_hbs)=aux2_box_val                            !ROH 
                          filtro=.FALSE.
                          DEALLOCATE(aux2_box_ind,aux2_box_val)                          !ROH 
                          lim_hbs=gg                                                     !ROH 
                       END IF                                                            !ROH 
                       aux_box_ind(gg)=acc_A(jj)                                         !ROH 
                       aux_box_val(gg)=dist2_h_acc                                        !ROH 
                    END IF                                                               !ROH 
                 END DO                                                                  !ROH 
                                                                                         !ROH 
                 IF (gg>0) THEN                                                          !ROH 
                    ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                  !ROH 
                    filtro(1:gg)=.TRUE.                                                  !ROH 
                    DO jj=1,gg                                                           !ROH 
                       ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                       !ROH 
                       filtro(ll)=.FALSE.                                                !ROH 
                       hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                              !ROH 
                       hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                              !ROH 
                    END DO                                                               !ROH 
                 END IF                                                                  !ROH 
                                                                                         !ROH 
                 num_hbs_b(hh)=gg                                                        !ROH 
                                                                                         !ROH 
              END DO                                                                     !ROH 
           END DO                                                                        !ROH 
                                                                                         !ROH 
        END IF


     CASE (3) ! R(o,o)-Ang(o,o,h)                                                                        !ROO_ANG 
        !!!! Source: A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)                          !ROO_ANG 
        DO ii=1,nA_don                                                                                   !ROO_ANG 
           don=don_A(ii)+1                                                                               !ROO_ANG 
           pos_don=coors1(don,:)                                                                         !ROO_ANG 
           DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                                           !ROO_ANG 
              don_h=don_H_A(hh)+1                                                                        !ROO_ANG 
              vect_don_h=coors1(don_h,:)-pos_don                                                         !ROO_ANG
              dist2_don_h=dot_product(vect_don_h,vect_don_h)                                             !ROO_ANG
              gg=0                                                                                       !ROO_ANG 
              DO jj=1,nB_acc                                                                             !ROO_ANG 
                 acc=acc_B(jj)+1                                                                         !ROO_ANG 
                 vect_don_acc=coors2(acc,:)-pos_don(:)                                                   !ROO_ANG 
                 IF (pbc_opt) CALL PBC(vect_don_acc,box1,ortho1)                                         !ROO_ANG 
                 dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)                                    !ROO_ANG
                 IF (dist2_don_acc<roo2_param) THEN                                                      !ROO_ANG 
                    aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))       !ROO_ANG
                    IF (aux_cos>cos_angooh_param) THEN                                                   !ROO_ANG
                       gg=gg+1                                                                           !ROO_ANG 
                       IF (gg>lim_hbs) THEN                                                              !ROO_ANG 
                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                          !ROO_ANG 
                          aux2_box_ind=aux_box_ind                                                       !ROO_ANG 
                          aux2_box_val=aux_box_val                                                       !ROO_ANG 
                          DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                     !ROO_ANG 
                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                           !ROO_ANG 
                          aux_box_ind(1:lim_hbs)=aux2_box_ind                                            !ROO_ANG 
                          aux_box_val(1:lim_hbs)=aux2_box_val                                            !ROO_ANG 
                          filtro=.FALSE.                                                                 !ROO_ANG
                          DEALLOCATE(aux2_box_ind,aux2_box_val)                                          !ROO_ANG 
                          lim_hbs=gg                                                                     !ROO_ANG 
                       END IF                                                                            !ROO_ANG 
                       aux_box_ind(gg)=acc_B(jj)                                                         !ROO_ANG 
                       aux_box_val(gg)=dist2_don_acc                                                     !ROO_ANG 
                    END IF                                                                               !ROO_ANG 
                 END IF                                                                                  !ROO_ANG
              END DO                                                                                     !ROO_ANG 
                                                                                                         !ROO_ANG 
              IF (gg>0) THEN                                                                             !ROO_ANG 
                 ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                                     !ROO_ANG 
                 filtro(1:gg)=.TRUE.                                                                     !ROO_ANG 
                 DO jj=1,gg                                                                              !ROO_ANG 
                    ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                          !ROO_ANG 
                    filtro(ll)=.FALSE.                                                                   !ROO_ANG 
                    hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                                 !ROO_ANG 
                    hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                                 !ROO_ANG 
                 END DO                                                                                  !ROO_ANG 
              END IF                                                                                     !ROO_ANG 
                                                                                                         !ROO_ANG 
              num_hbs_a(hh)=gg                                                                           !ROO_ANG 
                                                                                                         !ROO_ANG 
           END DO                                                                                        !ROO_ANG 
        END DO                                                                                           !ROO_ANG 
                                                                                                         !ROO_ANG 
                                                                                                         !ROO_ANG 
        IF (diff_set) THEN                                                                               !ROO_ANG 
                                                                                                         !ROO_ANG 
           DO ii=1,nB_don                                                                                !ROO_ANG 
              don=don_B(ii)+1                                                                            !ROO_ANG 
              pos_don=coors2(don,:)                                                                      !ROO_ANG 
              DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)                                                        !ROO_ANG 
                 don_h=don_H_B(hh)+1                                                                     !ROO_ANG 
                 vect_don_h=coors2(don_h,:)-pos_don                                                      !ROO_ANG
                 dist2_don_h=dot_product(vect_don_h,vect_don_h)                                          !ROO_ANG
                 gg=0                                                                                    !ROO_ANG 
                 DO jj=1,nA_acc                                                                          !ROO_ANG 
                    acc=acc_A(jj)+1                                                                      !ROO_ANG 
                    vect_don_acc=coors1(acc,:)-pos_don(:)                                                !ROO_ANG 
                    IF (pbc_opt) CALL PBC(vect_don_acc,box1,ortho1)                                      !ROO_ANG 
                    dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)                                 !ROO_ANG
                    IF (dist2_don_acc<roo2_param) THEN                                                     !ROO_ANG
                       aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))    !ROO_ANG
                       IF (aux_cos>cos_angooh_param) THEN                                                !ROO_ANG
                          gg=gg+1                                                                        !ROO_ANG 
                          IF (gg>lim_hbs) THEN                                                           !ROO_ANG 
                             ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                       !ROO_ANG 
                             aux2_box_ind=aux_box_ind                                                    !ROO_ANG 
                             aux2_box_val=aux_box_val                                                    !ROO_ANG 
                             DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                  !ROO_ANG 
                             ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                        !ROO_ANG 
                             aux_box_ind(1:lim_hbs)=aux2_box_ind                                         !ROO_ANG 
                             aux_box_val(1:lim_hbs)=aux2_box_val                                         !ROO_ANG 
                             filtro=.FALSE.                                                              !ROO_ANG
                             DEALLOCATE(aux2_box_ind,aux2_box_val)                                       !ROO_ANG 
                             lim_hbs=gg                                                                  !ROO_ANG 
                          END IF                                                                         !ROO_ANG 
                          aux_box_ind(gg)=acc_A(jj)                                                      !ROO_ANG 
                          aux_box_val(gg)=dist2_don_acc                                                  !ROO_ANG
                       END IF                                                                            !ROO_ANG 
                    END IF                                                                               !ROO_ANG
                 END DO                                                                                  !ROO_ANG 
                                                                                                         !ROO_ANG 
                 IF (gg>0) THEN                                                                          !ROO_ANG 
                    ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                                  !ROO_ANG 
                    filtro(1:gg)=.TRUE.                                                                  !ROO_ANG 
                    DO jj=1,gg                                                                           !ROO_ANG 
                       ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                       !ROO_ANG 
                       filtro(ll)=.FALSE.                                                                !ROO_ANG 
                       hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                                              !ROO_ANG 
                       hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                                              !ROO_ANG 
                    END DO                                                                               !ROO_ANG 
                 END IF                                                                                  !ROO_ANG 
                                                                                                         !ROO_ANG 
                 num_hbs_b(hh)=gg                                                                        !ROO_ANG 
                                                                                                         !ROO_ANG 
              END DO                                                                                     !ROO_ANG 
           END DO                                                                                        !ROO_ANG 
                                                                                                         !ROO_ANG 
        END IF                                                                                           !ROO_ANG


     !!!CASE (4) ! Donor-Acceptor-Number
     !!!   !!!! Source: A. D. Hammerich, V. J. Buch. J. Chem. Phys. 128, 111101 (2008)
     !!!   print*, "NOT IMPLEMENTED YET"
     !!! 
     !!!CASE (5) ! Topological
     !!!   !!!! Source: R. H. Henchman and S. J. Irudayam. J. Phys. Chem. B. 114, 16792-16810 (2010)
     !!! 
     !!!   DO ii=1,nA_don                                                                                 
     !!!      DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                                                  
     !!!         don_h=don_H_A(hh)+1                                                                             
     !!!         pos_h=coors1(don_h,:)                                                                       
     !!!         gg=0
     !!!         candidato=0
     !!!         val_out=100.0d0 !! Should be NaN
     !!!         DO jj=1,nB_acc                                                                           
     !!!            acc=acc_B(jj)+1                                                                          
     !!!            vect_h_acc=coors2(acc,:)-pos_don(:)                                                    
     !!!            IF (pbc_opt) CALL PBC(vect_h_acc,box1,ortho1)                                          
     !!!            dist_h_acc=dot_product(vect_h_acc,vect_h_acc)
     !!!            IF (dist_don_acc<val_out) THEN
     !!!               candidato=jj
     !!!            END IF
     !!!         END DO
     !!!         
     !!! 
     !!! 
     !!! 
     !!!                  gg=gg+1                                                                            
     !!!                  IF (gg>lim_hbs) THEN                                                               
     !!!                     ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                           
     !!!                     aux2_box_ind=aux_box_ind                                                        
     !!!                     aux2_box_val=aux_box_val                                                        
     !!!                     DEALLOCATE(aux_box_ind,aux_box_val)                                             
     !!!                     ALLOCATE(aux_box_ind(gg),aux_box_val(gg))                                       
     !!!                     aux_box_ind(1:lim_hbs)=aux2_box_ind                                             
     !!!                     aux_box_val(1:lim_hbs)=aux2_box_val                                             
     !!!                     DEALLOCATE(aux2_box_ind,aux2_box_val)                                           
     !!!                     lim_hbs=gg                                                                      
     !!!                  END IF                                                                             
     !!!                  aux_box_ind(gg)=acc_B(jj)                                                          
     !!!                  aux_box_val(gg)=val_out                                                            
     !!!               END IF                                                                                
     !!!            END IF                                                                                  
     !!!         END DO                                                                                      
     !!!                                                                                                     
     !!!         IF (gg>0) THEN                                                                              
     !!!            ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                                      
     !!!            ALLOCATE(filtro(gg))                                                                     
     !!!            filtro=.TRUE.                                                                            
     !!!            DO jj=1,gg                                                                               
     !!!               ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                           
     !!!               filtro(ll)=.FALSE.                                                                    
     !!!               hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                                  
     !!!               hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                                  
     !!!            END DO                                                                                   
     !!!            DEALLOCATE(filtro)                                                                       
     !!!         END IF                                                                                      
     !!!                                                                                                     
     !!!         num_hbs_a(hh)=gg                                                                            
     !!!                                                                                                     
     !!!      END DO                                                                                         
     !!!   END DO                                                                                            
     !!!                                                                                                     
     !!!                                                                                                     
     !!!   IF (diff_set) THEN                                                                                
     !!!                                                                                                     
     !!!      DO ii=1,num_don_B                                                                              
     !!!         don=don_B(ii,1)+1                                                                           
     !!!         pos_don=coors2(don,:)                                                                       
     !!!         DO hh=H_s_B(ii)+1,H_s_B(ii+1)                                                               
     !!!            don_h=H_B(hh)+1                                                                          
     !!!            pos_h=coors2(don_h,:)                                                                    
     !!!            vect_don_h=coors2(don_h,:)-pos_don_h                                                    
     !!!            gg=0                                                                                     
     !!!            DO jj=1,num_acc_A                                                                        
     !!!               acc=acc_A(jj)+1                                                                       
     !!!               vect_don_acc=coors1(acc,:)-pos_don(:)                                                 
     !!!               IF (pbc_opt) CALL PBC(vect_don_acc,box1,ortho1)                                       
     !!!               dist_don_acc=sqrt(dot_product(vect_don_acc,vect_don_acc))                            
     !!!               IF (dist_don_acc<roo_param) THEN                                                      
     !!!                  dist_don_h_acc=sqrt(dot_product(vect_don_h_acc,vect_don_h_acc))                   
     !!!                  aux_cos=dot_product(vect_don_h_acc,vect_don_acc)/(dist_don_acc*dist_don_h_acc)    
     !!!                  IF (aux_cos>cos_angooh_param) THEN                                                
     !!!                     gg=gg+1                                                                         
     !!!                     IF (gg>lim_hbs) THEN                                                            
     !!!                        ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                        
     !!!                        aux2_box_ind=aux_box_ind                                                     
     !!!                        aux2_box_val=aux_box_val                                                     
     !!!                        DEALLOCATE(aux_box_ind,aux_box_val)                                          
     !!!                        ALLOCATE(aux_box_ind(gg),aux_box_val(gg))                                    
     !!!                        aux_box_ind(1:lim_hbs)=aux2_box_ind                                          
     !!!                        aux_box_val(1:lim_hbs)=aux2_box_val                                          
     !!!                        DEALLOCATE(aux2_box_ind,aux2_box_val)                                        
     !!!                        lim_hbs=gg                                                                   
     !!!                     END IF                                                                          
     !!!                     aux_box_ind(gg)=acc_A(jj)                                                       
     !!!                     aux_box_val(gg)=val_out                                                         
     !!!                  END IF                                                                             
     !!!               END IF                                                                               
     !!!            END DO                                                                                   
     !!!                                                                                                     
     !!!            IF (gg>0) THEN                                                                           
     !!!               ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                                   
     !!!               ALLOCATE(filtro(gg))                                                                  
     !!!               filtro=.TRUE.                                                                         
     !!!               DO jj=1,gg                                                                            
     !!!                  ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                        
     !!!                  filtro(ll)=.FALSE.                                                                 
     !!!                  hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                                               
     !!!                  hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                                               
     !!!               END DO                                                                                
     !!!               DEALLOCATE(filtro)                                                                    
     !!!            END IF                                                                                   
     !!!                                                                                                     
     !!!            num_hbs_b(hh)=gg                                                                         
     !!!                                                                                                     
     !!!         END DO                                                                                      
     !!!      END DO                                                                                         
     !!!                                                                                                     
     !!!   END IF                                                                                           
     !!! 
     !!!CASE (6) ! Donor-Number-Ang(o,o,h)
     !!!   print*, "NOT IMPLEMENTED YET"
     !!! 
     !!!CASE (7) ! Nearest-Neighbour
     !!!   print*, "NOT IMPLEMENTED YET"

     CASE DEFAULT
        PRINT*, 'Error: Hbond definition unknown'
     
     END SELECT

     print*,SUM(num_hbs_a(:),DIM=1),SUM(num_hbs_b(:),DIM=1)



  ELSE
     PRINT*,'NOT IMPLEMENTED YET'
  END IF


  DEALLOCATE(aux_box_ind,aux_box_val)
  DEALLOCATE(hbs_a_ind,hbs_a_val)
  DEALLOCATE(hbs_b_ind,hbs_b_val)
  DEALLOCATE(num_hbs_a,num_hbs_b)
  DEALLOCATE(filtro)


END SUBROUTINE GET_HBONDS


SUBROUTINE PERPENDICULAR_NORMED_VECT (vect1,vect2,vect_out)

  DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::vect1,vect2
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::vect_out

  CALL PRODUCT_VECT(vect1,vect2,vect_out)
  CALL NORMALIZE_VECT(vect_out)

END SUBROUTINE PERPENDICULAR_NORMED_VECT


SUBROUTINE PRODUCT_VECT(a,b,normal)

  DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::a,b
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::normal
  
  normal(1)=a(2)*b(3)-a(3)*b(2)
  normal(2)=-a(1)*b(3)+a(3)*b(1)
  normal(3)=a(1)*b(2)-a(2)*b(1)

END SUBROUTINE PRODUCT_VECT

SUBROUTINE NORMALIZE_VECT (a)

  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::a
  DOUBLE PRECISION::norm

  norm=sqrt(dot_product(a,a))
  a=a/norm

END SUBROUTINE NORMALIZE_VECT


END MODULE HBONDS
