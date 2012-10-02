MODULE GLOB

CONTAINS

SUBROUTINE PBC(vector,box,ortho)
 
  IMPLICIT NONE
 
  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::vector
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
  INTEGER,INTENT(IN)::ortho
  INTEGER::i
  DOUBLE PRECISION::x,L,Lhalf
 
  IF (ortho) THEN
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

SUBROUTINE DISTANCE (diff_system,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,matrix)

IMPLICIT NONE

INTEGER,INTENT(IN)::diff_system,pbc_opt,ortho1
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

IF (diff_system) THEN
   IF (pbc_opt) THEN
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=1,n2
            aj=llist2(j)
            vect=(vect_aux-coors2(aj,:))
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
            vect=(vect_aux-coors2(aj,:))
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
         end do
      end do
   END IF
ELSE
   IF (pbc_opt) THEN
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=1,n2
            aj=llist2(j)
            if (aj>ai) then
               vect=(vect_aux-coors1(aj,:))
               CALL PBC (vect,box1,ortho1)
               val_aux=sqrt(dot_product(vect,vect))
               matrix(i,j)=val_aux
               matrix(j,i)=val_aux
            end if
         end do
      end do
   ELSE
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=1,n2
            aj=llist2(j)
            if (aj>ai) then
               vect=(vect_aux-coors1(aj,:))
               val_aux=sqrt(dot_product(vect,vect))
               matrix(i,j)=val_aux
               matrix(j,i)=val_aux
            end if
         end do
      end do
   END IF
END IF

DEALLOCATE(vect,vect_aux)
DEALLOCATE(llist1,llist2)

END SUBROUTINE DISTANCE


SUBROUTINE neighbs_ranking (diff_system,pbc_opt,limit,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,neighb_list) !before: neighb_dist,neighb_uvect

  IMPLICIT NONE

  INTEGER,INTENT(IN)::diff_system,pbc_opt,ortho1
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
  CALL distance (diff_system,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,dist_matrix)

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
        
  DEALLOCATE(dist_matrix,dist_aux,filter,neight_aux)

END SUBROUTINE neighbs_ranking


SUBROUTINE neighbs_dist (diff_system,pbc_opt,limit,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,contact_map,num_neighbs,dist_matrix) !before: neighb_dist,neighb_uvect

  IMPLICIT NONE

  INTEGER,INTENT(IN)::diff_system,pbc_opt,ortho1
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
  CALL DISTANCE (diff_system,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,dist_matrix)

  contact_map=0
  num_neighbs=0

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

END SUBROUTINE neighbs_dist

SUBROUTINE translate_list (sort,list,filter,distances,dim_out,n_list,trans_inds)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::n_list,dim_out,sort
  INTEGER,DIMENSION(n_list),INTENT(IN)::list,filter
  DOUBLE PRECISION,DIMENSION(n_list),INTENT(IN)::distances
  INTEGER,DIMENSION(dim_out),INTENT(OUT)::trans_inds

  LOGICAL,DIMENSION(:),ALLOCATABLE::ifilter
  INTEGER::ii,gg

  IF (sort) THEN

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
        IF (filter(ii)) THEN
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
DOUBLE PRECISION::sk_param,roh_param,roo_param,cos_angooh_param


!!Output
INTEGER,DIMENSION(:,:),ALLOCATABLE::hbonds_out

CONTAINS

SUBROUTINE FREE_MEMORY ()

  IF (ALLOCATED(hbonds_out)) DEALLOCATE(hbonds_out)

END SUBROUTINE FREE_MEMORY



END MODULE HBONDS
