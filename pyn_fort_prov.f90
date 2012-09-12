MODULE AUX

INTEGER,DIMENSION(:),ALLOCATABLE::T_ind,T_tau,T_start
INTEGER,DIMENSION(:,:),ALLOCATABLE:: labels

TYPE array_pointer
INTEGER,DIMENSION(:),POINTER::p1
END TYPE array_pointer

TYPE darray_pointer
DOUBLE PRECISION,DIMENSION(:),POINTER::d1
END TYPE darray_pointer

DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::histo_x,histo_y

CONTAINS

  SUBROUTINE free_mem ()
    IF (ALLOCATED(histo_x)) DEALLOCATE(histo_x)
    IF (ALLOCATED(histo_y)) DEALLOCATE(histo_y)
  END SUBROUTINE free_mem

SUBROUTINE traj2net(len_str,traj_full,ranges,num_parts,num_frames,dimensions,N_nodes,Ktot)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::len_str,num_parts,num_frames,dimensions
  INTEGER,DIMENSION(num_parts,num_frames,dimensions),INTENT(IN)::traj_full
  INTEGER,DIMENSION(dimensions,2),INTENT(IN)::ranges
  INTEGER,INTENT(OUT)::N_nodes,Ktot

  INTEGER::index
  INTEGER::ii,jj,gg,kk,ll,hh,hhh,aa,bb
  INTEGER::llevamos,cajon,contador
  INTEGER::x,x2
  INTEGER,DIMENSION(:,:),ALLOCATABLE::tray
  INTEGER,DIMENSION(:,:),ALLOCATABLE::indice
  LOGICAL,DIMENSION(:,:),ALLOCATABLE::ocupado
  INTEGER,DIMENSION(:),ALLOCATABLE::deantes,aux
  CHARACTER*40::f1,f2,f3

  INTEGER,DIMENSION(:),ALLOCATABLE::SL,W,K_out
  TYPE(array_pointer),DIMENSION(:),POINTER::C,WK_out
  INTEGER::Kmax
  INTEGER,DIMENSION(:),POINTER::aux_puntero,aux_puntero2
  LOGICAL::switch


  ALLOCATE(tray(num_parts,num_frames))

  IF ((dimensions>999999).or.(len_str>999999)) THEN
     print*, 'ERROR in fortran traj2net'
     stop
  END IF

  DO cajon=1,dimensions
     
     IF (cajon==1) THEN
        llevamos=1
        tray=1
     ELSE
        CALL SYSTEM ('mv trad_aux.aux trad_aux_old.aux')
     END IF
     
     ALLOCATE(ocupado(llevamos,ranges(cajon,1):ranges(cajon,2)),indice(llevamos,ranges(cajon,1):ranges(cajon,2)))
     ocupado=.false.

     DO ii=1,num_parts
        DO jj=1,num_frames
           bb=tray(ii,jj)
           aa=traj_full(ii,jj,cajon)
           ocupado(bb,aa)=.true.
        END DO
     END DO

     contador=0

     WRITE(f1,'(I6)') cajon
     WRITE(f2,'(I6)') len_str
     f3="(I,"//TRIM(ADJUSTL(f1))//"I"//TRIM(ADJUSTL(f2))//")"

     OPEN(21,FILE="trad_aux.aux",status="REPLACE",ACTION="WRITE")
     IF (cajon==1) THEN
        DO ii=1,llevamos
           DO jj=ranges(1,1),ranges(1,2)
              IF (ocupado(ii,jj)==.true.) THEN
                 contador=contador+1
                 indice(ii,jj)=contador
                 WRITE(21,f3) contador,jj
              END IF
           END DO
        END DO
     ELSE
        OPEN(61,FILE="trad_aux_old.aux",status="OLD",ACTION="READ")
        ALLOCATE(deantes(cajon-1))         
        DO ii=1,llevamos
           READ(61,*) aa,deantes(:)
           DO jj=ranges(cajon,1),ranges(cajon,2)
              IF (ocupado(ii,jj)==.true.) THEN
                 contador=contador+1
                 WRITE(21,f3) contador,deantes(:),jj
                 indice(ii,jj)=contador
              END IF
           END DO
        END DO
        DEALLOCATE(deantes)
        CLOSE(61)
     END IF
     
     CLOSE(21)


     DEALLOCATE(ocupado)

     DO ii=1,num_parts
        DO jj=1,num_frames
           bb=tray(ii,jj)
           aa=traj_full(ii,jj,cajon)
           tray(ii,jj)=indice(bb,aa)
        END DO
     END DO

     DEALLOCATE(indice)
     llevamos=contador
   
     IF (cajon/=1) CALL SYSTEM('rm trad_aux_old.aux')
     
     print*,'>>',cajon,llevamos
  END DO


  N_nodes=llevamos
  
  print*,'DENTRO'
  ALLOCATE(C(N_nodes),K_out(N_nodes))
  ALLOCATE(WK_out(N_nodes),SL(N_nodes),W(N_nodes))
  ALLOCATE(aux_puntero(25000),aux_puntero2(25000))

  DO ii=1,N_nodes
     ALLOCATE(C(ii)%p1(1),WK_out(ii)%p1(1))
     C(ii)%p1(:)=0
     WK_out(ii)%p1(:)=0
  END DO

  SL=0
  W=0
  K_out=0
  kk=0
  x=0
  x2=0
  aux_puntero(:)=0
  aux_puntero2(:)=0
  contador=0

  DO index=1,num_parts

     x=tray(index,1)
     contador=contador+1
     x2=tray(index,2)
     contador=contador+1

     DO ii=3,num_frames

        W(x)=W(x)+1
      
        switch=.true.
        DO gg=1,K_out(x)
           IF(C(x)%p1(gg)==x2) THEN
              WK_out(x)%p1(gg)=WK_out(x)%p1(gg)+1
              switch=.false.
              EXIT
           END IF
        END DO

        IF (switch==.true.) THEN
           IF (K_out(x)>0) THEN
              gg=K_out(x)
              IF (gg>25000) THEN
                 print*,'tenemos un problema'
                 stop
              END IF
              aux_puntero(1:gg)=C(x)%p1(:)
              aux_puntero2(1:gg)=WK_out(x)%p1(:)
              DEALLOCATE(C(x)%p1,WK_out(x)%p1)
              ALLOCATE(C(x)%p1(gg+1),WK_out(x)%p1(gg+1))
              C(x)%p1(1:gg)=aux_puntero(:)
              WK_out(x)%p1(1:gg)=aux_puntero2(:)
              C(x)%p1(gg+1)=x2
              WK_out(x)%p1(gg+1)=1
              K_out(x)=K_out(x)+1
           ELSE
              WK_out(x)%p1(1)=1
              C(x)%p1(1)=x2
              K_out(x)=1
           END IF
        END IF

        x=x2  !!bin
        
        x2=tray(index,ii)
        contador=contador+1
        
     END DO

     W(x)=W(x)+1
   
     switch=.true.
   
     DO gg=1,K_out(x)
        IF(C(x)%p1(gg)==x2) THEN
           WK_out(x)%p1(gg)=WK_out(x)%p1(gg)+1
           switch=.false.
           EXIT
        END IF
     END DO
   
     IF (switch==.true.) THEN
        IF (K_out(x)>0) THEN
           gg=K_out(x)
           IF (gg>25000) THEN
              print*,'tenemos un problema'
              stop
           END IF
           aux_puntero(1:gg)=C(x)%p1(:)
           aux_puntero2(1:gg)=WK_out(x)%p1(:)
           DEALLOCATE(C(x)%p1,WK_out(x)%p1)
           ALLOCATE(C(x)%p1(gg+1),WK_out(x)%p1(gg+1))
           C(x)%p1(1:gg)=aux_puntero(:)
           WK_out(x)%p1(1:gg)=aux_puntero2(:)
           C(x)%p1(gg+1)=x2
           WK_out(x)%p1(gg+1)=1
           K_out(x)=K_out(x)+1
        ELSE
           WK_out(x)%p1(1)=1
           C(x)%p1(1)=x2
           K_out(x)=1
        END IF
     END IF
 
     x=x2
     contador=contador+1

     W(x)=W(x)+1

  END DO

  Kmax=maxval(K_out(:))
  Ktot=sum(K_out(:))

  ALLOCATE(T_start(N_nodes+1),T_ind(Ktot),T_tau(Ktot))
  T_start=0
  T_ind=0
  T_tau=0

  ll=0
  DO ii=1,N_nodes

     hh=0
     switch=.false.
     DO gg=1,K_out(ii)
        IF (C(ii)%p1(gg)==ii) THEN
           switch=.true.
           exit
        END IF
     END DO

     T_start(ii)=ll
     IF (switch==.true.) THEN
        ll=ll+1
        T_ind(ll)=ii
        hh=WK_out(ii)%p1(gg)
        T_tau(ll)=hh
        hh=hh+hhh
        DO jj=1,K_out(ii)
           IF (jj/=gg) THEN
              ll=ll+1
              T_ind(ll)=C(ii)%p1(jj)
              hhh=WK_out(ii)%p1(jj)
              T_tau(ll)=hhh
              hh=hh+hhh
           END IF
        END DO
     ELSE
        DO jj=1,K_out(ii)
           ll=ll+1
           T_ind(ll)=C(ii)%p1(jj)
           hhh=WK_out(ii)%p1(jj)
           hh=hh+hhh
           T_tau(ll)=hhh
        END DO
     END IF

  END DO

  T_start(N_nodes+1)=ll

  DEALLOCATE(tray)

end subroutine traj2net


subroutine rao_stat_1(tw,traj_dists,limits,frames,len_lims,traj_bins)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::frames,len_lims,tw
  DOUBLE PRECISION,DIMENSION(frames),INTENT(IN)::traj_dists
  DOUBLE PRECISION,DIMENSION(len_lims),INTENT(IN)::limits
  INTEGER,DIMENSION(frames,len_lims+1),INTENT(OUT)::traj_bins

  INTEGER::ii,jj,kk,gg

  traj_bins=0
     
  DO jj=1,frames
     kk=len_lims+1
     DO gg=1,len_lims
        IF (traj_dists(jj)<limits(gg)) THEN
           kk=gg
           exit
        END IF
     END DO
     DO gg=jj-tw,jj+tw
        IF ((0<gg).and.(gg<=frames)) THEN
           traj_bins(gg,kk)=traj_bins(gg,kk)+1
        END IF
     END DO
  END DO

end subroutine rao_stat_1


subroutine ganna (opt_range,opt,ibins,immn,immx,idelta_x,traj,ksi,tw,len_traj,traj_out)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::opt_prec,opt_range,opt,bins,tw
  DOUBLE PRECISION,INTENT(IN)::delta_x,mmn,mmx,prec
  DOUBLE PRECISION,INTENT(IN)::ksi
  DOUBLE PRECISION,DIMENSION(len_traj),INTENT(IN)::traj
  DOUBLE PRECISION,DIMENSION(len_traj-2*tw),INTENT(OUT)::traj_out





  INTEGER::ii,jj,gg,represent,Ltw
  INTEGER::num_rep
  INTEGER,DIMENSION(:),ALLOCATABLE::list,aux_list
  DOUBLE PRECISION::dsm,dist
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::cumul_test,segmento
  LOGICAL::switch
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter

  !!! For histogram:

  bins=ibins
  max=imax
  min=imin
  delta_x=idelta_x

  IF (opt_range==0) THEN
     IF (opt==1) THEN
        bins=CEILING((max-min)/delta_x)
        sobra=(bins*delta_x-(max-min))/2.0d0
        bins=bins+1
        min=min-sobra
        max=max+sobra
     ELSE
        delta_x=(max-min)/(bins*1.0d0)
        sobra=delta_x/2.0d0
        min=min-sobra
        max=max+sobra
        bins=bins+1
     END IF
  ELSE
     IF (opt==1) THEN
        bins=CEILING((max-min)/delta_x)
     ELSE
        delta_x=(max-min)/(bins*1.0d0)
     END IF
  END IF

  ALLOCATE(cumul(bins))
  cumul=0.0d0

  !!!


  Ltw=(2*tw+1)
  traj_out=0
  dsm=sqrt(2.0d0/(1.0d0*Ltw))*ksi !Kolmogorov-Smirnov

  ALLOCATE(cumul_test(num_bins),segmento(Ltw))

  num_rep=1
  ALLOCATE(list(num_rep),filter(num_rep))
  list(1)=num_rep
  traj_out(1)=num_rep

  DO ii=1,len_traj-2*tw
     kk=ii+tw
     segmento(:)=traj(ii:(ii+Ltw))
     CALL histograma(1,opt_prec,opt_range,opt,1,segmento,bins,mmn,mmx,delta_x,prec,Ltw)
     

  DO ii=2,len_traj
     cumul_test=cumul_traj(ii,:)
     switch=.true.
     filter=.true.
     DO jj=ii-1,1,-1
        gg=traj_out(jj)
        represent=list(gg)
        IF (filter(gg)==.true.) THEN
           dist=MAXVAL(abs(cumul_traj(represent,:)-cumul_test(:)),DIM=1)
           IF (dist<=dsm) THEN
              traj_out(ii)=gg
              switch=.false.
              EXIT
           END IF
           filter(gg)=.false.
           IF (COUNT(filter,DIM=1)==0) THEN
              EXIT
           END IF
        END IF
     END DO
     IF (switch==.true.) THEN
        ALLOCATE(aux_list(num_rep))
        aux_list=list
        DEALLOCATE(list,filter)
        num_rep=num_rep+1
        ALLOCATE(list(num_rep),filter(num_rep))
        list(1:(num_rep-1))=aux_list
        DEALLOCATE(aux_list)
        list(num_rep)=ii
        traj_out(ii)=num_rep
     END IF
  END DO

END subroutine ganna

  SUBROUTINE histograma (opt_norm,opt_prec,opt_range,opt,opt_cumul,idatos,ibins,imin_x,imax_x,idelta_x,iprec,l)
    
    implicit none
    INTEGER,INTENT(IN)::opt_norm,opt_prec,opt_range,opt,l,ibins,opt_cumul
    DOUBLE PRECISION,DIMENSION(l),INTENT(IN)::idatos
    DOUBLE PRECISION,INTENT(IN)::idelta_x,imax_x,imin_x,iprec

    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::frecuencias
    INTEGER::i,j,k,prec,aux,bins,tt
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::datos
    DOUBLE PRECISION::max,min,delta_x,total,sobra
    
    bins=ibins

    ALLOCATE(datos(l))
    datos=0.0d0

    !! Por un problema que hay con python 2.6 corrijo la precision
    IF (opt_prec==1) THEN
       prec=nint(1.0/iprec)
       datos=0.0d0
       DO i=1,l
          aux=nint(idatos(i)*prec)
          datos(i)=(aux*1.0d0)/(prec*1.0d0)
       END DO
       max=(nint(imax_x*prec)*1.0d0)/(prec*1.0d0)
       min=(nint(imin_x*prec)*1.0d0)/(prec*1.0d0)
       delta_x=idelta_x
    ELSE
       datos=idatos
       max=imax_x
       min=imin_x
       delta_x=idelta_x
    END IF
    
    IF (opt_range==0) THEN
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
          sobra=(bins*delta_x-(max-min))/2.0d0
          bins=bins+1
          min=min-sobra
          max=max+sobra
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
          sobra=delta_x/2.0d0
          min=min-sobra
          max=max+sobra
          bins=bins+1
       END IF
    ELSE
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
       END IF
    END IF
    
    ALLOCATE(frecuencias(bins))
    frecuencias=0.0d0
    
    DO k=1,l
       tt=FLOOR((datos(k)-min)/delta_x) 
       frecuencias(tt)=frecuencias(tt)+1.0d0
    END DO
    

    IF (opt_norm==1) THEN
       total=SUM(frecuencias)*delta_x
       frecuencias=frecuencias/total
    END IF
    
    !! Output

    CALL free_mem()
    ALLOCATE(histo_y(bins),histo_x(bins))
    histo_y=frecuencias
    IF (opt_cumul==1) THEN
       DO i=1,bins-1
          histo_y(i+1)=histo_y(i+1)+histo_y(i)
       END DO
    END IF

    DO i=1,bins
       histo_x(i)=min+(i*1.0d0-0.50d0)*delta_x
    END DO
    
    
    DEALLOCATE(datos,frecuencias)

  END SUBROUTINE histograma
  


END MODULE AUX
