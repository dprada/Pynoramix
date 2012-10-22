MODULE stats

    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::histo_x,histo_y
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::histo_z

  CONTAINS

  SUBROUTINE free_mem ()
    IF (ALLOCATED(histo_x)) DEALLOCATE(histo_x)
    IF (ALLOCATED(histo_y)) DEALLOCATE(histo_y)
    IF (ALLOCATED(histo_z)) DEALLOCATE(histo_z)
  END SUBROUTINE free_mem

  SUBROUTINE average (idatos,l,aa,sigma)
    
    implicit none
    INTEGER,INTENT(IN)::l
    DOUBLE PRECISION,DIMENSION(l),INTENT(IN)::idatos
    DOUBLE PRECISION,INTENT(OUT)::aa,sigma

    INTEGER::i,j
    DOUBLE PRECISION::aux_ave,aux_ave2

    aux_ave=0.0d0
    aux_ave2=0.0d0

    DO i=1,l
       aux_ave=aux_ave+idatos(i)
       aux_ave2=aux_ave2+idatos(i)**2
    END DO

    aux_ave=aux_ave/(l*1.0d0)
    aux_ave2=aux_ave2/(l*1.0d0)

    aa=aux_ave
    sigma=sqrt(aux_ave2-aux_ave**2)

  END SUBROUTINE average
  


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
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
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
       tt=CEILING((datos(k)-min)/delta_x) 
       if (tt==0) tt=1
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
       histo_y=histo_y*delta_x
       DO i=1,bins-1
          histo_y(i+1)=histo_y(i+1)+histo_y(i)
       END DO
    END IF

    DO i=1,bins
       histo_x(i)=min+(i*1.0d0-0.50d0)*delta_x
    END DO
    
    
    DEALLOCATE(datos,frecuencias)

  END SUBROUTINE histograma
  
  SUBROUTINE histograma_2d (opt_norm,opt_prec,opt_range,opt,idatos,ibins,imin_x,imax_x,idelta_x,iprec,l)
    
    implicit none
    INTEGER,INTENT(IN)::opt_norm,opt_prec,opt_range,opt,l
    INTEGER,DIMENSION(2),INTENT(IN)::ibins
    DOUBLE PRECISION,DIMENSION(l,2),INTENT(IN)::idatos
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN)::idelta_x,imax_x,imin_x
    DOUBLE PRECISION,INTENT(IN)::iprec

    INTEGER::tt,qq
    DOUBLE PRECISION,DIMENSION(2)::max,min
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::frecuencias
    INTEGER::i,j,k,prec,aux
    INTEGER,DIMENSION(2)::bins
    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::datos
    DOUBLE PRECISION::omax(2),omin(2),delta_x(2)
    DOUBLE PRECISION::total,sobra
    
    bins=ibins

    ALLOCATE(datos(l,2))
    datos=0.0d0

    !! Por un problema que hay con python 2.6 corrijo la precision
    IF (opt_prec==1) THEN
       prec=nint(1.0/iprec)
       datos=0.0d0
       DO i=1,l
          aux=nint(idatos(i,1)*prec)
          datos(i,1)=(aux*1.0d0)/(prec*1.0d0)
          aux=nint(idatos(i,2)*prec)
          datos(i,2)=(aux*1.0d0)/(prec*1.0d0)
       END DO
       omax=(nint(imax_x*prec)*1.0d0)/(prec*1.0d0)
       omin=(nint(imin_x*prec)*1.0d0)/(prec*1.0d0)
       delta_x=idelta_x
    ELSE
       datos=idatos
       omax=imax_x
       omin=imin_x
       delta_x=idelta_x
    END IF
    
    IF (opt_range==0) THEN
       IF (opt==1) THEN
          DO i=1,2
             bins(i)=CEILING((omax(i)-omin(i))/delta_x(i))
          END DO
       ELSE
          DO i=1,2
             delta_x(i)=(omax(i)-omin(i))/(bins(i)*1.0d0)
          END DO
       END IF
    ELSE
       IF (opt==1) THEN
          DO i=1,2
             bins(i)=CEILING((omax(i)-omin(i))/delta_x(i))
          END DO
       ELSE
          DO i=1,2
             delta_x(i)=(omax(i)-omin(i))/(bins(i)*1.0d0)
          END DO
       END IF
    END IF
    
    ALLOCATE(frecuencias(bins(1),bins(2)))
    frecuencias=0.0d0
    
    DO k=1,l
       tt=FLOOR((datos(k,1)-omin(1))/delta_x(1)) 
       qq=FLOOR((datos(k,2)-omin(2))/delta_x(2))
       IF (tt==0) tt=1
       IF (qq==0) qq=1
       frecuencias(tt,qq)=frecuencias(tt,qq)+1
    END DO
    
    IF (opt_norm==1) THEN
       frecuencias=frecuencias/(l*delta_x(1)*delta_x(2))
    END IF
        
    !! Output
    CALL free_mem()
    ALLOCATE(histo_x(bins(1)),histo_y(bins(2)),histo_z(bins(1),bins(2)))
    histo_z=frecuencias
    DO i=1,bins(1)
       histo_x(i)=omin(1)+(i*1.0d0-0.50d0)*delta_x(1)
    END DO
    DO i=1,bins(2)
       histo_y(i)=omin(2)+(i*1.0d0-0.50d0)*delta_x(2)
    END DO
    
    DEALLOCATE(datos,frecuencias)

  END SUBROUTINE histograma_2d


  SUBROUTINE binning (opt_range,opt_delta_x,idatos,ibins,imin_x,imax_x,idelta_x,l,tray_bins)
    
    implicit none
    INTEGER,INTENT(IN)::opt_range,opt_delta_x,l,ibins
    DOUBLE PRECISION,DIMENSION(l),INTENT(IN)::idatos
    DOUBLE PRECISION,INTENT(IN)::idelta_x,imax_x,imin_x
    INTEGER,DIMENSION(l),INTENT(OUT)::tray_bins

    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::frecuencias
    INTEGER::i,j,k,aux,bins,tt
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::datos
    DOUBLE PRECISION::max,min,delta_x,total,sobra
    
    tray_bins=0
    bins=ibins


    ALLOCATE(datos(l))
    datos=0.0d0

    IF (opt_range==0) THEN
       max=MAXVAL(idatos(:),DIM=1)
       min=MINVAL(idatos(:),DIM=1)
       IF (opt_delta_x==1) THEN
          bins=CEILING((max-min)/delta_x)
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
       END IF
    ELSE
       IF (opt_delta_x==1) THEN
          bins=CEILING((max-min)/delta_x)
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
       END IF
    END IF
    
    DO k=1,l
       tt=CEILING((idatos(k)-min)/delta_x) 
       IF (tt==0) tt=1
       tray_bins(k)=tt-1
    END DO
    
    
    !! Output
    CALL free_mem()
    ALLOCATE(histo_x(bins))

    DO i=1,bins
       histo_x(i)=min+(i*1.0d0-0.50d0)*delta_x
    END DO
    
    
    DEALLOCATE(datos)

  END SUBROUTINE binning
  

  SUBROUTINE binning_x (opt_range,opt,ibins,imin_x,imax_x,idelta_x,delta_x)
    
    implicit none
    INTEGER,INTENT(IN)::opt_range,opt,ibins
    DOUBLE PRECISION,INTENT(IN)::idelta_x,imax_x,imin_x

    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::frecuencias
    INTEGER::i,j,k,aux,bins
    DOUBLE PRECISION::max,min,total,sobra
    DOUBLE PRECISION,INTENT(OUT)::delta_x
    
    bins=ibins

    !! Por un problema que hay con python 2.6 corrijo la precision

    max=imax_x
    min=imin_x
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
    
    !! Output
    CALL free_mem()
    ALLOCATE(histo_x(bins))

    DO i=1,bins
       histo_x(i)=min+(i*1.0d0-0.50d0)*delta_x
    END DO
    
  END SUBROUTINE binning_x

END MODULE stats

! f2py -c -m pyn_fort_math pyn_fort_math.f90

