MODULE GLOB


CONTAINS

SUBROUTINE ind_wat_limit_4_nosim (aux,hbs,hbdists,num_hbs,num_wats,num_atoms,mss)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_wats,num_hbs,num_atoms
  INTEGER,DIMENSION(num_atoms,2),INTENT(IN)::aux
  INTEGER,DIMENSION(num_hbs,3),INTENT(IN)::hbs
  DOUBLE PRECISION,DIMENSION(num_hbs),INTENT(IN)::hbdists
  INTEGER,DIMENSION(num_wats,17),INTENT(OUT)::mss

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dist_first_shell
  INTEGER,DIMENSION(:,:),ALLOCATABLE::first_shell
  INTEGER,DIMENSION(:),ALLOCATABLE::bonds_o
  INTEGER,DIMENSION(:),ALLOCATABLE::microstate

  INTEGER::ii,jj,kk,bb,gg
  INTEGER::ind_oh,hi,ind_o

  mss=0
  ALLOCATE(dist_first_shell(num_wats,4),first_shell(num_wats,4),bonds_o(num_wats))
  ALLOCATE(microstate(17))

  first_shell=0
  dist_first_shell=0.0d0
  bonds_o=0

  DO ii=1,num_hbs
     ind_oh=aux(hbs(ii,1)+1,1)+1
     hi=aux(hbs(ii,2)+1,2)
     ind_o=aux(hbs(ii,3)+1,1)+1
     IF (first_shell(ind_oh,hi)==0) THEN
        first_shell(ind_oh,hi)=ind_o
        dist_first_shell(ind_oh,hi)=hbdists(ii)
     ELSE
        IF (dist_first_shell(ind_oh,hi)>hbdists(ii)) THEN
           first_shell(ind_oh,hi)=ind_o
           dist_first_shell(ind_oh,hi)=hbdists(ii)
        END IF
     END IF
     IF (bonds_o(ind_o)==0) THEN
        first_shell(ind_o,3)=ind_oh
        dist_first_shell(ind_oh,3)=hbdists(ii)
        bonds_o(ind_o)=1
     ELSE
        IF (bonds_o(ind_o)==1) THEN
           IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
              first_shell(ind_o,4)=ind_oh
              dist_first_shell(ind_o,4)=hbdists(ii)
           ELSE
              first_shell(ind_o,4)=first_shell(ind_o,3)
              dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
              first_shell(ind_o,3)=ind_oh
              dist_first_shell(ind_o,3)=hbdists(ii)
           END IF
           bonds_o(ind_o)=2
        ELSE
           IF (bonds_o(ind_o)==2) THEN
              IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                 first_shell(ind_o,4)=first_shell(ind_o,3)
                 dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                 first_shell(ind_o,3)=ind_oh
                 dist_first_shell(ind_o,3)=hbdists(ii)
              ELSE
                 IF (dist_first_shell(ind_o,4)>hbdists(ii)) THEN
                    first_shell(ind_o,4)=ind_oh
                    dist_first_shell(ind_o,4)=hbdists(ii)
                 END IF
              END IF
           END IF
        END IF
     END IF
  END DO

  DO kk=1,num_wats

     microstate=0
     microstate(1)=kk
     
     microstate(2:5)=(/ first_shell(kk,:) /)

     gg=6
     DO ii=2,3
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg:gg+1)=(/ first_shell(bb,1:2) /)
           IF ((first_shell(bb,3)==kk).or.(first_shell(bb,4)==kk)) THEN
              IF (first_shell(bb,3)==kk) THEN
                 microstate(gg+2)=first_shell(bb,4)
              ELSE
                 microstate(gg+2)=first_shell(bb,3)
              END IF
              IF ((first_shell(bb,3)==kk).and.(first_shell(bb,4)==kk)) THEN
                 !print*,'problema ssn1' !!! Comprobar esto
                 microstate(gg+2)=first_shell(bb,3)
              END IF
           ELSE
              microstate(gg+2)=first_shell(bb,3)
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO

     DO ii=4,5
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg+1:gg+2)=(/ first_shell(bb,3:4) /)
           IF ((first_shell(bb,1)==kk).or.(first_shell(bb,2)==kk)) THEN
              IF (first_shell(bb,1)==kk) THEN
                 microstate(gg)=first_shell(bb,2)
              ELSE
                 microstate(gg)=first_shell(bb,1)
              END IF
              IF ((first_shell(bb,1)==kk).and.(first_shell(bb,2)==kk)) THEN
                 ! print*,'pproblema ssn2' !!! Comprobar esto
                 ! print*,a,'|',first_shell(b,:)
                 microstate(gg)=first_shell(bb,1)
                 ! print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
                 ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)
              END IF
           ELSE
              microstate(gg)=first_shell(bb,1) !!!! esto esta mal, hay que elegir entre el 1 y el 2
              !print*,'siiiiiiiiiiiiiiiiii'
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO
     
     !print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
     ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)

     !Quito indices de mols

     !!filtro=.true.
     !!aux=microstate
     !!mss_ind_wat=microstate
     !!DO i=1,17
     !!   IF (aux(i)<=0) filtro(i)=.false.
     !!END DO
     !! 
     !!shell_w(a,:)=0
     !!g=0
     !! 
     !!DO j=1,17
     !!   IF (filtro(j)==.true.) THEN
     !!      b=aux(j)
     !!      g=g+1
     !!      shell_w(a,g)=b
     !!      DO i=1,17
     !!         IF (filtro(i)==.true.) THEN
     !!            IF (aux(i)==b) THEN
     !!               microstate(i)=j
     !!               filtro(i)=.false.
     !!            END IF
     !!         END IF
     !!      END DO
     !!   END IF
     !!END DO

     mss(kk,:)=microstate

  END DO


  DEALLOCATE(dist_first_shell,first_shell,bonds_o)
  DEALLOCATE(microstate)


END SUBROUTINE ind_wat_limit_4_nosim

END MODULE GLOB
