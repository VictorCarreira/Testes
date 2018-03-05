PROGRAM teste_dists

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    !Programa de teste e comparação de duas diferentes distâncias                !
    !Orientador: Cosme Ferreira da Ponte Neto                                    !
    !Aluno: Victor Ribeiro Carreira                                              !
    !Categoria: classificador                                                    !
    !Subrotina Teste                                                             !
    !Para usar compilação com flags utilize:                                     !
    !gfortran -fbounds-check -fbacktrace -Wall -Wextra -pedantic                 !
    !"pasta/subpasta/nomedopragrama.f95" -o nomedoexecutável                     !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

                          !***********TABELA DE VARIÁVEIS***********!
                          !lito1: nuvem de pontos de treinamento    !
                          !lito2: o que vai ser classificado        !
                          !       (ponto ou nuvem)                  !
                          !np: número de pontos da/das nuvem/nuvens !
                          !r: raio do circulo                       !
                          !a: eixo menor da elipse                  !
                          !b: eixo maior da elipse                  !
                          !theta: angulo entre os eixos             !
                          !*min: valor minimo de uma variavel       !
                          !*max: valor maximo de uma variavel       !
                          !-----------------------------------------!

  IMPLICIT NONE

  INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
  INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

  INTEGER(KIND=SP):: i, ic1, ic2, ic3
  INTEGER(KIND=SP)::np1, np2, np3
  INTEGER(KIND=SP), PARAMETER:: semente=17

  REAL(KIND=DP)::inicio,final,x, px, py, fx
  REAL(KIND=DP)::mahaI, eucliI, mahaII, eucliII, L, LXmax, LXmin, LYmax, LYmin
  REAL(KIND=DP), PARAMETER::pi=3.141592653, rmax=3.0, rmin=0.0, thetamin=0.0, &
  thetamax=5.0*pi, a=15.0, b=3.0

  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:)::lito1,lito2,lito3,lito11,lito22,lito33
  
  CALL CPU_TIME(inicio)

  np1=500
  np2=500
  np3=500

  ALLOCATE(lito1(np1,2), lito2(np2,2), lito3(np3,2))

  lito1=0d0  !Zerando a matriz
  lito2=0d0  !Zerando a matriz
  lito3=0d0  !Zerando a matriz


   OPEN(1,FILE='teste_nuvem_T.txt')
   OPEN(2,FILE='teste_nuvem_D.txt')
   OPEN(3,FILE='teste_nuvem_C.txt')
   OPEN(4,FILE='AnaliseIII.txt')
   

   !Criando uma semente para o gerador de números aleatórios
   DO i=1,semente
    x=RAND()
   END DO  


  !Corrigindo o vício do sorteio aleatório
   L=2.0*rmax ! lado de um quadrado (válido somente para circunferências)
  !O zero está no centro da elipse
   LXmax=a !Lado máximo de um retângulo (válido somente para uma elipse)
   LXmin=-a !Lado mínimo de um retângulo (válido somente para uma elipse) 
   LYmax=b !Lado máximo de um retângulo (válido somente para uma elipse)
   LYmin=-b !Lado mínimo de um retângulo (válido somente para uma elipse) 


  ! Escrevendo os lados do quadrado 1 sorteio
  !DO i=1,np1
  ! lito1(i,1)=(L-rmin)*RAND()+rmin-rmax 
  ! lito1(i,2)=(L-rmin)*RAND()+rmin-rmax
  !END DO


  !Escrevendo os lados do retângulo 1 e sorteando os dados
  DO i=1,np1
   lito1(i,1)=(LXmax-LXmin)*RAND()+LXmin 
   lito1(i,2)=(LYmax-LYmin)*RAND()+LYmin
  END DO

!Criando um contador para a circunferência/Elipse 1
  ic1=0
  DO i=1,np1
   !IF (DSQRT(lito1(i,1)**2+lito1(i,2)**2) <= L/2.0) THEN ! válido para a circunferências. 
    px=lito1(i,1)
    py=lito1(i,2)
    fx=b*DSQRT(1.0-px**2.0/a**2.0)
    IF (py <= fx .AND. py >= -fx) THEN  !Válido para a elipse
     ic1=ic1+1!Contar quantos foram considerados
   END IF
  END DO 

!Escrevendo os lados do quadrado 2
DO i=1,np2
  lito2(i,1)=(L-rmin)*RAND()+rmin-rmax 
  lito2(i,2)=(L-rmin)*RAND()+rmin-rmax
END DO
!Criando um contador para a circunferência 2 
 ic2=0
 DO i=1,np2
  IF (DSQRT(lito2(i,1)**2+lito2(i,2)**2) <= L/2.0) THEN
   ic2=ic2+1!Contar quantos foram considerados
  END IF
 END DO 



!Escrevendo os lados do quadrado 3
DO i=1,np3
 lito3(i,1)=(L-rmin)*RAND()+rmin-rmax 
 lito3(i,2)=(L-rmin)*RAND()+rmin-rmax
END DO
!Criando um contador para a circunferência 3
 ic3=0
 DO i=1,np3
  IF (DSQRT(lito3(i,1)**2+lito3(i,2)**2) <= L/2.0) THEN
   ic3=ic3+1!Contar quantos foram considerados
  END IF
 END DO 

  

 ALLOCATE(lito11(ic1,2), lito22(ic2,2), lito33(ic3,2))
 lito11=0d0
 lito22=0d0
 lito33=0d0

 !Escrevendo os dados da nuvem 1 em uma elipse ou circunferência (basta descomentar um IF e comentar o outro)
  ic1=0
   DO i=1,np1
   !IF (DSQRT(lito1(i,1)**2+lito1(i,2)**2) <= L/2.0) THEN ! Circunferência
   px=lito1(i,1)
   py=lito1(i,2)
   fx=b*DSQRT(1.0-px**2.0/a**2.0)
   IF (py <= fx .AND. py >= -fx) THEN !Elipse
     ic1=ic1+1!Contar quantos foram considerados
     lito11(ic1,1)=lito1(i,1)+10.0! desloca o centro em 1 unidades
     lito11(ic1,2)=lito1(i,2)+15.0! descloca o centro em 15 unidades
     WRITE(1,FMT=*)lito11(ic1,1),lito11(ic1,2)
    END IF
   END DO 

!Escrevendo os dados da nuvem 2 em uma circunferência
 ic2=0
 DO i=1,np2
  IF (DSQRT(lito2(i,1)**2+lito2(i,2)**2) <= L/2.0) THEN
   ic2=ic2+1!Contar quantos foram considerados
   lito22(ic2,1)=lito2(i,1)+10.0! desloca o centro em 10 unidades
   lito22(ic2,2)=lito2(i,2)-5.0! descloca o centro em 5 unidades
   WRITE(3,FMT=*)lito22(ic2,1),lito22(ic2,2)
  END IF
 END DO 

!Escrevendo os dados da nuvem 3 em uma circunferência
 ic3=0
 DO i=1,np3
  IF (DSQRT(lito3(i,1)**2+lito3(i,2)**2) <= L/2.0) THEN
   ic3=ic3+1!Contar quantos foram considerados
   lito33(ic3,1)=lito3(i,1)+30.0! desloca o centro em 10 unidades
   lito33(ic3,2)=lito3(i,2)+15.0! descloca o centro em 5 unidades
   WRITE(2,FMT=*)lito33(ic3,1),lito33(ic3,2)
  END IF
 END DO 


 !Redefinindo o número de pontos/dados que entrarão na rotina
  np1=ic1
  np2=ic2
  np3=ic3


  CALL euclideana(lito11,lito22,eucliI)
  CALL euclideana(lito11,lito33,eucliII)
  CALL mahalanobeana(lito11,np1,lito22,np2,2,mahaI)
  CALL mahalanobeana(lito11,np1,lito33,np3,2,mahaII)

  PRINT*,'dist de euclides (1-2)=', eucliI
  PRINT*, 'dist de maha (1-2)=', mahaI
  PRINT*,'dist de euclides (1-3)=', eucliII
  PRINT*, 'dist de maha (1-3)=', mahaII

  CALL saida

  CLOSE(1)
  CLOSE(2)
  CLOSE(3)
  CLOSE(4)

  CALL CPU_TIME(final)
  PRINT*,'tempo de máquina=',final-inicio
  !*******************************************************************************************!
  CONTAINS

  SUBROUTINE saida

    ! Formatos dos arquivos de saida
     21 FORMAT(A22,2x,E12.3)
     22 FORMAT(A20,2x,I10)
    !Gravando os parâmetros do modelo
    WRITE(4,FMT=*)'-------Parâmetros do modelo-------'
    WRITE(4,FMT=21)'rmax=',rmax
    WRITE(4,FMT=21)'rmin=',rmin
    WRITE(4,FMT=21)'thetamin=',thetamin
    WRITE(4,FMT=21)'thetamax=',thetamax
    WRITE(4,FMT=21)'a=',a
    WRITE(4,FMT=21)'b=',b
    WRITE(4,FMT=22)'semente randômica=',semente
    WRITE(4,FMT=*)'---------------Dados--------------'
    WRITE(4,FMT=22)'cluster 1=',np1
    WRITE(4,FMT=22)'cluster 2=',np2
    WRITE(4,FMT=22)'cluster 3=',np3
    WRITE(4,FMT=*)'------------Distâncias------------'
    WRITE(4,FMT=21)'Euclideana(I-II)=', eucliI
    WRITE(4,FMT=21)'Mahalanobeana(I-II)=', mahaI
    WRITE(4,FMT=21)'Euclideana(I-III)=', eucliII
    WRITE(4,FMT=21)'Mahalanobeana(I-III)=', mahaII
  END SUBROUTINE saida 

  SUBROUTINE euclideana(lito1,lito2,eucli)

   IMPLICIT NONE
   INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
   INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

   REAL(KIND=DP), DIMENSION(:,:), INTENT(IN)::lito1, lito2
   REAL(KIND=DP), INTENT(OUT):: eucli
   REAL(KIND=DP)::media1, media2

   INTEGER(KIND=SP):: k

    eucli=0d0

    IF(SIZE(lito1(1,:)) /= SIZE(lito2(1,:)))THEN
      PRINT*,'WARNING! THE PROPERTIES NUMBER´S OF lito1 AND lito2 MUST BE THE SAME.'
      STOP
      RETURN
    END IF

    DO k=1,SIZE(lito1(1,:))  ! Inicia o laço da primeira até a última propriedade que é dado pelo size de lito
     media1=0d0 !zera as variáveis
     media2=0d0 !zera as variáveis para a cálculo do centróide 2
     media1=SUM(lito1(:,k))/SIZE(lito1(:,k)) !calcula as médias para as k propriedades
     media2=SUM(lito2(:,k))/SIZE(lito2(:,k)) !calcula as médias para k propriedades para uma segunda nuvem de pontos
     !eucli= eucli + (lito2(1,k)-media1)**2 ! Cálculo da medida de semelhança de euclides para um conjunto de pontos e um centróide
     eucli= eucli + (media2-media1)**2
    END DO ! Final do laço das k propriedades
    eucli=SQRT(eucli)

  
  END SUBROUTINE euclideana

  SUBROUTINE mahalanobeana(g11,np1,g22,np2,ndim,dist)

  !  	subrotina que calcula a distância de mahalanobis entre
  !  	dois agrupamentos de elementos com dimensão ndim

     IMPLICIT NONE
      INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
      INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

      INTEGER(KIND=SP), INTENT(IN):: np1, np2
      INTEGER(KIND=SP), INTENT(IN):: ndim
      REAL(KIND=DP),INTENT(OUT):: dist

      INTEGER(KIND=SP):: i,j,k
      REAL(KIND=DP),ALLOCATABLE,DIMENSION(:)::soma, xm1, xm2
      REAL(KIND=DP),ALLOCATABLE, DIMENSION(:,:):: g1, g2, g1T, g2T, cov1, cov2, &
      covag, g11, g22, md, mdT, alfa, d2

      ALLOCATE(soma(ndim),xm1(ndim),xm2(ndim))

      ALLOCATE(g1(np1,ndim),g2(np2,ndim),g1T(ndim,np1),g2T(ndim,np2),&
      cov1(ndim,ndim),cov2(ndim,ndim),covag(ndim,ndim),md(ndim,1),&
      mdT(1,ndim),alfa(1,ndim),d2(1,1))

      g1=g11
      g2=g22

  !  	grupo 1

    DO j=1,ndim
      soma(j)=0d0
      DO i=1,np1
        soma(j)=soma(j)+g1(i,j)
      END DO
    END DO

    DO i=1,ndim
      xm1(i)=soma(i)/dfloat(np1)
    END DO

  !  	grupo 2

    DO j=1,ndim
      soma(j)=0d0
      DO i=1,np2
        soma(j)=soma(j)+g2(i,j)
      END DO
    END DO

    DO i=1,ndim
      xm2(i)=soma(i)/dfloat(np2)
    END DO

  !  	vetor das diferenças - será escrito sobre a matrizes g1 e g2


    DO j=1,ndim
      DO i=1,np1
        g1(i,j)=g1(i,j)-xm1(j)
      END DO
    END DO

    DO  j=1,ndim
      DO i=1,np2
        g2(i,j)=g2(i,j)-xm2(j)
      END DO
    END DO

  !      --------GRUPO 1 ---------------------
  !  	criando a matriz transposta g1T
  !  	-------------- -------------------
    DO i=1,np1    !107 ! número de equações
      DO j=1,ndim   !2
        g1T(j,i)=g1(i,j)
      END DO
    END DO
  !  ----------------------------------------------------
  !  	 - multiplicação de matrizes
  !  	   multiplicação de g1T por g1

    DO k=1,ndim
      DO j=1,ndim
        cov1(j,k)=0.d0
        DO i=1,np1
          cov1(j,k)=cov1(j,k)+g1T(j,i)*g1(i,k)
        END DO
      END DO
    END DO

    write(6,*) '======covariância 1 ======'
    write(6,*) cov1(1,1),cov1(1,2)
    write(6,*) cov1(2,1),cov1(2,2)

    DO i=1,ndim
      DO j=1,ndim
        cov1(i,j)=cov1(i,j)/dfloat(np1)
      END DO
    END DO



  !      --------GRUPO 2 ---------------------
  !  	criando a matriz transposta g2T

    DO i=1,np2
      DO j=1,ndim
        g2T(j,i)=g2(i,j)
      END DO
    END DO

  !  ---------------------------------------------------
  !  	 - multiplicação de matrizes
  !  	   multiplicação de g2T por g2

     DO k=1,ndim
       DO j=1,ndim
         cov2(j,k)=0.d0
         DO i=1,np2
           cov2(j,k)=cov2(j,k)+g2T(j,i)*g2(i,k)
         END DO
       END DO
     END DO

     DO  i=1,ndim
       DO j=1,ndim
         cov2(i,j)=cov2(i,j)/dfloat(np2)
       END DO
     END DO

      WRITE(6,*) '======covariância 2 ======'
      WRITE(6,*) cov2(1,1),cov2(1,2)
      WRITE(6,*) cov2(2,1),cov2(2,2)


  !  	-------- covariância agrupada------

     DO i=1,ndim
       DO j=1,ndim
         covag(i,j)=dfloat(np1)*cov1(i,j)/(dfloat(np1+np2))+ &
         dfloat(np2)*cov2(i,j)/(dfloat(np1+np2))
       END DO
     END DO

      WRITE(6,*) '======covariância agrupada ======'
      WRITE(6,*) covag(1,1),covag(1,2)
      WRITE(6,*) covag(2,1),covag(2,2)

  !  	inversao da matriz covag - usando subrotina



     CALL INVERT(covag,ndim)
    ! Teste: dexando a matriz de cov =1
   !covag=0d0
    ! DO i=1,ndim
    !   covag(i,i)=1.0
     !END DO

      WRITE(6,*) '====== inv covariância agrupada ======'
      WRITE(6,*) covag(1,1),covag(1,2)
      WRITE(6,*) covag(2,1),covag(2,2)

  !  	diferenicas médias

     DO i=1,ndim
       md(i,1)=xm1(i)-xm2(i)
     END DO

  !  	write(6,*) '====== diferencias medias ======'
  !  	write(6,*) md(1,1)
  !  	write(6,*) md(2,1)

  !  	criando a matriz transposta mdT
  !  	---------------------------

    DO i=1,ndim
      DO j=1,1
        mdT(j,i)=md(i,j)
      END DO
    END DO

  !  ----------------------------------------------------
  !  	multiplicação de mdT por cov^-1
  !  	 - multiplicação de matrizes


    DO k=1,ndim
      DO j=1,1
        alfa(j,k)=0.d0
        DO i=1,ndim
          alfa(j,k)=alfa(j,k)+mdT(j,i)*covag(i,k)
        END DO
      END DO
    END DO

  !  ----------------------------------------------------
  !  	multiplicação de alfa por md
  !  	 - multiplicação de matrizes

    DO k=1,1
      DO j=1,1
        d2(j,k)=0.d0
        DO i=1,ndim  !2	!
          d2(j,k)=d2(j,k)+alfa(j,i)*md(i,k)
        END DO
      END DO
    END DO

    dist=dsqrt(d2(1,1))


  END SUBROUTINE mahalanobeana


  !--------------------------------------------------------------------------


     SUBROUTINE INVERT(A,i)
        integer i,im,j,k,l
        real*8 A(i,i),B(i)

         IM=I-1

         DO 5 K=1,I
           DO 2 J=1,IM
             2 B(J)=A(1,J+1)/A(1,1)
             B(I)=1.d0/A(1,1)
             DO 4 L=1,IM
               DO 3 J=1,IM
                 3 A(L,J)=A(L+1,J+1)-A(L+1,1)*B(J)
                 4 A(L,I)=-A(L+1,1)*B(I)
                 DO 5 J=1,I
                   5 A(I,J)=B(J)

     END SUBROUTINE INVERT




END PROGRAM teste_dists
