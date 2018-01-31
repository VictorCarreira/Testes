PROGRAM teste_dists

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    !Criação da subrotina do cálculo de distância euclideana                     !
    !Orientador: Cosme Ferreira da Ponte Neto                                    !
    !Aluno: Victor Ribeiro Carreira                                              !
    !Este programa validar a subrotina mahalanobeana                             !
    !Categoria: classificador                                                    !
    !Subrotina mahalanobeana                                                     !
    !Para usar compilação com flags utilize:                                     !
    !gfortran -fbounds-check -fbacktrace -Wall -Wextra -pedantic                 !
    !"pasta/subpasta/nomedopragrama.f95" -o nomedoexecutável                     !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


                          !lito1 nuvem de pontos
                          !lito2 um ponto no espaço
                          !np; número de pontos da nuvem

  IMPLICIT NONE

  INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
  INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

  INTEGER(KIND=SP):: i,j
  INTEGER(KIND=DP), PARAMETER::np=359

  REAL(KIND=DP)::r,theta, a, b
  REAL(KIND=DP)::maha, eucli
  REAL(KIND=DP), PARAMETER::pi=3.141592653, rmax=5.0, rmin=0.0, thetamin=0.0, &
  thetamax=2.0*pi, amax=0.1, amin=0.0, bmax=0.1,bmin=0.0

  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:)::lito1, lito2

  ALLOCATE(lito1(np,2), lito2(1,2))


  lito1=0d0  !Zerando a matriz
  lito2=0d0  !Zerando a matriz

  !Defindo as coordenadas do ponto no espaço (Centro do sistema de coordenadas)
  lito2(1,1)=100d0 ! atribuindo o valor zero a cooredentada x
  lito2(1,2)=100d0 ! atribuindo o valor zero a cooredentada y


   OPEN(1,FILE='teste_nuvem.txt')
   OPEN(2,FILE='teste_ponto.txt')


   ! Criando a nuvem de pontos no espaço em lito1 ELISE
   DO i=1,np
     a=(amax-amin)*RAND()+amin !Equação padrão para sortear números aleatórios em uma dist regular
     b=(bmax-bmin)*RAND()+bmin
     theta=(thetamax-thetamin)*RAND()+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
     lito1(i,1)=a*COS(theta) + 100.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
     lito1(i,2)=b*SIN(theta) + 100.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
     WRITE(1,FMT=*)lito1(i,1),lito1(i,2)
   END DO

   ! Criando a nuvem de pontos no espaço em lito1
   !DO i=1,np
    ! r=(rmax-rmin)*RAND()+rmin !Equação padrão para sortear números aleatórios em uma dist regular
     !theta=(thetamax-thetamin)*RAND()+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
     !lito1(i,1)=r*COS(theta) + 100.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
     !lito1(i,2)=r*SIN(theta) + 100.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
     !WRITE(1,FMT=*)lito1(i,1),lito1(i,2)
   !END DO

  ! Criando uma nuvem homogênia de pontos
  !DO j=1,359
    !r=(rmax-rmin)*1+rmin !Equação padrão para sortear números aleatórios em uma dist regular
    !  r=5
      !DO i=1,np
      !theta=(thetamax-thetamin)*i+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
    !   theta=DFLOAT(j)*pi/180.0
    !   lito1(j,1)=r*COS(theta) + 100.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
    !   lito1(j,2)=r*SIN(theta) + 100.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
      !END DO
    ! WRITE(1,FMT=*)lito1(j,1),lito1(j,2)
   !END DO

    WRITE(2,FMT=*)lito2(1,1),lito2(1,2)

    CLOSE(1)
    CLOSE(2)


  21 FORMAT(15(F4.2,2x))

  CALL euclideana(lito1,lito2,eucli)
  CALL mahalanobeana(lito1,np,lito2,1,2,maha)

  PRINT*,'dist de euclides=', eucli
  PRINT*, 'dist de maha=', maha




  !*******************************************************************************************!
  CONTAINS

  SUBROUTINE euclideana(lito1,lito2,eucli)

   IMPLICIT NONE
   INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
   INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

   REAL(KIND=DP), DIMENSION(:,:), INTENT(IN)::lito1, lito2
   REAL(KIND=DP), INTENT(OUT):: eucli
   REAL(KIND=DP)::media1

   INTEGER(KIND=SP):: k

    eucli=0d0

    IF(SIZE(lito1(1,:)) /= SIZE(lito2(1,:)))THEN
      PRINT*,'WARNING! THE PROPERTIES NUMBER´S OF lito1 AND lito2 MUST BE THE SAME.'
      STOP
      RETURN
    END IF

    DO k=1,SIZE(lito1(1,:))  ! Inicia o laço da primeira até a última propriedade que é dado pelo size de lito
     media1=0d0 !zera as variáveis
     media1=SUM(lito1(:,k))/SIZE(lito1(:,k)) !calcula as médias para as k propriedades
     eucli= eucli + (lito2(1,k)-media1)**2 ! Cálculo da medida de semelhança de euclides
    END DO ! Final do laço das k propriedades

    eucli=SQRT(eucli)
  END SUBROUTINE euclideana

  SUBROUTINE mahalanobeana(g11,np1,g22,np2,ndim,dist)

  !  	subrotina que calcula a distância de mahalanobis entre
  !  	dois agrupamentos de elementos com dimensão ndim

     IMPLICIT NONE
      INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
      INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

      INTEGER(KIND=DP), INTENT(IN):: np1
      INTEGER(KIND=SP), INTENT(IN):: np2, ndim
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
