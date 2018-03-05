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

  INTEGER(KIND=SP):: i,j, ic1, ic2
  INTEGER(KIND=SP)::np1, np2
  INTEGER(KIND=SP), PARAMETER:: semente=5

  REAL(KIND=DP)::r,theta, a, b, inicio,final,x
  REAL(KIND=DP)::maha, eucli, L, Lmax, Lmin
  REAL(KIND=DP), PARAMETER::pi=3.141592653, rmax=3.0, rmin=0.0, thetamin=0.0, &
  thetamax=5.0*pi, amax=15.0, amin=3.0, bmax=5,bmin=0.0

  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:)::lito1, lito2, lito11, lito22
  
  CALL CPU_TIME(inicio)

  np1=5000
  np2=500

  ALLOCATE(lito1(np1,2), lito2(np2,2))

  lito1=0d0  !Zerando a matriz
  lito2=0d0  !Zerando a matriz

  !Defindo as coordenadas do ponto no espaço (Centro do sistema de coordenadas)
  !lito2(1,1)=100d0 ! atribuindo o valor zero a cooredentada x
  !lito2(1,2)=100d0 ! atribuindo o valor zero a cooredentada y


   OPEN(1,FILE='teste_nuvem_T.txt')
   OPEN(2,FILE='teste_ponto.txt')
   OPEN(3,FILE='teste_nuvem_C.txt')
   OPEN(4,FILE='Parametros.txt')
   
   ! Criando a nuvem de pontos, em elipse, no espaço em lito1 (conjunto de treinamento)
   !DO i=1,np1
   !  a=(amax-amin)*RAND()+amin !Equação padrão para sortear números aleatórios em uma dist regular
   !  b=(bmax-bmin)*RAND()+bmin
   !  theta=(thetamax-thetamin)*RAND()+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
   !  lito1(i,1)=a*COS(theta) + 1.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
   !  lito1(i,2)=b*SIN(theta) + 15.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
   !  WRITE(1,FMT=*)lito1(i,1),lito1(i,2)
   !END DO

   ! Criando a nuvem de pontos no espaço em lito1, em círculo, no conjunto lito1 (treinamento)
   

   !Criando uma semente para o gerador de números aleatórios
   DO i=1,semente
    x=RAND()
   END DO  

   !DO i=1,np1
   ! r=(rmax-rmin)*RAND()+rmin !Equação padrão para sortear números aleatórios em uma dist regular
   ! theta=(thetamax-thetamin)*RAND()+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
   ! lito1(i,1)=r*COS(theta) + 1.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
   ! lito1(i,2)=r*SIN(theta) + 15.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
   ! WRITE(1,FMT=*)lito1(i,1),lito1(i,2)
   !END DO

  !Corrigindo o vício do sorteio aleatório
   L=2.0*rmax ! lado de um quadrado (válido somente para circunferências)
   Lmax=2.0*amax !Lado máximo de um retângulo (válido somente para uma elipse)
   Lmin=2.0*amin !Lado mínimo de um retângulo (válido somente para uma elipse) 
 
  ! circunferência 1
  !DO i=1,np1
   !lito1(i,1)=(L-rmin)*RAND()+rmin-rmax 
   !lito1(i,2)=(L-rmin)*RAND()+rmin-rmax
  !END DO

  !Elispe 1
  DO i=1,np1
   lito1(i,1)=(Lmax-amax)*RAND()+amin-amax 
   lito1(i,2)=(Lmin-amin)*RAND()+amin-amax
  END DO

 
  ic1=0
  DO i=1,np1
   !IF (DSQRT(lito1(i,1)**2+lito1(i,2)**2) <= L/2.0) THEN ! válido para a circunferências. Seleciona os dados de acordo com a hipotenusa
   IF (DSQRT(lito1(i,1)**2.0*lito1(i,2)**2-lito1(i,2)**2*RANGE(lito1(i,2))/lito1(i,1)**2.0) <= Lmax/2.0 .AND. & 
   DSQRT(lito1(i,1)**2.0*lito1(i,2)**2-lito1(i,2)**2*RANGE(lito1(i,2))/lito1(i,1)**2.0) >= Lmin/2.0) THEN  !Válido para a elipse
    ic1=ic1+1!Contar quantos foram considerados
  !  lito1(ic1,1)=lito1(i,1)+1.0! desloca o centro em 1 unidades
  !  lito1(ic1,2)=lito1(i,2)+15.0! descloca o centro em 15 unidades
  !  WRITE(1,FMT=*)lito1(i,1),lito1(i,2)
   END IF
  END DO 

!Circunferência 2
  DO i=1,np2
  lito2(i,1)=(L-rmin)*RAND()+rmin-rmax 
  lito2(i,2)=(L-rmin)*RAND()+rmin-rmax
 END DO

 ic2=0
 DO i=1,np2
  IF (DSQRT(lito2(i,1)**2+lito2(i,2)**2) <= L/2.0) THEN
   ic2=ic2+1!Contar quantos foram considerados
  ! lito2(ic2,1)=lito2(i,1)+10.0! desloca o centro em 10 unidades
  ! lito2(ic2,2)=lito2(i,2)+5.0! descloca o centro em 5 unidades
  ! WRITE(3,FMT=*)lito2(i,1),lito2(i,2)
  END IF
 END DO 
  

 ALLOCATE(lito11(ic1,2), lito22(ic2,2))
 lito11=0d0
 lito22=0d0

 ic1=0
 DO i=1,np1
  !IF (DSQRT(lito1(i,1)**2+lito1(i,2)**2) <= L/2.0) THEN ! Circunferência
  IF (DSQRT(lito1(i,1)**2.0*lito1(i,2)**2-lito1(i,2)**2*RANGE(lito1(i,2))/lito1(i,1)**2.0) <= Lmax/2.0 .AND. & 
  DSQRT(lito1(i,1)**2.0*lito1(i,2)**2-lito1(i,2)**2*RANGE(lito1(i,2))/lito1(i,1)**2.0) <= Lmin/2.0) THEN
   ic1=ic1+1!Contar quantos foram considerados
   lito11(ic1,1)=lito1(i,1)+10.0! desloca o centro em 1 unidades
   lito11(ic1,2)=lito1(i,2)+21.0! descloca o centro em 15 unidades
    WRITE(1,FMT=*)lito11(ic1,1),lito11(ic1,2)
  END IF
 END DO 


 ic2=0
 DO i=1,np2
  IF (DSQRT(lito2(i,1)**2+lito2(i,2)**2) <= L/2.0) THEN
   ic2=ic2+1!Contar quantos foram considerados
   lito22(ic2,1)=lito2(i,1)+1.0! desloca o centro em 10 unidades
   lito22(ic2,2)=lito2(i,2)+5.0! descloca o centro em 5 unidades
   WRITE(3,FMT=*)lito22(ic2,1),lito22(ic2,2)
  END IF
 END DO 


 !Redefinindo o número de pontos/dados que entrarão na rotina
  np1=ic1
  np2=ic2






  ! Criando uma nuvem homogênia de pontos
  !DO j=1,359
  !  r=(rmax-rmin)*RAND()+rmin !Equação padrão para sortear números aleatórios em uma dist regular
  !    r=5
  !    DO i=1,np
  !     theta=(thetamax-thetamin)*i+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
  !     theta=DFLOAT(j)*pi/180.0
  !     lito1(j,1)=r*COS(theta) + 100.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
  !     lito1(j,2)=r*SIN(theta) + 100.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
  !    END DO
  ! WRITE(1,FMT=*)lito1(j,1),lito1(j,2)
  !END DO

   !Transformando o lito2 em uma nuvem de pontos elipsoidal
   !DO i=1,np
   !  a=(amax-amin)*RAND()+amin !Equação padrão para sortear números aleatórios em uma dist regular
   !  b=(bmax-bmin)*RAND()+bmin
   !  theta=(thetamax-thetamin)*RAND()+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
   !  lito2(i,1)=a*COS(theta) + 10.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
   !  lito2(i,2)=b*SIN(theta) + 5.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
   !  WRITE(3,FMT=*)lito2(i,1),lito2(i,2)
   !END DO


     !Transformando o lito2 em uma nuvem de pontos circular
     !PRINT*,'np2=',np2
     !DO i=1,np2
     ! r=(rmax-rmin)*RAND()+rmin !Equação padrão para sortear números aleatórios em uma dist regular
     ! theta=(thetamax-thetamin)*RAND()+thetamin !Equação padrão para sortear números aleatórios em uma dist regular
     ! lito2(i,1)=r*COS(theta) + 10.0 !Preenche a coordenada x e desloca o centro para 100 unidades em x
     ! lito2(i,2)=r*SIN(theta) + 5.0 !Preenche a coordenada y e desloca o centro em 100 unidade em y
     !WRITE(3,FMT=*)lito2(i,1),lito2(i,2)
    !END DO
   


   !Gravando os arquivos de saída (lito2 é um ponto)
    !WRITE(2,FMT=*)lito2(1,1),lito2(1,2)

   !Gravando os arquivos de saída (lito2 é uma nuvem)
    !WRITE(3,FMT=*)lito2(1,1),lito2(1,2)
    !Gravando os parâmetros do modelo
    WRITE(4,FMT=*)'Parâmetros do modelo'
    WRITE(4,FMT=*)
    WRITE(4,FMT=21)'rmax=',rmax
    WRITE(4,FMT=21)'rmin=',rmin
    WRITE(4,FMT=21)'thetamin=',thetamin
    WRITE(4,FMT=21)'thetamax=',thetamax
    WRITE(4,FMT=21)'amax=',amax
    WRITE(4,FMT=21)'amin=',amin
    WRITE(4,FMT=21)'bmax=', bmax
    WRITE(4,FMT=21)'bmin=',bmin
    WRITE(4,FMT=*)'----------------------'
    WRITE(4,FMT=22)'cluster 1=',np1
    WRITE(4,FMT=22)'cluster 2=',np2


  ! Formatos dos arquivos de saida
  21 FORMAT(A9,2x,E12.2)
  22 FORMAT(A20,2x,I10)

  CALL euclideana(lito11,lito22,eucli)
  CALL mahalanobeana(lito11,np1,lito22,np2,2,maha)

  PRINT*,'dist de euclides=', eucli
  PRINT*, 'dist de maha=', maha

  CLOSE(1)
  CLOSE(2)
  CLOSE(3)
  CLOSE(4)

  CALL CPU_TIME(final)
  PRINT*,'tempo de máquina=',final-inicio
  !*******************************************************************************************!
  CONTAINS

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
