# 1Dconvection-diffusion
      
      PROGRAM TEAMWORK
*     分别表示工质速度和网格数
      REAL VELOCITY
      INTEGER N
      WRITE(*,*) '请输入工质速度'
      READ (*,*) VELOCITY
      N=10
      CALL   CENTRAL_DIF(N,VELOCITY)
C      CALL   用上风格式的解的函数
C      CALL   用混合格式解的函数
      CALL   SOLUTION(N,VELOCITY)
      PAUSE
      END

C     中心差分格式求解函数
      SUBROUTINE CENTRAL_DIF(N,VELOCITY)
      REAL VELOCITY
      INTEGER N
      REAL LEN,DEN,GAMA,F,D,DELTAX,FAI0,FAIL
*     标准形式里的各系数aP,aE,aW,C
      REAL AP(10),AE(10),AW(10),AC(10),X(10)
*     后面用来循环初值用的整数 
      INTEGER I
C      ALLOCATE(AP(N))
C      ALLOCATE(AE(N))
C      ALLOCATE(AW(N))
C      ALLOCATE(AC(N))
C      ALLOCATE(X(N))
      FAI0=1
      FAIL=0
      I=2
      LEN=1
      DEN=1
      GAMA=0.1
      DELTAX=LEN/N
      D=GAMA/DELTAX
      F=DEN*(VELOCITY)
      
*     左边界网格的系数
      AP(1)=3*(D)+F/2
      AE(1)=D-F/2
      AW(1)=0
      AC(1)=(F+2*(D))*FAI0
      
*     右边界网格的系数
      AP(10)=3*(D)-F/2
      AE(10)=0
      AW(10)=D+F/2
      AC(10)=(F+2*(D))*FAIL
      
*     中间网格的系数
      DO 10 ,I=2,N-1,1
          AP(I)=2*D
          AE(I)=D-F/2
          AW(I)=D+F/2
          AC(I)=0
10    CONTINUE      
      CALL TDMA(AP,AE,AW,AC,N,X)
      OPEN(UNIT=1,FILE='CENTRAL_DIF.txt')
      DO K=1,N,1
      WRITE(1,*)  X(K)
      END DO
      RETURN
      END


C     迎风格式求解函数


C     混合格式求解函数
 
 
C     解析解求解函数
      SUBROUTINE SOLUTION(N,VELOCITY)
      INTEGER N
      REAL VELOCITY
      REAL LEN,DEN,GAMA,DELTAX,FAI0,FAIL,X(10)
      LEN=1
      DEN=1
      GAMA=0.1
      DELTAX=LEN/N
      FAI0=1
      FAIL=0
      DO 30,I=1,N,1
          X(I)=(FAIL-FAI0)
     1*(EXP(DEN* VELOCITY*(I-0.5)* DELTAX/GAMA)-1)
     1/(EXP(DEN* VELOCITY* LEN/GAMA)-1)+FAI0 
30    CONTINUE          
      OPEN(UNIT=1,FILE='SOLUTION.txt')
      DO K=1,N,1
      WRITE(1,* )  X(K)
      END DO
      RETURN 
      END
 
 
C     TDMA解法
      SUBROUTINE TDMA(AP,AE,AW,AC,N,X)
      INTEGER N
      INTEGER I
      REAL AP(10),AE(10),AW(10),AC(10),X(10)
      AE(1)=AE(1)/(AP(1))
      AC(1)=AC(1)/(AP(1))
      DO 100, I=2,N,1
          AE(I)=(AE(I))/((AP(I))-(AW(I))* (AE(I-1)))
          AC(I)=((AC(I))+(AW(I))* (AC(I-1)))/((AP(I))-(AW(I))* (AE(I-1)))
100    CONTINUE
      X(10)=(AC(10))
      DO 20, I=N-1,1,-1
          X(I)=AE(I)* X(I+1)+AC(I)
20    CONTINUE
      END 
          
      
      
