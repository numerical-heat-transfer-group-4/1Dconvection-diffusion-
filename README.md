# 1Dconvection-diffusion
      
      PROGRAM TEAMWORK
*     分别表示工质速度和网格数
      REAL VELOCITY
      INTEGER N=10
      WRITE(*,*) '请输入工质速度'
      READ (*,*)VELOCITY
      CALL   CENTRRAL_DIF
      CALL   用上风格式的解的函数
      CALL   用混合格式解的函数
      CALL   用解析解解的函数
      PAUSE
      END

C     中心差分格式求解函数
      FUNCTION CENTRAL_DIF（N,VELOCITY）
      REAL LEN=1,DEN=1,GAMA=0.1
*     后面用来循环初值用的整数 
      INTEGER I=2
*     分别表征对流强度和扩散强度的系数,以及网格长度
      REAL F,D，DELTAX
*     标准形式里的各系数aP,aE,aW,C
      REAL,ALLOCATABLE::AP(:),AE(:),AW(:),AC(:)，X(:)
      ALLOCATE(AP(N))
      ALLOCATE(AE(N))
      ALLOCATE(AW(N))
      ALLOCATE(AC(N))
      ALLOCATE(X(N))
      
      DELTAX=LEN/N
      D=GAMA/DELTAX
      F=DEN*(VELOCITY)
      
*     左边界网格的系数
      AP(1)=3*(D)+F/2
      AE(1)=D-F/2
      AW(1)=0
      AC(1)=(F+2*(D))
      
*     右边界网格的系数
      AP(N)=3*(D)-F/2
      AE(N)=0
      AW(N)=D+F/2
      AC=(F+2*(D))
      
*     中间网格的系数
      
      DO 10,WHILE I<N
      AP(I)=2* D
      AE(I)=D-F/2
      AW(I)=D+F/2
      AC(I)=0
10    CONTINUE
      CALL TDMA(AP,AE,AW,AC,N,X)
      OPEN(UNIT=1,FILE='CENTRAL_DIF.txt'）
      DO K=1,N,1
      WIRTE(1,* )  X (K)
      END DO
      RETURN
      END


C     迎风格式求解函数


C     混合格式求解函数
 
 
C     解析解求解函数
 
 
 
C     TDMA解法
      FUCTION TDMA(AP,AE,AW,AC,N,X)
      AE(1)=AE(1)/AP(1)
      AC(1)=AC(1)/AP(1)
      DO I=2，N
          AE(I)=AE(I)/(AP(I)-AW(I)* AE(I-1))
          AC(I)=(AC(I)+AW(I)* AC(I-1))/(AP(I)-AW(I)* A(I-1))
      END DO
      X(N)=AC(N)
      DO I=N-1,1,-1
          X(I)=AE(I)* X(I+1)+AC(I)
      END DO
      END 
          
      
      
