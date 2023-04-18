*     Subroutine dsden compute the raw double smooth density estimate
*     Called from R
*
      SUBROUTINE dsden(x, nobs, wj, xout, yout, nout)
      DOUBLE PRECISION x(nobs), wj(nobs), xout(nout), yout(nout)
      DOUBLE PRECISION Gaussn, wjx, d, Kxix0, Kxix, Kxjx
      INTEGER*4 nobs, nout, I, J, K
* Compute wj
	  DO 20 j=1,nobs
         wjx = 0.d0
		 DO 10 I=1,nobs
		    wjx = wjx + Gaussn((x(i) - x(j))/dsqrt(2d0))
 10     CONTINUE
        wj(j) = DBLE(nobs)/wjx
 20   CONTINUE
* Loop through output
      DO 130 K=1,nout
        d = 0.d0
        DO 120 I=1,nobs
          Kxix0 = Gaussn((x(i)-xout(k))/dsqrt(3d0))
          DO 110 J=1,nobs
            Kxix = Gaussn((x(I)-x(J))/dsqrt(3d0))
            Kxjx = Gaussn((x(J)-xout(K))/dsqrt(3d0))
            d = d + Kxix0*wj(J)*Kxix*Kxjx/dsqrt(3d0)**3
 110      CONTINUE
 120    CONTINUE
        yout(k) = d
 130  CONTINUE
      RETURN
      END
**********************************************************************
*
*     Function Gaussn                            Called by: dsden
*
*     Calculate gaussian kernel easily, without the divisor
*
*     Arguments
*     ---------
*     dist the relative distance
*
**********************************************************************
      DOUBLE PRECISION FUNCTION Gaussn(dist)
*
*     subroutine arguments
*
      DOUBLE PRECISION dist
*
      Gaussn = dexp(-0.5D0*dist**2)
      RETURN
      END