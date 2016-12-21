! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
! ************************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
! ************************************************************************

      subroutine sgrnd(seed)
      implicit integer(a-z)

! Period parameters
      parameter(N     =  624)
      dimension mt(0:N-1)
! the array for the state vector
      common /block/mti,mt
      save   /block/

!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]

      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue

      end subroutine sgrnd

! ************************************************************************

      double precision function grnd()
      implicit integer(a-z)

! Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
! constant vector a
      parameter(UMASK = -2147483648)
! most significant w-r bits
      parameter(LMASK =  2147483647)
! least significant r bits
! Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)

      dimension mt(0:N-1)
! the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
! mti==N+1 means mt[N] is not initialized

      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01

! mag01(x) = x * MATA for x=0,1

      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)

      if(mti.ge.N) then
! generate N words at one time
        if(mti.eq.N+1) then
! if sgrnd() has not been called,
          call sgrnd(4357)
! a default initial seed is used
        endif

        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif

      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))

      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif

      end function grnd

! ************************************************************************

      subroutine gaussrnd(rr)
      implicit integer(i-n)
      implicit double precision(a-h,o-y)
      implicit complex*16(z)
      
 10   x1 = 2.0D0*grnd() - 1.0D0;
      x2 = 2.0D0*grnd() - 1.0D0;
      w = x1**2 + x2**2;

      if (w .ge. 1.0D0) then
         goto 10
      endif

      w = dsqrt((-2.0D0*dlog(w))/w);
      rr = x1*w

      return
      end subroutine gaussrnd

