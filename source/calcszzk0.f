      subroutine calcszzk0(n, xk, szz, npts, mdeg, 
     &                    beta, s0, s2, eps, ierr)
c-----------------------------------------------------------------------
c     Extrapolate the charge-charge structure factor S_zz(k) to k -> 0
c     and return the second-moment coefficient
c
c         S2 = lim_{k->0} S_zz(k) / k**2
c
c     together with the dielectric constant from the (Stillinger-Lovett)
c     second-moment condition, in Gaussian units,
c
c         1 - 1/eps = 4*pi*beta * S2   ->   eps = 1/(1 - 4*pi*beta*S2).
c
c     S_zz(k) is even in k and, for a charge-neutral solvent, behaves as
c
c         S_zz(k) = S2*k**2 + S4*k**4 + O(k**6),
c
c     so a polynomial in u = k**2 is least-squares fitted over the NPTS
c     lowest-k samples:
c
c         S_zz = a0 + a1*u + a2*u**2 + ...      (mdeg = highest power).
c
c     a0 is returned as S0; it is the extrapolated S_zz(0) and should
c     come out ~0.  Its departure from 0 is a useful diagnostic of the
c     long-range / neutrality handling of h(r).  a1 is S2.
c
c     Arguments
c       n     (in)  number of k samples
c       xk    (in)  k grid, length n, ascending, xk(i) > 0 (k=0 may be
c                   present as xk(1)=0; it is handled correctly)
c       szz   (in)  S_zz(k), length n
c       npts  (in)  number of lowest-k points used in the fit
c       mdeg  (in)  polynomial degree in u=k**2  (1, 2 or 3)
c       beta  (in)  1/(kB T) in the same units as the charges; if
c                   beta <= 0, eps is not evaluated (returned 0)
c       s0    (out) extrapolated S_zz(0)              [diagnostic, ~0]
c       s2    (out) lim S_zz/k**2
c       eps   (out) dielectric constant (0 if beta<=0 or singular)
c       ierr  (out) 0 = ok; nonzero = error (see code)
c
c     Notes
c       * Pass k and S_zz on the RISM transform grid as-is; any DST
c         convention (k_i = i*dk, (i-1)*dk, (i-0.5)*dk ...) is fine
c         because the actual k values are supplied in xk.
c       * Choose NPTS so the window stays in the linear-in-u regime;
c         check stability of S2 against NPTS and MDEG, and that S0 ~ 0.
c       * Unweighted fit.  Emphasis on small k is obtained by keeping
c         NPTS small rather than by weighting (which amplifies noise).
c       * For SI units replace 4*pi by 1/eps0 in the eps formula.
c-----------------------------------------------------------------------
      implicit none
c     ---- arguments ----
      integer    n, npts, mdeg, ierr
      real*8     xk(n), szz(n), beta, s0, s2, eps
c     ---- locals ----
      integer    mxc
      parameter (mxc = 4)
      integer    ncoef, i, p, q
      real*8     amat(mxc,mxc), bvec(mxc), cvec(mxc), px(mxc)
      real*8     u, umax, x, pi, fourpb, denom
      parameter (pi = 3.14159265358979324d0)
c
      ierr = 0
      s0   = 0.0d0
      s2   = 0.0d0
      eps  = 0.0d0
c
c     ---- argument checks ----
      if (mdeg .lt. 1 .or. mdeg .gt. mxc-1) then
         ierr = 3
         return
      end if
      ncoef = mdeg + 1
      if (npts .lt. ncoef) then
         ierr = 1
         return
      end if
      if (npts .gt. n) then
         ierr = 2
         return
      end if
c
c     ---- scale u = k**2 by its largest value in the window ----
c     ---- (improves conditioning of the normal equations)    ----
      umax = xk(npts)*xk(npts)
      if (umax .le. 0.0d0) then
         ierr = 4
         return
      end if
c
c     ---- assemble normal equations for the fit in x = u/umax ----
      do p = 1, ncoef
         bvec(p) = 0.0d0
         do q = 1, ncoef
            amat(p,q) = 0.0d0
         end do
      end do
      do i = 1, npts
         u  = xk(i)*xk(i)
         x  = u/umax
         px(1) = 1.0d0
         do p = 2, ncoef
            px(p) = px(p-1)*x
         end do
         do p = 1, ncoef
            bvec(p) = bvec(p) + szz(i)*px(p)
            do q = 1, ncoef
               amat(p,q) = amat(p,q) + px(p)*px(q)
            end do
         end do
      end do
c
c     ---- solve the (ncoef x ncoef) normal-equation system ----
      call gesolv(amat, bvec, cvec, ncoef, mxc, ierr)
      if (ierr .ne. 0) then
         ierr = 10 + ierr
         return
      end if
c
c     ---- unscale: coefficient of u**j is cvec(j+1)/umax**j ----
      s0 = cvec(1)
      s2 = cvec(2)/umax
c
c     ---- dielectric constant (optional) ----
      if (beta .gt. 0.0d0) then
         fourpb = 4.0d0*pi*beta
         denom  = 1.0d0 - fourpb*s2
         if (abs(denom) .gt. 1.0d-12) then
            eps = 1.0d0/denom
         else
            ierr = 5
         end if
      end if
c
      return
      end
c
c
      subroutine gesolv(a, b, x, n, lda, ierr)
c-----------------------------------------------------------------------
c     Solve A x = b for a small dense system (n <= 4) by Gaussian
c     elimination with partial pivoting.  Works on a local copy so the
c     caller's A and b are left unchanged.
c-----------------------------------------------------------------------
      implicit none
      integer  n, lda, ierr
      real*8   a(lda,*), b(*), x(*)
      integer  mxc
      parameter (mxc = 4)
      real*8   aa(mxc,mxc), bb(mxc), piv, fac, tmp
      integer  i, j, k, ip
c
      ierr = 0
      if (n .gt. mxc) then
         ierr = 1
         return
      end if
c     ---- local copy ----
      do j = 1, n
         bb(j) = b(j)
         do i = 1, n
            aa(i,j) = a(i,j)
         end do
      end do
c     ---- forward elimination with partial pivoting ----
      do k = 1, n-1
         ip  = k
         piv = abs(aa(k,k))
         do i = k+1, n
            if (abs(aa(i,k)) .gt. piv) then
               piv = abs(aa(i,k))
               ip  = i
            end if
         end do
         if (piv .le. 0.0d0) then
            ierr = 2
            return
         end if
         if (ip .ne. k) then
            do j = k, n
               tmp      = aa(k,j)
               aa(k,j)  = aa(ip,j)
               aa(ip,j) = tmp
            end do
            tmp    = bb(k)
            bb(k)  = bb(ip)
            bb(ip) = tmp
         end if
         do i = k+1, n
            fac = aa(i,k)/aa(k,k)
            do j = k+1, n
               aa(i,j) = aa(i,j) - fac*aa(k,j)
            end do
            bb(i) = bb(i) - fac*bb(k)
         end do
      end do
      if (abs(aa(n,n)) .le. 0.0d0) then
         ierr = 2
         return
      end if
c     ---- back substitution ----
      x(n) = bb(n)/aa(n,n)
      do i = n-1, 1, -1
         tmp = bb(i)
         do j = i+1, n
            tmp = tmp - aa(i,j)*x(j)
         end do
         x(i) = tmp/aa(i,i)
      end do
      return
      end


c-----------------------------------------------------------------------
c     Pick the number of lowest-k points for the k->0 fit in szzk0,
c     from the r-grid (ngrid, rdelta) and a k-space window ceiling kmax.
c
c       dk            = pi/(ngrid*rdelta)        (DST conjugate spacing)
c       points < kmax = kmax/dk = kmax*ngrid*rdelta/pi
c
c       ngrid  (in) number of r (= k) grid points
c       rdelta (in) r grid spacing [Angstrom]
c       kmax   (in) window ceiling [Angstrom^-1]
c                   ~0.3 typical; ~0.2 or less for strong-dipole /
c                   electrolyte solvents (narrower linear-in-u region)
c       mdeg   (in) polynomial degree in u=k**2 used by szzk0 (1..3)
c-----------------------------------------------------------------------
      integer function nptsk0(ngrid, rdelta, kmax, mdeg)
      implicit none
      integer  ngrid, mdeg, np, nmin
      real*8   rdelta, kmax, pi, dk
      parameter (pi = 3.14159265358979324d0)
c
      dk = pi/(dble(ngrid)*rdelta)
      np = nint(kmax/dk)
c     ---- need at least ncoef points; keep a few extra for stability --
      nmin = max(mdeg+1, 5)
      if (np .lt. nmin)  np = nmin
      if (np .gt. ngrid) np = ngrid
      nptsk0 = np
      return
      end
