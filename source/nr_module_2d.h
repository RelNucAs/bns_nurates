//...Module to implement Newton Rapson method for arbitrary dimension
//...The module contains 'ludcmp' which performs the decomposition of a matrix
//...into an LU product; 'lubskb' which solves a linear system and 'mnewt' which
//...implements Newton Rapson

	module nr_module

	use system_module

	implicit none

	integer, parameter :: m = 2 !Dimension of the system
	contains

	subroutine ludcmp(a,indx,d)
	use nrutil
	implicit none
	real, dimension(:,:), intent(inout) :: a
	integer, dimension(:), intent(out)  :: indx
	real, intent(out)		    :: d
!...Given and NxN input matrix a, this routine replaces it by the LU decomposition of
!...a rowwise permutation of itself. On output, a is arranged as in equation (2.3.14);
!...indx is and output vector of lenght N that records the row permutation effected by the
!...partial pivoting; d is output as +/-1 depending on whether the number of row interchanges
!...was even or odd, respectively. This routine is used in combination with lubksb to solve 
!...linear equations or invert a matrix.

	real, dimension(size(a,1)) :: vv !vv stores the implicit scaling of each row	
	real, parameter	:: tiny=1.0e-20 !a small number
	integer	:: j,n,imax
	
	n = assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
	d = 1.0	      !No row interchanges yet
	vv = maxval(abs(a),dim=2)  !Loop over rows to get the implicit scaling information
	if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp') !There is a row of zeros
	vv=1/vv		!Save the scaling
	do j=1,n
	  imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))  !Find the pivot row
	  if (j /= imax) then   !Do we need to interchange rows?
	    call swap(a(imax,:),a(j,:))  !Yes do so...  	
	    d=-d		!and change the parity of d
	    vv(imax)=vv(j)       !also interchange the scale factor
	  end if
	  indx(j)=imax
	  if (a(j,j) == 0.0) a(j,j)=tiny
!.......If the pivot element is zero the matrix is singular. For some applications on singular
!.......matrices, it is desirable to substitute tiny for zero
	  a(j+1:n,j) = a(j+1:n,j)/a(j,j)  !divide by the pivot element
	  a(j+1:n,j+1:n) = a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!........Reduce the remaining submatrix	 
	end do

	end subroutine ludcmp

	subroutine lubksb(a,indx,b)

	use nrutil, only: assert_eq

	implicit none

	real, dimension(:,:), intent(in)  :: a
	integer, dimension(:), intent(in) :: indx
	real, dimension(:), intent(inout) :: b
!...Solves the set of N linear equations A X = B. Here the NxN matrix a is input,
!...not as the original matrix A, but rather as its LU decomposition, determined by 
!...the routine ludcmp. indx is input as the permutation vector of lenght N returned by ludcmp.
!...b is input as the right-hand-side vector B, also of lenght N, and returns
!...with the solution vector X. a and indx are not modified by this routine and can be 
!...left in place for successive calls with different right-hand sides b. This routine 
!...takes into account the possibility that b will begin with many zero elements, so it is efficient
!...for use in matrix inversion.

	integer :: i, n, ii, ll
	real :: summ
	n = assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
	ii=0	
	do i=1,n
	  ll=indx(i)
	  summ=b(ll)
	  b(ll)=b(i)
	  if (ii /= 0) then
	    summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
	  else if (summ /= 0.0) then
	    ii=i
	  end if
	  b(i) = summ  
	end do
	do i=n,1,-1
	  b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
	end do
	end subroutine lubksb

void mnewt_mod(const int ntrail, 
	subroutine mnewt(ntrial,x,den,t,yle,ylmu,tolx,convergence,nnr)
	use nrutil
	implicit none

	integer, intent(in)		  :: ntrial      !Maximum number of trials
	real, intent(in)		  :: tolx	     !Tolerance
	double precision          :: den, t, yle, ylmu !Variables known in the system
	real, dimension(m), intent(inout) :: x
	integer, intent(out)		  :: convergence,nnr

!...Given an initial guess for a root in m dimensions, take ntrial Newton-Raphson steps
!...to improve the root. Stop if the root converges in either summed absolute variable
!...increments tolx or summed absolute function values tolf.
	
	integer	  :: i,j,jmax
	integer, dimension(size(x)) :: indx
	real :: d, tolf, err_new, err_old
	real, dimension(size(x))	  :: fvec, p, tmp, relat_error
	real, dimension(size(x),size(x))  :: fjac

!.......Set the tolerance on relative erros
	!tolf = 1.e-7
        tolf = 1.e-5
!.......Initialize the convergence flag 
!.......1: convergence reached, 2: convergence not reached
	convergence = 1  
!.......Initialize your equations with the input guess
        call system(den,t,yle,ylmu,x,fvec,fjac)
!.......Define relative errors
	relat_error(1) = fvec(1)/(x(1)+x(2))
	relat_error(2) = fvec(2)/(x(1)+x(2))
        err_new = sum(abs(relat_error))
        
        i=1
        do while(i<ntrial .and. err_new>tolf)
	  write(*,*)'NR trial: ',i
          err_old = err_new
	  p=-fvec                                !Right-hand side of linear equations
	  call ludcmp(fjac,indx,d)               !Solve linear equations using LU decomposition
	  call lubksb(fjac,indx,p)
	  tmp = x
          err_new = err_old + 1.
          jmax = 10
          j=1
	  do while(err_new>err_old .and. j<jmax)
           x = tmp+p/j
           call system(den,t,yle,ylmu,x,fvec,fjac)
	   relat_error(1) = fvec(1)/(x(1)+x(2))
	   relat_error(2) = fvec(2)/(x(1)+x(2))
           err_new = sum(abs(relat_error))
           j=j+1
          end do
          if(j+1 > jmax) then
            write(6,*)"Warning: j exceeded jmax"
          end if
          i=i+1
         end do

         if(i==ntrial) then
          call nrerror('mnewt exceeded maximum iterations')
          nnr = i
          convergence = 2
          return
         end if         

         if (sum(abs(p)) <= tolx .and. sum(abs(relat_error)) <= tolf) then 
          nnr = i
          return
         else
          call nrerror('cannot reach the precision')
          nnr = i
          convergence = 2
         end if
	end subroutine mnewt


	subroutine mnewt_cold(ntrial,x,den,tolx,convergence,nnr)
	use nrutil
	implicit none

	integer, intent(in)		  :: ntrial      !Maximum number of trials
	real, intent(in)		  :: tolx	     !Tolerance
	double precision          :: den         !Variables known in the system
	real, dimension(m), intent(inout) :: x
	integer, intent(out)		  :: convergence,nnr

!...Given an initial guess for a root in m dimensions, take ntrial Newton-Raphson steps
!...to improve the root. Stop if the root converges in both summed absolute variable
!...increments tolx or summed absolute function values tolf.
	
	integer	  :: i,j,jmax
	integer, dimension(size(x)) :: indx
	real :: d, tolf, err_new, err_old
	real, dimension(size(x))	  :: fvec, p, tmp, error
	real, dimension(size(x),size(x))  :: fjac

!.......Set the tolerance on relative erros
	tolf = 1.e-6

!.......Initialize the convergence flag 
!.......1: convergence reached, 2: convergence not reached
	convergence = 1  
!.......Initialize your equations with the input guess
    call system_cold_equilibrium(den,x,fvec,fjac)
!.......Define relative errors
	error(1) = fvec(1)
	error(2) = fvec(2)
    err_new = sum(abs(error))
        
    i=1
    do while(i<ntrial .and. err_new>tolf)
	  !write(*,*)'NR trial: ',i
      err_old = err_new
	  p=-fvec                                !Right-hand side of linear equations
	  call ludcmp(fjac,indx,d)               !Solve linear equations using LU decomposition
	  call lubksb(fjac,indx,p)
	  tmp = x
      err_new = err_old + 1.
      jmax = 10
      j=1
	  do while(err_new>err_old .and. j<jmax)
        x = tmp+p/j
        call system_cold_equilibrium(den,x,fvec,fjac)
	    error(1) = fvec(1)
	    error(2) = fvec(2)
        err_new = sum(abs(error))
        j=j+1
      end do
      if(j+1 > jmax) then
        write(6,*)"Warning: j exceeded jmax"
      end if
      i=i+1
    end do

    if(i==ntrial) then
      call nrerror('mnewt exceeded maximum iterations')
      nnr = i
      convergence = 2
      return
    end if         

    if (sum(abs(p)) <= tolx .and. sum(abs(error)) <= tolf) then 
      nnr = i
      return
    else
      call nrerror('cannot reach the precision')
      nnr = i
      convergence = 2
    end if

	end subroutine mnewt_cold


	end module nr_module
