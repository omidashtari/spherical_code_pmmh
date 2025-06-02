module mod_IterativeSolvers_test
	implicit none
	
	contains
	subroutine GMRES(A, m, n, itmax, tol, b, ubest, gmres_iters)
		implicit none
		
		! program inputs
		integer, intent(in) :: n                                 ! size of the linear system
		integer, intent(in) :: m                                 ! size of the largest Krylov subspace
		integer, intent(in) :: itmax                             ! maximum number of iterations
		double precision, dimension(n, n), intent(in) :: A       ! coefficient matrix
		double precision, dimension(n), intent(in) :: b          ! RHS vector
		double precision, intent(in) :: tol                      ! convergence criterion

		! program outputs
		double precision, dimension(n), intent(out) :: ubest     ! best x in Ax = b
		integer, intent(out) :: gmres_iters                      ! number of GMRES iterations

		! internal variables
		double precision, dimension(n) :: x
		double precision, dimension(m+1, m) :: h, h_
		double precision, dimension(n, m+1) :: v
		double precision :: res
		double precision, dimension(n) :: w, z
		double precision, dimension(m + 1) :: y, p
		double precision, dimension(4*m+1) :: work
		integer :: its, rank, i, j
		integer, dimension(m) :: piv
		double precision, save :: beta
		logical :: done

		! Initialize parameters
		its = 0
		x = 0.0    ! initial guess is zero vector
		v = 0.0

		outer_loop: do while (.true.)
			! Initialize beta and the Krylov basis
			beta = sqrt(dot_product(x,x))

			if (beta /= 0.0) then
				w = matmul(A,x)
			else
				w = 0.0
			end if

			w = b - w

			beta = sqrt(dot_product(w,w))

			v(:, 1) = w / beta

			! Initialize Hessenberg matrix h to zero
			h = 0.0

			inner_loop: do j = 1, m
				its = its + 1
				z = v(:, j)

				w = matmul(A,z)

				! Orthogonalize w against previous Krylov vectors
				do i = 1, j
					h(i, j) = dot_product(w, v(:, i))
					w = w - h(i, j) * v(:, i)
				end do

				h(j+1, j) = sqrt(dot_product(w,w))
				v(:, j+1) = w / h(j+1, j)

				! Solve the least squares problem using LAPACK's dgelsy
				p(1) = beta
				p(2:j+1) = 0.0
				h_(1:j+1, 1:j) = h(1:j+1, 1:j)

				call dgelsy(j+1, j, 1, h_, m+1, p, m+1, piv, m, rank, work, 4*m+1, i)

				if (i /= 0) then
					stop 'gmresm: dgelsy failed'
				end if

				y = p

				! Calculate the residual
				p(1:j+1) = - matmul(h(1:j+1, 1:j), y(1:j))
				p(1) = p(1) + beta

				res = sqrt(dot_product(p(1:j+1),p(1:j+1)))

				if (res > 1.0e5) then
					print*, "Residual too big, stopping simulation..."
					stop
				end if

				print*, 'gmresm: iteration = ', its, ' residual = ', res

				done = (res <= tol .or. its == itmax)

				if (done .or. j == m) then
					z = matmul(v(:, 1:j), y(1:j))
					x = x + z
					
					if (done) then
						exit outer_loop  ! Exit the outer loop as well
					end if

					exit inner_loop  ! Exit only the inner loop if needed
				end if
				
			end do inner_loop
		end do outer_loop

		gmres_iters = its
		ubest = x

	end subroutine GMRES
end module mod_IterativeSolvers_test
