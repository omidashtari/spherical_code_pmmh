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
		double precision, dimension(n) :: w, z, r
		double precision, dimension(m + 1) :: y, p
		double precision, dimension(4*m+1) :: work
		integer :: its, evals, rank, i, j
		integer, dimension(m) :: piv
		double precision, save :: beta
		double precision :: tolr, normr, normb
		logical :: done

		open(unit = 1, file = "./GMRES_residual.txt", status = "replace", action = "write")
		open(unit = 2, file = "./GMRES_evals.txt", status = "replace", action = "write")
		
		! Initialize parameters
		its = 0
		evals = 0
		x = 0.0    ! initial guess is zero vector
		v = 0.0
		normb = sqrt(dot_product(b, b))
		tolr = tol * normb

		r = b - matmul(A, x)
		evals = evals + 1
				
		normr = sqrt(dot_product(r, r))
		write(1, '(E12.5)') normr/normb
		write(2, '(I0)') evals
				
		if (normr <= tolr) then
			print*, "Initial guess is good enough!"
			stop
		end if
						
		outer_loop: do while (.true.)
			! Initialize beta and the Krylov basis
			beta = sqrt(dot_product(x,x))

			if (beta /= 0.0) then
				w = matmul(A, x)
				evals = evals + 1
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

				w = matmul(A, z)
				evals = evals + 1

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

				print*, 'gmresm: iteration = ', its, ' evals = ', evals, ' residual = ', res/normb
				write(1, '(E12.5)') res/normb
				write(2, '(I0)') evals

				done = (res <= tolr .or. its == itmax)

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

	    close(1)
	    close(2)
	end subroutine GMRES
	
	
	
	
	subroutine IDRs(A, b, s, n, tol, itmax, x0, ubest, idrs_iters)
		implicit none
		
		! program inputs ------------------------------------------------------------------------------------------------
		double precision, dimension(n, n), intent(in) :: A       ! coefficient matrix
		double precision, dimension(n), intent(in) :: b          ! RHS vector
		double precision, dimension(n), intent(in) :: x0         ! initial guess
		integer, intent(in) :: s                                 ! dimension of the shadow space
		integer, intent(in) :: n                                 ! size of the linear system
		integer, intent(in) :: itmax                             ! maximum number of iterations
		double precision, intent(in) :: tol                      ! convergence criterion (relative to the norm of b)

		! program outputs
		double precision, dimension(n), intent(out) :: ubest     ! best x for Ax = b
		integer, intent(out) :: idrs_iters                       ! number of IDR(s) iterations (A*x evaluations)
		
		! internal variables
		double precision, dimension(n) :: x                      ! solution iterates
		double precision, dimension(n) :: r                      ! residual vector
		double precision :: normr                                ! norm of the residual vector
		double precision :: normb                                ! norm of the RHS vector
		double precision :: tolr                                 ! absolute convergence criterion
		double precision, dimension(s, n) :: P                   ! Rows are orthonormal vectors spanning the shadow space
		double precision, dimension(n, s) :: dR, dX
		double precision, dimension(s, s) :: M, W
		double precision :: om
		integer :: its, evals, i, j, k, oldest
		double precision, dimension(n) :: random_vec, v, q, t
		double precision, dimension(s) :: c_, m_, dm_
		integer :: seed(33)                                      ! seed for random but reproducible vectors
		integer :: IPIV(s), INFO                                 ! working variables for x = inv(A) * b using LAPACK 
		
		open(unit = 1, file = "./IDRs_residual.txt", status = "replace", action = "write")
		open(unit = 2, file = "./IDRs_evals.txt", status = "replace", action = "write")
		
		! consistency checks --------------------------------------------------------------------------------------------
		if (s >= n) then
			print*, "The dimension of the shadow space must be smaller than the system!"
			stop
		end if
		
		! compute initial residual --------------------------------------------------------------------------------------
		its = 0
		evals = 0
		
		x = x0
		r = b - matmul(A, x)
		evals = evals + 1
		
		normr = sqrt(dot_product(r, r))
		normb = sqrt(dot_product(b, b))
		tolr = tol * normb
		
		write(1, '(E12.5)') normr/normb
		write(2, '(I0)') evals
		
		if (normr <= tolr) then
			print*, "Initial guess is good enough!"
			stop
		end if
		
		! construct the shadow space ------------------------------------------------------------------------------------
		do i = 1, 33
			seed(i) = i
		end do
    	call random_seed(put = seed)
    
		P(1,:) = r / normr
		if (s > 1) then
			do i = 2, s
				call random_number(random_vec)
				
				do j = 1, i-1
					random_vec = random_vec - dot_product(P(j,:), random_vec) * P(j,:)
				end do
				P(i,:) = random_vec / sqrt(dot_product(random_vec, random_vec))
			end do			
		end if
		
		do i = 1, s
			if (dot_product(P(i,:), P(i,:)) - 1.0 > 1.0d-12) then
				print*, "Problem with making a normal basis."
				stop
			end if
			
			if (i > 1) then
				do j = 1, i-1
					if (dot_product(P(i,:), P(j,:)) > 1.0d-12) then
						print*, "Problem with making an orthogonal basis."
						stop
					end if
				end do
			end if
		end do
		
		! produce start vectors -----------------------------------------------------------------------------------------
		dR = 0
		dX = 0
		
		do k = 1, s
			v = matmul(A, r)
			evals = evals + 1
			
			om = dot_product(v, r) / dot_product(v, v)
			
	        dX(:, k) = om * r
        	dR(:, k) = -om * v
        	
	        x = x + dX(:, k)
	        r = r + dR(:, k)
	        
	        M(:, k) = matmul(P, dR(:,k))
		end do
		
		! Main iteration loop, build G-spaces ---------------------------------------------------------------------------
		oldest = 1
    	m_ = matmul(P, r)
    	
    	IDRs_loop: do while (normr > tolr .and. its < itmax)
    		do k = 0, s
	    		c_ = m_
	    		W = M
    			call dgesv(s, 1, W, s, IPIV, c_, s, INFO)    ! directly solve the s-by-s system M * c_ = m_
    		
				if (INFO /= 0) then
					stop 'IDR(s): dgesv failed to solve an s-by-s linear system.'
				end if
			
				q = -matmul(dR, c_)
				v = r + q
    		
				if (k == 0) then
					t = matmul(A, v)
					evals = evals + 1
					
					om = dot_product(t, v) / dot_product(t, t)
	                dR(:, oldest) = q - om * t
	                dX(:, oldest) = -matmul(dX, c_) + om * v
				else
	                dX(:, oldest) = -matmul(dX, c_) + om * v
	                dR(:, oldest) = -matmul(A, dX(:, oldest))
	                evals = evals + 1
				end if
				
	            r = r + dR(:, oldest);
	            x = x + dX(:, oldest);
    	        
    	        normr = sqrt(dot_product(r, r))
    	        
    	        dm_ = matmul(P, dR(:, oldest))
    	        M(:, oldest) = dm_
    	        m_ = m_ + dm_
    	        
    	        oldest = oldest + 1
    	        if (oldest > s) then
    	        	oldest = 1
    	        end if
    		end do
    		
    		its = its + 1
    		print*, 'IDR(s): iteration = ', its, ' evals = ', evals, ' residual = ', normr/normb    		
			write(1, '(E12.5)') normr/normb
			write(2, '(I0)') evals
    	end do IDRs_loop
    	
	    close(1)
	    close(2)
	end subroutine IDRs
	
	
	
	
	subroutine BiCGSTAB(A, b, n, tol, itmax, x0, ubest, bicgstab_iters)
		implicit none
		
		! program inputs ------------------------------------------------------------------------------------------------
		double precision, dimension(n, n), intent(in) :: A       ! coefficient matrix
		double precision, dimension(n), intent(in) :: b          ! RHS vector
		double precision, dimension(n), intent(in) :: x0         ! initial guess
		integer, intent(in) :: n                                 ! size of the linear system
		integer, intent(in) :: itmax                             ! maximum number of iterations
		double precision, intent(in) :: tol                      ! convergence criterion (relative to the norm of b)

		! program outputs
		double precision, dimension(n), intent(out) :: ubest     ! best x for Ax = b
		integer, intent(out) :: bicgstab_iters                   ! number of BiCGSTAB iterations (A*x evaluations)

		! internal variables
		double precision, dimension(n) :: r                      ! residual vector
		double precision, dimension(n) :: x                      ! solution iterates
		double precision, dimension(n) :: rr, r_, p, t
		double precision, dimension(n) :: Ap, At
		double precision :: normr, normb, tolr
		double precision :: alpha, beta, xi
		integer :: its, evals
		
		open(unit = 1, file = "./BiCGSTAB_residual.txt", status = "replace", action = "write")
		open(unit = 2, file = "./BiCGSTAB_evals.txt", status = "replace", action = "write")
		
		! compute initial residual --------------------------------------------------------------------------------------
		its = 0
		evals = 0
		x = x0
		
		r = b - matmul(A, x)
		evals = evals + 1
		
		normr = sqrt(dot_product(r, r))
		normb = sqrt(dot_product(b, b))
		tolr = tol * normb
		
		write(1, '(E12.5)') normr/normb
		write(2, '(I0)') evals
			
		if (normr <= tolr) then
			print*, "Initial guess is good enough!"
			stop
		end if
		
		! BiCGSTAB loop -------------------------------------------------------------------------------------------------
		rr = r
		p = 0
		beta = 0
		xi = 0
		Ap = 0
		
		BiCGSTAB_loop: do while (normr > tolr .and. its < itmax)
			p = r + beta * (p - xi * AP)
			
			Ap = matmul(A, p)
			evals = evals + 1
			
			alpha = dot_product(rr, r) / dot_product(rr, Ap)
			t = r - alpha * Ap
			
			At = matmul(A, t)
			evals = evals + 1
			
			xi = dot_product(At, t) / dot_product(At, At)
			x = x + alpha * p + xi * t
			
			r_ = r
			r = t - xi * At
			normr = sqrt(dot_product(r, r))
			
			beta = (alpha / xi) * dot_product(rr, r) / dot_product(rr,r_)
			
			its = its + 1
			print*, 'BiCGSTAB: iteration = ', its, ' evals = ', evals, ' residual = ', normr/normb
			write(1, '(E12.5)') normr/normb
			write(2, '(I0)') evals
		end do BiCGSTAB_loop
		
	    close(1)
	    close(2)
	end subroutine BiCGSTAB
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	subroutine print_vector(vec, n, q)
	! prints the first q entries of the n-dimensional vector vec
		implicit none
		double precision, dimension(n), intent(in) :: vec
		integer, intent(in) :: n, q
		integer :: k

		do k = 1, q
			print*, ' ', vec(k)
		end do
	end subroutine print_vector
  
  
end module mod_IterativeSolvers_test
