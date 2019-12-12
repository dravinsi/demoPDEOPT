function GCR(x,A,b,P,tol,maxit,outstring)
	# Implements the Generalized Conjugate Residual method
	# ---------------------------------
	# Input: x = initial guess
	#        A = system matrix 
	#	 b = RHS
	# 	 P = (left) preconditioner
	#	 tol = tolerance
	# 	 maxit = maximum iterations
	#	 outstring = string identifier for output
	# Output: x = solution
	#	  iter = number of iterations
	# ---------------------------------
	res::Array{Float64,1} = b - A*x
	#nb = norm(b)
	pnb = norm( P \ b)
	#relres = 1
	Prelres = 1 # for Preconditioned residual
	search = P \ res # initialize
	Hn = 0. 	# Initialize
	H = zeros(length(x)) # Initialize so the GC does not delete it
	Hd = zeros(length(x)) 
	it::Int64 = 1
	while( (Prelres>tol) && (iter < maxit) )
		Asearch = A*search
		if( it > 1 )
			y = (H'*Asearch) ./ Hn
			Asearch = Asearch - H*y
			H = [H Asearch] 
			nAsearch = norm(Asearch)^2
			Hn = [Hn; nAsearch]
			search = search - Hd*y
			Hd = [Hd search]
		else
			H = Asearch
			nAsearch = norm(Asearch)^2
			Hn = nAsearch
			Hd = search
		end
 		c = Asearch'*res/nAsearch
		x = x + c*search
		res = res - c*Asearch
		nres = norm(res)
		search = P \ res
		Prelres = norm(search)/pnb # Preconditined relative tolerance
		#relres =nres/nb  # Relative tolerance
		it = it + 1
	end
	if(Prelres<tol)
		println(outstring,"GCR converged in ", iter," iterations.")
	else
		println(outstring,"GCR terminated without convergence in ", iter, " iterations.")
	end
	return x, it
end
