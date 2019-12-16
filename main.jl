# Demo-code for PDE-OPT optimal control written by Ivo Dravins
# Given as is without any gurantees.
using LinearAlgebra, SparseArrays, IterativeSolvers 
using AlgebraicMultigrid, JLD2, FileIO 
#using Makie
# Load auxilliary functions
include("aux.jl")
#####################################
#        DEFINE PRECOND STRUCTS        
#####################################
struct LT_precond
	# This is the structure that holds the
	# information composing the preconditioner
	# |  A     -K^T |     | 1    0  | | A  -K^T |
	# |(K+B)  (K+B) |  =  | 0  (K+B)| | 1    1  |
	Block1 # = P.K + P.B*P.S_e
	Block2 # = P.K' + P.S_e*P.A
	K_tilde # Stifness
	S_e # Epsilon scaling 
	p1 # AMG1
	p2 # AMG2 
	tol_cg # Tolerance
end

function Base.:\(P::LT_precond, b::Array{Float64,1})
	# Defines backslash operator for preconditioner
	# --------------------------
	#P_exact = [P.A -(P.K') ; (P.K+P.B) (P.K+P.B) ] 
	#return P_exact\b
	# We break it down in parts: P_LT = P_1 P_2
	# --------------------------
	n_block,m_block = size(P.K_tilde)
	x = zeros(2*n_block) # allocate memory
	# Solve P1 Accouting for epsilon-scaling
	b_top = b[1:n_block]
	#b_bot = (P.K + P.B*P.S_e) \ b[n_block+1:2*n_block] # Direct solve reference	
	b_bot = cg(P.Block1, b[n_block+1:2*n_block], Pl = P.p1 , tol = P.tol_cg ) # AMG precondtioned CG
	b_bot = P.S_e*b_bot # Recover substituion
	# Solve P2 Accouting for epsilon-scaling
	#x[1:n_block] = (P.K' + P.S_e*P.A) \ (P.S_e*( b_top + P.K_tilde'*b_bot)) # Direct solve	reference
	x[1:n_block] = cg(P.Block2,  (P.S_e*( b_top + P.K_tilde'*b_bot)) , Pl = P.p2 , tol = P.tol_cg ) # AMG precondtioned CG
	x[1+n_block:2*n_block] = -x[1:n_block]+b_bot
	return x
end
#####################################
#        END PRECOND STRUCTS       
#####################################
function Generate_activeset(sol::Array{Float64,1},alpha::Float64,beta::Function,a::Function,b::Function,ymin::Function,ymax::Function,c::Float64,epsilon::Float64,dof_map::Array{Float64,2})
	# Input: Solution to KKT-system: sol=(y,p,lambda)
	#	 L2-regularizaation parameter: alpha
	#	 L1-regularization parameter: beta
	#	 Control box constraints: [a,b]
	#	 State box-constraints: [ymin,ymax]
	# 	 Postive constant: c ( Assumes c=c_1=c_2 in article notation)
	#	 Moreau-Yosida regularization: epsilon
	#	 DoF = DOF handler
	#
	# Output: Resulting active sets and conditional matricies/vectors
	n_block::Int64 = Int64(length(sol)/3)
	#These two can be changed but these values seem stabe with c = 1/alpha, they correspond to c_1 and c_2 in the article. 
	sigma1::Float64 = c #
	sigma2::Float64 = c #
	# Start with no active points:
	activepoints_control::Array{Int64,1} =  Array{Int64}(zeros(n_block))
	activepoints_state::Array{Int64,1} =  Array{Int64}(zeros(n_block))
	#Initialize scalar memory
	ytemp::Float64 = 0
	ptemp::Float64 = 0
	ltemp::Float64 = 0
	xcoord::Float64 = 0
	ycoord::Float64 = 0
	h = dof_map[2,:][1]-dof_map[1,:][1] # Stepsize
	# locate memory
	C_y_val::Array{Float64,1} = zeros(n_block)
	C_p_val::Array{Float64,1} = zeros(n_block)  
	# RHS vectors
	c_rhs_y::Array{Float64,1} = zeros(n_block)  
	c_rhs_c::Array{Float64,1} = zeros(n_block) 
	# Scalars
	cval_p::Float64 = 0
	cval_y::Float64 = 0
	#Iterate over all points in the domain:
	for i =1:n_block
		xcoord,ycoord=dof_map[i,:]
		ytemp = sol[i] # State
		ptemp = sol[i+n_block] # Adjoint
		ltemp = sol[i+2*n_block] # Multiplier
		cval_p = 0
        	cval_y = 0
		#------------------------------------------------
		#Check condition Lambda_1+
		if( (1/alpha)*(ptemp-ltemp) + sigma1*(ltemp - beta(xcoord,ycoord)) > 0 )
			#Check condition Lambda_2+
			if(  sigma2*((1/alpha)*(ptemp-ltemp) - b(xcoord,ycoord)) + sigma1*(ltemp - beta(xcoord,ycoord)) > 0.)
				activepoints_control[i] = 2 #condition Lambda_2+
				cval_p = 0	
				c_rhs_c[i] = b(xcoord,ycoord)
			else		
				activepoints_control[i] = 1 #condition Lambda_1+
				cval_p = 1	
				c_rhs_c[i] = - beta(xcoord,ycoord)/alpha
			end
		#Check condition Lambda_1-
		elseif( (1/alpha)*(ptemp-ltemp) + sigma1*(ltemp + beta(xcoord,ycoord)) < 0 )
			#Check condition Lambda_2-
			if(  sigma2*((1/alpha)*(ptemp-ltemp) - a(xcoord,ycoord)) + sigma1*(ltemp + beta(xcoord,ycoord)) < 0 )
				activepoints_control[i] = -2
				cval_p = 0	
				c_rhs_c[i] = a(xcoord,ycoord)
			else
				activepoints_control[i] = -1
				cval_p =  1	
				c_rhs_c[i] = beta(xcoord,ycoord)/alpha
			end
		end
		#Check Moreau Yosida condition +
        	if( ytemp - ymax(xcoord,ycoord)  > 0. )
			activepoints_state[i] = 1
			cval_y = cval_y + 1
			c_rhs_y[i] = c_rhs_y[i] + (1. /epsilon)*ymax(xcoord,ycoord)
		end 
		#Check Moreau Yosida condition -
		if( ytemp - ymin(xcoord,ycoord)  < 0. )
			activepoints_state[i] = -1
			cval_y = cval_y + 1
			c_rhs_y[i] = c_rhs_y[i] + (1. /epsilon)*ymin(xcoord,ycoord)
		end
		#------------------------------------------------
		# Save terminal values to vector:
		C_y_val[i] = cval_y
		C_p_val[i] = cval_p 
	end
	# Consttruct diagonal matricies:
	C_y = spdiagm(0 => C_y_val)
	C_p = spdiagm(0 => C_p_val)
	return C_y::SparseMatrixCSC{Float64,Int64}, C_p::SparseMatrixCSC{Float64,Int64}, c_rhs_y::Array{Float64,1}, c_rhs_c::Array{Float64,1}, activepoints_control::Array{Int64,1}, activepoints_state::Array{Int64,1}
end

function GCR(x::Array{Float64,1},A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},P,tol::Float64,maxit::Int64,outstring::String)
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
	#nb::Float64 = norm(b)
	pnb::Float64 = norm( P \ b)
	#relres::Float64 = 1.
	Prelres::Float64 = 1. # for Preconditioned residual
	search::Array{Float64,1} = P \ res # initialize
	Asearch::Array{Float64,1} = zeros(length(x)) 
	it::Int64 = 0
	# Initialize variables which are not type-stable
	Hn = zeros(length(x)) 	
	H = zeros(length(x)) 
	Hd = zeros(length(x)) 
	while( (Prelres>tol) && (it < maxit) )
		Asearch = A*search
		it = it + 1
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
	end
	if(Prelres<tol)
		println(outstring,"GCR converged in ", it," iterations.")
	else
		println(outstring,"GCR terminated without convergence in ", it, " iterations.")
	end
	return x, it
end

#########################################
######## BEGIN REDUCED FUNCTIONS
#########################################
function recover_lagrange(solred::Array{Float64,1},alpha::Float64,beta::Function,a::Function,b::Function,activepoints_control::Array{Int64,1},dof_map::Array{Float64,2})
	# Input: Reduced solution solred = (y,p)^T and active points corresponding Lambda_1+- and Lambda_2+-
	# 	 and problem parameters alpha,beta,a,b
	# Output: Recovers and returns associated lambda values
	n_block::Int64 = length(activepoints_control)
	lagrangeval::Array{Float64,1} = zeros(n_block) 	#Allocate memory for the Lagrange values
	for i=1:n_block
		xcoord::Float64,ycoord::Float64 = dof_map[i,:]
		ptemp::Float64 = solred[i+n_block] # Adjoint
		if(activepoints_control[i]>0)
			if(activepoints_control[i]>1)
			# Both lambda1+ and lambda2+ conditions active
				lagrangeval[i] = ptemp - alpha*b(xcoord,ycoord)
			else
			# Only lambda1+ active
				lagrangeval[i] = beta(xcoord,ycoord)
			end
		elseif(activepoints_control[i]<0)
			if(activepoints_control[i]<-1)
			# Both lambda1- and lambda2- conditions active			
				lagrangeval[i] = ptemp - alpha*a(xcoord,ycoord)
			else
			# Only lambda1- active
				lagrangeval[i] = -beta(xcoord,ycoord)
			end
		else
			# No conditions active
			lagrangeval[i] = ptemp
		end
	end
	return lagrangeval::Array{Float64,1}
end

function Reducednewton(xnew::Array{Float64,1},n::Int64,alpha::Float64,beta::Function,a::Function,b::Function,ymin::Function,ymax::Function,c::Float64,epsilon::Float64,func_desired::Function,lintol,CGtol)
	# Input: inital guess x0
	# 	 n to identify which FEM files to load
	#	 L2-regularizaation parameter: alpha
	#	 L1-regularization parameter: beta
	#	 Control box constraints: [a,b]
	#	 State box-constraints: [ymin,ymax]
	#
	# Output: Performs Netwon iterations and returns results
	# Read in files:
	K::SparseMatrixCSC{Float64,Int64},M::SparseMatrixCSC{Float64,Int64},negator::SparseMatrixCSC{Float64,Int64},dof_map::Array{Float64,2} = readin(n)
	y_d = desired_state(func_desired,dof_map)
	h = dof_map[2,:][1]-dof_map[1,:][1]
	println("h = ", h)
	# Construct auxillary matricies:
	n_block::Int64, m_block::Int64 = size(K)
	I_block::SparseMatrixCSC{Float64,Int64} = sparse(I, n_block, m_block)
	boundaryzeros::SparseMatrixCSC{Float64,Int64} = I_block + negator # This is a unity matrix which has zeros on the rows where we have Dirichlet conditions
	boundaryones::SparseMatrixCSC{Float64,Int64} = I_block - boundaryzeros # This is a diagonal matrix which has ones on the rows where we have Dirichlet conditions, 0 elsewhere
	O_block = spzeros(n_block,m_block)
	println("Block dimension = ",n_block)
	println("dim of system = ",2*n_block)
	# Impose boundary conditions
	M = boundaryzeros * M # This zeros the rows of M
	K = boundaryzeros * K # This zeros the rows of K
	Kc::SparseMatrixCSC{Float64,Int64} = K * boundaryones # K compensator to fix RHS after symmetrization
	K = sqrt(alpha)*K * boundaryzeros + boundaryones
	# Initialize iteration paramters
	flag::Bool = true
	itt::Int64 = 1
	maxit::Int64 = 40
	tol::Float64 = 1e-6 # 2-norm change of soltuion convergence, currently inactive
	tolpoints::Float64 = 1 # Number of points changing sets, currently active
	errornew::Float64 = 1e6 # Initialize by something ridiclous
	xnewred::Array{Float64,1} = xnew[1:Int64(2*n_block)]
	s1c::Float64 = 0
	s2c::Float64 = 0
	s3c::Float64 = 0
	totaliter::Int64 = 0
	xnewoldred = zeros(2*n_block)
	# Initial matrix xonstruction
 	PI_y::SparseMatrixCSC{Float64,Int64}, PI_p::SparseMatrixCSC{Float64,Int64}, c_rhs_y::Array{Float64,1}, c_rhs_c::Array{Float64,1}, activepoints_control::Array{Int64,1}, activepoints_state::Array{Int64,1} = Generate_activeset(xnew,alpha,beta,a,b,ymin,ymax,c,epsilon,dof_map)
	while(flag)
		println("Running iteration ", itt )
		xold::Array{Float64,1} = xnew # Update
		errorold::Float64 = errornew
		epsilonsqrt::SparseMatrixCSC{Float64,Int64} = diag_sqrt((I_block+(1/epsilon)*boundaryzeros*PI_y))
		epsilonsqrtinv::SparseMatrixCSC{Float64,Int64} = diag_sqrt_inv( (I_block+(1/epsilon)*PI_y) ) 
		K_tilde = K*epsilonsqrtinv 
		B = (PI_p*M)
		system::SparseMatrixCSC{Float64,Int64} = [M -K_tilde' ; K_tilde B]
		rhs::Array{Float64,1} = [ epsilonsqrtinv*(M*(c_rhs_y+y_d)) ;  sqrt(alpha)*(M*c_rhs_c-Kc*y_d) + boundaryones*y_d  ] 
		# Set up algebraic multigrid for this step
		ml1 = ruge_stuben(K+B*epsilonsqrt)
		p1  = aspreconditioner(ml1)
		ml2 = ruge_stuben(K'+epsilonsqrt*M) 
		p2 = aspreconditioner(ml2)
		P = LT_precond(K+B*epsilonsqrt,K'+epsilonsqrt*M,K_tilde,epsilonsqrt,p1,p2,CGtol) # Construct preconditioner
		if(itt == 1)
			xnewred, currentiter = GCR(zeros(2*n_block),system, rhs,P,lintol,50,"    GCR: ") # zero initial guess in first step
		else
			xnewred, currentiter = GCR(xnewoldred,system, rhs,P,lintol,50,"    GCR: ") # Previous sol as initial guess in later steps
		end
		totaliter = totaliter + currentiter # Linear iteration counter
		xnewoldred = xnewred
		xnewred[1:n_block] = epsilonsqrtinv*xnewred[1:n_block]
		xnewred[n_block+1:2*n_block] = -sqrt(alpha)*xnewred[n_block+1:2*n_block]  # Resubstitue p = -sqrt(alpha)*p_hat
		lval::Array{Float64,1} = recover_lagrange(xnewred,alpha,beta,a,b,activepoints_control,dof_map) # recover Lagrangian
		xnew = [xnewred ; lval]
		activepoints_control_old = activepoints_control
		activepoints_state_old = activepoints_state
		# Generate sets and matricies for next
 		PI_y, PI_p, c_rhs_y, c_rhs_c, activepoints_control, activepoints_state = Generate_activeset(xnew,alpha,beta,a,b,ymin,ymax,c,epsilon,dof_map)
		println("-------------")
		# Check termination criteria
		s1c =  norm(activepoints_control-activepoints_control_old,1)
		s2c =  norm(activepoints_state-activepoints_state_old,1)
		println(" Active set control change ", s1c )
		println(" Active set state change ", s2c )
		x2norm = norm( xnew - xold )
		maxpointchange = maximum( abs.(xnew - xold)) 
		println( "||xnew-xold|| = " , x2norm )
		# Check stopping criteria
		if( s1c < tolpoints && s2c < tolpoints)
			println(" Set convergence in ", itt, " iterations.")
			println("Average linear iterations = ", totaliter/itt )
			println("Total linear iterations = ", totaliter )
			flag=false
			return xnew::Array{Float64,1}
		end
		if( x2norm < -tol) # Currently inactive, remove - to activiate
			println("Solution convergence in ", itt, " iterations.")
			println("Average linear iterations = ", totaliter/itt )
			println("Total linear iterations = ", totaliter )
			flag=false
			return xnew::Array{Float64,1}
		end
		if( maxpointchange < -h^2) # Currently inactive, remove - to activiate
			println("Max point movement convergence in ", itt, " iterations.")
			println("Average linear iterations = ", totaliter/itt )
			println("Total linear iterations = ", totaliter )
			flag=false
			return xnew::Array{Float64,1}
		end
		if(itt>maxit)
			println("Maxit reached, terminating...")
			println("Average linear iterations = ", totaliter/itt )
			println("Total linear iterations = ", totaliter )
			flag=false
			return xnew::Array{Float64,1}
		end	
		itt = itt + 1 # Update iteration counter
		#--------------------------------------------------------
	end
end
#---------------------------------------#
#        END REDUCED FUNCTIONS          #
#---------------------------------------#

##############################
#       CONTROL PANEL        #
##############################
# Initial coarse grid
n = 32 
alpha= 1e-6 #L2-regulization 
c =  1. /alpha # value assigned to both c_1 and c_2 (in article notation)
epsilon = alpha^(1/4) # Moreau-Yosida reg.
lintol = 1e-6 # Tolerance for linear solver (GCR)
CGtol = 1e-6 # Tolerance for the block solvers in the preconditioner (CG)

function beta(x::Float64,y::Float64)
	# L1 regularization (can be a function)
	# Here x and y denotes spacial coordinates. 
	return 1e-4
end

a(x::Float64,y::Float64) = -30 #lower limit control function
b(x::Float64,y::Float64) = 30  #upper limit control function

ymin(x::Float64,y::Float64) = -0.2 # lower limit state
ymax(x::Float64,y::Float64) = 0.3 # upper limit state 

function desiredstate(x::Float64,y::Float64)
	# Defines the desired state, here x and y denotes spacial coordinates. 
	# Note: Boundary conditions will also be set by the desired state. 
	# ------
	return sin(2*pi*x)*sin(2*pi*y)*exp(2*x)/6. # PROBLEM 1
	#return -exp( abs(x-0.5)+abs(y-0.5) ) # PROBLEM 2
	#return abs(sin(2*pi*x)*sin(2*pi*y)) # PROBLEM 3
end
##############################
#          SOLVING           #
##############################
# Initial guess on coarsest grid (constant ones)
dim = (n+1)^2 
x0 = ones(3*dim)*1.

# First solve
sol =  Reducednewton(x0,n,alpha,beta,a,b,ymin,ymax,c,epsilon,desiredstate,lintol,CGtol)

# Interpolate solution and solve again
sol = interpolate_mesh(sol,n,n*2)
n = n*2
sol = Reducednewton(sol,n,alpha,beta,a,b,ymin,ymax,c,epsilon,desiredstate,lintol,CGtol)

# Interpolate solution and solve again
x0 = interpolate_mesh(sol,n,n*2)
n = n*2
sol = Reducednewton(x0,n,alpha,beta,a,b,ymin,ymax,c,epsilon,desiredstate,lintol,CGtol)

# Interpolate solution and solve again
#x0 = interpolate_mesh(sol,n,n*2)
#n = n*2
#sol = Reducednewton(x0,n,alpha,beta,a,b,ymin,ymax,c,epsilon,desiredstate,lintol,CGtol)

# Interpolate solution and solve again
#x0 = interpolate_mesh(sol,n,n*2)
#n = n*2
#sol = Reducednewton(x0,n,alpha,beta,a,b,ymin,ymax,c,epsilon,desiredstate,lintol,CGtol)



# Interpolate solution and solve again, we time the solve on the finest grid
x0 = interpolate_mesh(sol,n,n*2)
n = n*2
start =  time()
sol = Reducednewton(x0,n,alpha,beta,a,b,ymin,ymax,c,epsilon,desiredstate,lintol,CGtol)
elapsed = time() - start
println("Time to solve on the fine grid: ", round(elapsed,digits=2), " seconds.")
##############################
#          Plotting          #
##############################
# We give two simple scatter plot functions using different libraries
# They should be enough to help any user create their preferred plots

# Uncomment to plot with Pyplot:-----------------------
#include("Pyplot_Plot.jl")
#plot_solution(n,sol,desiredstate,alpha)

# Comment: This can be quite sluggish but at least it visualizes the result.

# Uncomment to plot with Makie:-----------------------
#include("Makie_Plot.jl")
#state,control=plot_solution(n, sol,desiredstate,alpha)
#update_cam!(control, FRect(Vec3f0(0),Vec3f0(1)) )

# Comment: Control which plot to view using display(state), display(control)
#	   update camera with update_cam!(control, FRect(Vec3f0(0),Vec3f0(1)) )
#-----------------------------------------------------
println("Program ended")
