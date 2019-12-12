

#####################################
#     Define Auxillary functions   
#####################################

function readin(n::Int64)
	println("Begin read-in and formatting for n = ",n,"...")
	@load "data"*string(n)*".jld2" K M negator dof_map 
	println("Finsihes read-in and formatting for n = ",n)
	println("--------------------------------------------")
	return K::SparseMatrixCSC{Float64,Int64}, M, negator, dof_map
end

function interpolate_mesh(sol_coarse::Array{Float64,1},n_coarse::Int64,n_fine::Int64)
	# Input: coarse solution: sol_coarse = (y,p,lambda)
	#	 coarse DoF-map: dof_coarse ordered by x-coordinate first, y-coordinate second with stepsize h
	#	 Fine   DoF-map: dof_fine   ordered by x-coordinate first, y-coordinate second with stepsize h/2
	#
	# Assumptions: both discretizations are uniform square grids in the unit square with seq. numbering
	#
	# Output: Linear interpolation of the coarse solution on the fine mesh. 
	dof_fine = load("data"*string(n_fine)*".jld2", "dof_map")
	dof_coarse = load("data"*string(n_coarse)*".jld2", "dof_map")
	#Find stepsizes:
	h_coarse = dof_coarse[2,:][1] - dof_coarse[1,:][1]
	h_fine = h_coarse/2

	# Find number of nodes per size (assuming unit square)
	side_coarse = Int64(1/h_coarse + 1)
	side_fine =  Int64(1/h_fine + 1)

	# With which we can calculate the number of points
	n_block_coarse = side_coarse^2
	n_block_fine = side_fine^2

	# This allows us to allocate memory for the interpolated solution
	y_fine = zeros(side_fine^2) # State
	p_fine = zeros(side_fine^2) # adjoint
	l_fine = zeros(side_fine^2) # multiplier

	# Now loop over the coarse DoF-handler
	# This loop will transfer the known values from the coarse
	# grid to the fine grid and interpolate the coarse points between them

	# to map the coarse grid to the fine we need to keep track of the indices
	indexadjust=0 # This variable will ensure dof_coarse[i]=dof_fine[i+indexadjust]
	for i = 1:n_block_coarse
		#println("i = ",i)
		# Current x,y values:
		x_coarse = dof_coarse[i,:][1]
		y_coarse = dof_coarse[i,:][2]
		# Extract current values from coarse solution
		#statetemp = sol_coarse[i]
		#adjointtemp = sol_coarse[i+n_block_coarse]
		#lagrangetemp = sol_coarse[i+2*n_block_coarse]


		# In each step we will transfer the known value and
		# construct at most 3 points in the interpolated solution
		# If possible we interpolate to the right, up and right/up diagonally
		# ----------
		# Thus, in each step we need to check wether we are in the right conrner or top of the coarse grid
		# if neither: we send 4 values into the coarse grid and update indexajust+=1
		# if at the right: we only interpolate up and update indexajust += side_fine as we move to the next row
		# If we are at the top and not at the right: we interpolate right and update indexadjust +=1

		# Check if we are at the top (in course grid), i.e. y=1, we use h_fine  
		if( y_coarse + h_fine > 1)
			# Now we are at the top, check if we are also on the right i.e. x=1
			if( x_coarse + h_fine > 1)
				# This is the final point we need only transfer the final value.
				y_fine[i+indexadjust] = sol_coarse[i]
				p_fine[i+indexadjust] = sol_coarse[i+n_block_coarse]
				l_fine[i+indexadjust] = sol_coarse[i+2*n_block_coarse]
			else
				# We are at the top but not at the far right so we need oo interpolate 
				# one step to the right and set indexadjust += 1
				
				# Transfer known value:
				y_fine[i+indexadjust] = sol_coarse[i]
				p_fine[i+indexadjust] = sol_coarse[i+n_block_coarse]
				l_fine[i+indexadjust] = sol_coarse[i+2*n_block_coarse]

				# Interpolate one step to the right:
				y_fine[i+indexadjust+1] = (sol_coarse[i] + sol_coarse[i+1])/2
				p_fine[i+indexadjust+1] = (sol_coarse[i+n_block_coarse] + sol_coarse[i+1+n_block_coarse])/2
				l_fine[i+indexadjust+1] = (sol_coarse[i+2*n_block_coarse] + sol_coarse[i+1+2*n_block_coarse])/2

				#Update indexadjust +=1
				indexadjust = 1 + indexadjust
			end

		# If not at the top check if we are at the right:
		elseif( x_coarse + h_fine > 1 )
			# Transfer known value:
			y_fine[i+indexadjust] = sol_coarse[i]
			p_fine[i+indexadjust] = sol_coarse[i+n_block_coarse]
			l_fine[i+indexadjust] = sol_coarse[i+2*n_block_coarse]

			# Interpolate up
			y_fine[i+indexadjust+side_fine] = (sol_coarse[i]+sol_coarse[i+side_coarse])/2 
			p_fine[i+indexadjust+side_fine] = (sol_coarse[i+n_block_coarse]+sol_coarse[i+side_coarse+n_block_coarse])/2 
			l_fine[i+indexadjust+side_fine] = (sol_coarse[i+2*n_block_coarse]+sol_coarse[i+side_coarse+2*n_block_coarse])/2 

			#Update indexadjust += n_block_fine
			indexadjust = indexadjust + side_fine	

		# We are neither at the nor the right, i.e. this is a "normal" point 
		#and we should transfer 4 values to the fine grid
		else
			# Transfer known value:
			y_fine[i+indexadjust] = sol_coarse[i]
			p_fine[i+indexadjust] = sol_coarse[i+n_block_coarse]
			l_fine[i+indexadjust] = sol_coarse[i+2*n_block_coarse]

			# Interpolate up:
			y_fine[i+indexadjust+side_fine] = (sol_coarse[i]+sol_coarse[i+side_coarse])/2 
			p_fine[i+indexadjust+side_fine] = (sol_coarse[i+n_block_coarse]+sol_coarse[i+side_coarse+n_block_coarse])/2 
			l_fine[i+indexadjust+side_fine] = (sol_coarse[i+2*n_block_coarse]+sol_coarse[i+side_coarse+2*n_block_coarse])/2 

			# Interpolate one step to the right:
			y_fine[i+indexadjust+1] = (sol_coarse[i] + sol_coarse[i+1])/2
			p_fine[i+indexadjust+1] = (sol_coarse[i+n_block_coarse] + sol_coarse[i+1+n_block_coarse])/2
			l_fine[i+indexadjust+1] = (sol_coarse[i+2*n_block_coarse] + sol_coarse[i+1+2*n_block_coarse])/2

			# Interpolate diagonally up-right:
			y_fine[i+indexadjust+side_fine+1] = ( sol_coarse[i] + sol_coarse[i+1] + sol_coarse[i+side_coarse] + sol_coarse[i+side_coarse+1] )/4
			p_fine[i+indexadjust+side_fine+1] = ( sol_coarse[i+n_block_coarse] + sol_coarse[i+1+n_block_coarse] + sol_coarse[i+side_coarse+n_block_coarse] + sol_coarse[i+side_coarse+1+n_block_coarse] )/4
			l_fine[i+indexadjust+side_fine+1] = ( sol_coarse[i+2*n_block_coarse] + sol_coarse[i+1+2*n_block_coarse] + sol_coarse[i+side_coarse+2*n_block_coarse] + sol_coarse[i+side_coarse+1+2*n_block_coarse] )/4
			#Update indexadjust +=1
			indexadjust = 1 + indexadjust

		end

	end
	return [y_fine; p_fine ; l_fine]
end

function diag_sqrt_inv(A)
	# Returns inv(sqrt(A)) for a diagonal matrix A
	n_block,m_block = size(A)
	A_sqrt_inv = zeros(n_block)
	for i=1:n_block
		A_sqrt_inv[i]=1/sqrt(A[i,i])
	end
	return spdiagm(0 => A_sqrt_inv)
end


function diag_sqrt(A)
	# Returns (sqrt(A)) for a diagonal matrix A
	n_block,m_block = size(A)
	A_sqrt = zeros(n_block)
	for i=1:n_block
		A_sqrt[i]=sqrt(A[i,i])
	end
	return spdiagm(0 => A_sqrt)
end

function desired_state(func::Function,dof_map::Array{Float64,2})
	# Evaluates desired state on Dofs and returns vector
	n_block::Int64 = size(dof_map)[1]
	yd::Array{Float64,1} = zeros(n_block)
	for i=1:n_block
		xtemp,ytemp=dof_map[i,:]
		yd[i] = func(xtemp,ytemp)
	end
	return yd
end

#########################################
########           END           ########
#########################################










