using PyPlot

function plot_solution(n::Int64,sol::Array{Float64,1},desiredstate::Function,alpha::Float64)
	# Simple scatter-plot of the nodes 
	# Input: solution sol=(y,p,lambda)
	#        Desired state: y_d
	#	 L2-regularization alpha
	
	n_block::Int64 = Int64(length(sol)/3)

	dof_map = load("data"*string(n)*".jld2", "dof_map")
	xval::Array{Float64,1} = zeros(n_block)
	yval::Array{Float64,1} = zeros(n_block)

	stateval::Array{Float64,1} = zeros(n_block)
	controlval::Array{Float64,1} = zeros(n_block)
	adjointval::Array{Float64,1} = zeros(n_block)
	lagrangeval::Array{Float64,1} = zeros(n_block)
	y_d::Array{Float64,1} = zeros(n_block)

	# Separating all variables into vectors
	for i=1:n_block
		xval[i] = dof_map[i,:][1]
		yval[i] = dof_map[i,:][2]
		y_d[i] = desiredstate(xval[i],yval[i])
		stateval[i] = sol[i]
		adjointval[i] = sol[i+n_block]
		lagrangeval[i] = sol[i+2*n_block]
		controlval[i] = (1. /alpha)*(adjointval[i] - lagrangeval[i] )
	end

	# Plot state
	figure("1")
	scatter3D(xval, yval, stateval, s=0.1, c="k", label="State")
	scatter3D(xval, yval, y_d, s=0.1, c="b", alpha=0.1, label="Desired state")
	title("State")
	xlabel("x_1")
	ylabel("x_2")
	zlabel("y")
	legend(loc="best")

	# Plot Control
	figure("2")
	scatter3D(xval, yval, controlval, s=0.1, c="k")
	title("Control")
	xlabel("x_1")
	ylabel("x_2")
	zlabel("u")
end
