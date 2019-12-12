using Makie, AbstractPlotting

function plot_solution(n::Int64,sol::Array{Float64,1},desiredstate::Function,alpha::Float64)
	# Simple scatter-plot of the nodes 
	# Input: solution sol=(y,p,lambda)
	#        Desired state: y_d
	#	 L2-regularization alpha

	n_block::Int64 = Int64(length(sol)/3)
	n::Int64 = Int64( sqrt(n_block)-1 ) 

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
	scene1 = Scene()
	scatter!(scene1 ,xval, yval, stateval, markersize=0.005, color = stateval ,colormap = :viridis)
	scale!(scene1, 1, 1, 1/(maximum(abs.(y_d))*2 ) )  
	stateaxis  = scene1[Axis]
	stateaxis[:names, :axisnames] = ("x_1", "x_2", "y")

	# Plot Control
	scene2 = Scene()
	s = scatter!(scene2, xval, yval, controlval, markersize=0.005, color= controlval, colormap=:viridis )#color= :Black)
	scale!(scene2, 1, 1, 1/(maximum(abs.(controlval))*2 ) )  
	controlaxis  = scene2[Axis]
	controlaxis[:names, :axisnames] = ("x_1", "x_2", "u")
	display(scene2)

	return scene1,scene2
end
