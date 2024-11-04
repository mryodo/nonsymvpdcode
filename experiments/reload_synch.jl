using LinearAlgebra, Arpack
using SparseArrays
using StatsBase
using DifferentialEquations
using Interpolations
using Peaks
using Polynomials
using Statistics

using DataInterpolations
using BenchmarkTools

using Plots, ColorSchemes, LaTeXStrings, StatsPlots
pgfplotsx()
theme(:dao)
Plots.scalefontsizes(2.25)

rwth = Dict(
      "blue" => [colorant"#00549F", colorant"#407FB7", colorant"#8EBAE5", colorant"#C7DDF2",  colorant"#E8F1FA",],
      "black" => [ colorant"#000000", colorant"#646567", colorant"#9C9E9F", colorant"#CFD1D2", colorant"#ECEDED",],
      "magenta" => [ colorant"#E30066", colorant"#E96088", colorant"#F19EB1", colorant"#F9D2DA", colorant"#FDEEF0",],
      "yellow" => [colorant"#FFED00", colorant"#FFF055", colorant"#FFF59B", colorant"#FFFAD1", colorant"#FFFDEE",],
      "petrol" => [colorant"#006165", colorant"#2D7F83", colorant"#7DA4A7", colorant"#BFD0D1", colorant"#E6ECEC",],
      "turq" => [colorant"#0098A1", colorant"#00B1B7", colorant"#89CCCF", colorant"#CAE7E7", colorant"#EBF6F6", ],
      "green" => [ colorant"#57AB27", colorant"#8DC060", colorant"#B8D698", colorant"#DDEBCE", colorant"#F2F7EC",],
      "maigreen" => [ colorant"#BDCD00", colorant"#D0D95C", colorant"#E0E69A", colorant"#F0F3D0", colorant"#F9FAED",],
      "orange" => [colorant"#F6A800", colorant"#FABE50", colorant"#FDD48F", colorant"#FEEAC9", colorant"#FFF7EA", ],
      "red" => [colorant"#CC071E", colorant"#D85C41", colorant"#E69679", colorant"#F3CDBB", colorant"#FAEBE3",],
      "bord" => [colorant"#A11035", colorant"#B65256", colorant"#E5C5C0", colorant"#F0F3D0", colorant"#F5E8E5", ],
      "violet" => [colorant"#612158", colorant"#834E75", colorant"#A8859E", colorant"#D2C0CD", colorant"#EDE5EA",],
      "lilac" => [colorant"#7A6FAC", colorant"#9B91C1", colorant"#BCB5D7", colorant"#DEDAEB", colorant"#F2F0F7"]
)
custom_cmap = cgrad( [ rwth["blue"][4], rwth["blue"][1] ], scale = :log10 )
custom_cmap2 = cgrad( [ rwth["magenta"][4], rwth["magenta"][1] ], scale = :log10 )
custom_cmap3 = cgrad( [ rwth["petrol"][4], rwth["petrol"][1] ], scale = :log10 ) 

##########################################################
#                      VpD solvers                       #
##########################################################


function f(du, u, p, t)
      Δω, μ1, μ2=p;
      du[1]=u[2];
      du[2]=(1-u[1]^2)*u[2]-(1-Δω)*u[1]-μ1*(u[2]-u[4]);
      du[3]=u[4];
      du[4]=(1-u[3]^2)*u[4]-(1+Δω)*u[3]-μ2*(u[4]-u[2]);
end

function vpdSolve(problem::ODEProblem, interp::Bool, mult::Int64)    
      sol=solve(problem, Tsit5(); reltol=1e-12, abstol=1e-12);
      x=reduce(hcat, sol.u)[1, :];
      dx=reduce(hcat, sol.u)[2, :];
      y=reduce(hcat, sol.u)[3, :];
      dy=reduce(hcat, sol.u)[4, :];
      t=sol.t;
      tspan=problem.tspan;
      Δt = minimum( diff( t ) );
      if Δt < 1e-6
          Δt = 1e-6
      end
      num = Int(floor(( tspan[2] - tspan[1] ) / Δt));
      if interp 
          t_new=range(tspan[1], tspan[2], num);
          #itp=LinearInterpolation(t, x); 
          #x=itp(t_new);
          #itp=LinearInterpolation(t, dx); dx=itp(t_new);
          #itp=LinearInterpolation(t, y); y=itp(t_new);
          #itp=LinearInterpolation(t, dy); dy=itp(t_new);
          itp = CubicSpline(x, t); x = itp( t_new )
          itp = CubicSpline(dx, t); dx = itp( t_new )
          itp = CubicSpline(y, t); y = itp( t_new )
          itp = CubicSpline(dy, t); dy = itp( t_new )
          t=t_new;
          return t, x, dx, y, dy
      else
          return t, x, dx, y, dy
      end
end


function cleanMinima( x, y, pksx, pksy; window = 20 )
      pksx = pksx[ x[pksx] .< 0  ]
      pksy = pksy[ y[pksy] .< 0  ]
      tmpx = x[pksx][ end - 2*window : end ]
      indx1 = pksx[ end - 2*window : end ][ ( tmpx .- mean(tmpx) ) .> 0 ]
      indx2 = pksx[ end - 2*window : end ][ ( tmpx .- mean(tmpx) ) .< 0 ]

      tmpy = y[pksy][ end - 2*window : end ]
      indy1 = pksy[ end - 2*window : end ][ ( tmpy .- mean(tmpy) ) .> 0 ]
      indy2 = pksy[ end - 2*window : end ][ ( tmpy .- mean(tmpy) ) .< 0 ]

      
      return abs( mean( x[indx1] ) -  mean( x[indx2] ) ) / abs( maximum([ mean( x[indx1] ) , mean( x[indx2] ) ])) > 0.001 ? ( mean( x[indx1] ) < mean( x[ indx2 ] ) ? indx1 : indx2  ) : pksx[ end - window + 1 : end ] , 
            abs( mean( y[indy1] ) -  mean( y[indy2] ) ) / abs( maximum([ mean( y[indy1] ) , mean( y[indy2] ) ])) > 0.001 ? ( mean( y[indy1] ) < mean( y[ indy2 ] ) ? indy1 : indy2  ) : pksy[ end - window + 1 : end ] 
end

stable( x;  ) = std( x ) / mean( x )

function checkWeakSynch( μ1, μ2 ; Δω = 0.2, n = 100, window = 20 )
      T=2*π;
      tspan=(0.0, n*T);
      p=(Δω, μ1, μ2);

      u0 = 3.0 * rand( 4, 1 );
      problem=ODEProblem(f, u0, tspan, p);
      t, x, dx, y, dy=vpdSolve(problem, false, 10);
      pksx, _ = findminima(x)
      pksy, _ = findminima(y) 
      px, py = cleanMinima( x, y, pksx, pksy; window = window )

      timesx, timesy = t[ px ], t[ py ] 
      freqx, freqy = 1 ./ diff( timesx ), 1 ./ diff( timesy ) 

      sameCheck =  mean( abs.( freqx[ end -  minimum([ size(freqx, 1) , size(freqy, 1) ] ) + 1: end] - freqy[ end -  minimum([ size(freqx, 1) , size(freqy, 1) ] ) + 1: end] ) )

      if sameCheck < 5*1e-4
            cors = [ cor( x[ px[1] : end - ind ], y[ px[1] + ind : end  ]  ) for ind in 1:20:(px[2]-px[1]) ]
            return stable( freqx; ), stable( freqy;  ), sameCheck, (t[2]-t[1]) * argmax(cors), maximum(cors)
      else
            return stable( freqx; ), stable( freqy;  ), sameCheck, NaN, NaN
      end
end



T=2*π;
n=100;
tspan=(0.0, n*T);
Δω=0.6; μ1=0.2; μ2=1.2; p=(Δω, μ1, μ2);

u0=[3; 0; 3; 0];
problem=ODEProblem(f, u0, tspan, p);
t, x, dx, y, dy=vpdSolve(problem, true, 10);

pNum = Int( round( size( x, 1 ) / 100 ) )

begin
      plot()
      
      plot!(
             t[end-7*pNum:10:end], x[end-7*pNum:10:end], lw = 3, color = rwth["blue"][1], labels = L"x(t)" )
      plot!(
             t[end-7*pNum:10:end], y[end-7*pNum:10:end], lw = 3, color = rwth["magenta"][1], labels = L"y(t)" )

      plot!( size = (1000, 300) )
end


pksx, _ = findminima(x)
pksy, _ = findminima(y) 

begin
      plot()


      plot!(
            t[end-7*pNum:10:end], x[end-7*pNum:10:end], lw = 3, color = rwth["blue"][1], labels = "", alpha = 0.2 )
       plot!(
            t[end-7*pNum:10:end], y[end-7*pNum:10:end], lw = 3, color = rwth["magenta"][1], labels = "", alpha = 0.2 )


      scatter!(  t[ pksx ][end - 6:end] , x[ pksx ][end - 6:end],  color = rwth["blue"][1], labels = "", marker = 6, alpha = 1.0  )
      scatter!(  t[ pksy ][end - 6:end] , y[ pksy ][end - 6:end], color = rwth["magenta"][1], labels = "", marker = 6, alpha = 1.0  )


      plot!( size = (1000, 300) )
end


window = 20
x[pksx][ end - 2*window : end ]

pksy = pksy[ y[pksy] .< 0  ]

tmpy = y[pksy][ end - 2*window : end ]
indy1 = pksy[ end - 2*window : end ][ ( tmpy .- mean(tmpy) ) .> 0 ]
indy2 = pksy[ end - 2*window : end ][ ( tmpy .- mean(tmpy) ) .< 0 ]



px, py = cleanMinima( x, y, pksx, pksy; window = 20 )


timesx, timesy = t[ px ], t[ py ] 
freqx, freqy = 1 ./ diff( timesx ), 1 ./ diff( timesy ) 

sameCheck =  mean( abs.( freqx[ end -  minimum([ size(freqx, 1) , size(freqy, 1) ] ) + 1: end] - freqy[ end -  minimum([ size(freqx, 1) , size(freqy, 1) ] ) + 1: end] ) )


begin
      plot()


      plot!(
            t[end-7*pNum:10:end], x[end-7*pNum:10:end], lw = 3, color = rwth["blue"][1], labels = "", alpha = 0.2 )
       plot!(
            t[end-7*pNum:10:end], y[end-7*pNum:10:end], lw = 3, color = rwth["magenta"][1], labels = "", alpha = 0.2 )


      scatter!(  t[ px ] , x[ px ],  color = rwth["blue"][1], labels = "", marker = 6, alpha = 1.0  )
      scatter!(  t[ py ] , y[ py ], color = rwth["magenta"][1], labels = "", marker = 6, alpha = 1.0  )


      plot!( size = (1000, 300) )
end







Δω = 0.6
m1 = 0.01:0.03:2*Δω
m2 = 0.01:0.03:2*Δω

sx, sy, szero, ϕs, crs = zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) )

for i in eachindex( m1 )
      for j in eachindex( m2 )
            μ1 = m1[ i ]
            μ2 = m2[ j ]
            sx[ i , j ], sy[ i, j ], szero[ i, j ], ϕs[ i, j ], crs[ i, j ] = checkWeakSynch( μ1, μ2; Δω = Δω )
      end
      println("progress:  ",  i / size( m1, 1 ) , "  ||  ")
end






begin
      plot( layout = grid(2, 3), margin=5Plots.mm )

      heatmap!(m1, m2, sx, sp = 1, cmap = custom_cmap, aspect_ratio = :equal ) 
      heatmap!(m1, m2, sy, sp = 2, cmap = custom_cmap, aspect_ratio = :equal )
      heatmap!(m1, m2, szero, sp = 3, cmap = custom_cmap, aspect_ratio = :equal )
      heatmap!(m1, m2, ϕs, sp = 4, cmap = custom_cmap2, aspect_ratio = :equal ) 
      heatmap!(m1, m2, crs, sp = 5, cmap = custom_cmap3, aspect_ratio = :equal )


      xlims!( m1[1], m1[end] )
      ylims!( m2[1], m2[end] )

      
      xlabel!( L"\mu_2", sp = 4 )
      xlabel!( L"\mu_2", sp = 5 )
      xlabel!( L"\mu_2", sp = 3 )
      ylabel!( L"\mu_1", sp = 1 )
      title!( L"\textrm{std}(\Omega_x) / \textrm{mean}(\Omega_x)", sp = 1 )
      title!( L"\textrm{std}(\Omega_y) / \textrm{mean}(\Omega_y)", sp = 2 )
      title!( L"\textrm{mean}(|\Omega_x - \Omega_y|)", sp = 3 )
      title!( L"\Delta\varphi", sp = 4 )
      title!( L"\max\mathbf{Corr}_{x,y}(\tau)", sp = 5 )
      plot!(sp = 6, framestyle = :none)
      plot!( size = (1000, 600) )
end

savefig("whale_new_06.pdf")