using LinearAlgebra, Arpack
using SparseArrays
using StatsBase
using DifferentialEquations
using Interpolations
using Peaks
using Polynomials
using Statistics

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
          itp=LinearInterpolation(t, x); 
          x=itp(t_new);
          itp=LinearInterpolation(t, dx); dx=itp(t_new);
          itp=LinearInterpolation(t, y); y=itp(t_new);
          itp=LinearInterpolation(t, dy); dy=itp(t_new);
          t=t_new;
          return t, x, dx, y, dy
      else
          return t, x, dx, y, dy
      end
end


T=2*π;
n=100;
tspan=(0.0, n*T);
Δω=0.2; μ1=0.21; μ2=0.21; p=(Δω, μ1, μ2);

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

# same length is not guaranteed
timesx, timesy = t[ pksx ], t[ pksy ] 
freqx, freqy = 1 ./ diff( timesx ), 1 ./ diff( timesx ) 


phases = ( size( timesx, 1 ) < size( timesy, 1 ) ) ?  [  minimum( abs.( timesy .- timesx[i] ) ) for i in eachindex( timesx )
] :  [ minimum( abs.( timesx[ maximum([i - 10,0]) : minimum([i + 10, end]) ] .- timesy[i] ) ) for i in eachindex( timesy )
]  

stable( x; window = 20 ) = std( x[ end - window : end ] ) / mean( x[ end - window : end ] )


function checkSynch( μ1, μ2 ; Δω = 0.2, n = 100, window = 20 )
      T=2*π;
      tspan=(0.0, n*T);
      p=(Δω, μ1, μ2);

      u0=[3; 0; 3; 0];
      problem=ODEProblem(f, u0, tspan, p);
      t, x, dx, y, dy=vpdSolve(problem, false, 10);
      pksx, _ = findminima(x)
      pksy, _ = findminima(y) 

      timesx, timesy = t[ pksx ], t[ pksy ] 
      freqx, freqy = 1 ./ diff( timesx ), 1 ./ diff( timesx ) 

      phases = ( size( timesx, 1 ) < size( timesy, 1 ) ) ?  [  minimum( abs.( timesy .- timesx[i] ) ) for i in eachindex( timesx )
      ] :  [ minimum( abs.( timesx .- timesy[i] ) ) for i in eachindex( timesy )
      ]  
      

      return stable( freqx; window = window ), stable( freqy; window = window ), stable( phases; window = window )

end


m1 = 0.01:0.005:0.4
m2 = 0.01:0.005:0.4

sx, sy, sϕ = zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) )

for i in eachindex( m1 )
      for j in eachindex( m2 )
            μ1 = m1[ i ]
            μ2 = m2[ j ]
            sx[ i , j ], sy[ i, j ], sϕ[ i, j ] = checkSynch( μ1, μ2 )
      end
end

Δω = 0.2
thr = m2 .^ 2 ./ ( 2 * Δω .- m2 )
thr2 = m1 .^ 2 ./ ( 2 * Δω .- m1 )

custom_cmap = cgrad( [ rwth["blue"][4], rwth["blue"][1] ], scale = :log10 )
begin
      plot( layout = grid(1, 3), margin=5Plots.mm )

      heatmap!(m1, m2, sx, sp = 1, cmap = custom_cmap, aspect_ratio = :equal ) 
      heatmap!(m1, m2, sy, sp = 2, cmap = custom_cmap, aspect_ratio = :equal )
      heatmap!(m1, m2, sϕ, sp = 3, cmap = custom_cmap, aspect_ratio = :equal )

      xlims!( m1[1], m1[end] )
      ylims!( m2[1], m2[end] )

      #=plot!( thr[20:39], m2[20:39], color = rwth["magenta"][1], ls = :dash, lw = 1, labels = "",
            sp = 3
      )
      plot!( m1[20:39], thr2[20:39], color = rwth["red"][1], ls = :dash, lw = 1, labels = "",
            sp = 3
      )=#

      xlabel!( L"\mu_1" )
      ylabel!( L"\mu_2", sp = 1 )
      title!( L"\textrm{std}(\Omega_x) / \textrm{mean}(\Omega_x)", sp = 1 )
      title!( L"\textrm{std}(\Omega_y) / \textrm{mean}(\Omega_y)", sp = 2 )
      title!( L"\textrm{std}(\Delta\varphi) / \textrm{mean}(\Delta\varphi)", sp = 3 )
      plot!( size = (1000, 300) )
end

savefig("synch_heatmap.tex")
savefig("synch_heatmap.pdf")

#@btime checkSynch( 0.25, 0.25 )