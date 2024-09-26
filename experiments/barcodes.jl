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
          #itp=LinearInterpolation(t, x); 
          #x=itp(t_new);
          #itp=LinearInterpolation(t, dx); dx=itp(t_new);
          #itp=LinearInterpolation(t, y); y=itp(t_new);
          #itp=LinearInterpolation(t, dy); dy=itp(t_new);
          itp = cubic_spline_interpolation(t, x); x = itp( t_new )
          itp = cubic_spline_interpolation(t, dx); dx = itp( t_new )
          itp = cubic_spline_interpolation(t, y); y = itp( t_new )
          itp = cubic_spline_interpolation(t, dy); dy = itp( t_new )
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

      #u0=[3; 0; 3; 0];
      u0 = 3.0 * rand( 4, 1 );
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


Δω = 0.2
m1 = 0.01:0.005:2*Δω
m2 = 0.01:0.005:2*Δω

sx, sy, sϕ = zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) )

for i in eachindex( m1 )
      for j in eachindex( m2 )
            μ1 = m1[ i ]
            μ2 = m2[ j ]
            sx[ i , j ], sy[ i, j ], sϕ[ i, j ] = checkSynch( μ1, μ2; Δω = Δω )
      end
end





#Δω = 0.2
#thr = m2 .^ 2 ./ ( 2 * Δω .- m2 )
#thr2 = m1 .^ 2 ./ ( 2 * Δω .- m1 )

custom_cmap = cgrad( [ rwth["blue"][4], rwth["blue"][1] ], scale = :log10 )
begin
      plot( layout = grid(1, 3), margin=5Plots.mm )

      heatmap!(m1, m2, log10.( abs.( log10.(sx) ) ), sp = 1, cmap = custom_cmap, aspect_ratio = :equal ) 
      heatmap!(m1, m2, log10.( abs.(  log10.(sy) ) ), sp = 2, cmap = custom_cmap, aspect_ratio = :equal )
      heatmap!(m1, m2, log10.( abs.(- log10.(sϕ) ) ), sp = 3, cmap = custom_cmap, aspect_ratio = :equal )

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

savefig("04.png")

#savefig("synch_heatmap.tex")
#savefig("synch_heatmap.pdf")

#@btime checkSynch( 0.25, 0.25 )

#using JLD
#@save "synch.jld" "freqx" sx "freqy" sy "sphi" sϕ


heatmap(m1[40:end], m2[40:end], sϕ[40:end, 40:end], cmap = custom_cmap, aspect_ratio = :equal ) 






Δω = 0.6
m1 = 0.01:0.01:2*Δω
m2 = 0.01:0.01:2*Δω

sx, sy, sϕ = zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) ), zeros( size( m1, 1), size( m2, 1 ) )

for i in eachindex( m1 )
      for j in eachindex( m2 )
            μ1 = m1[ i ]
            μ2 = m2[ j ]
            sx[ i , j ], sy[ i, j ], sϕ[ i, j ] = checkSynch( μ1, μ2; Δω = Δω )
      end
end



phis_06 = sϕ 
sx_06 = sx
sy_06 = sy
m1_06 = m1 
m2_06 = m2


phis_08 = sϕ 
sx_08 = sx
sy_08 = sy
m1_08 = m1 
m2_08 = m2

phis_04 = sϕ 
sx_04 = sx
sy_04 = sy
m1_04 = m1 
m2_04 = m2



res = [ findall( sϕ[:, i ] .< 1e-1)[1] for i in axes(sϕ , 2 ) ]

m1_01 = m1[1:end]
m2_01 = m2[res][1:end]

m1_02 = m1[1:end-4]
m2_02 = m2[res][1:end-4]

m1_04 = m1[1:end-11]
m2_04 = m2[res][1:end-11]

m1_08 = m1[1:39]
m2_08 = m2[res][1:39]


begin
      plot()

      plot!( m1_01, m2_01, marker = 6, color = rwth["magenta"][1], lw = 2, labels = L"\Delta\omega = 0.1" )
      plot!( m1_02, m2_02, marker = 6, color = rwth["blue"][1], lw = 2, labels = L"\Delta\omega = 0.2" )
      plot!( m1_04, m2_04, marker = 6, color = rwth["petrol"][1], lw = 2, labels = L"\Delta\omega = 0.4" )
      #m2_08[4:7] .= 0.61
      plot!( m1_08, m2_08, marker = 6, color = rwth["lilac"][1], lw = 2, labels = L"\Delta\omega = 0.8" )

      #=
      plot!(
            [ 0, 0.2 ], [ 0.2, 0 ], color = rwth["magenta"][1], ls = :dash, 
            labels = L"\Delta\omega = 0.1 \textrm{, kuramoto}", dash_pattern = "on 0.8cm off 0.2cm"
      )
      plot!(
            [ 0, 0.4 ], [ 0.4, 0 ], color = rwth["blue"][1], ls = :dash, 
            labels = L"\Delta\omega = 0.2 \textrm{, kuramoto}", dash_pattern = "on 0.8cm off 0.2cm"
      )
      plot!(
            [ 0, 0.6 ], [ 0.6, 0 ], color = rwth["petrol"][1], ls = :dash, 
            labels = L"\Delta\omega = 0.3 \textrm{, kuramoto}", dash_pattern = "on 0.8cm off 0.2cm"
      )
      =#
      xlabel!(L"\mu_1")
      ylabel!(L"\mu_2")

      plot!( size = (550, 500 ), legend = :outerright,  )
end

savefig("synch_boundary_2.tex")
savefig("synch_boundary_2.pdf")






begin
      plot( )

      plot!( m1[1:end-17], m2[res][1:end-17] , #.* ( 2 * Δω .- m1[1:end-14]),# ./ ( 0.1 .- 0.3 * m1[1:end-8] ), 
            lw = 4, color = rwth["blue"][1],  )
      #ylims!( ( 0.95, 1.05 ) )
      plot!()
end


Polynomials.fit(
      m1[1:end-14], m2[res][1:end-14] .* ( 2 * Δω .- m1[1:end-14]),  1
)







Δω = 0.6
n = 100;
T=2*π;
tspan=(0.0, n*T);
u0 = 3.0 * rand( 4, 1 );
begin
      plot( layout = grid(3, 1) )

      μ1, μ2 = 0.08, 0.7
      p=(Δω, μ1, μ2);
      problem=ODEProblem(f, u0, tspan, p);
      t, x, dx, y, dy=vpdSolve(problem, false, 10);
      pksx, _ = findminima(x)
      pksy, _ = findminima(y) 

      plot!(
            t[ pksx[end-7]:end ], x[ pksx[end-7]:end ],
            lw = 4, color = rwth["blue"][1], labels = L"x(t)",
            sp = 1
      )
      plot!(
            t[ pksx[end-7]:end ], y[ pksx[end-7]:end ],
            lw = 4, color = rwth["magenta"][1], labels = L"y(t)",
            sp = 1
      )
      title!(L"\mu_1 = %$(μ1), \; \mu_2 = %$(μ2)", sp = 1 )

      μ1, μ2 = 0.16, 0.7
      p=(Δω, μ1, μ2);
      problem=ODEProblem(f, u0, tspan, p);
      t, x, dx, y, dy=vpdSolve(problem, false, 10);
      pksx, _ = findminima(x)
      pksy, _ = findminima(y) 

      plot!(
            t[ pksx[end-7]:end ], x[ pksx[end-7]:end ],
            lw = 4, color = rwth["blue"][1], labels = L"x(t)",
            sp = 2
      )
      plot!(
            t[ pksx[end-7]:end ], y[ pksx[end-7]:end ],
            lw = 4, color = rwth["magenta"][1], labels = L"y(t)",
            sp = 2
      )
      title!(L"\mu_1 = %$(μ1), \; \mu_2 = %$(μ2)", sp = 2 )
      
      μ1, μ2 = 0.32, 0.7
      p=(Δω, μ1, μ2);
      problem=ODEProblem(f, u0, tspan, p);
      t, x, dx, y, dy=vpdSolve(problem, false, 10);
      pksx, _ = findminima(x)
      pksy, _ = findminima(y) 

      plot!(
            t[ pksx[end-7]:end ], x[ pksx[end-7]:end ],
            lw = 4, color = rwth["blue"][1], labels = L"x(t)",
            sp = 3
      )
      plot!(
            t[ pksx[end-7]:end ], y[ pksx[end-7]:end ],
            lw = 4, color = rwth["magenta"][1], labels = L"y(t)",
            sp = 3
      )
      title!(L"\mu_1 = %$(μ1), \; \mu_2 = %$(μ2)", sp = 3 )


      plot!( size = (1000, 750) )
end


savefig("examples_trio.tex")
savefig("examples_trio.pdf")



μ1, μ2 = 0.16, 0.7
p=(Δω, μ1, μ2);
problem=ODEProblem(f, u0, tspan, p);
t, x, dx, y, dy=vpdSolve(problem, false, 10);
pksx, _ = findminima(x)
pksy, _ = findminima(y) 

begin
      plot()

      plot!( t[pksy[end-7]:end], y[pksy[end-7]:end], color = rwth["magenta"][1], alpha = 0.5, lw = 4, label=L"y(t)"
      )
      scatter!(
            t[pksy[end-7:end]],  y[pksy[end-7:end]], m = 6, color = rwth["magenta"][1], label=""
      )

      plot!( size = (1000, 300) )
end

y[pksx[end-7:end]]


