using LinearAlgebra, Arpack
using SparseArrays
using StatsBase
using DifferentialEquations
using Interpolations
using Peaks
using Polynomials
using Statistics

using BenchmarkTools
using Printf

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
      print("!")
      x=reduce(hcat, sol.u)[1, :];
      dx=reduce(hcat, sol.u)[2, :];
      y=reduce(hcat, sol.u)[3, :];
      dy=reduce(hcat, sol.u)[4, :];
      t=sol.t;
      tspan=problem.tspan;
      Δt = minimum( diff( t ) );
      if Δt < 1e-4
          Δt = 1e-4
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
          print("?")
          return t, x, dx, y, dy
      else
          return t, x, dx, y, dy
      end
end


function getLimCycleNaive(t, x, dx, y, dy)
    inds_x=findminima(x)[1];
    γ_x=x[inds_x[end-1]:inds_x[end]];
    γ_dx=dx[inds_x[end-1]:inds_x[end]];

    inds_y=findminima(y)[1];
    γ_y=y[inds_y[end-1]:inds_y[end]];
    γ_dy=dy[inds_y[end-1]:inds_y[end]];
    return γ_x, γ_dx, γ_y, γ_dy
end


function getPhiC(mu::Float64, Tstar::Float64; dw = 0.2)
    u0=[3; 0; 3; 0];
     #mu=0.25; inherited from input
    p=(dw, mu);
    T=2*π;
    n=150;
    tspan=(0.0, n*T);  
    problem=ODEProblem(f, u0, tspan, p);
    t, x, dx, y, dy=vpdSolve(problem, true, 10);
    t_ind=findfirst(x->x==true, t.>Tstar);

    h=t[2]-t[1];
    p_approx=Int(round(T/h));
    step=2;
    phi=Int(round(p_approx/4));
    shifts=[maximum((x[t_ind+phi:end]-y[t_ind:end-phi])) for phi in 1:step:p_approx-1];
    C, phi_ind=findmin(shifts);
    phi=(1:step:p_approx-1)[phi_ind]*h
    phi=phi % (2*π);

    γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
    γ=(γ_x, γ_dx, γ_y, γ_dy);
    pNum=size(γ_x, 1);
    Ax = 0.5*[ maximum( x[i : i + pNum-1 ] ) - minimum( x[i : i + pNum-1 ] ) for i in size(x, 1) - 2*pNum + 1 : 10 : size(x, 1) - pNum + 1 ]
    Ay = 0.5*[ maximum( y[i : i + pNum-1 ] ) - minimum( y[i : i + pNum-1 ] ) for i in size(x, 1) - 2*pNum + 1 : 10 : size(y, 1) - pNum + 1] 

    return phi, C - (Ax[end]-Ay[end]), abs(abs(C) - abs(Ax[end]-Ay[end]))
end

function rectangle(w, h, x, y)
    Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
end

function rectplot( w, h, x, y; lw = 2, ls = :dash, color = :white, sp = 1, labels = "" )
    plot!( [x, x], [ y, y + h ] , lw = lw, ls = ls, color = color , sp = sp, labels = labels)
    plot!( [x, x+w], [y,  y  ], lw = lw, ls = ls, color = color , sp = sp , labels = labels)
    plot!( [x+w, x+w], [y, y + h ] , lw = lw, ls = ls, color = color , sp = sp, labels = labels)
    plot!( [x, x+w], [ y+h, y + h ] , lw = lw, ls = ls, color = color , sp = sp, labels = labels)
end

function getShifts( step, pNum, x, y )
    inds = 1:step:pNum-1
    res = zeros( size( inds, 1 ) )
    for i in eachindex( inds )
        res[ i ] = maximum( x[ inds[i]:end ] - y[1:end - inds[i] + 1 ] )
    end
    return res
end

stable( x; window = 20 ) = std( x[ end - window : end ] ) / mean( x[ end - window : end ] )



Δω=0.2;
m1 = [ 0.01:0.02:2*Δω; 2*Δω:0.05:1.0 ; 1.1:0.1:4.5]
m2 = [ 0.01:0.02:2*Δω; 2*Δω:0.05:1.0 ; 1.1:0.1:4.5]

Cs = zeros( size(m1, 1), size(m2, 1) )
ΔA = zeros( size(m1, 1), size(m2, 1) )
phis = zeros( size(m1, 1), size(m2, 1) )


for i in eachindex( m1 )
    @time begin
        for j in eachindex( m2 )
            u0=[3; 0; 3; 0] + 0.01 * rand(4, 1);
             μ1=m1[i]; μ2=m2[j]; p=(Δω, μ1, μ2);
            T=2*π;
            n=120;
            tspan=(0.0, n*T);  
            problem=ODEProblem(f, u0, tspan, p);
            t, x, dx, y, dy=vpdSolve(problem, true, 10);
            pksx, _ = findminima(x)
            pksy, _ = findminima(y) 
            print(".")

            timesx, timesy = t[ pksx ], t[ pksy ] 
            #freqx, freqy = 1 ./ diff( timesx ), 1 ./ diff( timesx ) 
      
            phases = ( size( timesx, 1 ) < size( timesy, 1 ) ) ?  [  minimum( abs.( timesy .- timesx[i] ) ) for i in eachindex( timesx )
            ] :  [ minimum( abs.( timesx .- timesy[i] ) ) for i in eachindex( timesy )
            ]  
            print(":")
            if ( stable( phases ) <= 0.01 ) || ( max( μ1, μ2 ) > 2 * Δω )
                Tstar = 100 * T
                #t_ind=findfirst(x->x==true, t.>Tstar);
                t_ind = Int( round( 100/n* size( t, 1 ) ) )

                h=t[2]-t[1];
                #γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
                #print("!")
                #γ=(γ_x, γ_dx, γ_y, γ_dy);
                #pNum=Int( round( median( [ size(γ_x, 1), size(γ_dx, 1), size(γ_y, 1), size(γ_dy, 1) ] #) ) )
                pNum = Int( round( 0.5 *( pksx[end] - pksx[end-1] + pksy[end] - pksy[end-1] ) ) )
                print(pNum)
                step=5;
                #phi=Int(round(pNum/4));
                #shifts=[maximum((x[t_ind+phi:end]-y[t_ind:end-phi])) for phi in 1:step:pNum-1];
                shifts = getShifts( step, pNum, x[t_ind:end], y[t_ind:end] )
                print(";")
                C, phi_ind=findmin(shifts);
                print("!")
                phi=(1:step:pNum-1)[phi_ind]*h
                phi=phi % (2*π);

                Ax = 0.5*[ maximum( x[i : i + pNum-1 ] ) - minimum( x[i : i + pNum-1 ] ) for i in size(x, 1) - 2*pNum + 1 : 10 : size(x, 1) - pNum + 1 ]
                Ay = 0.5*[ maximum( y[i : i + pNum-1 ] ) - minimum( y[i : i + pNum-1 ] ) for i in size(x, 1) - 2*pNum + 1 : 10 : size(y, 1) - pNum + 1] 
                print("?")
                
                ΔA[ i, j ] = Ax[end] - Ay[end] 
                Cs[ i, j ] = C 
                phis[ i, j ] = phi 
            else
                ΔA[ i, j ] = NaN
                Cs[ i, j ] = NaN 
                phis[ i, j ] = NaN 
            end
            #phis[i], Cs[i], Cs_test[i]=getPhiC(mus[i], Tstar; dw = dw);
            @printf "coupling: %f / %f   |||  phase dif:  %f   :::  vert shift: %f\n " m1[i] m2[j] phis[i,j] Cs[i, j]-ΔA[i,j]
        end
    end
end



using JLD
@save "cphi_02_2.jld" "freqdiff" Δω "vert" Cs "ampl" ΔA "phiss" phis
 



begin
    l = @layout [  a{0.5w} [ b{0.66w} c; d e ] ]
    plot( layout = l, margin=5Plots.mm )
    
    custom_cmap = cgrad( [ rwth["black"][1], rwth["violet"][1], rwth["bord"][1], rwth["magenta"][2], rwth["orange"][2] ], ) #[0, 0.05, 0.1, 0.12, 0.13] * 5/3 )
    heatmap!( (m1), (m2), (Cs-ΔA), #clim = (0, 0.15), 
        sp = 1,  cmap = custom_cmap, )
    rectplot( m1[30] - m1[1], m1[30] - m1[1], m1[1], m1[1]; color = :white, )
    rectplot(m1[end] - m1[10], m1[end] - m1[10], m1[10], m1[10] ; color = rwth["yellow"][1], )
    annotate!(
        m1[30]+0.1, m1[30]+0.1, text( L"\textbf{(A)}", 12, :white ),
        sp = 1
    )
    annotate!(
        m1[10]+0.15, m1[end]-0.15, text( L"\textbf{(B)}", 12, rwth["yellow"][1] ),
        sp = 1
    )
    annotate!(
        0.5, 2.9, text( L"\textbf{(C)}", 12, rwth["yellow"][1] ),
        sp = 1
    )
    plot!( [0.1, 0.35], [ 2.6, 2.9 ], lw = 1, sp = 1, labels="", color = rwth["yellow"][1] )
    annotate!(
        3.5, 0.45, text( L"\textbf{(D)}^\top", 12, rwth["yellow"][1] ),
        sp = 1
    )
    plot!( [3.6, 3.8], [ 0.4, 0.1 ], lw = 1, sp = 1, labels="", color = rwth["yellow"][1] )
    #plot!(
    #        rectangle( m1[30] - m1[1], m1[30] - m1[1], m1[1], m1[1] ), fillcolor = :false, linecolor = :white, ls = :dash, sp = 1, labels=""
    #)    
    #plot!(
    #        rectangle( m1[end] - m1[10], m1[end] - m1[10], m1[10], m1[10] ), fillcolor = :false, linecolor = rwth["yellow"][1], ls = :dash, sp = 1, labels=""
    #)

    heatmap!( (m1)[1:30], (m2)[1:30], (Cs-ΔA)[1:30, 1:30], #clim = (0, 0.15), 
        sp = 2,  cmap = custom_cmap,     colorbar_tickfontsize = 10,)
    heatmap!( (m2)[10:end], (m1)[10:end], (Cs-ΔA)[10:end, 10:end], #clim = (0, 0.15), 
        sp = 4,  cmap = custom_cmap,    colorbar_tickfontsize = 10, )

    heatmap!(  (m1)[1:10], (m2)[10:end], (Cs-ΔA)[1:10, 10:end]', #clim = (0, 0.15), 
    sp = 3,  cmap = custom_cmap,     colorbar_tickfontsize = 10,)
    heatmap!( (m2)[1:10], (m1)[10:end],  (Cs-ΔA)[10:end, 1:10], #clim = (0, 0.15), 
    sp = 5,  cmap = custom_cmap,     colorbar_tickfontsize = 10,)
   
    tmp = (Cs-ΔA)[10:end, 10:end]
    inds = [ nanargmax( tmp[i, 5:end] ) for i in 1:size(tmp , 1) ]
    plot!(    m1[ 15 .+ inds ][1:end-31],  m2[10:end-31],      lw = 2, color = rwth["black"][3], sp = 4, labels = "" )

    xlims!( m1[1], m1[end], sp = 1)
    ylims!( m1[1], m1[end], sp = 1)
    xlims!( m1[10], m1[end], sp = 4)
    ylims!( m1[10], m1[end], sp = 4)

    ylabel!( L"\mu_1", sp = 1 )
    xlabel!( L"\mu_2", sp = 1 ) 
    xlabel!( L"\mu_2", sp = 2, xguidefontsize=14)
    xlabel!( L"\mu_2", sp = 3, xguidefontsize=14 )
    xlabel!( L"\mu_2", sp = 4, xguidefontsize=14 )
    xlabel!( L"\mu_1", sp = 5, xguidefontsize=14 )

    ylabel!( L"\mu_1", sp = 2, yguidefontsize=14)
    ylabel!( L"\mu_1", sp = 3, yguidefontsize=14 )
    ylabel!( L"\mu_1", sp = 4, yguidefontsize=14 )
    ylabel!( L"\mu_2", sp = 5, yguidefontsize=14 )

    xticks!( [0.01, 0.2], [L"0", L"0.2"], sp = 3)
    xticks!( [0.01, 0.2], [L"0", L"0.2"], sp = 5)

    title!( L"C - \Delta A, \textrm{full}", sp = 1 )
    title!( L"\textbf{(A), }\textrm{small values}", sp = 2, titlefontsize = 16  )
    title!( L"\textbf{(B), }\mu_i > \Delta\omega", sp = 4, titlefontsize = 16 )
    title!( L"\textbf{(C), }\mu_1 <  \Delta\omega", sp = 3, titlefontsize = 16 )
    title!( L"\textbf{(D), }\mu_2 <  \Delta\omega", sp = 5, titlefontsize = 16 )
    plot!( size = (1000, 500) )
end

savefig("vertshift_heatmap.tex")
savefig("vertshift_heatmap.pdf")



function getInds( δ )
#δ = 0.5
    top = findlast( m1 * ( 1 + δ ) .<= maximum(m2) ) 
    bot = findfirst( m1 * ( 1 + δ ) .>= minimum(m2) ) 
    inds = [ argmin( abs.( m2 .- (1 + δ) * m1[i] ) ) for i in bot:top ]
    return bot:top, inds 
end



r1, r2 = getInds( 0.0 )
[ phis[ r1[i], r2[i] ] for i in eachindex(r1) ]

begin
    plot( layout = grid( 1, 2), margins = 10Plots.mm )

    custom_cmap = cgrad( [ rwth["blue"][4], rwth["blue"][1] ] )
    heatmap!(m1, m2, phis, sp = 1, cmap = custom_cmap, aspect_ratio=:equal) 

    r1, r2 = getInds( 0.0 )
    plot!(
        [ (Cs-ΔA)[ r1[i], r2[i] ] for i in eachindex(r1) ], [ phis[ r1[i], r2[i] ] for i in eachindex(r1) ],  
        sp = 2, lw = 3, color = rwth["blue"][1], labels = L"\delta = 0"    
   )
   plot!(
        m1[r1], m2[r2], lw = 1, color=rwth["black"][2], labels="", ls = :dash
   )
    annotate!( 3.5, 3.5+0.2, text(L"\delta = 0", 12, color=rwth["black"][2], rotation = 180 * atan( 1 + 0.0) / π ) )

   r1, r2 = getInds( 0.2 )
    plot!(
        [ (Cs-ΔA)[ r1[i], r2[i] ] for i in eachindex(r1) ], [ phis[ r1[i], r2[i] ] for i in eachindex(r1) ],  
        sp = 2, lw = 2, color = rwth["magenta"][1], alpha =0.5, labels = L"\delta = 0.1"
    )
    plot!(
        m1[r1], m2[r2], lw = 1, color=rwth["black"][2], labels="", ls = :dash
   )
   annotate!( 2.8, 2.8 * (1 + 0.2)+0.2, text(L"\delta = 0.2", 12, color=rwth["black"][2], rotation = 180 * atan( 1 + 0.2) / π ) )
    
   r1, r2 = getInds( -0.2 )
    plot!(
        [ (Cs-ΔA)[ r1[i], r2[i] ] for i in eachindex(r1) ], [ phis[ r1[i], r2[i] ] for i in eachindex(r1) ],  
        sp = 2, lw = 2, color = rwth["magenta"][1], ls = :dash, labels = L"\delta = -0.1"
    )
    plot!(
        m1[r1], m2[r2], lw = 1, color=rwth["black"][2], labels="", ls = :dash
   )
   annotate!( 2.8, 2.8 * (1 - 0.2)+0.2, text(L"\delta = -0.2", 12, color=rwth["black"][2], rotation = 180 * atan( 1 - 0.2) / π ) )

    r1, r2 = getInds( 0.4 )
    plot!(
        [ (Cs-ΔA)[ r1[i], r2[i] ] for i in eachindex(r1) ], [ phis[ r1[i], r2[i] ] for i in eachindex(r1) ],  
        sp = 2, lw = 2, color = rwth["petrol"][1], alpha = 0.5, labels = L"\delta = 0.2"
    )
    plot!(
        m1[r1], m2[r2], lw = 1, color=rwth["black"][2], labels="", ls = :dash
   )
   annotate!( 2.5, 2.5 * (1 + 0.4)+0.2, text(L"\delta = 0.4", 12, color=rwth["black"][2], rotation = 180 * atan( 1 + 0.4) / π ) )

    r1, r2 = getInds( -0.4 )
    plot!(
        [ (Cs-ΔA)[ r1[i], r2[i] ] for i in eachindex(r1) ], [ phis[ r1[i], r2[i] ] for i in eachindex(r1) ],  
        sp = 2, lw = 2, color = rwth["petrol"][1], ls = :dash, labels = L"\delta = -0.2"
    )
    plot!(
        m1[r1], m2[r2], lw = 1, color=rwth["black"][2], labels="", ls = :dash
   )
   annotate!( 2.5, 2.5 * (1 - 0.4)+0.2, text(L"\delta = -0.4", 12, color=rwth["black"][2], rotation = 180 * atan( 1 - 0.4) / π ) )

    xlims!(m1[1], m1[end], sp = 1, aspect_ratio=:equal)
    ylims!(m1[1], m1[end], sp = 1)

    xlabel!( L"\mu_2", sp = 1)
    xlabel!( L"C(\mu_1, \mu_2), \textrm{vertical shift}", sp = 2)
    ylabel!( L"\mu_1", sp = 1)
    ylabel!( L"\varphi(\mu_1, \mu_2), \textrm{phase diff}", sp = 2)
    title!( L"\varphi(\mu_1, \mu_2)" , sp = 1 )
    title!( L"\frac{\mu_1}{\mu_2} = 1 + \delta", sp = 2)
    plot!( legend = :topright )
    plot!( size = (800, 400 ) )
end



savefig("cphi_heatmap.tex")
savefig("cphi_heatmap.pdf")






























begin
    plot( layout = grid(1, 3),)# margin=5Plots.mm )
    custom_cmap = cgrad( [ rwth["blue"][4], rwth["blue"][1] ] )
    heatmap!(m1, m2, phis, sp = 1, cmap = custom_cmap, aspect_ratio = :equal ) 
   # custom_cmap = cgrad( [ rwth["magenta"][1], rwth["blue"][1] ] )
   custom_cmap = cgrad( [ rwth["black"][1], rwth["violet"][1], rwth["bord"][1], rwth["magenta"][2], rwth["orange"][2] ], ) #[0, 0.05, 0.1, 0.12, 0.13] * 5/3 )
    heatmap!( (m1), (m2), (Cs-ΔA), #clim = (0, 0.15), 
        sp = 2,  cmap = custom_cmap, aspect_ratio = :equal )
    custom_cmap = cgrad( [ rwth["petrol"][4], rwth["petrol"][1] ] )
    heatmap!(m1, m2, ΔA, sp = 3, cmap = custom_cmap, aspect_ratio = :equal )

    #xlims!( m1[1], m1[end] )
    #ylims!( m2[1], m2[end] )

    #=plot!( thr[20:39], m2[20:39], color = rwth["magenta"][1], ls = :dash, lw = 1, labels = "",
          sp = 3
    )
    plot!( m1[20:39], thr2[20:39], color = rwth["red"][1], ls = :dash, lw = 1, labels = "",
          sp = 3
    )=#

    xlabel!( L"\mu_1" )
    ylabel!( L"\mu_2", sp = 1 )
    title!( L"\Delta\varphi( \mu_1, \mu_2)", sp = 1 )
    title!( L"C - \Delta A", sp = 2 )
    title!( L"\Delta A(\mu_1, \mu_2)", sp = 3 )
    plot!( size = (1200, 300) )
end

custom_cmap = cgrad( [ rwth["black"][1], rwth["violet"][1], rwth["bord"][1], rwth["magenta"][2], rwth["orange"][2] ], [0, 0.05, 0.1, 0.12, 0.13] * 5/3 )

custom_cmap = cgrad( [ rwth["magenta"][4], rwth["magenta"][1] ] )
heatmap(m1, m2, (Cs-ΔA), cmap = custom_cmap, aspect_ratio = :equal )

gr()
begin
    heatmap(  m2[22:end], m1[20:end], ( Cs-ΔA )[20:end, 22:end], size = ( 500, 500 ), clim = (0, 0.13) )
    heatmap!(  m2[1:22], m1[20:end], ( Cs-ΔA )[20:end, 1:22], size = ( 500, 500 ), clim = (0, 0.13) )
    heatmap!(  m2[22:end], m1[1:20], ( Cs-ΔA )[1:20, 22:end], size = ( 500, 500 ), clim = (0, 0.13) ) 
end

p = heatmap(data)
scatter!(p, Tuple.(findall(isnan, data)), label="NaNs", markercolor=:red)
scatter!(p, Tuple.(findall(isinf, data)), label="Infs", markercolor=:black)



heatmap(
    [1, 2,], [1, 2, 4], [ 1 2 3; 4 5 6;]
)