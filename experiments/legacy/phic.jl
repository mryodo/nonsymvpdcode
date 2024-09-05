using LinearAlgebra, Arpack
using BenchmarkTools
#math
using DifferentialEquations
using Interpolations
using Peaks
using Polynomials

#plots
using Plots
pgfplotsx()
theme(:mute)
using ColorSchemes
cols=ColorSchemes.Spectral_11;
Plots.scalefontsizes(1.5);
using LaTeXStrings
using DelimitedFiles

using Printf

include("vpdModule.jl")


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

dw = 0.2
T=2*π;
#getPhiC(0.25, 50*T)
#getPhiC(1.5, 50*T)
#getPhiC(6.0, 50*T)


function several( ; dw = 0.2, N = 100  )

      mus=[range( dw * 1.1, 1.75, 50); range(1.8, 6, 20)];
      Tstar=N*T;

      Cs=zeros(size(mus, 1));
      Cs_test=zeros(size(mus, 1));
      phis=zeros(size(mus, 1));

      
      for i in 1:size(mus,1)
            @time phis[i], Cs[i], Cs_test[i]=getPhiC(mus[i], Tstar; dw = dw);
            @printf "coupling: %f   |||  phase dif:  %f   :::  vert shift: %f\n " mus[i] phis[i] Cs[i]
      end
      return phis, Cs
end

#phis, Cs = several( dw = 0.2 )
#phis03, Cs03 = several( dw = 0.3 )
#phis04, Cs04 = several( dw = 0.4 )
#phis001, Cs001 = several( ; dw = 0.01 )
#phis005, Cs005 = several( ; dw = 0.05 )
#phis05, Cs05 = several( ; dw = 0.5 )
#phis095, Cs095 = several( ; dw = 0.95 )

#HERE LIES THE PLOTTER, ONLY AESTHETICS


#=
writedlm( "phi_dw001.out", phis001 )
writedlm( "phi_dw005.out", phis005 )
writedlm( "phi_dw01.out", phis01 )
writedlm( "phi_dw05.out", phis05 )
writedlm( "phi_dw09.out", phis09 )
writedlm( "phi_dw095.out", phis095 )
writedlm( "C_dw001.out", Cs001 )
writedlm( "C_dw005.out", Cs005 )
writedlm( "C_dw01.out", Cs01 )
writedlm( "C_dw05.out", Cs05 )
writedlm( "C_dw09.out", Cs09 )
writedlm( "C_dw095.out", Cs095 )
=#

Cs = readdlm( "C_dw02.out" )
Cs01 = readdlm( "C_dw01.out" )
Cs001 = readdlm( "C_dw001.out" )
Cs005 = readdlm( "C_dw005.out" )
Cs03 = readdlm( "C_dw03.out" )
Cs04 = readdlm( "C_dw04.out" )
Cs05 = readdlm( "C_dw05.out" )
Cs09 = readdlm( "C_dw09.out" )
Cs095 = readdlm( "C_dw095.out" )
phis = readdlm( "phi_dw02.out" )
phis01 = readdlm( "phi_dw01.out" )
phis001 = readdlm( "phi_dw001.out" )
phis005 = readdlm( "phi_dw005.out" )
phis03 = readdlm( "phi_dw03.out" )
phis04 = readdlm( "phi_dw04.out" )
phis05 = readdlm( "phi_dw05.out" )
phis09 = readdlm( "phi_dw09.out" )
phis095 = readdlm( "phi_dw095.out" )

ind=findmax(Cs[3:end])[2];
plot_font = "Computer Modern";
grad=cgrad(:Spectral_5, 1000, categorical=true)

ψ=range(π/12+π/2, 0, 1000);
rad=0.25

dw = 0.2
mus=[range( dw * 1.1, 1.75, 50); range(1.8, 6, 20)];




begin
      plot()
      #scatter!(Cs[4:end], phis[4:end], xscale=:log10, yscale=:log10, color=:firebrick, labels="")
      scatter!(Cs[4:end], phis[4:end], xscale=:log10, yscale=:log10, color=grad[Int.(round.(( mus[4:end] .- minimum(mus[4:end]) )/(maximum(mus[4:end]) - minimum(mus[4:end]))*998)).+1], labels="", colorbar=true)

     
      #=
      annotate!(
            0.14, 0.5, text(L"\mathbf{\Delta\omega = 0.4}", :left, 12, :black)
      )
      annotate!(
            0.135, 0.35, text(L"\mathbf{\Delta\omega = 0.3}", :left, 12, :blue)
      )
      annotate!(
            0.025, 0.875, text(L"\mathbf{\Delta\omega = 0.2}", :left, 12, :red)
      )
      =#
      plot!(Cs[4:end], phis[4:end], xscale=:log10, yscale=:log10, alpha=0.5, lw=4, color=:firebrick, linestyle=:solid, labels=L"{\Delta\omega = 0.2}", )
       plot!(
            Cs005[4:end], phis005[4:end], xscale=:log10, yscale=:log10, color=:orange, alpha = 0.5, lw=3,  line=:dash, labels=L"{\Delta\omega = 0.05}"
      )
       plot!(
            Cs01[4:end], phis01[4:end], xscale=:log10, yscale=:log10, color=:green, alpha = 0.5, lw=3,  line=:dash, labels=L"{\Delta\omega = 0.1}"
      )

      
       plot!(
            Cs03[4:end], phis03[4:end], xscale=:log10, yscale=:log10, color=:blue, alpha=0.5,  lw=3, line=:dash, labels=L"{\Delta\omega = 0.3}"
      )

      plot!(
            Cs04[4:end], phis04[4:end], xscale=:log10, yscale=:log10, color=:black, alpha = 0.5, lw=3,  line=:dash, labels=L"{\Delta\omega = 0.4}"
      )
     
     
      plot!(
            Cs05[4:end], phis05[4:end], xscale=:log10, yscale=:log10, color=:brown, alpha = 0.5, lw=3,  line=:dash, labels=L"{\Delta\omega = 0.5}"
      )
      plot!(
            Cs09[11:end], phis09[11:end], xscale=:log10, yscale=:log10, color=:purple, alpha = 0.5, lw=3,  line=:dash, labels=L"{\Delta\omega = 0.9}"
      )
      plot!(
            Cs095[11:end], phis095[11:end], xscale=:log10, yscale=:log10, color=:pink, alpha = 0.5, lw=3,  line=:dash, labels=L"{\Delta\omega = 0.95}"
      )
      scatter!([Cs[ind-1]], [phis[ind-1]], color=:black, markersize=7, labels="")
      annotate!(Cs[ind-1]*0.35*3, 1.0*11/9*phis[ind-1]/4, 
                  text(L"\mu^* \approx 0.97", :left, 11
                  )
                  )
      annotate!(Cs[ind-1]*0.35*3, 1.0*phis[ind-1]/4, 
                  text(L"\mathbf{C(\mu^*)} \approx 0.1", :left, 11
                  )
                  )
      annotate!(Cs[ind-1]*0.35*3, 9/11*phis[ind-1]/4, 
                  text(L"\mathbf{\varphi(\mu^*)}\approx 0.23", :left, 11
                  )
                  )

      plot!([ Cs[ind-1]*1.1, 1.0*11/9*phis[ind-1]/2.5 ], [ phis[ind-1], Cs[ind-1]*0.35*2.75 ], color=:black, lw=1.5, labels="")

      annotate!(1.2*Cs[3], 0.9*phis[3], text(L"\mathbf{\mu \approx \Delta\omega}", :right, 10))
      #annotate!(0.9*Cs[end], 1.1*phis[end], text(L"\mathbf{\mu \to \infty}"*("\n")*L"(C(\mu), \varphi(\mu) \to (0, 0)", :right, 10))

      annotate!(0.7*Cs[end]/3, 1.1*phis[end]/3, text(L"\mathbf{\mu \to \infty}", :right, 13))
      annotate!(0.9*Cs[end]/3, 0.975*phis[end]/3, text(L"C(\mu), \varphi(\mu) \to (0, 0)", :right, 13))
      #C(\mu), \varphi(\mu) \to (0, 0)"

      plot!(10 .^ (rad*cos.(ψ) .- 1.75)/2, 10 .^ (rad*sin.(ψ) .- 0.6), color=:black, lw=2, labels="")
      scatter!( [10 .^ (rad*cos.(ψ[end]) .- 1.75)]/2, [10 .^ (rad*sin.(ψ[end]) .- 0.6)], color=:black, markershape=:dtriangle, labels="")
      annotate!(10 .^ (rad*cos.(ψ[1]) .- 1.75)/2, 1.1*10 .^ (rad*sin.(ψ[1]) .- 0.6), text(L"\mathbf{\mu}\textrm{\; increases}", :left, :bottom, 14))

      

      wid=0.0006
      x1=10^(-1.1)
      for i in 1:size(grad)[1]
            plot!([x1; x1*(10^wid)], [0.04/4; 0.045/4], color=grad[i], lw=1.5, labels="")
            x1=x1*(10^wid)
      end
      annotate!(0.15, 0.05/4, text(L"\textrm{coupling, \; }\mu", :center, 10))
      annotate!(10^(-1.1), 0.0375/4, text(L"\Delta\omega", :center, 10))
      annotate!(x1/(10^(5*wid)), 0.0375/4, text(L"6", :right, 10))


     # xlims!((0.0075, 0.3))
      xlabel!(L"\textrm{V}\textrm{ertical \; mismatch, }C(\mu)")
      ylabel!(L"\textrm{P}\textrm{hase \; difference, }\varphi(\mu)")
      plot!(size=(600, 600), framestyle=:axes, legend = :left)

end

savefig("figure17.png")

savefig("figure6_semi.png")
savefig("figure6.svg")

savefig("Cphi.png")
savefig("Cphi.svg")








begin
    plot()
    #scatter!(Cs[4:end], phis[4:end], xscale=:log10, yscale=:log10, color=:firebrick, labels="")
    scatter!(Cs_test[4:end], phis[4:end], xscale=:log10, yscale=:log10, color=grad[Int.(round.(( mus[4:end] .- minimum(mus[4:end]) )/(maximum(mus[4:end]) - minimum(mus[4:end]))*998)).+1], labels="", colorbar=true)
    
    plot!(Cs_test[4:end], phis[4:end], xscale=:log10, yscale=:log10, alpha=0.2, lw=4, color=:firebrick, linestyle=:solid, labels="", )
    scatter!([Cs_test[ind-1]], [phis[ind-1]], color=:black, markersize=7, labels="")
    annotate!(Cs_test[ind-1]*1.1, 1.1*phis[ind-1], 
                text(L"\mathbf{C(\mu^*)} \approx 0.1", :left, 11
                )
                )
    annotate!(Cs_test[ind-1]*1.1, 0.90*phis[ind-1], 
                text(L"\mathbf{\varphi(\mu^*)}\approx 0.23", :left, 11
                )
                )
    
    annotate!(1.6*Cs_test[3], 0.8*phis[3], text(L"\mathbf{\mu \approx \Delta\omega}", :right, 10))
    #annotate!(0.9*Cs[end], 1.1*phis[end], text(L"\mathbf{\mu \to \infty}"*("\n")*L"(C(\mu), \varphi(\mu) \to (0, 0)", :right, 10))
    
    annotate!(0.7*Cs_test[end], 1.1*phis[end], text(L"\mathbf{\mu \to \infty}", :right, 13))
    annotate!(0.9*Cs_test[end], 0.975*phis[end], text(L"C(\mu), \varphi(\mu) \to (0, 0)", :right, 13))
    #C(\mu), \varphi(\mu) \to (0, 0)"
    
    plot!(10 .^ (rad*cos.(ψ) .- 1.75), 10 .^ (rad*sin.(ψ) .- 0.6), color=:black, lw=2, labels="")
    scatter!( [10 .^ (rad*cos.(ψ[end]) .- 1.75)], [10 .^ (rad*sin.(ψ[end]) .- 0.6)], color=:black, markershape=:dtriangle, labels="")
    annotate!(10 .^ (rad*cos.(ψ[1]) .- 1.75), 1.1*10 .^ (rad*sin.(ψ[1]) .- 0.6), text(L"\mathbf{\mu}\textrm{\; increases}", :left, :bottom, 14))
    
    
    wid=0.0006
    x1=10^(-1.1)
    for i in 1:size(grad)[1]
    plot!([x1; x1*(10^wid)], [0.04; 0.045], color=grad[i], lw=1.5, labels="")
    x1=x1*(10^wid)
    end
    annotate!(0.15, 0.05, text(L"\textrm{coupling, \; }\mu", :center, 10))
    annotate!(10^(-1.1), 0.0375, text(L"0.2", :center, 10))
    annotate!(x1/(10^(5*wid)), 0.0375, text(L"6.0", :right, 10))
    
    
    xlims!((0.0075, 0.3))
    xlabel!(L"\textrm{V}\textrm{ertical \; displacement, }C(\mu)")
    ylabel!(L"\textrm{P}\textrm{hase \; difference, }\varphi(\mu)")
    plot!(size=(500, 500))
    
end




mus02=[range( 0.2 * 1.1, 1.75, 50); range(1.8, 6, 20)];
mus03=[range( 0.3 * 1.1, 1.75, 50); range(1.8, 6, 20)];
mus04=[range( 0.4 * 1.1, 1.75, 50); range(1.8, 6, 20)];

phis = readdlm( "phi_dw02.out" )
phis03 = readdlm( "phi_dw03.out" )
phis04 = readdlm( "phi_dw04.out" )
Cs = readdlm( "C_dw02.out" )
Cs03 = readdlm( "C_dw03.out" )
Cs04 = readdlm( "C_dw04.out" )

begin
      plot()

      plot!( mus02[2:end], Cs[2:end], lw = 3, color = cols[4], labels = L"\Delta\omega = 0.2" )
      plot!( mus03[3:end], Cs03[3:end], lw = 3, color = cols[9], labels = L"\Delta\omega = 0.3" )
      plot!( mus04[3:end], Cs04[3:end], lw = 3, color = cols[10], labels = L"\Delta\omega = 0.4" )
      
      xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu")
      ylabel!(L"\textrm{V}\textrm{ertical \; displacement, }C(\mu)")
      plot!(size=(900, 350), legend = :topright)
end

savefig("figure16.png")

Ns = [75, 100, 125, 150, 175, 200]

function severalT()
      #Ns = [75, 100, 125, 150, 175, 200]
      Ns = [3, 5, 7, 9, 11 , 13]

      Cs=zeros(size(Ns, 1));
      Cs_test=zeros(size(Ns, 1));
      phis=zeros(size(Ns, 1));

      dw = 0.2
      mu = 0.5
      
      for i in 1:size(Ns,1)
            N = Ns[i]
            Tstar=N*T;
            @time phis[i], Cs[i], Cs_test[i]=getPhiC(mu, Tstar; dw = dw);
            @printf "TSTAR: %f   |||  phase dif:  %f   :::  vert shift: %f\n " Tstar phis[i] Cs[i]
      end
      return phis, Cs, Cs_test
end

phis, Cs, Cs_test = severalT()


begin
      plot() 

      plot!( Ns * T, Cs,  lw = 3, mark=:circle, color= cols[4], yscale = :log10 )

      xlabel!(L"\textrm{T}\textrm{ime \; of \; calculation,} \, T^{*}")
      ylabel!(L"\textrm{V}\textrm{ertical \; displacement, }C(T^{*})")
      plot!(size=(900, 350), legend = :topright)
end
