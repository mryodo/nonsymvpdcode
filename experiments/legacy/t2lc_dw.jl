using LinearAlgebra, Arpack
using BenchmarkTools
#math
using DifferentialEquations
using Interpolations
using Peaks
using Polynomials
using Statistics
#plots
using Plots
pgfplotsx()
theme(:mute)
using ColorSchemes
cols=ColorSchemes.Spectral_11;
Plots.scalefontsizes(1.5);

using Polynomials


using Printf
using LaTeXStrings

include("vpdModule.jl")


function getTau(dists_x2, t2, D_0, thr, pNum, multiplier)
      tail=Int(round(size(dists_x2, 1)/10));
      tail=5*pNum;
      m=mean(dists_x2[end-tail:end]);
      thr0=m+multiplier*(maximum(dists_x2[end-tail:end])-m);
      
      indx=findfirst(x->x<thr, dists_x2);
      if indx<pNum
          @printf "SHIT SHIT SHIT \n"
      end
      x=t2[indx+1-(pNum ÷ 2)*min(indx ÷ pNum, multiplier):indx];
      y=log10.( dists_x2[indx+1-(pNum ÷ 2)*min(indx ÷ pNum, multiplier):indx]);
      ks2=coef(x, y);
      return ks2, -1/ks2[2]/2/π, indx
end
  
function getDistsSmooth(dists_x, t, pNum)
      dists_x_sm=[mean((dists_x[maximum([1, i - pNum ÷ 2]) : minimum([i+ pNum ÷ 2, size(dists_x, 1)])])) for i in 1:size(dists_x,1)];
     
      t_sm=t;
      return dists_x_sm, t_sm
end
  
function coef(x, y)
      n=size(x, 1);
      xy=sum(x .* y); sx=sum(x); sy=sum(y);
      sx2=sum(x .* x);
      noma=sy*sx2-sx*xy;
      nomb=n*xy-sum(x)*sum(y);
      denom=n*sx2-sum(x)*sum(x);
      return [noma/denom; nomb/denom]
      
end
  


T=2*π;
n=100;
tspan=(0.0, n*T);
Δω=0.2; μ=0.6; p=(Δω, μ);

function form_mu( dw )
      return [range(dw*1.22, dw*2.0, 32); range(dw*2.05, dw*4.0, 20); range(dw*4.05, dw * 6.0, 15); range(sqrt(dw *6.05), sqrt(5.0), 20) .^ 2] 
end




num=4;
dws = [0.1, 0.2, 0.5, 0.75 ];
mus= [ form_mu( dws[i] ) for i in 1:num ]
share002=zeros(size(mus[1], 1), num);

#D_0 = 0.5
for jj in 3:3
    global mus, share002
      for i in 1:size(mus[1], 1)
            T=2*π;
            n=100;
            tspan=(0.0, n*T);
            Δω=dws[jj]; μ=mus[jj][i]; p=(Δω, μ);
            
            u0=[3; 0; 3; 0];
            problem=ODEProblem(f, u0, tspan, p);
            mult=10;
            t, x, dx, y, dy=vpdSolve(problem, true, mult);
            γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
            γ=(γ_x, γ_dx, γ_y, γ_dy);

            #D_0=0.1; thr=0.01;

            pNum=size(γ_x, 1);
            #j=2;
            D_0 = 0.5
            if jj == 3
                  j = 3.25
            else 
                  j = 3.5
            end
            ϕ=π*(j)/7;
            k=tan(ϕ);
            indx=findmin(abs.(γ_dx ./ γ_x .- k) )[2];
            (2* indx > pNum+12) ? (indx-= Int(round(pNum/2))) : 0;
      
            u0=getNewDot(D_0, γ_x[indx], γ_dx[indx], γ_y[indx], γ_dy[indx], p);
            if μ<2
                  n2=10;
                  mult2=mult;
            else
                  n2=30; 
                  mult2=mult;
            end
            tspan2=(0.0, n2*T);
            problem2=ODEProblem(f, u0, tspan2, p);
            t2, x, dx, y, dy=vpdSolve(problem2, true, mult2);
      
            γ_x2, γ_dx2, γ_y2, γ_dy2=getLimCycleNaive(t2, x, dx, y, dy);
            pNum2=size(γ_x2, 1);
            
            if μ<2
                  till=8;
            else
                  till=16;
            end
      
            radius=20;
            dists_x, dists_y=getDists(x[1:Int(round(till*size(x, 1)/n2))],
                  dx[1:Int(round(till*size(x, 1)/n2))],
                  y[1:Int(round(till*size(x, 1)/n2))],
                  dy[1:Int(round(till*size(x, 1)/n2))],
                  γ, radius
            );
      
            dists_x3, t32=getDistsSmooth(dists_x, t2[1:Int(round(till*size(x, 1)/n2))], pNum2);
            dists_x3 /= D_0; 
            dists_x3, t32=getDistsSmooth(dists_x3, t32, pNum2);
            thr= 10*1e-3; indx=findfirst(x->x<thr, dists_x3); share002[i, jj]=indx/pNum2;
            #thr= 50*1e-3; indx=findfirst(x->x<thr, dists_x3); share002_005[i, jj]=indx/pNum2;
            #thr= 20*1e-3; indx=findfirst(x->x<thr, dists_x3); share002_002[i, jj]=indx/pNum2;
            @printf "coupling: %f (%d) \n" round(mus[jj][i]; digits=3) jj;
        end
    
end

using DelimitedFiles
writedlm("t2lc_dws_2.out", share002)


share002 = readdlm("t2lc_dws.out")

pgfplotsx()
begin
      plot()
      plot!(mus[1], share002[:, 1], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.1", color=cols[4] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus[2], share002[:, 2], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.2", color=cols[10] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus[3], share002[:, 3], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.5", color=cols[9] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus[4], share002[:, 4], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.9", color=cols[1] , legend_position=:topright, legend_columns=1, legendfont=font(16))

      xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu", sp=1)
      ylabel!(L"\textrm{T}\textrm{ime to  the  LC,} \, \tau_{LC}/(2\pi), \; \phi=\frac{\pi}{2}", sp=1)

      plot!(mus[1][end-20:end], 
        share002[end-20:end, 1]/dws[1], 
        lw=3, alpha=1.0,
        labels="", 
        marker=(:circle, 3),
        inset=(1, bbox(0.075, 0.4, 0.5, 0.65, :bottom, :left)),
        subplot=2,
        xtickfont = font(10),
        ytickfont = font(10),
        yticks = [],
        framestyle = :box, color=cols[4],
      )
      plot!(mus[2][end-20:end], share002[end-20:end, 2]/dws[2], 
        lw=3, alpha=1.0, labels="",  marker=(:circle, 3), color=cols[10], sp=2
      )
      plot!(mus[3][end-30:end], 1.25*share002[end-30:end, 3]/dws[3] .- 2.5, 
        lw=3, alpha=1.0, labels="",  marker=(:circle, 3), color=cols[9], sp=2
      ) 
      plot!(mus[4][end-35:end], share002[end-35:end, 4]/dws[4] .- 7, 
        lw=3, alpha=1.0, labels="",  marker=(:circle, 3), color=cols[1], sp=2
      )  
      ylabel!(L"\tau_{LC}/(2\pi \Delta\omega)" , sp = 2, font = 8)


      plot!(size=(900, 350), left_margin = 8Plots.mm, bottom_margin=8Plots.mm)
      #ylims!(0.9, 6, sp=1)
end
savefig("figure20new.png")

begin
      plot()
      plot!(mus[1]/dws[1] .- 1, share002[:, 1]/dws[1], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.1", color=cols[4] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus[2]/dws[2] .- 1, share002[:, 2]/dws[2], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.2", color=cols[10] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus[3]/dws[3] .- 1, share002[:, 3]/dws[3], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.5", color=cols[9] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus[4]/dws[4] .- 1, share002[:, 4]/dws[4], lw=3, alpha=1.0, labels=L"\Delta\omega = 0.9", color=cols[1] , legend_position=:topright, legend_columns=1, legendfont=font(16))

      
      xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu")
      ylabel!(L"\textrm{T}\textrm{ime to  the  LC,} \, \tau_{LC}/(2\pi), \; \phi=\frac{4\pi}{7}")
      plot!(size=(900, 350), left_margin = 8Plots.mm, bottom_margin=8Plots.mm, xscale =:log10)
      #ylims!(0.9, 6, sp=1)
end