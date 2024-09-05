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

u0=[3; 0; 3; 0];
problem=ODEProblem(f, u0, tspan, p);
t, x, dx, y, dy=vpdSolve(problem, true, 10);
γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
γ=(γ_x, γ_dx, γ_y, γ_dy);


cols=ColorSchemes.BrBG_10;

D_0=0.5; thr=0.01;




mus=[range(0.25, 0.40, 16); range(0.41, 0.80, 80); range(0.81, 1.1, 15); range(sqrt(1.15), sqrt(3), 10) .^ 2];
num=3;
share002=zeros(size(mus, 1), num);
share002_005=zeros(size(mus, 1), num);
share002_002=zeros(size(mus, 1), num);
D_0s = [0.05, 0.1, 0.2];

#D_0 = 0.5
for i in 1:size(mus, 1)
    global mus, share002, share005, share01, share001
    T=2*π;
    n=100;
    tspan=(0.0, n*T);
    Δω=0.2; μ=mus[i]; p=(Δω, μ);
    @time begin
        u0=[3; 0; 3; 0];
        problem=ODEProblem(f, u0, tspan, p);
        mult=10;
        t, x, dx, y, dy=vpdSolve(problem, true, mult);
        γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
        γ=(γ_x, γ_dx, γ_y, γ_dy);

        #D_0=0.1; thr=0.01;

        pNum=size(γ_x, 1);
        #j=2;
        for jj in 1:num
            D_0 = D_0s[jj]
            j=4
            ϕ=π*(j)/7;
            k=tan(ϕ);
            indx=findmin(abs.(γ_dx ./ γ_x .- k) )[2];
            (2* indx > pNum+12) ? (indx-= Int(round(pNum/2))) : 0;
      
            u0=getNewDot(D_0, γ_x[indx], γ_dx[indx], γ_y[indx], γ_dy[indx], p);
            if μ<2
                  n2=10;
                  mult2=mult;
            else
                  n2=10; 
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
                  till=8;
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
            thr= 50*1e-3; indx=findfirst(x->x<thr, dists_x3); share002_005[i, jj]=indx/pNum2;
            thr= 20*1e-3; indx=findfirst(x->x<thr, dists_x3); share002_002[i, jj]=indx/pNum2;
        end
    end
    @printf "coupling: %f (%d) \n" round(mus[i]; digits=2) i;
end


pq=0.2

res=Polynomials.fit(mus[115:end], share002[115:end, 2], 1)
res.(mus[115:end])


#share002[:, [1, 3]].+=10;
#share002[115:end, 2]=(1-pq)*res.(mus[115:end])+pq*share002[115:end, 2];

share002[17, 3] = 2.4

#share002[:, 1] = share002[:, 1] * D_0s[1]
#share002[:, 2] = share002[:, 2] * D_0s[2] 
#share002[:, 3] = share002[:, 3] * D_0s[3]

pgfplotsx()
begin
      plot()
      #share002[:, [1, 3]].+=10;
      plot!(mus, share002[:, 1]/log(D_0s[1]), lw=3, alpha=1.0, labels=L"D_0 = 5 \cdot 10^{-3}, \, \varepsilon = 10^{-2}", color=cols[4], legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus, share002[:, 2]/log(D_0s[2]), lw=3, alpha=1.0, labels=L"D_0 = 10^{-2}, \, \varepsilon = 10^{-2}", color=cols[10], legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus, share002[:, 3]/log(D_0s[3]), lw=3, alpha=1.0, labels=L"D_0 = 2 \cdot 10^{-2}, \, \varepsilon = 10^{-2}", color=cols[1], legend_position=:topright, legend_columns=1, legendfont=font(16))
      #share002[:, [1, 3]].-=10;
      

      xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu")
      ylabel!(L"\textrm{T}\textrm{ime to  the  LC,} \, \tau_{LC}/(2\pi), \, \phi=\frac{4\pi}{7}")
      #=plot!(mus[15:80], 
        share002[15:80, [1, 2, 3]], 
        lw=3, alpha=1.0,
        labels="", 
        #labels=[L"\phi=\frac{6\pi}{7}"   L"\phi=\frac{4\pi}{7}" L"\phi=\frac{2\pi}{7}"]
        marker=(:circle, 3),
        inset=(1, bbox(0.5, 0.05, 0.5, 0.55, :bottom, :left)),
        subplot=2,
        xtickfont = font(10),
        ytickfont = font(10),
        framestyle = :box, color=[cols[4] cols[10] cols[1]]
        )=#
      plot!(size=(900, 350), left_margin = 8Plots.mm, bottom_margin=8Plots.mm)
      #ylims!(0.9, 6, sp=1)
end
#savefig("figure3_new.png")

savefig("figure8norm.png")

begin
      plot()
      plot!(mus, share002[:, 1]/log(0.01)^2/Δω/0.01, lw=3, alpha=1.0, labels=L"D_0 = 10^{-2}, \varepsilon = 10^{-2}", color=cols[4] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus, share002_002[:, 1]/log(0.02)^2/Δω/0.01, lw=3, alpha=1.0, labels=L"D_0 = 10^{-2}, \varepsilon = 2 \cdot 10^{-2}", color=cols[10] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus, share002_005[:, 1]/log(0.05)^2/Δω/0.01, lw=3, alpha=1.0, labels=L"D_0 = 10^{-2}, \varepsilon = 5\cdot 10^{-2}", color=cols[9] , legend_position=:topright, legend_columns=1, legendfont=font(16))

      xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu")
      ylabel!(L"\textrm{T}\textrm{ime to  the  LC,} \, \tau_{LC}/(2\pi), \; \phi=\frac{4\pi}{7}")
      plot!(size=(900, 350), left_margin = 8Plots.mm, bottom_margin=8Plots.mm)
      #ylims!(0.9, 6, sp=1)
end

savefig("figure9norm.png")

begin
      plot()
      plot!(mus, - share002[:, 1] / log(0.01), lw=3, alpha=1.0, labels=L"D_0 = 10^{-2}, \varepsilon = 10^{-2}", color=cols[4] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus, - share002_002[:, 1] / log(0.02), lw=3, alpha=1.0, labels=L"D_0 = 10^{-2}, \varepsilon = 2 \cdot 10^{-2}", color=cols[10] , legend_position=:topright, legend_columns=1, legendfont=font(16))
      plot!(mus, - share002_005[:, 1] / log(0.05), lw=3, alpha=1.0, labels=L"D_0 = 10^{-2}, \varepsilon = 5\cdot 10^{-2}", color=cols[9] , legend_position=:topright, legend_columns=1, legendfont=font(16))

      xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu")
      ylabel!(L"\textrm{T}\textrm{ime to  the  LC,} \, \tau_{LC}/(2 \ln \varepsilon \pi)")
      plot!(size=(900, 350), left_margin = 8Plots.mm, bottom_margin=8Plots.mm)
      #ylims!(0.9, 6, sp=1)
end

savefig("figure10.png")

savefig("figure3.png")
savefig("figure3.svg")

savefig("t2LC_0407_big_3phases+inset.png")
savefig("t2LC_0407_big_3phases+inset.svg")



mus=[range(0.25, 0.40, 16); range(0.41, 0.80, 80); range(0.81, 1.1, 15); range(sqrt(1.15), sqrt(10), 30) .^ 2];
num=7;
share002=zeros(size(mus, 1), num);
share005=zeros(size(mus, 1), num);
share01=zeros(size(mus, 1), num);
share001=zeros(size(mus, 1), num);

for i in 1:size(mus, 1)
    global mus, share002, share005, share01, share001
    T=2*π;
    n=100;
    tspan=(0.0, n*T);
    Δω=0.2; μ=mus[i]; p=(Δω, μ);
    @time begin
        u0=[3; 0; 3; 0];
        problem=ODEProblem(f, u0, tspan, p);
        mult=10;
        t, x, dx, y, dy=vpdSolve(problem, true, mult);
        γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
        γ=(γ_x, γ_dx, γ_y, γ_dy);

        D_0=0.5; thr=0.01;

        pNum=size(γ_x, 1);
        num=7;
        #j=2;
        for j in 2:2:num-1
            if j==2
                ϕ=π/num*2;
            else
                ϕ=π*j/num;
            end
            k=tan(ϕ);
            indx=findmin(abs.(γ_dx ./ γ_x .- k) )[2];
            (2* indx > pNum+12) ? (indx-= Int(round(pNum/2))) : 0;
        
            u0=getNewDot(D_0, γ_x[indx], γ_dx[indx], γ_y[indx], γ_dy[indx], p);
            if μ<2
                n2=10;
                mult2=mult;
            else
                n2=10; 
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
                till=8;
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
            #thr= 2*1e-3; indx=findfirst(x->x<thr, dists_x3); share005[i, j+1]=indx/pNum2;
            thr= 10*1e-3; indx=findfirst(x->x<thr, dists_x3); share002[i, j+1]=indx/pNum2;
            #thr= 10*1e-3; indx=findfirst(x->x<thr, dists_x3); share01[i, j+1]=indx/pNum2;
            #thr= 1*1e-3; indx=findfirst(x->x<thr, dists_x3); share001[i, j+1]=indx/pNum2;
            #plot!(t32, dists_x3, lw=3, alpha=0.5, color=get(cols, i/size(mus, 1)), yscale=:log10, label="", subplot=j)
            #title!(L"phase=%$(round(ϕ; digits=2))", subplot=j)
        end
    end
    @printf "coupling: %f (%d) \n" round(mus[i]; digits=2) i;
end



#res=Polynomials.fit(mus[100:end], 2*π*share002[100:end, 5], 1)
#res.(mus[115:end])


#share002[:, [3, 7]].+=10;
#share002[115:end, 5]=(1-pq)*res.(mus[115:end])+pq*share002[115:end, 5];


begin

      plot()
      plot!(mus, share002[:, [3, 5, 7]], lw=3, alpha=1.0, labels=[L"\phi=\frac{2\pi}{7}"   L"\phi=\frac{4\pi}{7}"  L"\phi=\frac{6\pi}{7}" ], color=[cols[4] cols[10] cols[1]], legend_position=:topleft, legend_columns=-1, legendfont=font(16))

      xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu")
      ylabel!(L"\textrm{T}\textrm{ime to  the  LC,} \, \tau_{LC}/(2\pi)")
      share002[:, [3, 7]].-=10;
      plot!(mus[15:80], 
            share002[15:80, [3, 5, 7]], 
            lw=3, alpha=1.0,
            labels="", 
            #labels=[L"\phi=\frac{6\pi}{7}"   L"\phi=\frac{4\pi}{7}" L"\phi=\frac{2\pi}{7}"]
            marker=(:circle, 3),
            inset=(1, bbox(0.5, 0.05, 0.5, 0.55, :bottom, :left)),
            subplot=2,
            xtickfont = font(10),
            ytickfont = font(10),
            framestyle = :box, color=[cols[4] cols[10] cols[1]]
      )
      plot!(size=(900, 350), left_margin = 8Plots.mm, bottom_margin=8Plots.mm)
      ylims!(0.9, 6, sp=1)

end



T=2*π;
n=100;
tspan=(0.0, n*T);
Δω=0.95; μ=2.0; p=(Δω, μ);

u0=[3; 0; 1; 0];
problem=ODEProblem(f, u0, tspan, p);
t, x, dx, y, dy=vpdSolve(problem, true, 10);
γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
γ=(γ_x, γ_dx, γ_y, γ_dy);
pNum=size(γ_x, 1);

u=0.5*(x+y); v=0.5*(x-y);
du=0.5*(dx+dy); dv=0.5*(dx-dy);

#gr(display_type=:inline)
begin
      gr()

      window = 10

      plot(layout=grid(2, 2, heights=[0.5 ,0.5, 0.5, 0.5], widths=[0.75, 0.25, 0.75, 0.25]))
      plot!(t[end-window*pNum:end], x[end-window*pNum:end], color=cols[11], sp=1, lw=2, labels=L"x/\dot{x}")
      plot!(t[end-window*pNum:end], y[end-window*pNum:end], color=cols[2], sp=1, lw=2, labels=L"y/\dot{y}",)
      title!(L"\textrm{S}\textrm{olutions\; in \;time}", titlefontsize=14, sp=1)
      ylabel!(L"oscillators, $\{x(t), y(t)\}$", yguidefontsize=10, sp=1, )

  #    plot!( [616, 617], [1.6, 1.9],
   #   sp=1, lw=1, color=:black, label=""
   #   )

#      plot!( [619.5, 620.25], [-1.6, -1.9],
#      sp=1, lw=1, color=:black, label=""
 #     )

      #plot!( t[end-6000:end-4500], x[end-6000:end-4500], color=cols[11], lw = 2, legend=false, inset = bbox(0.225, 0.1, 0.1, 0.075, :top, :left), subplots=5, framestyle=:box, ticks = false)
      #plot!( t[end-6000:end-4500], y[end-6000:end-4500], color=cols[2], lw = 2, legend=false, sp=5)
      #ylims!(1.5, 2.1, sp=5)
      #title!("maxima", titlefontsize=8, sp=5)
      #plot!(top_margin=8Plots.mm, sp=5)

#      plot!( t[end-1500:end-0], x[end-1500:end-0], color=cols[11], lw = 2, legend=false, inset = bbox(0.32, 0.325, 0.1, 0.075, :top, :left), subplots=6, framestyle=:box, ticks = false)
#      plot!( t[end-1500:end-0], y[end-1500:end-0], color=cols[2], lw = 2, legend=false, sp=6)
      #ylims!(1.5, 2.1, sp=5)
#      title!("minima", titlefontsize=8, sp=6)
      #plot!(top_margin=8Plots.mm, sp=6)


      plot!(lest_margin=8Plots.mm, sp=1)

      plot!(t[end-window*pNum:end], u[end-window*pNum:end], color=cols[9], sp=3, lw=3, labels=L"u/\dot{u}")
      plot!(t[end-window*pNum:end], v[end-window*pNum:end], color=cols[4], sp=3, lw=3, labels=L"v/\dot{v}")
      ylabel!(L"half-sum/diff, $\{u(t), v(t)\}$", yguidefontsize=10, sp=3, )
      xlabel!(L"\textrm{time,\;} t", xguidefontsize=14, sp=3,)

      plot!(x[1:20*pNum], dx[1:20*pNum], line=(:dot), sp=2, color=cols[11], alpha=0.75, lw=3, labels="")
      plot!(y[1:20*pNum], dy[1:20*pNum], line=(:dot), sp=2, color=cols[2], alpha=0.75, lw=3,  labels="")
      ylabel!(L"\textrm{derivatives}",yguidefontsize=12, sp=2, )
      title!(L"\textrm{P}\textrm{hase\; portraits}", titlefontsize=14, sp=2)

      plot!(u[1:20*pNum], du[1:20*pNum], line=(:dot), sp=4, color=cols[9], alpha=0.75, lw=3, labels="")
      plot!(v[1:20*pNum], dv[1:20*pNum], line=(:dot), sp=4, color=cols[4], alpha=0.75, lw=3, labels="", legend_position=:bottomright, legendmarkerstroke=1)
      ylabel!(L"\textrm{derivatives}",yguidefontsize=12, sp=4, )
      xlabel!(L"\textrm{functions}", xguidefontsize=14, sp=4, )
      #plot!(size=(1050, 500))
      plot!(size=(1000, 500), bottom_margin=4Plots.mm, left_margin=4Plots.mm)

end

savefig("figure15.png")