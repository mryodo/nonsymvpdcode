
function f(du, u, p, t)
    Δω, μ=p;
    du[1]=u[2];
    du[2]=(1-u[1]^2)*u[2]-(1-Δω)*u[1]-μ*(u[2]-u[4]);
    du[3]=u[4];
    du[4]=(1-u[3]^2)*u[4]-(1+Δω)*u[3]-μ*(u[4]-u[2]);
end


function getDD(x::Float64, dx::Float64, y::Float64, dy::Float64, p::Tuple{Float64, Float64})
    Δω, μ=p;
    ddx=(1-x^2)*dx-(1-Δω)*x-μ*(dx-dy);
    ddy=(1-y^2)*dy-(1+Δω)*y-μ*(dy-dx);
    return ddx, ddy
end

function getNewDot(eps::Float64, x::Float64, dx::Float64, y::Float64, dy::Float64, p::Tuple{Float64, Float64})
    ddx, ddy=getDD(x, dx, y, dy, p);
    dir_x=[-ddx; dx]; dir_y=[-ddy; dy];
    dir_x/=norm(dir_x, 2); dir_y/=norm(dir_y, 2);
    (dir_x[2]<0) && (dir_x*=-1);
    (dir_y[2]<0) && (dir_y*=-1); 
    return [x+eps*dir_x[1]; dx+eps*dir_x[2]; y+eps*dir_y[1]; dy+eps*dir_y[2]]
end


@inline dsquare(x,y) = (x[1]-y[1])^2 + (x[2]-y[2])^2;

function vpdSolve(problem::ODEProblem, interp::Bool, mult::Int64)    
    sol=solve(problem, Tsit5(); reltol=1e-12, abstol=1e-12);
    x=reduce(hcat, sol.u)[1, :];
    dx=reduce(hcat, sol.u)[2, :];
    y=reduce(hcat, sol.u)[3, :];
    dy=reduce(hcat, sol.u)[4, :];
    t=sol.t;
    tspan=problem.tspan;
    if interp 
        t_new=range(tspan[1], tspan[2], mult*size(t, 1));
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

function vpdSolve2(problem::ODEProblem, interp::Bool, mult::Int64)    
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



function getLimCycleNaive(t, x, dx, y, dy)
    inds_x=findminima(x)[1];
    γ_x=x[inds_x[end-1]:inds_x[end]];
    γ_dx=dx[inds_x[end-1]:inds_x[end]];

    inds_y=findminima(y)[1];
    γ_y=y[inds_y[end-1]:inds_y[end]];
    γ_dy=dy[inds_y[end-1]:inds_y[end]];
    return γ_x, γ_dx, γ_y, γ_dy
end

function minimumdistance(x,y)
    jatom = 0
    dmin = +Inf
    @inbounds for j in 1:size(y, 1)
      d = dsquare(x,y[j, :])
      if d < dmin
        jatom = j
        dmin = d
      end
    end
    return sqrt(dmin), jatom
end

function getDists(x, dx, y, dy, γ, radius::Int)
    γ_x, γ_dx, γ_y, γ_dy=γ;
    dists_x=zeros(size(x, 1));
    dists_y=zeros(size(y, 1));

    minval, ind=minimumdistance([x[1] dx[1]], [γ_x γ_dx]);
    dists_x[1]=minval;
    gx3=[γ_x; γ_x; γ_x];
    gdx3=[γ_dx; γ_dx; γ_dx];
    shift=size(γ_x, 1);
    for j in 2:size(x, 1)
        minval, ind_new=minimumdistance([x[j] dx[j]], [gx3[shift+ind-radius:shift+ind+radius] gdx3[shift+ind-radius:shift+ind+radius]]);
        dists_x[j]=minval;
        ind=ind-radius+ind_new;
        (ind > shift) ? (ind=ind-shift) : 0; 
    end

    minval, ind=minimumdistance([y[1] dy[1]], [γ_y γ_dy]);
    dists_y[1]=minval;
    gy3=[γ_y; γ_y; γ_y];
    gdy3=[γ_dy; γ_dy; γ_dy];
    shift=size(γ_y, 1);
    for j in 2:size(y, 1)
        minval, ind_new=minimumdistance([y[j] dy[j]], [gy3[shift+ind-radius:shift+ind+radius] gdy3[shift+ind-radius:shift+ind+radius]]);
        dists_y[j]=minval;
        ind=ind-radius+ind_new;
        (ind > shift) ? (ind=ind-shift) : 0; 
    end

    return dists_x, dists_y
end
