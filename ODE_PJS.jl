#Define analytic function
function ode_analy(x,a,fa,D,omega)
    so=sqrt(omega)
    return (D/so).*sinh.(so.*(x.-a))+fa.*cosh.(so.*(x.-a))
end
#********************************************************
#Set domain
#********************************************************
a=1.0; b=2.0
#********************************************************
#Set boundary conditions
#********************************************************
fa=1.0; D=10.0; omega=10.0
n=[11:10:101;]
nN=length(n)
dx=zeros(nN,1)
err=zeros(nN,1)
#********************************************************
#Cycle through discretisations
#********************************************************
for j in 1:nN
    dx[j]=(b-a)/(n[j]-1)
    xi=collect(range(a,b,n[j]))
    #Calculate numerical solution
    fnum=zeros(n[j],1)
    fnum[1]=fa
    for i in 1:n[j]-1
        if i==1
            fnum[i+1]=0.5*fa*(2.0+omega*dx[j]^2)+D*dx[j]
        else
            fnum[i+1]=fnum[i]*(2.0+omega*dx[j]^2)-fnum[i-1]
        end
    end
    #Calculate analytic solution
    fan=ode_analy(xi,a,fa,D,omega)
    err[j]=maximum(abs.(fnum.-fan))
end