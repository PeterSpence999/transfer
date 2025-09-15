using Colors, Images, FileIO, ColorSchemes, Plots
using Base.Threads
gr()
#******************************************************************
function fn_mandel(cnum,max_iter)
    zn=0+0im
    for cnt in 1:max_iter
        zn1=zn^2+cnum
        magz=abs(zn1)
        if magz>2.0
            return cnt
        end
        zn=zn1
    end
    return max_iter
end
#******************************************************************
max_iter=511
nx,ny=2001, 1001
#xs,xe=0.375, 0.4
#ys,ye=0.25, 0.29
xs,xe=-2.0, 2.0
ys,ye=0.0,2.0
dx=(xe-xs)/(nx-1)
dy=(ye-ys)/(ny-1)
x=range(xs,xe,length=nx)
y=range(ys,ye,length=ny)
mandelmat=Array{Int}(undef,ny,nx) #matrix of iterations
@threads for j in 1:ny
    for i in 1:nx
        cnum=x[i]+y[j]*im
        mandelmat[j,i]=fn_mandel(cnum,max_iter)
    end
end
#mandelmat_r=reverse(mandelmat[2:end,:],dims=1)
yout=vcat(-1.0.*reverse(y[2:end,1]),y); #reflect y
finalmat=vcat(reverse(mandelmat[2:end,:],dims=1),mandelmat) #add reflected mandelbrot set
#heatmap(mandelmat, c=:viridis, clims=(0,max_iter))
#heatmap(x,y,mandelmat_r, c=:plasma, xlabel= 'x', ylabel='y')
#mmin=log10(minimum(finalmat))
#mmax=log10(maximum(finalmat))
#clevels=10.0.^range(mmin,mmax,length=max_iter);
p=heatmap(x,yout,finalmat, c=:coolwarm)
#contour!(p,x,yout,finalmat,levels=clevels, c=:white, linewidth=1.0)
display(p)
savefig(p,"mandelbrot1a.png")
#******************************************************************
# Map iterations → colors
#img = [get(ColorSchemes.coolwarm, sqrt(mandelmat[j,i] / max_iter)) for j in 1:ny, i in 1:nx]
img = [get(ColorSchemes.coolwarm, sqrt(finalmat[j,i] / max_iter)) for j in 1:length(yout), i in 1:nx]
rgb_img = colorview(RGB, permutedims(img, (1, 2)))
save("mandelbrot.png", rgb_img)
a4_width, a4_height=3508, 2480
# Resize the image to A4 landscape dimensions
resized_img = imresize(rgb_img, (a4_width, a4_height))
save("mandelbrot_a4_landscape.png", resized_img)
#******************************************************************
println("End")
