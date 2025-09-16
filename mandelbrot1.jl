using Colors, Images, FileIO, ColorSchemes, Plots
using Base.Threads
pyplot()
#******************************************************************
function fn_mandel(cnum,max_iter,iter_pow)
    zn=0.0+0.0im
    for cnt in 1:max_iter
        zn1=zn^iter_pow+cnum
        magz=abs(zn1)
        if magz>2.0
            return cnt
        end
        zn=zn1
    end
    return max_iter
end
#******************************************************************
println("Enter iteration power: ")
iter_pow=parse(Float64,readline())
println("Enter colour scale power: ")
col_pow=parse(Float64,readline())
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
stime=time()
@threads for j in 1:ny
    for i in 1:nx
        cnum=x[i]+y[j]*im
        mandelmat[j,i]=fn_mandel(cnum,max_iter,iter_pow)
    end
end
elapsed=time()-stime
println("Time taken=$elapsed seconds")
yout=vcat(-1.0.*reverse(y[2:end,1]),y); #reflect y to create full range
finalmat=(vcat(reverse(mandelmat[2:end,:],dims=1),mandelmat)./max_iter).^col_pow #add y reflected mandelbrot set and scale
p=heatmap(x,yout,finalmat, c=:viridis, aspect_ratio=:equal)
display(p)
savefig(p,"mandelbrot1a.png")
#******************************************************************
#img = [get(ColorSchemes.coolwarm, finalmat[j,i]) for j in 1:length(yout), i in 1:nx]
img = [get(ColorSchemes.viridis, finalmat[j,i]) for j in 1:length(yout), i in 1:nx]
rgb_img = colorview(RGB, permutedims(img, (1, 2)))
save("mandelbrot.png", rgb_img)
a4_width, a4_height=3508, 2480
# Resize the image to A4 landscape dimensions
resized_img = imresize(rgb_img, (a4_width, a4_height))
save("mandelbrot_a4_landscape.png", resized_img)
#******************************************************************
println("End")
