#mandelbrot_module.jl
module MOD_mandel_PJS
export fn_mandel
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
end #end of module
