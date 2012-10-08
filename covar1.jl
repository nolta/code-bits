# Cholesky decomposition of covariance matrix R, where
# R_{ij} = r(|j-i|) = \sum_{m=1}^k a(m) b(m)^|j-i| is real for i,j = 1,2,...
# and |b(m)|<1 for m = 1,...,k.

function covchol(a,b,n)
    if length(a)!=length(b)
        error("covchol: a,b must be vectors of same length")
    end
    k=length(a)
    d=zeros(Complex128,k,n)
    d[:,1]=a./sqrt(real(sum(a)))
    accum=(d[:,1].*b)*(d[:,1].*b)'
    # Proceed row by row to factor
    bb=b*b'
    for m=2:n
        v=a-sum(accum,2)
        d[:,m]=v./sqrt(real(sum(v)))
        accum=(accum+d[:,m]*d[:,m]').*bb
    end
    dsum=real(sum(d,1))
    dsuminv=1./dsum
    cholsav=(k,n,copy(b),d,dsum,dsuminv)
end


# Matrix-vector product U' x = y, where R = U' U is the covariance matrix
# described above and cholsav was produced by covchol above.  The number
# of unknowns n can be less than or equal to the size for which covchol
# was created.

function covprodut(x,cholsav)
    n=length(x)
    (k,nsav,b,d,dsum,dsuminv)=cholsav
    if n>nsav
        error("covprodu: length(x) greater than that provided in cholsav")
    end
    ss=zeros(Complex128,k)
    y=zeros(Float64,n)
    for i=1:n
        y[i]=x[i]*dsum[i]+real(sum(ss))
        ss=(ss+x[i]*conj(d[:,i])).*conj(b)
    end
    return y
end


# Solve linear system U' x = y, where R = U' U is the covariance matrix
# described above and cholsav was produced by covchol above.  The number
# of unknowns n can be less than or equal to the size for which covchol
# was created.

function covsolut(y,cholsav)
    n=length(y)
    (k,nsav,b,d,dsum,dsuminv)=cholsav
    if n>nsav
        error("covsolu: length(y) greater than that provided in cholsav")
    end
    ss=zeros(Complex128,k)
    x=zeros(Float64,n)
    for i=1:n
        x[i]=(y[i]-real(sum(ss)))*dsuminv[i]
        ss=(ss+x[i]*conj(d[:,i])).*conj(b)
    end
    return x
end


const a=[10.0+0.0im;
   0.194487+0.405512im;
   0.194487-0.405512im;
   -0.4358-0.0374477im;
   -0.4358+0.0374477im;
   0.4986+0.31128im;
   0.4986-0.31128im;
   0.385488-0.00129318im;
   0.385488+0.00129318im;
   -0.283494-0.291219im;
   -0.283494+0.291219im]
const b=[0.12+0.0im;
   -0.320372+0.797491im;
   -0.320372-0.797491im;
   0.720776+0.102379im;
   0.720776-0.102379im;
   0.370054-0.0357288im;
   0.370054+0.0357288im;
   -0.652465+0.506429im;
   -0.652465-0.506429im;
   0.696761+0.622623im;
   0.696761-0.622623im]
const n=1000000
const x=rand(n)

tic()
const cholsav=covchol(a,b,n)
toc()
tic()
const y=covprodut(x,cholsav)
toc()
tic()
const z=covsolut(y,cholsav)
println(norm(x-z)/norm(z))
toc()
