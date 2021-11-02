using Combinatorics

function compute_c(k::Int64,r::Int64,n::Int64,j::Int64,tilde=false)
#computes individual c_values, as elements of the Polynomial reconstruction C matrix
# k Stencil size, Order of the reconstruction polynomial + 1
# r Rightshift of the stencil or position of the reconstruction cell within the stencil from the left, r=0 left most cell
# n derivative in Question / row within C matrix n=0 => first row
# j counter / collumn within C matrix j=0 => first collumn
# tilde control reconstruction position left cell interface for tilde=false, right interface tilde=true 
c=0;
for l in j+1:k
    temp_f=(-1)^(k-l)*factorial(n+1)/factorial(l)/factorial(k-l);
    temp_s=0;
    S=r+1 -tilde .-vcat(0:l-1,l+1:k)
    I=Combinatorics.powerset(S,k-n-1,k-n-1);
    for set in I
        temp_s=temp_s+prod(set);
    end
    c=c+temp_f*temp_s;
end
return c
end

function compute_C(k::Int64,r::Int64,tilde=false)
# computes Polynomial reconstruction C matrices
# k Stencil size, Order of the reconstruction polynomial + 1
# r Rightshift of the stencil or position of the reconstruction cell within the stencil from the left, r=0 left most cell
# tilde control reconstruction position left cell interface for tilde=false, right interface tilde=true
    C=zeros(k,k);
    for n in 0:k-1, j in 0:k-1
        C[n+1,j+1]=compute_c(k,r,n,j,tilde);
    end
    return C;
end

function compute_C_constants(n_max::Int64,tilde=false)
    # computes the whole set of polynomial reconstruction C matrices up to degree n_max
    # tilde control reconstruction position left cell interface for tilde=false, right interface tilde=true
    Konstanten=Array{Array{Float64,2},2}(undef,n_max,n_max);
    for n in 1:n_max, r in 0:n-1
        Konstanten[n,r+1]=compute_C(n,r,tilde);
    end
    return Konstanten;
end

function compute_B(C_rk)
# computes Smoothness indicator matrix B_r^k, degree k, rightshift r
# r Rightshift of the stencil or position of the reconstruction cell within the stencil from the left, r=0 left most cell
    k=size(C_rk,1);  # k Stencil size, Order of the reconstruction polynomial + 1
    B=zeros(Float64,k,k);
    for j1 in 0:k-1, j2 in 0:k-1
        for l in 1:k-1, n1 in l:k-1, n2 in l:k-1
            B[j1+1,j2+1] += (-1)^(n1+n2-2*l)/(factorial(n1-l)*factorial(n2-l)*(n1+n2-2*l+1))*C_rk[n1+1,j1+1]*C_rk[n2+1,j2+1];
        end
    end
    return B
end

function compute_D(k::Int64,C::Array{Array{Float64,2},2})
# Computes WENO prototype weights for stencils of length k / optimal weights for smooth data reconstruction
# k Stencil size, Order of the reconstruction polynomial + 1
# C Array of reconstruction coefficients as computed by Compute_C_constants, at least to order 2*k-1
    M=zeros(k,k);
    for z=1:k,s=z:k
        M[z,s]=C[k,s][1,s-z+1];
    end
    D=zeros(k,1);
    for z=1:k
        D[z]=C[2*k-1,k][1,k-z+1]
    end
    return M\D;
end

function compute_D_constants(n_max::Int64,C)
# Computes all WENO Prototype constants up to order n_max
# input C is the WENO reconstruction constants as computed by compute_C_constants
    D=Array{Array{Float64,2},1}(undef,n_max);
    for k=1:n_max
        D[k]=compute_D(k,C);
    end
    return D
end

# wip marker, got till here

function compute_B_constants(n_max::Int64,C)
# computes the complete set of Smoothness indicator matrices up to degree n_max
# C the Reconstruction coefficients as computed by compute_C_constants
    B_constants=Array{Array{Float64,2},2}(undef,n_max,n_max);
    for n in 1:n_max, r in 0:n-1
        B_constants[n,r+1]=compute_B(C[n,r+1]);
    end
    return B_constants;
end

struct PR_constants
    n::Int64;
    C::Array{Array{Float64,2},2};
    C_tilde::Array{Array{Float64,2},2};
    B::Array{Array{Float64,2},2};
    D::Array{Array{Float64,2},1};
    D_tilde::Array{Array{Float64,2},1};
    function PR_constants(n::Int64)
        C=compute_C_constants(n);
        C_tilde=compute_C_constants(n,true);
        B=compute_B_constants(n,C);
        D=compute_D_constants((n+1)÷2,C);
        D_tilde=compute_D_constants((n+1)÷2,C_tilde);
        new(n,C,C_tilde,B,D,D_tilde);
    end
end

export PR_constants
