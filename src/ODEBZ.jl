module ODEBZ

using DifferentialEquations, ParameterizedFunctions

FLorenz = @ode_def LorenzExample begin
    dx = σ*(y-x)
    dy = x*(ρ-z) - y
    dz = x*y - β*z
  end σ ρ β


  f(x,y,d) = x*y/d
  NJ = @ode_def FuncTest begin
    dx = a*x - b*x*y
    dy = -c*y + f(x,y,d)
  end a b c d


y_tilde(x, z, v, α, k₁, k₂, k₆, kf, V₀, X₀ ,Y₀, Z₀, A, H) = 
((α*k₆*Z₀*V₀*z*v)/(k₁*H*X₀*x+k₂*A*H^2+kf))/Y₀

X₀(k₂, k₅,  A, H)= (k₂*A*H^2)/k₅
Y₀(k₂, k₅, A, H) = 4*X₀(k₂, k₅, A, H)  
Z₀(A, C, M) = (C*A)/(40*M)
V₀(A, C, H, M) = 4*A*H*C/M^2
T₀(k₂, A, C, H) = 1/(10*k₂*A*H*C)

y_tilde(x, z, v, α, k₁, k₂, k₅, k₆, kf, A, C, H, M) = y_tilde(x, z, v, α, k₁, k₂, k₆, kf,
                                                          V₀(A, C, H, M),
                                                          X₀(k₂, k₅,  A, H), 
                                                          Y₀(k₂, k₅, A, H), 
                                                          Z₀(A, C, M), 
                                                          A, H)

FBZ = @ode_def BZExample begin
    dx = T₀(k₂, A, C, H)*
        (-k₁*H*Y₀(k₂, k₅, A, H)*x*y_tilde(x, z, v, α, k₁, k₂, k₅, k₆, kf, A, C, H, M)+
        k₂*A*H^2*4*y_tilde(x, z, v, α, k₁, k₂, k₅, k₆, kf, A, C, H, M)
        -2*k₃*X₀(k₂, k₅,  A, H)*x^2 + 0.5*k₄*A^0.5*H^1.5*X₀(k₂, k₅,  A, H)^(-0.5)*
        (C-Z₀(A, C, M)*z)*x^(0.5)-0.5*k₅*Z₀(A, C, M)*x*z-kf*x)
    dz = T₀(k₂, A, C, H)*(k₄*A^0.5*H^(1.5)*X₀(k₂, k₅,  A, H)^(0.5)*(C/Z₀(A, C, M)-z)*
      x^(0.5) - k₅*X₀(k₂, k₅,  A, H)*x*z - α*k₆*V₀(A, C, H, M)*z*v - β*k₇*M*z - kf*z)
    dv = T₀(k₂, A, C, H)*(2*k₁*H*X₀(k₂, k₅,  A, H)*Y₀(k₂, k₅, A, H)/V₀(A, C, H, M)*x*y_tilde(x, z, v, α, k₁, k₂, k₅, k₆, kf, A, C, H, M)  
    +k₂*A*H^2*Y₀(k₂, k₅, A, H)/V₀(A, C, H, M)*y_tilde(x, z, v, α, k₁, k₂, k₅, k₆, kf, A, C, H, M) 
    +k₃*X₀(k₂, k₅,  A, H)^2/V₀(A, C, H, M)*x^2 - α*k₆*Z₀(A, C, M)*z*v - kf*v)
end α β k₁ k₂ k₃ k₄ k₅ k₆ k₇ kf A C H M 

p = [333.3; 0.2609; 4.0e06; 2.0; 3000; 55.2; 7000; 0.09; 0.23; 2.16e-3; 0.14; 0.001; 0.26; 0.3]
u0 = [1.1702; 6.3063; 0.3949]

u0LF = [0.0468627 ; 0.8987; 0.846515]
pLF = [666.7; 0.3478; 4.0e06; 2.0; 3000; 55.2; 7000; 0.09; 0.23; 3.9e-4; 0.1; 0.000833 ; 0.26; 0.25]


using DifferentialEquations
tspan = (0.0, 1000.0)
probHF = ODEProblem(ODEBZ.FBZ, u0, tspan, p)
probLF = ODEProblem(ODEBZ.FBZ, u0LF, tspan, pLF)

#solHF = solve(probHF, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10)
#solLF = solve(probHF, Rosenbrock23(); abstol = 1e-10, reltol = 1e-10)

end # module
