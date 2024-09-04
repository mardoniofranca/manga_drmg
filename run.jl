using ITensors, TimeEvoMPS
using LinearAlgebra
using DelimitedFiles
using ForwardDiff

include("lib.jl")

ITensors.space(::SiteType"S=1/2") = 2

let
    N=8;  @show N
    J=1.0; D=0; E=0
    p=Int(N/2+1)
    @show p
    L=N/2 ;  maxdim = 256;  h0=0
    theta_g    = [15,75,200]
    theta_list = [pi/12, (5/12)*pi, (10/9)*pi] #theta_list = [0]
    nsweeps = 10 # Int(floor(1.5*N))
    ek = nsweeps
    sites = siteinds("S=1/2",N;conserve_qns=true,conserve_sz=false)
    sweeps = Sweeps(nsweeps) # number of sweeps is 5
    maxdim!(sweeps,50, 100, 200, 400)
    cutoff!(sweeps,1E-10) # desired truncation errorhtop
  
    state = [isodd(n) ? "Up" : "Dn"  for n in 1:N]
    psi1 = randomMPS(sites,state,linkdims = ek) #tentar 40

    state = [isodd(n) ? "Dn" : "Up"  for n in 1:N]
    psi2 = randomMPS(sites,state,linkdims = ek) #tentar 40
    psi0=+(psi1,psi2)
    println(psi0)
   
    normalize!(psi0)
   
    index = 0
    
    index = 1 # Inicialize o index antes do loop





    for i in 1:1:360
        println(i)
        theta = deg2rad(i)

        println("DMRG_BASE")
        J_1 = J * cos(theta)
        J_2 = J * sin(theta)
        println("J_1, J_2: ", J_1, ", ", J_2)
        @show h0
        H = Set_Hamiltonian(N, J_1, J_2, D, E, h0, sites)
        energy_base, psi_base = dmrg(H, psi0; nsweeps, eigsolve_krylovdim=3)
        entropy = VN_entropy(psi_base, p)
        println("")
        println("Theta: ", theta_g[index], " | Energy: ", energy_base, " | Entropy: ", entropy)
        println("")
        data = string(i, ";", theta_g[index], ";", energy_base, ";", entropy)
        println(data)
        open("data/resultados.dat", "a") do io
            println(io, data) # Corrigido para escrever a linha de dados
        end
        
    end

end