using ITensors, TimeEvoMPS
using LinearAlgebra
using DelimitedFiles
using ForwardDiff

include("lib.jl")

ITensors.space(::SiteType"S=1/2") = 2

let
    N = 4
    J = 1.0
    D = 0
    E = 0
    p = Int(N/4 + 1)
    L = N/2
    maxdim = 256
    h0 = 0
    #theta_list = [pi/12, (5/12)*pi, (10/9)*pi]
    theta_list = [pi]
    nsweeps = Int(floor(2.0 * N))
    ek = nsweeps
    sites = siteinds("S=1/2", N; conserve_qns=true)
    sweeps = Sweeps(nsweeps)
    maxdim!(sweeps, 50, 100, 200, 400, 600, 800)
    cutoff!(sweeps, 1E-12)
    
    state = fill("Up", N)  # Iniciar com todos os spins para cima
    psi0 = randomMPS(sites, state, linkdims = ek)
    normalize!(psi0)
  
    for theta in theta_list

        J_1 = J * cos(theta)

        J_2 = J * sin(theta)

        print(J_1, J_2)

        println("DMRG_BASE")

        println("----")
        
        H = Set_Hamiltonian(N, J_1, J_2, D, E, h0, sites)
        
        energy_base, psi_base = dmrg(H, psi0; nsweeps, eigsolve_krylovdim=3*ek)

        entropy = VN_entropy(psi_base,p)
        print("")
        print("Energy:", energy_base,"|" ,"Entropy:", entropy)
        print("")
        println("----")
    end
end