function calc_serca!(du, u, p, t, a)

    @unpack Vnonjunct1, cpumps, k4, k3, J_SERCASR_multiplier, Vnonjunct2, Vnonjunct3, Vss, k1, k2, J_bulkSERCA_multiplier = p
    
    CaSR1, serca_a1, CaSR2, serca_a2, CaSR3, serca_a3, CaSR4, serca_ass, Cai1, Cai2, Cai3, Cass = u[1], u[38], u[2], u[39], u[3], u[40], u[4], u[41], u[5], u[6], u[7], u[9]
    
    J_SERCASR1 = J_SERCASR_multiplier * ((-k3) * ((CaSR1) ^ 2.0) * (cpumps - cpumps * serca_a1) + k4 * cpumps * serca_a1) * Vnonjunct1 * 2.0
    J_SERCASR2 = J_SERCASR_multiplier * ((-k3) * ((CaSR2) ^ 2.0) * (cpumps - cpumps * serca_a2) + k4 * cpumps * serca_a2) * Vnonjunct2 * 2.0
    J_SERCASR3 = J_SERCASR_multiplier * ((-k3) * ((CaSR3) ^ 2.0) * (cpumps - cpumps * serca_a3) + k4 * cpumps * serca_a3) * Vnonjunct3 * 2.0
    J_SERCASRss = J_SERCASR_multiplier * ((-k3) * ((CaSR4) ^ 2.0) * (cpumps - cpumps * serca_ass) + k4 * cpumps * serca_ass) * Vss * 2.0
    J_bulkSERCA1 = J_bulkSERCA_multiplier * (k1 * ((Cai1) ^ 2.0) * (cpumps - cpumps * serca_a1) - k2 * cpumps * serca_a1) * Vnonjunct1 * 2.0
    J_bulkSERCA2 = J_bulkSERCA_multiplier * (k1 * ((Cai2) ^ 2.0) * (cpumps - cpumps * serca_a2) - k2 * cpumps * serca_a2) * Vnonjunct2 * 2.0
    J_bulkSERCA3 = J_bulkSERCA_multiplier * (k1 * ((Cai3) ^ 2.0) * (cpumps - cpumps * serca_a3) - k2 * cpumps * serca_a3) * Vnonjunct3 * 2.0
    J_bulkSERCAss = J_bulkSERCA_multiplier * (k1 * ((Cass) ^ 2.0) * (cpumps - cpumps * serca_ass) - k2 * cpumps * serca_ass) * Vss * 2.0
    du[38] = d_serca_a1 = 0.5 * ((-J_SERCASR1) + J_bulkSERCA1) / Vnonjunct1 / cpumps
    du[39] = d_serca_a2 = 0.5 * ((-J_SERCASR2) + J_bulkSERCA2) / Vnonjunct2 / cpumps
    du[40] = d_serca_a3 = 0.5 * ((-J_SERCASR3) + J_bulkSERCA3) / Vnonjunct3 / cpumps
    du[41] = d_serca_ass = 0.5 * ((-J_SERCASRss) + J_bulkSERCAss) / Vss / cpumps
    
    @pack! a = J_SERCASR1, J_SERCASR2, J_SERCASR3, J_SERCASRss, J_bulkSERCA1, J_bulkSERCA2, J_bulkSERCA3, J_bulkSERCAss
    nothing
end
