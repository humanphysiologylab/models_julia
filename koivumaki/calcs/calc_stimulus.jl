function calc_stimulus!(du, u, p, t, a)

    @unpack STIM_PERIOD, STIM_LEVEL, STIM_DURATION, amplitude = p

    i_stim = (t % STIM_PERIOD ≤ STIM_DURATION) ? STIM_LEVEL * amplitude : 0
    
    @pack! a = i_stim
    nothing
end
