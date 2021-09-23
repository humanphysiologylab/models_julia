function calc_stimulus!(du, u, p, t, a)

    @unpack STIM_PERIOD, STIM_LEVEL, STIM_DURATION, amplitude = p

    i_stim = (t % STIM_PERIOD ≤ STIM_DURATION) ? STIM_LEVEL * amplitude : 0
    @show i_stim
    
    @pack! a = i_stim
    nothing
end


function calc_stimulus_v2!(du, u, p, t, a)

    @unpack STIM_LEVEL, amplitude = p

    i_stim = STIM_LEVEL * amplitude
    
    @pack! a = i_stim
    nothing
end