using DifferentialEquations


function reset_dt!(integrator)
    set_proposed_dt!(integrator, 1e-6)
end


function stimulate!(integrator)
    @unpack STIM_PERIOD, STIM_DURATION = integrator.p
    t = integrator.t
    if (t % STIM_PERIOD ≤ STIM_DURATION)
        integrator.p.STIM_LEVEL = 1.
    else
        integrator.p.STIM_LEVEL = 0.
    end
    # @show t, integrator.p.STIM_LEVEL
end


function positive_mod(x, y)
    ((x % y) + y) % y
end


function is_shift(t, shift, period)
    half_period = period / 2
    x = t - shift + half_period
    positive_mod(x, period) - half_period
end


function is_before_stimulation(u, t, integrator; margin=1e-6)
    period = integrator.p.STIM_PERIOD
    shift = period - margin
    condition = is_shift(t, shift, period)
    # if condition ≈ 0.
    #     @show "before", t
    # end
    condition
end


function is_after_stimulation(u, t, integrator; margin=1e-6)
    period = integrator.p.STIM_PERIOD
    shift = integrator.p.STIM_DURATION + margin
    is_shift(t, shift, period)
    condition = is_shift(t, shift, period)
    # if condition ≈ 0.
    #     @show "after", t
    # end
    condition
end


function change_dtmax!(integrator; dtmax=1e-5)
    integrator.opts.dtmax = dtmax
    # @show integrator.t, get_proposed_dt(integrator)
end


function do_nothing(x...)
    nothing
end


function create_cb_set()

    margin = 1e-6

    cb_before = ContinuousCallback((u, t, integrator) -> is_before_stimulation(u, t, integrator; margin=margin),
                                   i -> change_dtmax!(i; dtmax=5e-5),
                                   do_nothing;
                                   save_positions=(false, false))

    cb_after = ContinuousCallback((u, t, integrator) -> is_after_stimulation(u, t, integrator; margin=margin),
                                  i -> change_dtmax!(i; dtmax=1.),
                                  do_nothing;
                                  save_positions=(false, false))

    cb_stim = DiscreteCallback((u, t, integrator) -> true,
                               stimulate!;
                               save_positions=(false, false))

    cb_set = CallbackSet(cb_before, cb_after, cb_stim)
    # cb_set = CallbackSet(cb_before, cb_after)

    return cb_set
end
