"""
Generic, domain-agnostic helpers shared across the package.
"""
module Utils

export @timed_step

"""
    @timed_step timings key msg begin ... end

Log `msg`, run the block, record its duration in `timings[key]`, and fold the
duration into the log line: on a terminal the in-progress line is replaced with
`msg (1.234s)`; otherwise the timing is appended on a new line. Like `@elapsed`,
the block is spliced into the surrounding scope (no new scope), so variables it
defines remain visible to later steps.
"""
macro timed_step(timings, key, msg, ex)
    return quote
        local _msg = $(esc(msg))
        @info _msg
        local _t0 = time_ns()
        $(esc(ex))
        local _dt = (time_ns() - _t0) / 1e9
        $(esc(timings))[$(esc(key))] = _dt
        if stderr isa Base.TTY
            print(stderr, "\e[1A\e[2K")          # move up one line and clear it
            @info "$_msg ($(round(_dt; digits=3))s)"
        else
            @info "  ($(round(_dt; digits=3))s)"
        end
        _dt
    end
end

end # module Utils
