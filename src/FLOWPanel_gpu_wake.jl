#=##############################################################################
# DESCRIPTION
    Task 052 stage A: host/device seams for a CuArray-backed wake
    `FLOWVPM.ParticleField` inside `PanelParticleWake`.

    A device-resident particle field keeps the per-step N-body work (UJ radix
    FMM, SFS, euler convection, relaxation) on the GPU; everything that is
    host-only by construction stays on the host through the seams below:

      - particle maintenance (trim policies + `FLOWVPM.merge_particles!`):
        runs on a cached HOST MIRROR ParticleField — D2H of the live particle
        prefix, host maintenance, H2D of the surviving prefix, side-buffer
        (splitting state / filament edge graph) sync both ways. ~2 x 67 MB
        of transfer at production shape (~15-25 ms).
      - the per-step read-only host consumers (WakeHealthMonitor,
        WakeInventoryMonitor, wake VTK writer): read the same host mirror,
        refreshed on demand (D2H only).
      - `apply_freestream!`: whole-field broadcast instead of the scalar loop.
      - the SFS AfterUJ hook when the rectangular-direct seam handled the
        wake self-influence pass (radix already accumulated E_str; only the
        dynamic-procedure coefficient update remains — device broadcasts in
        FLOWVPM_subfilterscale_gpu.jl).

    Nothing here activates for a plain Matrix-backed field: every entry point
    guards on `pfield.particles isa Array`.

# AUTHORSHIP
  * Created by  : task 052 stage A (agent), Aug 2026
  * License     : GNU Public License
=###############################################################################

"Whether a FLOWVPM particle field is device (CuArray) backed."
_gpu_pfield_on_device(pfield::FLOWVPM.ParticleField) =
    !(pfield.particles isa Array)

################################################################################
# HOST MIRROR
################################################################################

# host-mirror ParticleFields keyed by device-field identity
const _GPU_PFIELD_MIRRORS = Dict{UInt,Any}()

"Drop cached host mirrors (e.g. when a wake is rebuilt)."
clear_gpu_pfield_mirrors!() = (empty!(_GPU_PFIELD_MIRRORS); nothing)

function _gpu_pfield_mirror(pfield::FLOWVPM.ParticleField)
    key = objectid(pfield)
    m = get(_GPU_PFIELD_MIRRORS, key, nothing)
    if m === nothing
        m = FLOWVPM.ParticleField(pfield.maxparticles, eltype(pfield.particles);
                kernel=pfield.kernel, transposed=pfield.transposed,
                formulation=pfield.formulation)
        _GPU_PFIELD_MIRRORS[key] = m
    end
    return m::FLOWVPM.ParticleField
end

# copy the host-side per-particle bookkeeping (splitting state + filament
# edge graph) between two fields; these live in plain host arrays on BOTH
# device- and host-backed fields, so plain copyto! suffices
function _gpu_copy_side_buffers!(dst::FLOWVPM.ParticleField,
                                 src::FLOWVPM.ParticleField)
    for (d, s) in ((dst.splitting_state, src.splitting_state),
                   (dst.filament_edge_graph, src.filament_edge_graph))
        for fname in fieldnames(typeof(s))
            a = getfield(s, fname)
            a isa AbstractArray && copyto!(getfield(d, fname), a)
        end
    end
    return nothing
end

"D2H: refresh the host mirror's live prefix from the device field."
function _gpu_sync_mirror_from_device!(mirror::FLOWVPM.ParticleField,
                                       pfield::FLOWVPM.ParticleField)
    np = pfield.np
    if np > 0
        # columns 1:np are a contiguous linear prefix; the 5-arg copyto! keeps
        # the fast CuArray↔Array method (a 2-arg SubArray copy falls back to
        # element-wise scalar indexing, disallowed on GPU arrays)
        copyto!(mirror.particles, 1, pfield.particles, 1,
                size(pfield.particles, 1) * np)
    end
    mirror.np = np
    mirror.t = pfield.t
    mirror.nt = pfield.nt
    _gpu_copy_side_buffers!(mirror, pfield)
    return mirror
end

"H2D: push the host mirror's live prefix back to the device field."
function _gpu_sync_device_from_mirror!(pfield::FLOWVPM.ParticleField,
                                       mirror::FLOWVPM.ParticleField)
    np = mirror.np
    if np > 0
        # same contiguous-prefix rationale as _gpu_sync_mirror_from_device!
        copyto!(pfield.particles, 1, mirror.particles, 1,
                size(mirror.particles, 1) * np)
    end
    pfield.np = np
    _gpu_copy_side_buffers!(pfield, mirror)
    return pfield
end

################################################################################
# HOST-ONLY CONSUMER SEAMS
################################################################################

"""
Device path of `apply_particle_maintenance!`: D2H -> host trim/merge -> H2D.
The trim policies and `FLOWVPM.merge_particles!` are host scalar-index code;
running them on the mirror and pushing the surviving prefix back is the
least-invasive correct seam (all per-particle state lives in the `particles`
matrix; side buffers are synced explicitly).
"""
function _apply_particle_maintenance_device!(pfield::FLOWVPM.ParticleField,
        maintenance, ctx)
    mirror = _gpu_pfield_mirror(pfield)
    _gpu_sync_mirror_from_device!(mirror, pfield)
    apply_particle_maintenance!(mirror, maintenance, ctx)
    _gpu_sync_device_from_mirror!(pfield, mirror)
    return nothing
end

"Broadcast (device-safe) freestream application to particle velocities."
function _gpu_apply_freestream_device!(pfield::FLOWVPM.ParticleField, uinf)
    np = pfield.np
    np == 0 && return nothing
    u = SVector{3,eltype(pfield.particles)}(uinf[1], uinf[2], uinf[3])
    view(pfield.particles, FLOWVPM.U_INDEX, 1:np) .+= u
    return nothing
end

"""
Host view of a particle field for read-only host consumers (wake monitors,
VTK writer). Returns the field itself when host-backed; otherwise refreshes
and returns the host mirror (D2H of the live prefix, no write-back).
"""
function _wake_monitor_host_pfield(pfield::FLOWVPM.ParticleField)
    _gpu_pfield_on_device(pfield) || return pfield
    mirror = _gpu_pfield_mirror(pfield)
    # 052c: nested timer isolates the D2H refresh cost wherever a host
    # consumer (monitor, VTK writer) triggers it.
    _step_timer_measure(:d2h_mirror_sync; nested=true) do
        _gpu_sync_mirror_from_device!(mirror, pfield)
    end
    return mirror
end

################################################################################
# SFS POST HOOK FOR SEAM-HANDLED WAKE SELF-INFLUENCE (outputs === nothing)
################################################################################

# When the rectangular seam fully handles a pass (influence! passes
# `outputs=nothing` to the post hook), the only self-paired SFS-enabled
# particle field it can have handled is a device-backed one (host SFS
# self-pairs make the seam fall back to fmm!). The radix lifecycle already
# ACCUMULATED E_str into SFS_INDEX during the seam's UJ call, mirroring the
# host Estr_fmm! convention, so only the dynamic-procedure AfterUJ remains.
function post_evaluate_influence!(pfield::FLOWVPM.ParticleField,
        source::FLOWVPM.ParticleField, backend::FastMultipoleBackend,
        outputs::Nothing; i_target::Int=1, i_source::Int=1)
    pfield === source || return nothing
    FLOWVPM.isSFSenabled(pfield.SFS) || return nothing
    _gpu_pfield_on_device(pfield) || error(
        "seam-handled SFS self-influence post hook reached with a host " *
        "particle field; the rectangular seam should have fallen back to fmm!")
    _step_timer_measure(:wake_sfs; nested=true) do
        pfield.SFS(pfield, FLOWVPM.AfterUJ())
    end
    return nothing
end

################################################################################
# TIMERS (task 052 stage-A instrumentation; env-gated, default off)
################################################################################

const _GPU_TIMERS = Ref{Union{Nothing,Bool}}(nothing)

"Whether FLOWPANEL_GPU_TIMERS is armed (read lazily, cached)."
function gpu_timers_enabled()
    if _GPU_TIMERS[] === nothing
        _GPU_TIMERS[] = lowercase(get(ENV, "FLOWPANEL_GPU_TIMERS", "0")) in
            ("1", "true", "on")
    end
    return _GPU_TIMERS[]::Bool
end

"Log an env-gated timing line (used by the influence seam and the solve)."
function _gpu_timer_log(label::AbstractString, seconds::Real)
    gpu_timers_enabled() || return nothing
    @info "gpu_timer $(label) $(round(seconds; sigdigits=5)) s"
    return nothing
end
