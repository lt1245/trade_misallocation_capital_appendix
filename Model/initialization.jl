# =============================================================================
# initialization.jl
#
# Activates the environment shipped with this folder (Project.toml /
# Manifest.toml) and downloads every recorded dependency, so that a fresh Julia
# installation can run the replication package without any manual setup.
#
# Included as the first executable statement of the two entry points:
#   figures_tables.jl  (front-loads its own plotting imports)
#   code_new.jl        (in turn covers main_script.jl and search_developed.jl,
#                       both of which open by including code_new.jl)
#
# Because the environment is anchored to this file's own directory via
# @__DIR__, the package no longer depends on the working directory Julia was
# launched from.
#
# Idempotent: repeat includes are no-ops, so the cost is paid once per session.
# =============================================================================

if !isdefined(Main, :REPLICATION_ENV_READY)
    using Pkg

    # Manifest.toml pins the exact versions used for the published results and
    # records julia_version = "1.12.5". The version guard lives here rather
    # than in a Project.toml [compat] section: adding [compat] changes the
    # recorded project_hash, which makes Pkg advise running Pkg.resolve() on
    # every run, and following that advice would destroy the pinned manifest.
    if VERSION < v"1.12"
        @warn """
        This replication package was built and tested on Julia 1.12.5, but you
        are running $VERSION. Instantiating the pinned Manifest.toml may fail;
        if it does, the fallback below will resolve different package versions
        and the results may not match the published tables. Installing Julia
        1.12.x is strongly recommended.
        """
    end

    Pkg.activate(@__DIR__)

    try
        Pkg.instantiate()
    catch err
        @warn """
        Could not reproduce the recorded Manifest.toml. Falling back to
        re-resolving Project.toml. Note this rewrites Manifest.toml in place,
        and resolved package versions may differ from those used for the
        published results, so check the output against the reference tables.
        """ exception = err
        Pkg.resolve()
        Pkg.instantiate()
    end

    @eval Main REPLICATION_ENV_READY = true
end
