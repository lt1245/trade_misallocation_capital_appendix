# =============================================================================
# initialization.jl
#
# Activates the environment shipped with this folder (Project.toml plus the
# manifest matching the running Julia) and downloads every recorded
# dependency, so that a fresh Julia installation can run the replication
# package without any manual setup.
#
# Included as the first executable statement of the two entry points:
#   figures_tables.jl  (front-loads its own plotting imports)
#   code_new.jl        (in turn covers main_script.jl and search_developed.jl,
#                       both of which open by including code_new.jl)
#
# Three pinned manifests ship with the package:
#   Manifest.toml        Julia 1.12.5 -- the environment behind the published
#                        results, and the reference for reproducing them
#   Manifest-v1.11.toml  Julia 1.11
#   Manifest-v1.10.toml  Julia 1.10 (LTS)
#
# Julia loads Manifest-vX.Y.toml in preference to Manifest.toml, so the right
# one is selected automatically by the running version and none of the three is
# ever rewritten. The two older manifests exist so the package installs and
# runs on the Julia versions a replicator is most likely to already have.
#
# They pin the same version of every registered package as Manifest.toml does.
# What differs is only what ships inside Julia itself: stdlibs such as
# LinearAlgebra, SparseArrays and Statistics, and the bundled binaries, most
# consequentially OpenBLAS -- 0.3.29 on 1.12.5, 0.3.27 on 1.11, 0.3.23 on
# 1.10. Nothing numerical is pinned to a different release, but a different
# BLAS can still move the last digits, so output from 1.10 or 1.11 should be
# checked against the reference tables. The only packages the older trees drop
# are JuliaSyntaxHighlighting (1.12-only) and, on 1.10, StyledStrings, both
# terminal text-styling stdlibs with no bearing on results.
#
# Verified on 1.11.9 and 1.10.12: all 22 project dependencies resolve,
# precompile (300 and 299 packages) and load.
#
# Because the environment is anchored to this file's own directory via
# @__DIR__, the package no longer depends on the working directory Julia was
# launched from.
#
# Idempotent: repeat includes are no-ops, so the cost is paid once per session.
# =============================================================================

"""
    _replication_is_download_failure(err) -> Bool

True when `err` was caused by not being able to fetch bytes over the network (a
blocked package server, a blocked artifact host, an offline machine) rather
than by an unsatisfiable set of package versions.

The two need opposite responses: re-resolving the environment is the remedy for
a version conflict, but it cannot conjure up a binary that could not be
downloaded, and attempting it there buries the real cause under a second,
misleading error.

The match runs against the printed exception, which includes nested causes,
because `Pkg.instantiate()` downloads in parallel and reports failures wrapped
in task and composite exceptions.
"""
function _replication_is_download_failure(err)
    hints = (
        "unable to automatically download/install artifact",
        "unable to automatically install",
        "requesterror",
        "connection was reset",
        "connection reset",
        "connection refused",
        "recv failure",
        "could not download",
        "could not resolve host",
        "operation timed out",
        "timeout was reached",
        "failed to clone",
        "certificate",
    )
    io = IOBuffer()
    try
        showerror(io, err)
    catch
        print(io, err)
    end
    text = lowercase(String(take!(io)))
    return any(hint -> occursin(hint, text), hints)
end

if !isdefined(Main, :REPLICATION_ENV_READY)
    using Pkg

    # The version guard lives here rather than in a Project.toml [compat]
    # section: adding [compat] changes the recorded project_hash, which makes
    # Pkg advise running Pkg.resolve() on every run, and following that advice
    # would destroy the pinned manifest.
    shipped_for_version = joinpath(@__DIR__, "Manifest-v$(VERSION.major).$(VERSION.minor).toml")

    if VERSION >= v"1.12"
        # Manifest.toml applies as recorded; nothing to warn about.
    elseif isfile(shipped_for_version)
        @warn """
        This replication package was built and tested on Julia 1.12.5, but you
        are running $VERSION. The environment pinned for this version,
        $(basename(shipped_for_version)), is selected automatically, so the
        package installs and runs as it stands.

        It is not quite the environment behind the published results. The
        registered packages pin to the same versions, but Julia ships its own
        stdlibs and bundled binaries, so LinearAlgebra, SparseArrays and above
        all OpenBLAS differ from those used for the published run. That can
        move the last digits of anything numerical. Check the output against
        the reference tables, and use Julia 1.12.x to reproduce the published
        environment exactly.
        """
    else
        @warn """
        This replication package was built and tested on Julia 1.12.5, and
        ships pinned environments for Julia 1.10 and 1.11 as well. You are
        running $VERSION, for which nothing is pinned, so Pkg resolves the
        dependencies from scratch and the versions it picks may not match the
        published results. Installing Julia 1.12.x is strongly recommended.

        Do NOT follow any Pkg advice to run Pkg.resolve() against
        Manifest.toml: it records stdlibs that exist only in Julia 1.12
        (JuliaSyntaxHighlighting v1.12.0), so resolving it under an older
        Julia fails with "Unsatisfiable requirements detected".
        """
    end

    Pkg.activate(@__DIR__)

    try
        Pkg.instantiate()
    catch err
        if _replication_is_download_failure(err)
            @error """
            Could not download the packages or binary artifacts this package
            depends on. That is a connectivity problem on this machine, not a
            problem with the pinned environment, so the package versions are
            left untouched and no fallback is attempted.

            Julia needs plain HTTPS (not just git) access to:
              pkg.julialang.org           package tarballs and binary artifacts
              *.githubusercontent.com     GitHub release assets, the artifact fallback
              julialang-s3.julialang.org  Julia and artifact downloads
              github.com                  package sources

            A firewall or TLS-intercepting proxy usually shows up as
            "Connection was reset". Being able to git clone from github.com is
            not enough: Julia has a git fallback for package source but none
            for binary artifacts (GR, Qt6, FFMPEG, JpegTurbo, ...), which can
            only arrive over HTTPS from the hosts above.

            If those hosts cannot be unblocked, copy the Julia depot from a
            machine where this package already runs -- the artifacts, packages,
            compiled and registries subfolders of its .julia home folder --
            and rerun. Setting JULIA_PKG_OFFLINE=true afterwards stops Pkg
            from reaching out at all.
            """
            rethrow()
        end

        # Whichever manifest the running Julia actually selected: the
        # version-specific one where it is shipped, Manifest.toml otherwise.
        # The fallback has to act on that file rather than on Manifest.toml
        # unconditionally, or it would discard the wrong environment. Resolved
        # here rather than before the try so that the successful path does no
        # extra work and gains no extra way to fail.
        manifest = Pkg.Types.Context().env.manifest_file
        pinned_copy = first(splitext(manifest)) * ".pinned.toml"

        @warn """
        Could not install $(basename(manifest)) under Julia $VERSION. Falling
        back to resolving Project.toml from scratch, which replaces that
        manifest. Resolved package versions may differ from those used for the
        published results, so check the output against the reference tables.
        The pinned manifest is preserved next to it as $(basename(pinned_copy));
        restore that file to get the shipped environment back.
        """ exception = (err, catch_backtrace())

        isfile(manifest) && cp(manifest, pinned_copy; force = true)
        # Emptied rather than deleted: an empty Manifest-vX.Y.toml still takes
        # precedence over Manifest.toml, so deleting it would silently hand the
        # run back to the 1.12 manifest that just failed.
        write(manifest, "")
        try
            Pkg.resolve()
            Pkg.instantiate()
        catch
            # Leave the shipped environment as we found it.
            isfile(pinned_copy) && cp(pinned_copy, manifest; force = true)
            rm(pinned_copy; force = true)
            rethrow()
        end
    end

    @eval Main REPLICATION_ENV_READY = true
end
