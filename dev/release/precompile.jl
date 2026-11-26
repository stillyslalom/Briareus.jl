# The precompile run of the release build (build.jl): what it runs is compiled into the
# system image, so a first call of each command needs no compilation. It runs every
# headless command of `briareus` on a starter design and on a copy with one beam moved, in
# a temporary folder, and opens and closes the window when BRIAREUS_RELEASE_WINDOW is not
# "false".

using BriareusGUI, BriareusGUI.Controllers, BriareusGUI.Briareus
const A = BriareusGUI

mktempdir() do dir
    sess = DesignSession(resolution = :coarse)
    a = joinpath(dir, "a.briar")
    save_session(a, sess)
    bs = sess.beams.beams[]
    set_beams!(sess.beams, [k == 1 ? remake(b; power = 0.9 * b.power) : b for (k, b) in enumerate(bs)])
    b = joinpath(dir, "b.briar")
    save_session(b, sess)
    out(name) = joinpath(dir, name)
    redirect_stdout(devnull) do
        A._main(["--help"])
        A._main(["talis", a, "--out", out("talis.csv"), "--json"])
        A._main(["review-1w", a, "--out", out("review"), "--json"])
        A._main(["roe", a, "--out", out("roe"), "--json"])
        A._main(["diff", a, b, "--out", out("diff"), "--json"])
        script = out("script.jl")
        write(script, "sess = DesignSession(resolution = :coarse); save_session(ARGS[1], sess)\n")
        A._main(["run", script, out("c.briar")])
    end
end

if get(ENV, "BRIAREUS_RELEASE_WINDOW", "true") != "false"
    # open the window on a starter design, let it draw a few frames, and close it
    t0 = time()
    A._TICK_HOOK[] = sh -> time() - t0 > 5 && A.request_close()
    try
        briareus(DesignSession(resolution = :coarse))
    finally
        A._TICK_HOOK[] = nothing
    end
end
