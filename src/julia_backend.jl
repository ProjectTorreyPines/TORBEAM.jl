# The pure-Julia backend: trace the beam, absorb, build the deposition profile
# and pack everything into TORBEAM's `BeamOutputs` layout.

"""
    julia_beam(inputs::BeamInputs, torbeam_params::TorbeamParams)

Run one launcher with the Julia backend and return `BeamOutputs` in the layout
of the Fortran library (cm, MW, kA; see `BeamOutputs`).
"""
function julia_beam(inputs::BeamInputs, torbeam_params::TorbeamParams)
    nprofv = torbeam_params.nprofv
    m = PlasmaModel(inputs)
    l = Launch(inputs)
    nmax = torbeam_params.npow == 0 ? 0 : torbeam_params.nmaxh
    b = trace_beam(m, l; rhostop=torbeam_params.rhostop, nmax, reltol=torbeam_params.xrtol, abstol=torbeam_params.xatol)

    # ray points every ~1 cm, as the Fortran stores them
    npts = max(2, min(ceil(Int, b.length / 0.01) + 1, NDAT))
    ss = range(0.0, b.length; length=npts)
    iend = length(ss)
    t1data = zeros(6 * iend)
    t1tdata = zeros(6 * iend)
    for (i, s) in enumerate(ss)
        bw = beam_widths(b, s)
        x = bw.x
        up = x + bw.wp * bw.ep
        lo = x - bw.wp * bw.ep
        le = x + bw.wh * bw.eh
        ri = x - bw.wh * bw.eh
        t1data[i] = 100 * hypot(x[1], x[2])
        t1data[iend+i] = 100 * x[3]
        t1data[2iend+i] = 100 * hypot(up[1], up[2])
        t1data[3iend+i] = 100 * up[3]
        t1data[4iend+i] = 100 * hypot(lo[1], lo[2])
        t1data[5iend+i] = 100 * lo[3]
        t1tdata[i] = 100 * x[1]
        t1tdata[iend+i] = 100 * x[2]
        t1tdata[2iend+i] = 100 * le[1]
        t1tdata[3iend+i] = 100 * le[2]
        t1tdata[4iend+i] = 100 * ri[1]
        t1tdata[5iend+i] = 100 * ri[2]
    end

    if torbeam_params.ncd == 1 && nmax > 0
        table = CurrentDriveTable(m, m.Zeff)
        efficiency = (u, s) -> cd_efficiency(table, state(m, u[1], u[2], u[3]), l.wave, u[4:6]; nmax)
    else
        efficiency = nothing
    end
    dep = deposition(b, m; efficiency)
    t2ndata = zeros(3 * NPNT)
    t2ndata[1:NPNT] = dep.ρ
    t2ndata[NPNT+1:2NPNT] = dep.dPdV ./ 1e6          # MW/m³
    # the Fortran reports the driven current in the direction of the plasma current:
    # our j∥ is along B, so multiply by sign(B0) sign(Ip), with sign(Ip) = sgnm
    # (COCOS 11: psi increases outwards for Ip > 0)
    sgn_j = sign(m.B0) * inputs.floatinbeam[34]
    t2ndata[2NPNT+1:3NPNT] = sgn_j .* dep.j ./ 1e6            # MA/m²
    Icd = sgn_j * sum(dep.Jbin) / (2π * m.R_axis) / 1e3      # kA: ∫ j dA ≈ ∫ j dV / (2π R)

    Pabs = (l.power - power(b, b.length)) / 1e6      # MW
    rhoresult = fill(-1.0, MAXRHR)
    if Pabs > 0
        k = argmax(dep.Pbin)
        rhoresult[1] = dep.ρ[k]
        # R, Z of the maximum absorption along the central ray
        Ps = [power(b, s) for s in ss]
        imax = argmax(-diff(Ps))
        x = b.sol(ss[imax])[1:3]
        rhoresult[2] = 100 * hypot(x[1], x[2])
        rhoresult[3] = 100 * x[3]
    end
    rhoresult[5] = 3 * nprofv
    rhoresult[13] = Icd
    rhoresult[14] = Pabs
    rhoresult[20] = b.exit == :absorbed || Pabs > 1e-6 * l.power / 1e6 ? 0.0 : (b.exit == :grid ? 1.0 : 2.0)

    # volume profile on nprofv points
    ρv = range(0.0, 1.0; length=nprofv)
    volprof = zeros(2 * MAXVOL)
    volprof[1:nprofv] = ρv
    volprof[nprofv+1:2nprofv] = cubic_resample(dep.ρ, dep.V, collect(ρv))
    t2data = zeros(3 * nprofv + 9)
    t2data[2nprofv+1:3nprofv] = volprof[nprofv+1:2nprofv]
    bw = beam_widths(b, b.length)
    vend = b.sol(b.length, Val{1})[1:3]
    vend ./= norm(vend)
    t2data[3nprofv+1:3nprofv+3] = vend
    t2data[3nprofv+4] = 100 * bw.wh
    t2data[3nprofv+5] = 100 * bw.wv

    return BeamOutputs(rhoresult, iend, t1data, t1tdata, 0, t2data, t2ndata, 0, 0, volprof)
end
