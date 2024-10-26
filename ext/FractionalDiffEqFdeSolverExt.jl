module FractionalDiffEqFdeSolverExt

using SciMLBase, BoundaryValueDiffEq, ODEInterface, RecursiveArrayTools, ConcreteStructs,
      Setfield
using PreallocationTools
import SciMLBase: AbstractDiffEqInterpolation, StandardBVProblem, __solve, _unwrap_val
import ODEInterface: OptionsODE, OPT_ATOL, OPT_RTOL, OPT_METHODCHOICE, OPT_DIAGNOSTICOUTPUT,
                     OPT_ERRORCONTROL, OPT_SINGULARTERM, OPT_MAXSTEPS, OPT_BVPCLASS,
                     OPT_SOLMETHOD, OPT_RHS_CALLMODE, OPT_COLLOCATIONPTS, OPT_ADDGRIDPOINTS,
                     OPT_MAXSUBINTERVALS, RHS_CALL_INSITU, evalSolution
import ODEInterface: Bvpm2, bvpm2_init, bvpm2_solve, bvpm2_destroy, bvpm2_get_x
import ODEInterface: bvpsol
import ODEInterface: colnew

import FastClosures: @closure
import ForwardDiff

#------
# BVPM2
#------
function SciMLBase.__solve(prob::BVProblem, alg::BVPM2; dt = 0.0, reltol = 1e-3, kwargs...)
    if !(prob.problem_type isa TwoPointBVProblem)
        throw(ArgumentError("`BVPM2` only supports `TwoPointBVProblem!`"))
    end

    t₀, t₁ = prob.tspan
    u0_ = __extract_u0(prob.u0, prob.p, t₀)
    u0_size = size(u0_)
    n = __initial_guess_length(prob.u0)

    n == -1 && dt ≤ 0 && throw(ArgumentError("`dt` must be positive."))

    mesh = __extract_mesh(prob.u0, t₀, t₁, ifelse(n == -1, dt, n - 1))
    n = length(mesh) - 1
    no_odes = length(u0_)

    if prob.f.bcresid_prototype !== nothing
        left_bc, right_bc = prob.f.bcresid_prototype.x
        left_bc_size, right_bc_size = size(left_bc), size(right_bc)
        no_left_bc = length(left_bc)
    else
        left_bc = prob.f.bc[1](u0_, prob.p) # Guaranteed to be out of place here
        no_left_bc = length(left_bc)
    end

    obj = Bvpm2()
    if prob.u0 isa Function
        guess_function = @closure (x, y) -> (y .= vec(__initial_guess(prob.u0, prob.p, x)))
        bvpm2_init(obj, no_odes, no_left_bc, mesh, guess_function,
            eltype(u0_)[], alg.max_num_subintervals)
    else
        u0 = __flatten_initial_guess(prob.u0)
        bvpm2_init(
            obj, no_odes, no_left_bc, mesh, u0, eltype(u0)[], alg.max_num_subintervals)
    end

    bvp2m_f = if isinplace(prob)
        @closure (t, u, du) -> prob.f(reshape(du, u0_size), reshape(u, u0_size), prob.p, t)
    else
        @closure (t, u, du) -> du .= vec(prob.f(reshape(u, u0_size), prob.p, t))
    end
    bvp2m_bc = if isinplace(prob)
        @closure (ya, yb, bca, bcb) -> begin
            prob.f.bc[1](reshape(bca, left_bc_size), reshape(ya, u0_size), prob.p)
            prob.f.bc[2](reshape(bcb, right_bc_size), reshape(yb, u0_size), prob.p)
            return nothing
        end
    else
        @closure (ya, yb, bca, bcb) -> begin
            bca .= vec(prob.f.bc[1](reshape(ya, u0_size), prob.p))
            bcb .= vec(prob.f.bc[2](reshape(yb, u0_size), prob.p))
            return nothing
        end
    end

    opt = OptionsODE(OPT_RTOL => reltol, OPT_METHODCHOICE => alg.method_choice,
        OPT_DIAGNOSTICOUTPUT => alg.diagnostic_output,
        OPT_SINGULARTERM => alg.singular_term, OPT_ERRORCONTROL => alg.error_control)

    sol, retcode, stats = bvpm2_solve(obj, bvp2m_f, bvp2m_bc, opt)
    retcode = retcode ≥ 0 ? ReturnCode.Success : ReturnCode.Failure
    destats = SciMLBase.DEStats(
        stats["no_rhs_calls"], 0, 0, 0, stats["no_jac_calls"], 0, 0, 0, 0, 0, 0, 0, 0)

    x_mesh = bvpm2_get_x(sol)
    evalsol = evalSolution(sol, x_mesh)
    ivpsol = SciMLBase.build_solution(prob, alg, x_mesh,
        map(x -> reshape(convert(Vector{eltype(evalsol)}, x), u0_size), eachcol(evalsol));
        retcode, stats = destats, original = (sol, retcode, stats))

    bvpm2_destroy(obj)
    bvpm2_destroy(sol)

    return ivpsol
end

#-------
# BVPSOL
#-------
function SciMLBase.__solve(prob::BVProblem, alg::BVPSOL; maxiters = 1000,
        reltol = 1e-3, dt = 0.0, verbose = true, kwargs...)
    if !(prob.problem_type isa TwoPointBVProblem)
        throw(ArgumentError("`BVPSOL` only supports `TwoPointBVProblem!`"))
    end
    if !__has_initial_guess(prob.u0)
        throw(ArgumentError("Initial Guess is required for `BVPSOL`"))
    end

    t₀, t₁ = prob.tspan
    u0_ = __extract_u0(prob.u0, prob.p, t₀)
    u0_size = size(u0_)
    n = __initial_guess_length(prob.u0)

    n == -1 && dt ≤ 0 && throw(ArgumentError("`dt` must be positive."))
    u0 = __flatten_initial_guess(prob.u0)
    mesh = __extract_mesh(prob.u0, t₀, t₁, ifelse(n == -1, dt, n - 1))
    if u0 === nothing
        # initial_guess function was provided
        u0 = mapreduce(@closure(t->vec(__initial_guess(prob.u0, prob.p, t))), hcat, mesh)
    end

    if prob.f.bcresid_prototype !== nothing
        left_bc, right_bc = prob.f.bcresid_prototype.x
        left_bc_size, right_bc_size = size(left_bc), size(right_bc)
        no_left_bc = length(left_bc)
    else
        left_bc = prob.f.bc[1](u0_, prob.p) # Guaranteed to be out of place here
        no_left_bc = length(left_bc)
    end

    opt = OptionsODE(
        OPT_RTOL => reltol, OPT_MAXSTEPS => maxiters, OPT_BVPCLASS => alg.bvpclass,
        OPT_SOLMETHOD => alg.sol_method, OPT_RHS_CALLMODE => RHS_CALL_INSITU)

    bvpsol_f = if isinplace(prob)
        @closure (t, u, du) -> prob.f(reshape(du, u0_size), reshape(u, u0_size), prob.p, t)
    else
        @closure (t, u, du) -> du .= vec(prob.f(reshape(u, u0_size), prob.p, t))
    end

    bvpsol_bc = if isinplace(prob)
        @closure (ya, yb, r) -> begin
            left_bc = reshape(@view(r[1:no_left_bc]), left_bc_size)
            right_bc = reshape(@view(r[(no_left_bc + 1):end]), right_bc_size)
            prob.f.bc[1](left_bc, reshape(ya, u0_size), prob.p)
            prob.f.bc[2](right_bc, reshape(yb, u0_size), prob.p)
            return nothing
        end
    else
        @closure (ya, yb, r) -> begin
            r[1:no_left_bc] .= vec(prob.f.bc[1](reshape(ya, u0_size), prob.p))
            r[(no_left_bc + 1):end] .= vec(prob.f.bc[2](reshape(yb, u0_size), prob.p))
            return nothing
        end
    end

    sol_t, sol_x, retcode, stats = bvpsol(bvpsol_f, bvpsol_bc, mesh, u0, alg.odesolver, opt)

    if verbose
        if retcode == -3
            @warn "Integrator failed to complete the trajectory"
        elseif retcode == -4
            @warn "Gauss Newton method failed to converge"
        elseif retcode == -5
            @warn "Given initial values inconsistent with separable linear bc"
        elseif retcode == -6
            @warn "Iterative refinement faild to converge for `sol_method=0` \
                   Termination since multiple shooting condition or \
                   condition of Jacobian is too bad for `sol_method=1`"
        elseif retcode == -8
            @warn "Condensing algorithm for linear block system fails, try `sol_method=1`"
        elseif retcode == -9
            @warn "Sparse linear solver failed"
        elseif retcode == -10
            @warn "Real or integer work-space exhausted"
        elseif retcode == -11
            @warn "Rank reduction failed - resulting rank is zero"
        end
    end

    ivpsol = SciMLBase.build_solution(prob, alg, sol_t,
        map(x -> reshape(convert(Vector{eltype(u0_)}, x), u0_size), eachcol(sol_x));
        retcode = retcode ≥ 0 ? ReturnCode.Success : ReturnCode.Failure,
        stats, original = (sol_t, sol_x, retcode, stats))

    return ivpsol
end

#-------
# COLNEW
#-------
function SciMLBase.__solve(prob::BVProblem, alg::COLNEW; maxiters = 1000,
        reltol = 1e-3, dt = 0.0, verbose = true, kwargs...)
    dt ≤ 0 && throw(ArgumentError("`dt` must be positive"))

    t₀, t₁ = prob.tspan
    u0_ = __extract_u0(prob.u0, prob.p, t₀)
    u0_size = size(u0_)
    n = __initial_guess_length(prob.u0)

    u0 = __flatten_initial_guess(prob.u0)
    mesh = __extract_mesh(prob.u0, t₀, t₁, ifelse(n == -1, dt, n - 1))
    if u0 === nothing
        # initial_guess function was provided
        u0 = mapreduce(@closure(t->vec(__initial_guess(prob.u0, prob.p, t))), hcat, mesh)
    end

    no_odes = length(u0_)

    # has_initial_guess = prob.u0 isa AbstractVector{<:AbstractArray}
    # dt ≤ 0 && throw(ArgumentError("dt must be positive"))
    # no_odes, n, u0 = if has_initial_guess
    #     length(first(prob.u0)), (length(prob.u0) - 1), reduce(hcat, prob.u0)
    # else
    #     length(prob.u0), Int(cld((prob.tspan[2] - prob.tspan[1]), dt)), prob.u0
    # end

    T = eltype(u0)
    # mesh = collect(range(prob.tspan[1], stop = prob.tspan[2], length = n + 1))
    orders = ones(Int, no_odes)
    _tspan = [prob.tspan[1], prob.tspan[2]]
    iip = SciMLBase.isinplace(prob)

    rhs = @closure (t, u, du) -> begin
        if iip
            prob.f(du, u, prob.p, t)
        else
            (du .= prob.f(u, prob.p, t))
        end
    end

    if prob.f.jac === nothing
        if iip
            jac = (df, u, p, t) -> begin
                _du = similar(u)
                prob.f(_du, u, p, t)
                _f = @closure (du, u) -> prob.f(du, u, p, t)
                ForwardDiff.jacobian!(df, _f, _du, u)
                return
            end
        else
            jac = (df, u, p, t) -> begin
                _du = prob.f(u, p, t)
                _f = @closure (du, u) -> (du .= prob.f(u, p, t))
                ForwardDiff.jacobian!(df, _f, _du, u)
                return
            end
        end
    else
        jac = prob.f.jac
    end
    Drhs = @closure (t, u, df) -> jac(df, u, prob.p, t)

    bcresid_prototype, _ = __get_bcresid_prototype(prob.problem_type, prob, u0)

    if prob.problem_type isa TwoPointBVProblem
        n_bc_a = length(first(bcresid_prototype))
        n_bc_b = length(last(bcresid_prototype))
        zeta = vcat(fill(first(prob.tspan), n_bc_a), fill(last(prob.tspan), n_bc_b))
        bc = @closure (i, z, resid) -> begin
            tmpa = copy(z)
            tmpb = copy(z)
            tmp_resid_a = zeros(T, n_bc_a)
            tmp_resid_b = zeros(T, n_bc_b)
            prob.f.bc[1](tmp_resid_a, tmpa, prob.p)
            prob.f.bc[2](tmp_resid_b, tmpb, prob.p)

            for j in 1:n_bc_a
                if i == j
                    resid[1] = tmp_resid_a[j]
                end
            end
            for j in 1:n_bc_b
                if i == (j + n_bc_a)
                    resid[1] = tmp_resid_b[j]
                end
            end
        end

        Dbc = @closure (i, z, dbc) -> begin
            for j in 1:n_bc_a
                if i == j
                    dbc[i] = 1.0
                end
            end
            for j in 1:n_bc_b
                if i == (j + n_bc_a)
                    dbc[i] = 1.0
                end
            end
        end
        fixed_points = nothing
    else
        zeta = sort(alg.zeta)
        bc = alg.bc_func
        Dbc = alg.dbc_func
        left_index = findlast(x -> x ≈ t₀, zeta) + 1
        right_index = findfirst(x -> x ≈ t₁, zeta) - 1
        fixed_points = alg.zeta[left_index:right_index]
    end

    if fixed_points === nothing
        opt = OptionsODE(
            OPT_BVPCLASS => alg.bvpclass, OPT_COLLOCATIONPTS => alg.collocationpts,
            OPT_MAXSTEPS => maxiters, OPT_DIAGNOSTICOUTPUT => alg.diagnostic_output,
            OPT_MAXSUBINTERVALS => alg.max_num_subintervals, OPT_RTOL => reltol)
    else
        opt = OptionsODE(
            OPT_BVPCLASS => alg.bvpclass, OPT_COLLOCATIONPTS => alg.collocationpts,
            OPT_MAXSTEPS => maxiters, OPT_DIAGNOSTICOUTPUT => alg.diagnostic_output,
            OPT_MAXSUBINTERVALS => alg.max_num_subintervals,
            OPT_RTOL => reltol, OPT_ADDGRIDPOINTS => fixed_points)
    end

    sol, retcode, stats = colnew(_tspan, orders, zeta, rhs, Drhs, bc, Dbc, nothing, opt)

    if verbose
        if retcode == 0
            @warn "Collocation matrix is singular"
        elseif retcode == -1
            @warn "The expected no. of subintervals exceeds storage(try to increase \
                   `OPT_MAXSUBINTERVALS`)"
        elseif retcode == -2
            @warn "The nonlinear iteration has not converged"
        elseif retcode == -3
            @warn "There is an input data error"
        end
    end

    evalsol = evalSolution(sol, mesh)
    destats = SciMLBase.DEStats(
        stats["no_rhs_calls"], 0, 0, 0, stats["no_jac_calls"], 0, 0, 0, 0, 0, 0, 0, 0)

    return DiffEqBase.build_solution(
        prob, alg, mesh, collect(Vector{eltype(evalsol)}, eachrow(evalsol));
        retcode = retcode > 0 ? ReturnCode.Success : ReturnCode.Failure,
        stats = destats, original = (sol, retcode, stats))
end

export BVPM2, BVPSOL, COLNEW

end
