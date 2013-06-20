module FM

#
# Single cell dynamics of the fuel moisture model including the Kalman predict/update mechanisms.
# 
#

type FMModel
    k::Int
    m_ext::Vector{Float64}
    P::Matrix{Float64}
    Tk::Vector{Float64}
    Trk::Float64
    r0::Float64
    rk::Float64
    S::Float64
    latlon::Tuple

    function FMModel(latlon, k, m0, P0, Tk)
        # construct initial state
        m_ext = zeros(2*k+3)
        m_ext[1:k] = m0

        # initialize with defaults
        new(k, m_ext, copy(P0), Tk * 3600, 14.0 * 3600, 5.0, 8.0, 2.5, latlon)
    end

end


function advance_model(c::FMModel, Ed::Float64, Ew::Float64,
                       r::Float64, dt::Float64, Q::Matrix{Float64})

    # acquire local variables
    Tk = c.Tk
    k = c.k
    m_ext = c.m_ext
    m = m_ext[1:k]
    dTk = m_ext[k+1:2*k]
    dE, dS, dTrk = m_ext[2*k+1:2*k+3]

    # add assimilation deltas to moisture equilibria
    Ed = max(0, Ed + dE)
    Ew = max(0, Ew + dE)

    # compute equilibria
    equi = copy(m)
    rlag = zeros(k)
    model_ids = zeros(Int32, k)
    
    # switch between models, first branch is rain model, rest are for no-rain conditions
    if r > c.r0
        equi[:] = c.S + dS
        model_ids[:] = 3
        rlag[:] = 1.0 / (c.Trk + dTrk) * (1.0 - exp(- (r - c.r0) / c.rk))
    else
        model_ids[:] = 4
        for i in [1:k]
            # drying model
            if equi[i] > Ed
                model_ids[i] = 1
                equi[i] = Ed
            # wetting model
            elseif equi[i] < Ew
                model_ids[i] = 2
                equi[i] = Ew
            end
        end
        rlag[:] = 1.0 ./ (Tk + dTk)
    end

    # select integration method according to actual change
    # in each fuel
    change = dt * rlag
    m_new = zeros(k)
    for i in 1:k
        if change[i] < 0.01
            m_new[i] = m[i] + (equi[i] - m[i]) * (1.0 - exp(-change[i]))
        else
            m_new[i] = m[i] + (equi[i] - m[i]) * change[i] * (1.0 - 0.5 * change[i])
        end
    end

    # compute Jacobian and update local covariance matrix
    J = zeros((2*k+3,2*k+3))
    for i in 1:k

        if change[i] < 0.01
            J[i,i] = model_ids[i] < 4 ? exp(-change[i]) : 1.0
            dmi_dchng = (equi[i] - m[i]) * exp(-change[i])
            dmi_dequi = (1.0 - exp(-change[i]))
        else
            J[i,i] = model_ids[i] < 4 ? 1.0 - change[i] * (1.0 - 0.5 * change[i]) : 1.0
            dmi_dchng = (equi[i] - m[i]) * (1.0 - change[i])
            dmi_dequi = change[i] * (1.0 - 0.5 * change[i])
        end

        if r < c.r0 # drying/wetting model or dead zone
            if model_ids[i] < 4 # if not dead zone
                J[i,k+i] = dmi_dchng * (-dt) * (Tk[i] + dTk[i])^(-2)
                J[i,2*k+1] = dmi_dequi
            end
        else # rain model
            J[i,2*k+2] = dmi_dequi
            J[i,2*k+3] = dmi_dchng * dt * (exp(-(r - c.r0) / c.rk) - 1.0) * (c.Trk + dTrk)^(-2)
        end

        J[k+i,k+i] = 1.0

    end
    J[2*k+1,2*k+1] = 1.0
    J[2*k+2,2*k+2] = 1.0
    J[2*k+3,2*k+3] = 1.0

    # update covariance (using computed Jacobian) and state
    c.P = J*c.P*J' + Q
    c.m_ext[1:k] = m_new
end


function kalman_update(c::FMModel, O::Vector{Float64},
                       V::Matrix{Float64}, fuel_types::Vector{Int})
    P = c.P
    Nf = size(fuel_types,1)
    H = zeros((Nf,2*c.k+3))

    # construct observation operator for the given fuel types
    ndx = 1
    for i in fuel_types
        H[ndx,i] = 1.0
        ndx += 1
    end
    
    # Kalman update matrices
    K = P * H' / (H * P * H' + V)

    # update the state
    c.m_ext = c.m_ext + K * (O - c.m_ext[fuel_types])
    c.P -= K * H * P

    # return the kalman update
    return K[fuel_types]

end


end