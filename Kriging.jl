module Kriging

#
#  The Kriging module provides two types of service:
#
#  * universal kriging with isotropic covariance or correlation
#  * trend surface model kriging
#
#

using Stations
import Stations.nearest_grid_point, Stations.obs_variance, Stations.obs_value

using Storage
import Storage.spush

    
function numerical_solve_bisect(e2, eps2, k)

    N = size(e2,1)
    tgt = N - k
    s2_eta_left = 0.0
    s2_eta_right = 0.1

    val_left = sum(e2 ./ eps2)
    val_right = sum(e2 ./ (eps2 + s2_eta_right))
    if val_left < tgt
        return -1.0
    end

    while val_right > tgt
        s2_eta_right *= 2.0
        val_right = sum(e2 ./ (eps2 + s2_eta_right))
    end

    # bisection implementation (initialized with s2_eta = 0)
    while val_left - val_right > 1e-6
        
        # compute new value at center of eta interval
        s2_eta = 0.5 * (s2_eta_left + s2_eta_right)
        val = sum(e2 ./ (eps2 + s2_eta))

        if val > tgt
            val_left, s2_eta_left = val, s2_eta
        else
            val_right, s2_eta_right = val, s2_eta
        end

    end

    return 0.5 * (s2_eta_left + s2_eta_right)

end


function trend_surface_model_kriging(obs_data, X, K, V)
    """
    Trend surface model kriging, which assumes spatially uncorrelated errors.

    WARNING: The variable X is clobbered.

    The kriging results in the matrix K, which contains the kriged observations
    and the matrix V, which contains the kriging variance.
    """
    Nobs = length(obs_data)
    Ncov_all = size(X,3)

    dsize = size(X)[1:2]
    y = zeros((Nobs,1))
    m_var = zeros(Nobs)

    Xobs = zeros(Nobs, Ncov_all)
    for (obs,i) in zip(obs_data, 1:Nobs)
    	p = nearest_grid_point(obs)
        y[i] = obs_value(obs)
        Xobs[i,:] = X[p[1], p[2], :]
        m_var[i] = obs_variance(obs)
    end

    # if we have covariates full of zeros (e.g. there is no rain in the entire domain at
    # the current simulation time), we must exclude them or a singular exception
    # will be thrown by '\'
    cov_ids = find(map(i -> sum(Xobs[:,i].^2) > 0, 1:Ncov_all))
    Ncov = length(cov_ids)
    X = X[:,:,cov_ids]
    Xobs = Xobs[:, cov_ids]

    # if there are less observations than covariates, remove
    # Ncov - Nobs covariates from end
    if Ncov > Nobs
        Ncov = Nobs
        cov_ids = cov_ids[1:Ncov]
        X = X[:,:,1:Ncov]
        Xobs = Xobs[:, 1:Ncov]
    end

    # quick pre-conditioning hack
    # rescale all columns of Xobs to have norm of first column
    norm_1 = sum(Xobs[:,1].^2)^0.5
    for i in 2:Ncov
        norm_i = sum(Xobs[:,i].^2)^0.5
        if norm_i > 0.0
            Xobs[:,i] *= norm_1 / norm_i
            X[:,:,i] *= norm_1 / norm_i
        end
    end

    # initialize iterative algorithm
    s2_eta_hat_old = 10.0
    s2_eta_hat = 0.0
    XtSX = nothing
    beta = nothing

    iters = 0
    subzeros = 0

    # the while loop contains protection against division by zero in case s2_eta_hat_old
    # is below zero (happens if the estimate of s2_eta_hat fails and this gets copied into s2_eta_hat_old)
    while abs( (s2_eta_hat_old - s2_eta_hat) / max(s2_eta_hat_old, 1e-8)) > 1e-2
    
        # shift current estimate to old var
        s2_eta_hat_old = s2_eta_hat

        # recompute the covariance matrix
        Sigma = diagm(m_var) + s2_eta_hat * eye(Nobs)
        XtSX = Xobs' * (Sigma \ Xobs)

        # QR solution method of least squares
        Sigma_1_2 = diagm(diag(Sigma).^-0.5)   # Sigma^(-1/2)
        yt = Sigma_1_2 * y
        Q, R = qr(Sigma_1_2 * Xobs)
        beta = R \ (Q' * yt)
        res = y - Xobs * beta

        # compute new estimate of variance of microscale variability
        s2_array = res.^2 - m_var
        for j in 1:Nobs
            s2_array[j] += dot(vec(Xobs[j,:]), vec(XtSX \ Xobs[j,:]'))
        end
        s2_eta_hat2 = sum(s2_array) / Nobs

        # solve equation by bisection
        s2_eta_hat = numerical_solve_bisect(res.^2, m_var, Ncov)
        if s2_eta_hat < 0.0
            s2_eta_hat = 0.0
            println("TSM: s2_eta_hat estimate below zero, culling to zero.")
        end

        subzeros = sum(s2_array .< 0)
        iters += 1
    end

    # compute the OLS fit of the covariates to the observations
    spush("kriging_xtx_cond", cond(XtSX))
    spush("kriging_errors", (y - Xobs * beta)')

    # printing construction that makes sure order of printed betas does not vary
    # across times even if there are zero covariates
    beta_ext = ones((Ncov_all,1)) * NaN
    beta_ext[cov_ids] = beta
    spush("kriging_beta", beta_ext)

    spush("kriging_sigma2_eta", s2_eta_hat)
    spush("kriging_iters", iters)
    spush("kriging_subzero_s2_estimates", subzeros)

    # compute kriging field and kriging variance 
    for i in 1:dsize[1]
        for j in 1:dsize[2]
            x_ij = squeeze(X[i,j,:], 1)'   # convert covariates at position i,j into a column vector
            K[i,j] = dot(vec(x_ij), vec(beta))
            V[i,j] = s2_eta_hat + dot(vec(x_ij), vec(XtSX \ x_ij))
        end
    end

end


end