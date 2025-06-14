module ClimateModelnd2
# nondimensional version
export model!,jacobian!

function model!(du, u, p, t)
	
    T_A, T_W1, T_W2, T_L1, T_L2, mu2, nu2, C_f, C_s, C_p = u
	# unpack the parameter array
	parms, funcs = p
	C0, nu0, T_fs, p_W, p_L, sigm, beta, epsilon, xi_R, xi_A, q, c_1, c_2, c_L, L, 
	 H_W1, H_W2, H_L1, H_L2, B_c, theta_v, A1, A2, alp_Ww, alp_Wc, alp_Lw, alp_Lc, B_alpha, 
	 G_C, G_CH4, G_W1, G_W2, eta_Cl, delta, k_nu, d5f, d5s, d5p, Q10, B5, B10, 
	 T_Wb, k_W1, k_W2, k_Wb, k_L, F_A, F_W = parms
	# unpack the functions array
	mu1, nu1, Qfunc = funcs

	mu = mu1(t) + mu2
	nu = nu1(t) + nu2
	(Qval, QWval, QLval) = Qfunc(t)
    eta_val = eta(T_A,mu,nu, G_C, G_CH4, delta, G_W2, G_W1, eta_Cl)
    FcW= F_C(T_W1, A1, A2)
    FcL= F_C(T_L1, A1, A2)
    cW1= c_W(T_W1, T_fs, c_1, c_2, B_c, L)
    cW2= c_W(T_W2, T_fs, c_1, c_2, B_c, L)
    cWL1= c_W(T_L1, 1.0, c_1, c_2, B_c, L)
    cWL2= c_W(T_L2, 1.0, c_1, c_2, B_c, L)
	kW = k_W(T_W1,T_fs,k_W1,k_W2,B_c)
    alp_W = alpha_x(T_W1, T_fs, alp_Ww, alp_Wc, B_alpha)
    alp_L = alpha_x(T_L1, 1.0, alp_Lw, alp_Lc, B_alpha)
	Q10factor = decay_equation(T_L2,Q10,B5,B10)
    df= d5f * Q10factor
    ds= d5s * Q10factor
    dp= d5p * Q10factor

	du[1] = (xi_A * Qval + F_A + eta_val * sigm*(p_W * T_W1^4 + p_L * T_L1^4) +
             p_W * FcW + p_L * FcL - epsilon * sigm * T_A^4)
	du[2] = ((1 - alp_W) * (1 - xi_A - xi_R) * QWval +

			 beta*epsilon * sigm * T_A^4 - sigm * T_W1^4 - FcW - kW*(T_W1-T_W2))/(H_W1 * cW1)
	du[3] = (-kW*(T_W2-T_W1) - k_Wb*(T_W2-T_Wb) + F_W)/(H_W2 * cW2)
	du[4] = ((1 - alp_L) * (1 - xi_A - xi_R) * QLval + beta*epsilon * sigm * T_A^4 - 
			sigm * T_L1^4 - FcL - k_L*(T_L1-T_L2))/(H_L1 * (c_L + theta_v * cWL1))
	du[5] = (-k_L*(T_L2 - T_L1))/(H_L2 * (c_L + theta_v * cWL2))
	val = (df*C_f + ds*C_s + dp*C_p) 
	du[6] = val * q * C0 + k_nu * nu2
	du[7] = val * (1 - q) * C0 - k_nu * nu2
    du[8] = -df * C_f
    du[9] = -ds * C_s
    du[10] = -dp* C_p
    return du
end

function jacobian!(jac, u, p, t)

    T_A, T_W1, T_W2, T_L1, T_L2, mu2, nu2, C_f, C_s, C_p = u
	# unpack the parameter array
	parms, funcs = p
	C0, nu0, T_fs, p_W, p_L, sigm, beta, epsilon, xi_R, xi_A, q, c_1, c_2, c_L, L, 
	 H_W1, H_W2, H_L1, H_L2, B_c, theta_v, A1, A2, alp_Ww, alp_Wc, alp_Lw, alp_Lc, B_alpha, 
	 G_C, G_CH4, G_W1, G_W2, eta_Cl, delta, k_nu, d5f, d5s, d5p, Q10, B5, B10, 
	 T_Wb, k_W1, k_W2, k_Wb, k_L, F_A, F_W = parms
	# unpack the functions array
	mu1, nu1, Qfunc = funcs

    mu = mu1(t) + mu2
	nu = nu1(t) + nu2
	(Qval, QWval, QLval) = Qfunc(t)

	# First reset all entries to zero
	for j in 1:10, i in 1:10
		jac[i,j] = 0.0
	end
	# compute value of eta
	etaval = eta(T_A, mu, nu, G_C, G_CH4, delta, G_W2, G_W1, eta_Cl)
	# compute the value of (1-ξ_A-ξ_R)*QL(t) and (1-ξ_A-ξ_R)*QW(t)
	xiQWval = (1.0 - xi_A - xi_R)*QWval
	xiQLval = (1.0 - xi_A - xi_R)*QLval
	# define a function for the derivative of F_C
	FCdiff(T) = A1 + A1^2*(T-1.0)/sqrt(A1^2*(T-1.0)^2+A2^2)
	# define a function for the derivative of alpha_x
	alphadiff(T,S,alp_w,alp_c) = 0.5*(alp_w - alp_c)*(sech((T-S)/B_alpha)^2)/B_alpha
	# define a function for the derivative of c_W
	function c_Wdiff(T,S)
		if T-S <= -B_c
			return 0.0
		elseif T - S >= B_c
			return 0.0
		else
			local dw = π/(2*B_c)
			local w = dw*(T - S)
			local x = (c_2 - c_1)*cos(w)
			return dw*(x - ((c_2 - c_1)*(T - S) + L)*sin(w)*dw*0.5)
		end
	end
	function k_Wdiff(T,S)
		if T-S <= -B_c
			return 0.0
		elseif T - S >= B_c
			return 0.0
		else
			local dw = π/(2*B_c)
			return 0.5*(k_W2 - k_W1)*cos(dw*(T-S))*dw
		end
	end
	function decaydiff(T,Q10,B5,B10,decayval)
		if T <= 1.0
			return 0.0
		elseif T < 1.0 + B5
			return decayval*log(Q10)/B10 + (0.5*π/B5)*sin((T-1.0)*π/B5)*Q10^((T-1.0-B5)/B10)
		else
			return decayval*log(Q10)/B10
		end
	end

	# First row of jacobian  T_A
	w = sigm*(p_W*T_W1^4 + p_L*T_L1^4)
	# derivative of Clausius-Clapeyron part
	CCdiff = CC(T_A,G_W1)*(G_W1/T_A - 1.0)/T_A
	jac[1,1] = (etaval - 1.0)*(-delta* G_W2)*CCdiff*w - 4*epsilon*sigm*T_A^3
	jac[1,2] = etaval *sigm*p_W*4*T_W1^3 + p_W*FCdiff(T_W1)
	jac[1,4] = etaval *sigm*p_L*4*T_L1^3 + p_L*FCdiff(T_L1)
	jac[1,6] = (etaval - 1.0)*(-G_C)*w
	jac[1,7] = (etaval - 1.0)*(-G_CH4)*w

	# next values used in rows 2 and 3
	kW = k_W(T_W1,T_fs,k_W1,k_W2,B_c)
	kWdiff = k_Wdiff(T_W1,T_fs)

	# Second row of jacobian  T_W1
	v = c_W(T_W1, T_fs, c_1, c_2, B_c, L)
	w = 1.0/(H_W1*v)
	jac[2,1] = w*beta*epsilon*sigm*4*T_A^3
	jac[2,2] = w*(-alphadiff(T_W1,T_fs,alp_Ww,alp_Wc) * xiQWval - 4*sigm*T_W1^3 - FCdiff(T_W1) - 
				  kW - (T_W1-T_W2)*kWdiff
				  - (c_Wdiff(T_W1,T_fs)/v)*((1-alpha_x(T_W1,T_fs,alp_Ww,alp_Wc,B_alpha))*xiQWval + 
				beta*epsilon*sigm*T_A^4 - sigm*T_W1^4 - F_C(T_W1,A1,A2) -kW*(T_W1-T_W2)))
	jac[2,3] = kW*w

	# Third row of jacobian  T_W2
	v = c_W(T_W2, T_fs, c_1, c_2, B_c, L)
	w = 1.0/(H_W2*v)
	jac[3,2] = w*(kW - kWdiff*(T_W2-T_W1))
	jac[3,3] = w*((-kW - k_Wb ) - (c_Wdiff(T_W2,T_fs)/v)*(-kW*(T_W2-T_W1)-k_Wb*(T_W2-T_Wb)+F_W))

	# Fourth row of jacobian   T_L1
	v = c_L + theta_v*c_W(T_L1, 1.0, c_1, c_2, B_c, L)
	w = 1.0/(H_L1*v)
	jac[4,1] = w*beta*epsilon*sigm*4*T_A^3
	jac[4,4] = w*(-alphadiff(T_L1,1.0,alp_Lw,alp_Lc) * xiQLval - 4*sigm*T_L1^3 - FCdiff(T_L1) - k_L
				- (theta_v*c_Wdiff(T_L1,1.0)/v)*((1-alpha_x(T_L1,1.0,alp_Lw,alp_Lc,B_alpha))*xiQLval +
				beta*epsilon*sigm*T_A^4 - sigm*T_L1^4 - F_C(T_L1,A1,A2) - k_L*(T_L1-T_L2)))
	jac[4,5] = k_L*w

	# Fifth row of the jacobian  T_L2
	v = c_L + theta_v*c_W(T_L2, 1.0, c_1, c_2, B_c, L)
	w = 1.0/(H_L2*v)
	jac[5,4] = k_L*w
	jac[5,5] = w*(-k_L - (theta_v*c_Wdiff(T_L2,1.0)/v)*(-k_L*(T_L2 - T_L1) ))

	# next values used in all remaining rows
	decayval = decay_equation(T_L2,Q10,B5,B10)
	decaydiffval = decaydiff(T_L2,Q10,B5,B10,decayval)

	# Sixth row of jacobian  mu2
	v = q*C0
	w = v*decayval
	jac[6,5] = v*(C_f*d5f + C_s*d5s + C_p*d5p)*decaydiffval
	jac[6,7] = k_nu
	jac[6,8] = d5f*w
	jac[6,9] = d5s*w
	jac[6,10] = d5p*w

	# Seventh row of jacobian  nu2
	v = (1.0-q)*C0
	w = v*decayval
	jac[7,5] = v*(C_f*d5f + C_s*d5s + C_p*d5p)*decaydiffval
	jac[7,7] = -k_nu
	jac[7,8] = d5f*w
	jac[7,9] = d5s*w
	jac[7,10] = d5p*w

	# Eighth row of jacobian  C_f
	jac[8,5] = -C_f*d5f*decaydiffval
	jac[8,8] = -d5f*decayval

	# Ninth row of jacobian  C_s
	jac[9,5] = -C_s*d5s*decaydiffval
	jac[9,9] = -d5s*decayval

	# Tenth row of jacobian  C_p
	jac[10,5] = -C_p*d5p*decaydiffval
	jac[10,10] = -d5p*decayval
    return 
end

function decay_equation(T_L, Q10, B5, B10)
	# change of carbon decay rates due to temperature
	if T_L <= 1.0
		val = 0.0
	else
		val = Q10^((T_L - 1.0 - B5) / B10)
		if T_L < 1.0 + B5
			val *= 0.5*(1 - cos((T_L-1.0)*π/B5))
		end
	end
	return val
end

function c_W(T, S, c_1, c_2, B_c, L)
	# The specific heat of water taking into account latent heat of freezing.
    if T - S <= -B_c
        return c_1
    elseif B_c <= T - S
        return c_2
    else
		x = pi*(T - S)/(2*B_c)
		return c_1 + 0.5*(c_2 - c_1) * (1.0 + sin(x)) + ((c_2 - c_1) * (T - S) + L) * (pi / (4*B_c)) * cos(x)
    end
end

function k_W(T, S, k_W1, k_W2, B_c)
	# The heat transfer coefficient between the water layers
    if T - S <= -B_c
        return k_W1
    elseif B_c <= T - S
        return k_W2
    else
		return 0.5*((k_W1 + k_W2) + (k_W2 - k_W1)*sin(pi*(T - S)/(2*B_c)))
    end
end

function alpha_x(T, S, alpha_xw, alpha_xc, B_alpha)
	# The albedo.
    return 0.5 * (alpha_xw + alpha_xc) + 0.5 * (alpha_xw - alpha_xc) * tanh((T - S) / B_alpha)
end

function CC(T_A,G_W1)
	# The Clausius-Clapeyron equation for absorption due to water vapour
	return (1.0 / T_A) * exp(G_W1 * (1.0 - 1.0 / T_A))
end

function eta(T_A, mu, nu, G_C, G_CH4, delta, G_W2, G_W1, eta_Cl)
	# absorption of longwave radiation
	return  1.0 - exp(-mu*G_C - nu*G_CH4 - delta*G_W2*CC(T_A,G_W1)) * (1 - eta_Cl)
end

function F_C(T, A1, A2)
	# sensible heat (evapotranspiration and thermals)
	val = A1 * (T - 1.0)
	return val + sqrt(val^2 + A2^2)
end

end
