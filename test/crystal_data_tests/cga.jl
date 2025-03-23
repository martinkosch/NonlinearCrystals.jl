using NonlinearCrystals
using Test
using Unitful

@test CGA.n_x_principal(default_lambda(CGA), default_temp(CGA)) > 1
@test CGA.n_y_principal(default_lambda(CGA), default_temp(CGA)) > 1
@test CGA.n_z_principal(default_lambda(CGA), default_temp(CGA)) > 1