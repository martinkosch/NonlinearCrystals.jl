using NonlinearCrystals
using Test
using Unitful

@test typeof(NonlinearCrystals.create_cga()) <: UnidirectionalCrystal

@test CGA.n_X_principal(default_lambda(CGA), default_temp(CGA)) > 1
@test CGA.n_Y_principal(default_lambda(CGA), default_temp(CGA)) > 1
@test CGA.n_Z_principal(default_lambda(CGA), default_temp(CGA)) > 1