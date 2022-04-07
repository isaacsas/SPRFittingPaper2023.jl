# unit conversions 

# # transforms RU to μM (is this really always right?)
# RU_to_muM(x) = x/149/26795*1000000 

# # transforms RU to (nm)⁻³
# RU_to_inv_cubic_nm(x) = x * (.6023/149/26795)

# transforms μM to (nm)⁻³ and vice versa
muM_to_inv_cubic_nm(x) = 6.023e-7 * x
inv_cubic_nm_to_muM(x) = x / 6.023e-7
