from models_2 import quadratic, quadratic_err, inverse_quadratic
calibration = lambda x : quadratic(x, 5.71510100787135e-06,0.9565998648050038,101.6381072097872)
calibration_error = lambda x, x_err : quadratic_err(x,x_err,5.71510100787135e-06,4.069501834853449e-06,0.9565998648050038,0.03874379703859934,101.6381072097872,89.36358988785692)
uncalibration = lambda true : inverse_quadratic(true,5.71510100787135e-06,0.9565998648050038,101.6381072097872)
