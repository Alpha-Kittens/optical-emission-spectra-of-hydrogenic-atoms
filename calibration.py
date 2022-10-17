from models_2 import quadratic, quadratic_err, inverse_quadratic
calibration = lambda x : quadratic(x, 1.055593492700244e-06,0.9985783005925553,14.076802380736353)
calibration_error = lambda x, x_err : quadratic_err(x,x_err,1.055593492700244e-06,4.2998702311272066e-08,0.9985783005925553,0.00041228204819659123,14.076802380736353,0.9581444594453459)
uncalibration = lambda true : inverse_quadratic(true,1.055593492700244e-06,0.9985783005925553,14.076802380736353)
