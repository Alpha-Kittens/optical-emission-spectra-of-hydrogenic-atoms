from models_2 import quadratic, quadratic_err
calibration = lambda x : quadratic(x, 1.084383924470547e-06,0.9982938422265092,14.753637496936658)
calibration_error = lambda x, x_err : quadratic_err(x,x_err,1.084383924470547e-06,3.026805352848886e-08,0.9982938422265092,0.00028897077306917746,14.753637496936658,0.665947001740797)
