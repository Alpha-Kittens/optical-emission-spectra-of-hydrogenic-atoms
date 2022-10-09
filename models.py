
def linear(x, m, b):
    """
    linear model y=mx+b

    Arguments: 
        * `x` : parameter value to evaluate the function at
        * `m` : slope of the line
        * `b` : y-intercept of the line

    Returns: 
        * `y` result of y=mx+b
    """
    return (m*x)+b

def quadratic(x, a, b, c):
    """
    quadratic model y=ax^2 + bx + c

    Arguments: 
        * `x` : parameter value to evaluate the function at
        * `a` : leading coefficient
        * `b` : 2nd term coefficient
        * `c` : y-interecept

    Returns: 
        * `y` result of y=ax^2 + bx + c
    """
    return (a*(x**2)) + (b*x) + c
