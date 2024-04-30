from scipy import optimize

def f(x, a):
    return x**2 - a*x + 4

cg = optimize.minimize(f, x0=0, args=(4), method="cg")
print("cg", cg, '\n', cg.x, cg.fun)
# print(round(cg.x[0],2))
